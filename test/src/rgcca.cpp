// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>

#include "../ipopt_options.h"

using namespace fdapde;
using namespace fdapde::rgcca;
using fdapde::test::almost_equal;

namespace {
constexpr int n_blocks = 4;
constexpr int n_comp = 3;
constexpr int n_obs = 200;
constexpr int n_locs = 101;
constexpr double lambda_weights = 1e-3;

enum class FunctionalDiscretization { MV, FEM, Splines };

Eigen::Matrix<double, Dynamic, Dynamic> load_market_matrix(const std::string& path) {
    Eigen::SparseMatrix<double> matrix;
    Eigen::loadMarket(matrix, path);
    return Eigen::Matrix<double, Dynamic, Dynamic>(matrix);
}

Eigen::Matrix<double, Dynamic, Dynamic> as_double_matrix(const Eigen::Matrix<bool, Dynamic, Dynamic>& matrix) {
    Eigen::Matrix<double, Dynamic, Dynamic> out(matrix.rows(), matrix.cols());
    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            out(i, j) = matrix(i, j) ? 1.0 : 0.0;
        }
    }
    return out;
}

void write_market_matrix(const Eigen::Matrix<double, Dynamic, Dynamic>& matrix, const std::string& path) {
    std::filesystem::create_directories(std::filesystem::path(path).parent_path());

    int nonzeros = 0;
    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            if (matrix(i, j) != 0.0) ++nonzeros;
        }
    }

    std::ofstream file(path);
    file << "%%MatrixMarket matrix coordinate real general\n";
    file << matrix.rows() << " " << matrix.cols() << " " << nonzeros << "\n";
    file << std::setprecision(17);
    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            if (matrix(i, j) != 0.0) file << i + 1 << " " << j + 1 << " " << matrix(i, j) << "\n";
        }
    }
}

bool update_functional_references() {
    const char* flag = std::getenv("FDAPDE_UPDATE_RGCCA_FUNCTIONAL_REFERENCES");
    bool update = flag != nullptr && std::string(flag) == "1";
    if (update) std::cout << "Updating reference" << std::endl;
    return update;
}

bool almost_equal_with_sign(const Eigen::VectorXd& actual, const Eigen::Matrix<double, Dynamic, Dynamic>& reference, int sign) {
    Eigen::Matrix<double, Dynamic, Dynamic> signed_actual = sign * actual;
    return almost_equal<double>(signed_actual, reference);
}

int sign_against_reference(const Eigen::VectorXd& actual, const Eigen::Matrix<double, Dynamic, Dynamic>& reference) {
    return actual.dot(reference.col(0)) < 0 ? -1 : 1;
}

void connect_reference_design(RGCCA<IndependentSampling>& rgcca) {
    rgcca.connect(0, 1);
    rgcca.connect(0, 2);
    rgcca.connect(0, 3);
    rgcca.connect(1, 3);
    rgcca.connect(2, 3);
}

Eigen::VectorXd centered_wave(const int n, const double frequency) {
    Eigen::VectorXd out(n);
    for (int i = 0; i < n; ++i)
        out[i] = std::sin(frequency * static_cast<double>(i));
    out.array() -= out.mean();
    out.normalize();
    return out;
}

Eigen::Matrix<double, Dynamic, Dynamic> two_column_block(
    const Eigen::VectorXd& first,
    const Eigen::VectorXd& second
) {
    Eigen::Matrix<double, Dynamic, Dynamic> out(first.size(), 2);
    out.col(0) = first;
    out.col(1) = second;
    return out;
}

void check_tau(
  const std::vector<Result>& results, const std::string& reference_path, int n_blocks, int n_comp) {
    Eigen::Matrix<double, Dynamic, Dynamic> actual_tau(n_blocks * n_comp, 1);
    for (int h = 0; h < n_comp; ++h) {
        for (int j = 0; j < n_blocks; ++j) {
            actual_tau(h * n_blocks + j, 0) = results[h].tau_values[j];
        }
    }
    EXPECT_TRUE(almost_equal<double>(actual_tau, load_market_matrix(reference_path)));
}

void write_component_significance_reference(
  const std::vector<Result>& results, const std::string& reference_path, int n_blocks, int n_comp) {
    Eigen::Matrix<double, Dynamic, Dynamic> significance(n_comp, 4);
    Eigen::Matrix<double, Dynamic, Dynamic> active_blocks(n_comp, n_blocks);

    for (int h = 0; h < n_comp; ++h) {
        significance(h, 0) = results[h].rho_tot;
        significance(h, 1) = results[h].rho_tot_p_value;
        significance(h, 2) = static_cast<double>(results[h].rho_tot_bootstrap_count);
        significance(h, 3) = results[h].component_significant ? 1.0 : 0.0;

        for (int j = 0; j < n_blocks; ++j)
            active_blocks(h, j) = results[h].active_blocks[j] ? 1.0 : 0.0;
    }

    write_market_matrix(significance, reference_path + "ref_component_significance.mtx");
    write_market_matrix(active_blocks, reference_path + "ref_active_blocks.mtx");
}

void check_market_matrix(const Eigen::Matrix<double, Dynamic, Dynamic>& actual, const std::string& reference_path) {
    EXPECT_TRUE(almost_equal<double>(actual, load_market_matrix(reference_path)));
}

void check_component_significance_reference(
  const std::vector<Result>& results, const std::string& reference_path, int n_blocks, int n_comp) {
    Eigen::Matrix<double, Dynamic, Dynamic> actual_significance(n_comp, 4);
    Eigen::Matrix<double, Dynamic, Dynamic> actual_active_blocks(n_comp, n_blocks);

    for (int h = 0; h < n_comp; ++h) {
        actual_significance(h, 0) = results[h].rho_tot;
        actual_significance(h, 1) = results[h].rho_tot_p_value;
        actual_significance(h, 2) = static_cast<double>(results[h].rho_tot_bootstrap_count);
        actual_significance(h, 3) = results[h].component_significant ? 1.0 : 0.0;

        EXPECT_TRUE(std::isfinite(results[h].rho_tot_raw));
        EXPECT_EQ(results[h].rho_tot_null_valid_count, results[h].rho_tot_bootstrap_count);
        EXPECT_TRUE(std::isfinite(results[h].rho_tot_null_mean));
        EXPECT_TRUE(std::isfinite(results[h].rho_tot_null_q95));
        EXPECT_TRUE(std::isfinite(results[h].rho_tot_null_max));
        EXPECT_LE(results[h].rho_tot_null_mean, results[h].rho_tot_null_max);
        EXPECT_LE(results[h].rho_tot_null_q95, results[h].rho_tot_null_max);

        for (int j = 0; j < n_blocks; ++j)
            actual_active_blocks(h, j) = results[h].active_blocks[j] ? 1.0 : 0.0;
    }

    check_market_matrix(actual_significance, reference_path + "ref_component_significance.mtx");
    check_market_matrix(actual_active_blocks, reference_path + "ref_active_blocks.mtx");
}

void check_block_component(
  RGCCA<IndependentSampling>::Block& block, const std::string& reference_path, int h) {
    const std::string name = block.name();
    SCOPED_TRACE(name + "_comp" + std::to_string(h + 1));

    const auto weights_ref =
      load_market_matrix(reference_path + "ref_weights_" + name + "_comp" + std::to_string(h + 1) + ".mtx");
    const int sign = sign_against_reference(block.weights_m().col(h), weights_ref);

    EXPECT_TRUE(almost_equal_with_sign(block.weights_m().col(h), weights_ref, sign));

    const auto weights_star_ref =
      load_market_matrix(reference_path + "ref_weights_star_" + name + "_comp" + std::to_string(h + 1) + ".mtx");
    EXPECT_TRUE(almost_equal_with_sign(block.weights_star_m().col(h), weights_star_ref, sign));

    const auto components_ref =
      load_market_matrix(reference_path + "ref_components_" + name + "_comp" + std::to_string(h + 1) + ".mtx");
    EXPECT_TRUE(almost_equal_with_sign(block.components_m().col(h), components_ref, sign));
}

void write_block_component_reference(
  RGCCA<IndependentSampling>::Block& block, const std::string& reference_path, int h) {
    const std::string name = block.name();

    write_market_matrix(
      block.weights_m().col(h), reference_path + "ref_weights_" + name + "_comp" + std::to_string(h + 1) + ".mtx");
    write_market_matrix(
      block.weights_star_m().col(h),
      reference_path + "ref_weights_star_" + name + "_comp" + std::to_string(h + 1) + ".mtx");
    write_market_matrix(
      block.components_m().col(h),
      reference_path + "ref_components_" + name + "_comp" + std::to_string(h + 1) + ".mtx");
}

Eigen::Matrix<double, Dynamic, Dynamic> read_rgcca_block(const std::string& data_path, int block_id) {
    return read_csv<double>(data_path + "X" + std::to_string(block_id) + ".csv").as_matrix();
}

void add_multivariate_blocks(RGCCA<IndependentSampling>& rgcca, const std::string& data_path) {
    for (int i = 1; i <= n_blocks; ++i) {
        Eigen::Matrix<double, Dynamic, Dynamic> X = read_rgcca_block(data_path, i);
        rgcca.add_multivariate_block("X" + std::to_string(i), std::move(X));
    }
}

void add_fem_functional_blocks(RGCCA<IndependentSampling>& rgcca, const std::string& data_path) {
    for (int i = 1; i <= n_blocks; ++i) {
        Eigen::Matrix<double, Dynamic, Dynamic> X = read_rgcca_block(data_path, i);

        Triangulation<1, 1> D = Triangulation<1, 1>::UnitInterval(n_locs);
        GeoFrame data(D);
        auto& level = data.insert_scalar_layer<POINT>("data", MESH_NODES);
        level.load_blk("X" + std::to_string(i), X.transpose());

        FeSpace Vh(D, P1<1>);
        TrialFunction f(Vh);
        TestFunction v(Vh);
        auto a = integral(D)(dot(grad(f), grad(v)));
        ZeroField<1> u;
        auto F = integral(D)(u * v);

        rgcca.add_functional_block(
          "X" + std::to_string(i), data, std::move(X), fe_normcovmax_elliptic(a, F));
    }
}

void add_spline_functional_blocks(RGCCA<IndependentSampling>& rgcca, const std::string& data_path) {
    for (int i = 1; i <= n_blocks; ++i) {
        Eigen::Matrix<double, Dynamic, Dynamic> X = read_rgcca_block(data_path, i);

        Triangulation<1, 1> D = Triangulation<1, 1>::UnitInterval(n_locs);
        GeoFrame data(D);
        auto& level = data.insert_scalar_layer<POINT>("data", MESH_NODES);
        level.load_blk("X" + std::to_string(i), X.transpose());

        BsSpace Bh(D, 3);
        TrialFunction f(Bh);
        TestFunction v(Bh);
        auto a = integral(D)(dxx(f) * dxx(v));
        ZeroField<1> u;
        auto F = integral(D)(u * v);

        rgcca.add_functional_block(
          "X" + std::to_string(i), data, std::move(X), bs_normcovmax_elliptic(a, F));
    }
}

void check_rgcca_against_cran(const std::string& reference_case, Mode mode) {
    const std::string data_path = "../data/models/rgcca/";
    const std::string reference_path = data_path + reference_case + "/";

    RGCCA<IndependentSampling>::Options options;
    options.mode = mode;

    RGCCA<IndependentSampling> rgcca(n_obs, options, n_comp);

    for (int i = 1; i <= n_blocks; ++i) {
        Eigen::Matrix<double, Dynamic, Dynamic> X = read_rgcca_block(data_path, i);
        rgcca.add_multivariate_block("X" + std::to_string(i), std::move(X));
    }

    connect_reference_design(rgcca);

    const auto results = rgcca.fit();
    check_tau(results, reference_path + "ref_tau.mtx", n_blocks, n_comp);

    for (const auto& block : rgcca.blocks()) {
        for (int h = 0; h < n_comp; ++h) {
            check_block_component(*block, reference_path, h);
        }
    }
}

void check_rgcca_against_first_run(
  const std::string& reference_case, FunctionalDiscretization discretization,
  WeightSignConstraint weight_sign_constraint = WeightSignConstraint::None) {
    const std::string data_path = "../data/models/rgcca/";
    const std::string reference_path = data_path + "functional/" + reference_case + "/";

    RGCCA<IndependentSampling>::Options options;
    options.mode = Mode::CovMax;
    options.weight_sign_constraint = weight_sign_constraint;

    RGCCA<IndependentSampling> rgcca(n_obs, options, n_comp);

    switch (discretization) {
        case FunctionalDiscretization::MV: add_multivariate_blocks(rgcca, data_path); break;
        case FunctionalDiscretization::FEM: {
            add_fem_functional_blocks(rgcca, data_path);
            rgcca.set_lambda_weights_all(lambda_weights);
            break;
        }
        case FunctionalDiscretization::Splines: {
            add_spline_functional_blocks(rgcca, data_path);
            rgcca.set_lambda_weights_all(lambda_weights);
            break;
        }
    }

    connect_reference_design(rgcca);

    const auto results = rgcca.fit();

    if (update_functional_references()) {
        for (const auto& block : rgcca.blocks()) {
            for (int h = 0; h < n_comp; ++h) {
                write_block_component_reference(*block, reference_path, h);
            }
        }
        return;
    }

    for (const auto& block : rgcca.blocks()) {
        for (int h = 0; h < n_comp; ++h) {
            check_block_component(*block, reference_path, h);
        }
    }
}

std::vector<double> bootstrap_lambda_grid() {
    return {1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2};
}

std::string comp_suffix(int h) {
    return "_comp" + std::to_string(h + 1);
}

const RGCCA<IndependentSampling>::BootstrapResult* find_bootstrap_result(
  const std::vector<RGCCA<IndependentSampling>::BootstrapResult>& bootstrap_results, int h) {
    for (const auto& result : bootstrap_results) {
        if (result.h == h) return &result;
    }
    return nullptr;
}

void check_compact_bootstrap_storage(const RGCCA<IndependentSampling>::BootstrapResult& component) {
    for (std::size_t i = 0; i < component.lambda_grid.size(); ++i) {
        const int B_used = component.B_used_by_lambda[i];
        EXPECT_EQ(component.corr_boot_by_lambda[i].cols(), B_used);
        EXPECT_EQ(component.candidate_ids_by_lambda[i].size(), static_cast<std::size_t>(B_used));
        for (const auto& weights : component.w_boot_by_lambda[i])
            EXPECT_EQ(weights.cols(), B_used);
    }
}

void write_bootstrap_references(
  const RGCCA<IndependentSampling>& rgcca, const std::vector<Result>& results, const std::string& reference_path) {
    const auto bootstrap_results = rgcca.bootstrap_selection_results();
    ASSERT_EQ(static_cast<int>(bootstrap_results.size()), n_comp);
    ASSERT_EQ(static_cast<int>(results.size()), n_comp);

    for (int h = 0; h < n_comp; ++h) {
        const auto* component = find_bootstrap_result(bootstrap_results, h);
        ASSERT_NE(component, nullptr);
        ASSERT_GE(component->lambda_opt_index, 0);
        check_compact_bootstrap_storage(*component);

        const std::string suffix = comp_suffix(h);
        Eigen::Matrix<double, Dynamic, Dynamic> lambda_opt(1, 1);
        lambda_opt(0, 0) = component->lambda_opt;
        write_market_matrix(lambda_opt, reference_path + "ref_lambda_opt" + suffix + ".mtx");

        write_market_matrix(as_double_matrix(results[h].C), reference_path + "ref_connections" + suffix + ".mtx");

        const int lambda_i = component->lambda_opt_index;
        for (int j = 0; j < n_blocks; ++j) {
            const std::string block_name = component->block_names[j];
            write_market_matrix(
              component->w_min_by_lambda[lambda_i][j],
              reference_path + "ref_wmin_" + block_name + suffix + ".mtx");

            const auto [ci_low, ci_high] = rgcca.bootstrap_weights_ci(h, j, rgcca.blocks()[j]->Psi_D());
            write_market_matrix(ci_low, reference_path + "ref_weights_ci_low_" + block_name + suffix + ".mtx");
            write_market_matrix(ci_high, reference_path + "ref_weights_ci_high_" + block_name + suffix + ".mtx");
        }

        write_market_matrix(component->corr_ci_low_by_lambda[lambda_i], reference_path + "ref_corr_ci_low" + suffix + ".mtx");
        write_market_matrix(component->corr_ci_high_by_lambda[lambda_i], reference_path + "ref_corr_ci_high" + suffix + ".mtx");
    }
}

void check_bootstrap_references(
  const RGCCA<IndependentSampling>& rgcca, const std::vector<Result>& results, const std::string& reference_path) {
    const auto bootstrap_results = rgcca.bootstrap_selection_results();
    ASSERT_EQ(static_cast<int>(bootstrap_results.size()), n_comp);
    ASSERT_EQ(static_cast<int>(results.size()), n_comp);

    for (int h = 0; h < n_comp; ++h) {
        const auto* component = find_bootstrap_result(bootstrap_results, h);
        ASSERT_NE(component, nullptr);
        ASSERT_GE(component->lambda_opt_index, 0);
        check_compact_bootstrap_storage(*component);

        const std::string suffix = comp_suffix(h);
        const auto lambda_ref = load_market_matrix(reference_path + "ref_lambda_opt" + suffix + ".mtx");
        ASSERT_EQ(lambda_ref.rows(), 1);
        ASSERT_EQ(lambda_ref.cols(), 1);
        EXPECT_DOUBLE_EQ(component->lambda_opt, lambda_ref(0, 0));

        check_market_matrix(as_double_matrix(results[h].C), reference_path + "ref_connections" + suffix + ".mtx");

        const int lambda_i = component->lambda_opt_index;
        for (int j = 0; j < n_blocks; ++j) {
            const std::string block_name = component->block_names[j];
            SCOPED_TRACE("bootstrap_" + block_name + suffix);

            check_market_matrix(
              component->w_min_by_lambda[lambda_i][j],
              reference_path + "ref_wmin_" + block_name + suffix + ".mtx");

            const auto [ci_low, ci_high] = rgcca.bootstrap_weights_ci(h, j, rgcca.blocks()[j]->Psi_D());
            check_market_matrix(ci_low, reference_path + "ref_weights_ci_low_" + block_name + suffix + ".mtx");
            check_market_matrix(ci_high, reference_path + "ref_weights_ci_high_" + block_name + suffix + ".mtx");
        }

        check_market_matrix(component->corr_ci_low_by_lambda[lambda_i], reference_path + "ref_corr_ci_low" + suffix + ".mtx");
        check_market_matrix(component->corr_ci_high_by_lambda[lambda_i], reference_path + "ref_corr_ci_high" + suffix + ".mtx");
    }
}

void write_multivariate_bootstrap_references(
  const RGCCA<IndependentSampling>& rgcca, const std::vector<Result>& results, const std::string& reference_path) {
    const int n_results = static_cast<int>(results.size());
    ASSERT_EQ(n_results, n_comp);

    for (const auto& block : rgcca.blocks()) {
        for (int h = 0; h < n_results; ++h) {
            write_block_component_reference(*block, reference_path, h);
        }
    }

    const auto bootstrap_results = rgcca.bootstrap_selection_results();
    ASSERT_EQ(static_cast<int>(bootstrap_results.size()), n_results);

    for (int h = 0; h < n_results; ++h) {
        const auto* component = find_bootstrap_result(bootstrap_results, h);
        ASSERT_NE(component, nullptr);
        ASSERT_EQ(component->lambda_opt_index, 0);
        check_compact_bootstrap_storage(*component);
        ASSERT_EQ(static_cast<int>(component->lambda_grid.size()), 1);
        EXPECT_TRUE(std::isnan(component->lambda_grid.front()));
        EXPECT_TRUE(std::isnan(component->lambda_opt));

        const std::string suffix = comp_suffix(h);
        write_market_matrix(as_double_matrix(results[h].C), reference_path + "ref_connections" + suffix + ".mtx");

        const int lambda_i = component->lambda_opt_index;
        for (int j = 0; j < n_blocks; ++j) {
            const std::string block_name = component->block_names[j];
            write_market_matrix(
              component->w_min_by_lambda[lambda_i][j],
              reference_path + "ref_wmin_" + block_name + suffix + ".mtx");

            const auto [ci_low, ci_high] = rgcca.bootstrap_weights_ci(h, j, rgcca.blocks()[j]->Psi_D());
            write_market_matrix(ci_low, reference_path + "ref_weights_ci_low_" + block_name + suffix + ".mtx");
            write_market_matrix(ci_high, reference_path + "ref_weights_ci_high_" + block_name + suffix + ".mtx");
        }

        write_market_matrix(component->corr_ci_low_by_lambda[lambda_i], reference_path + "ref_corr_ci_low" + suffix + ".mtx");
        write_market_matrix(component->corr_ci_high_by_lambda[lambda_i], reference_path + "ref_corr_ci_high" + suffix + ".mtx");
    }
}

void check_multivariate_bootstrap_references(
  const RGCCA<IndependentSampling>& rgcca, const std::vector<Result>& results, const std::string& reference_path) {
    const int n_results = static_cast<int>(results.size());
    ASSERT_EQ(n_results, n_comp);

    for (const auto& block : rgcca.blocks()) {
        for (int h = 0; h < n_results; ++h) {
            check_block_component(*block, reference_path, h);
        }
    }

    const auto bootstrap_results = rgcca.bootstrap_selection_results();
    ASSERT_EQ(static_cast<int>(bootstrap_results.size()), n_results);

    for (int h = 0; h < n_results; ++h) {
        const auto* component = find_bootstrap_result(bootstrap_results, h);
        ASSERT_NE(component, nullptr);
        ASSERT_EQ(component->lambda_opt_index, 0);
        check_compact_bootstrap_storage(*component);
        ASSERT_EQ(static_cast<int>(component->lambda_grid.size()), 1);
        EXPECT_TRUE(std::isnan(component->lambda_grid.front()));
        EXPECT_TRUE(std::isnan(component->lambda_opt));

        const std::string suffix = comp_suffix(h);
        check_market_matrix(as_double_matrix(results[h].C), reference_path + "ref_connections" + suffix + ".mtx");

        const int lambda_i = component->lambda_opt_index;
        for (int j = 0; j < n_blocks; ++j) {
            const std::string block_name = component->block_names[j];
            SCOPED_TRACE("mv_bootstrap_" + block_name + suffix);

            check_market_matrix(
              component->w_min_by_lambda[lambda_i][j],
              reference_path + "ref_wmin_" + block_name + suffix + ".mtx");

            const auto [ci_low, ci_high] = rgcca.bootstrap_weights_ci(h, j, rgcca.blocks()[j]->Psi_D());
            check_market_matrix(ci_low, reference_path + "ref_weights_ci_low_" + block_name + suffix + ".mtx");
            check_market_matrix(ci_high, reference_path + "ref_weights_ci_high_" + block_name + suffix + ".mtx");
        }

        check_market_matrix(component->corr_ci_low_by_lambda[lambda_i], reference_path + "ref_corr_ci_low" + suffix + ".mtx");
        check_market_matrix(component->corr_ci_high_by_lambda[lambda_i], reference_path + "ref_corr_ci_high" + suffix + ".mtx");
    }
}

void check_fem_bootstrap_rgcca_against_first_run() {
    const std::string data_path = "../data/models/rgcca/";
    const std::string reference_path = data_path + "functional/fem_bootstrap_cov/";

    RGCCA<IndependentSampling>::Options options;
    options.mode = Mode::CovMax;
    options.lambda_selection_weights = LambdaSelection::Automatic;
    options.block_deactivation = true;
    options.connection_deactivation = true;

    RGCCA<IndependentSampling> rgcca(n_obs, options, n_comp);
    add_fem_functional_blocks(rgcca, data_path);
    connect_reference_design(rgcca);
    rgcca.set_lambda_grid_weights(bootstrap_lambda_grid());

    RGCCA<IndependentSampling>::BootstrapConfig bootstrap_config;
    bootstrap_config.max_threads = 1;
    bootstrap_config.check_every = 5*12;
    bootstrap_config.patience = 1;
    rgcca.set_bootstrap_config(bootstrap_config);

    const auto results = rgcca.fit();

    if (update_functional_references()) {
        for (const auto& block : rgcca.blocks()) {
            for (int h = 0; h < n_comp; ++h) {
                write_block_component_reference(*block, reference_path, h);
            }
        }
        write_bootstrap_references(rgcca, results, reference_path);
        return;
    }

    for (const auto& block : rgcca.blocks()) {
        for (int h = 0; h < n_comp; ++h) {
            check_block_component(*block, reference_path, h);
        }
    }
    check_bootstrap_references(rgcca, results, reference_path);
}

void check_multivariate_model_selection_without_lambda_grid() {
    const std::string data_path = "../data/models/rgcca/";
    const std::string reference_path = data_path + "bootstrap/mv_deactivation_cov/";

    RGCCA<IndependentSampling>::Options options;
    options.mode = Mode::CovMax;
    options.block_deactivation = true;

    RGCCA<IndependentSampling> rgcca(n_obs, options, n_comp);
    add_multivariate_blocks(rgcca, data_path);
    connect_reference_design(rgcca);

    RGCCA<IndependentSampling>::BootstrapConfig bootstrap_config;
    bootstrap_config.max_threads = 1;
    bootstrap_config.check_every = 5*12;
    bootstrap_config.patience = 1;
    rgcca.set_bootstrap_config(bootstrap_config);

    const auto results = rgcca.fit();

    if (update_functional_references()) {
        write_multivariate_bootstrap_references(rgcca, results, reference_path);
        return;
    }

    check_multivariate_bootstrap_references(rgcca, results, reference_path);
}

void check_fem_cov_component_significance() {
    const std::string data_path = "../data/models/rgcca/";
    const std::string reference_path = data_path + "functional/fem_cov_component_significance/";

    RGCCA<IndependentSampling>::Options options;
    options.mode = Mode::CovMax;
    options.component_significance = true;

    RGCCA<IndependentSampling> rgcca(n_obs, options, n_comp);
    add_fem_functional_blocks(rgcca, data_path);
    rgcca.set_lambda_weights_all(lambda_weights);
    connect_reference_design(rgcca);

    RGCCA<IndependentSampling>::BootstrapConfig bootstrap_config;
    bootstrap_config.max_threads = 1;
    rgcca.set_bootstrap_config(bootstrap_config);

    const auto results = rgcca.fit();

    ASSERT_EQ(static_cast<int>(results.size()), n_comp);

    if (update_functional_references()) {
        for (const auto& block : rgcca.blocks()) {
            for (int h = 0; h < n_comp; ++h) {
                write_block_component_reference(*block, reference_path, h);
            }
        }
        write_component_significance_reference(results, reference_path, n_blocks, n_comp);
        return;
    }

    for (const auto& block : rgcca.blocks()) {
        for (int h = 0; h < n_comp; ++h) {
            check_block_component(*block, reference_path, h);
        }
    }
    check_component_significance_reference(results, reference_path, n_blocks, n_comp);
}
}   // namespace

TEST(rgcca, R_GCCA_cov) {
    check_rgcca_against_cran("cov", Mode::CovMax);
}

TEST(rgcca, R_GCCA_cor) {
    check_rgcca_against_cran("cor", Mode::CorMax);
}

TEST(rgcca, R_RGCCA) {
    check_rgcca_against_cran("rgcca", Mode::Regularized);
}

TEST(rgcca, connection_uncertainty_retains_threshold_overlap) {
    const double z = fdapde::internals::standard_normal_quantile(0.975);
    EXPECT_NEAR(z, 1.959963984540054, 1e-12);

    EXPECT_GE(fdapde::internals::wilson_score_upper_bound(110, 120, z), 0.95);
    EXPECT_LT(fdapde::internals::wilson_score_upper_bound(109, 120, z), 0.95);

    std::vector<double> overlapping(120, 0.04);
    std::fill(overlapping.begin() + 65, overlapping.end(), 0.06);
    EXPECT_GE(fdapde::internals::median_confidence_upper_bound(overlapping, z), 0.05);

    std::vector<double> clearly_weak(120, 0.04);
    EXPECT_LT(fdapde::internals::median_confidence_upper_bound(clearly_weak, z), 0.05);
}

TEST(rgcca, nonnegative_weight_solver_finds_active_boundary_optimum) {
    write_ipopt_options();

    Eigen::SparseMatrix<double> Psi(2, 2);
    Psi.setIdentity();

    Eigen::SparseMatrix<double> Omega(2, 2);
    Omega.insert(0, 0) = 2.0;
    Omega.insert(0, 1) = 1.0;
    Omega.insert(1, 0) = 1.0;
    Omega.insert(1, 1) = 2.0;
    Omega.makeCompressed();

    Eigen::Vector2d z;
    z << 1.0, -0.25;

    ::fdapde::internals::NonNegativeWeightSolver solver(Psi, Omega, false);
    const Eigen::VectorXd weights = solver.solve(z);

    EXPECT_GE(weights.minCoeff(), 0.0);
    EXPECT_NEAR(weights[0], 1.0 / std::sqrt(2.0), 1e-9);
    EXPECT_NEAR(weights[1], 0.0, 1e-12);
    EXPECT_NEAR(weights.dot(Omega * weights), 1.0, 1e-10);
    EXPECT_NEAR(z.dot(weights), 1.0 / std::sqrt(2.0), 1e-9);
}

TEST(rgcca, GCCA_NN_cov) {

    write_ipopt_options();

    check_rgcca_against_first_run("nn_cov",  FunctionalDiscretization::MV, WeightSignConstraint::NonNegative);
}

TEST(rgcca, F_GCCA_fem_cov) {
    check_rgcca_against_first_run("fem_cov", FunctionalDiscretization::FEM);
}

TEST(rgcca, F_GCCA_splines_cov) {
    check_rgcca_against_first_run("splines_cov", FunctionalDiscretization::Splines);
}

TEST(rgcca, F_GCCA_NN_fem_cov) {

    write_ipopt_options();

    check_rgcca_against_first_run(
      "fem_nn_cov", FunctionalDiscretization::FEM, WeightSignConstraint::NonNegative);
}

TEST(rgcca, F_GCCA_NN_splines_cov) {

    write_ipopt_options();

    check_rgcca_against_first_run(
      "splines_nn_cov", FunctionalDiscretization::Splines, WeightSignConstraint::NonNegative);
}

TEST(rgcca, F_GCCA_fem_bootstrap_cov) {
    check_fem_bootstrap_rgcca_against_first_run();
}

TEST(rgcca, GCCA_bootstrap_model_selection_without_lambda_grid) {
    check_multivariate_model_selection_without_lambda_grid();
}

TEST(rgcca, inactive_design_stops_remaining_components_without_significance) {
    constexpr int n = 4;
    constexpr int n_comp_local = 3;

    RGCCA<IndependentSampling>::Options options;
    options.mode = Mode::CovMax;
    options.component_significance = false;

    RGCCA<IndependentSampling> rgcca(n, options, n_comp_local);
    Eigen::Matrix<double, Dynamic, Dynamic> X1(n, 2);
    Eigen::Matrix<double, Dynamic, Dynamic> X2(n, 2);
    X1 << 1.0, 0.0,
          0.0, 1.0,
          1.0, 1.0,
          0.0, 0.0;
    X2 << 0.0, 1.0,
          1.0, 0.0,
          1.0, 1.0,
          0.0, 0.0;
    rgcca.add_multivariate_block("X1", std::move(X1));
    rgcca.add_multivariate_block("X2", std::move(X2));
    rgcca.connect(0, 1, false);

    int callback_count = 0;
    const auto results = rgcca.fit([&](auto&, const Result& result) {
        EXPECT_EQ(result.h, callback_count);
        ++callback_count;
    });

    EXPECT_EQ(callback_count, 1);
    ASSERT_EQ(static_cast<int>(results.size()), 1);
    EXPECT_EQ(rgcca.n_comp_effective(), 1);
    for (const auto& result : results) {
        EXPECT_FALSE(result.component_significant);
        EXPECT_DOUBLE_EQ(result.rho_tot, 0.0);
        EXPECT_DOUBLE_EQ(result.inner_ave, 0.0);
        EXPECT_DOUBLE_EQ(result.rho_tot_p_value, 1.0);
        for (bool active : result.active_blocks)
            EXPECT_FALSE(active);
    }
}

TEST(rgcca, component_diagnostics_track_correlation_and_deflated_variance) {
    constexpr int n = 80;
    constexpr int n_comp_local = 2;

    const Eigen::VectorXd signal = centered_wave(n, 0.17);
    Eigen::VectorXd second = centered_wave(n, 0.71);
    second -= signal * signal.dot(second);
    second.normalize();

    RGCCA<IndependentSampling>::Options options;
    options.mode = Mode::CovMax;

    RGCCA<IndependentSampling> rgcca(n, options, n_comp_local);
    rgcca.add_multivariate_block("X1", two_column_block(4.0 * signal, 2.0 * second));
    rgcca.add_multivariate_block("X2", two_column_block(3.0 * signal, 1.0 * second));
    rgcca.connect(0, 1);

    const auto results = rgcca.fit();

    ASSERT_EQ(static_cast<int>(results.size()), n_comp_local);
    for (int h = 0; h < n_comp_local; ++h) {
        const auto& result = results[h];
        EXPECT_TRUE(std::isfinite(result.rho_tot));
        EXPECT_TRUE(std::isfinite(result.rho_tot_raw));
        EXPECT_TRUE(std::isfinite(result.inner_ave));
        EXPECT_NEAR(result.rho_tot_raw, result.rho_tot, 1e-12);
        EXPECT_NEAR(result.inner_ave, result.rho_tot * result.rho_tot, 1e-12);
        EXPECT_GE(result.inner_ave, 0.0);
        EXPECT_LE(result.inner_ave, 1.0);

        for (int j = 0; j < 2; ++j) {
            ASSERT_GT(result.block_variance_initial[j], 0.0);
            EXPECT_LE(result.block_variance_after[j], result.block_variance_before[j] + 1e-12);
            EXPECT_NEAR(
                result.block_variance_explained[j],
                (result.block_variance_before[j] - result.block_variance_after[j]) /
                    result.block_variance_initial[j],
                1e-12
            );
            EXPECT_NEAR(
                result.block_variance_explained_cumulative[j],
                (result.block_variance_initial[j] - result.block_variance_after[j]) /
                    result.block_variance_initial[j],
                1e-12
            );
        }
    }

    for (int j = 0; j < 2; ++j) {
        EXPECT_NEAR(results[0].block_variance_before[j], results[0].block_variance_initial[j], 1e-12);
        EXPECT_NEAR(results[1].block_variance_before[j], results[0].block_variance_after[j], 1e-12);
        EXPECT_GE(
            results[1].block_variance_explained_cumulative[j] + 1e-12,
            results[0].block_variance_explained_cumulative[j]
        );
    }
}

TEST(rgcca, block_importance_detects_current_component_blocks) {
    constexpr int n = 80;

    const Eigen::VectorXd signal = centered_wave(n, 0.17);
    Eigen::VectorXd noise = centered_wave(n, 0.71);
    noise -= signal * signal.dot(noise);
    noise.normalize();

    RGCCA<IndependentSampling>::Options options;
    options.mode = Mode::CovMax;
    options.block_importance = true;

    RGCCA<IndependentSampling> rgcca(n, options, 1);
    rgcca.add_multivariate_block("signal_1", two_column_block(10.0 * signal, Eigen::VectorXd::Zero(n)));
    rgcca.add_multivariate_block("signal_2", two_column_block(10.0 * signal, Eigen::VectorXd::Zero(n)));
    rgcca.add_multivariate_block("noise", two_column_block(10.0 * noise, Eigen::VectorXd::Zero(n)));

    RGCCA<IndependentSampling>::BootstrapConfig bootstrap_config;
    bootstrap_config.max_threads = 1;
    bootstrap_config.block_importance_resamples = 19;
    bootstrap_config.block_importance_alpha = 0.1;
    rgcca.set_bootstrap_config(bootstrap_config);

    const auto results = rgcca.fit();

    ASSERT_EQ(static_cast<int>(results.size()), 1);
    EXPECT_EQ(results[0].block_importance_bootstrap_count, 19);
    ASSERT_EQ(static_cast<int>(results[0].block_importance.size()), 3);
    EXPECT_TRUE(std::isfinite(results[0].block_importance[0]));
    EXPECT_TRUE(std::isfinite(results[0].block_importance[1]));
    EXPECT_TRUE(std::isnan(results[0].block_importance[2]));
    EXPECT_TRUE(std::isnan(results[0].block_importance_p_values[2]));
    EXPECT_TRUE(results[0].block_importance_significant[0]);
    EXPECT_TRUE(results[0].block_importance_significant[1]);
    EXPECT_FALSE(results[0].block_importance_significant[2]);
}

TEST(rgcca, block_importance_tests_each_disconnected_design_component) {
    constexpr int n = 80;

    const Eigen::VectorXd signal_1 = centered_wave(n, 0.17);
    Eigen::VectorXd signal_2 = centered_wave(n, 0.71);
    signal_2 -= signal_1 * signal_1.dot(signal_2);
    signal_2.normalize();
    Eigen::VectorXd isolated = centered_wave(n, 1.31);
    isolated -= signal_1 * signal_1.dot(isolated);
    isolated -= signal_2 * signal_2.dot(isolated);
    isolated.normalize();

    RGCCA<IndependentSampling>::Options options;
    options.mode = Mode::CovMax;
    options.block_importance = true;

    RGCCA<IndependentSampling> rgcca(n, options, 1);
    rgcca.add_multivariate_block("pair_1a", two_column_block(10.0 * signal_1, Eigen::VectorXd::Zero(n)));
    rgcca.add_multivariate_block("pair_2a", two_column_block(8.0 * signal_2, Eigen::VectorXd::Zero(n)));
    rgcca.add_multivariate_block("pair_1b", two_column_block(10.0 * signal_1, Eigen::VectorXd::Zero(n)));
    rgcca.add_multivariate_block("pair_2b", two_column_block(8.0 * signal_2, Eigen::VectorXd::Zero(n)));
    rgcca.add_multivariate_block("isolated", two_column_block(6.0 * isolated, Eigen::VectorXd::Zero(n)));
    rgcca.connect(0, 2);
    rgcca.connect(1, 3);

    RGCCA<IndependentSampling>::BootstrapConfig bootstrap_config;
    bootstrap_config.max_threads = 1;
    bootstrap_config.block_importance_resamples = 19;
    bootstrap_config.block_importance_alpha = 0.1;
    rgcca.set_bootstrap_config(bootstrap_config);

    const auto results = rgcca.fit();

    ASSERT_EQ(static_cast<int>(results.size()), 1);
    ASSERT_EQ(static_cast<int>(results[0].block_importance.size()), 5);
    for (int j = 0; j < 4; ++j) {
        EXPECT_NEAR(results[0].block_importance[j], 1.0, 1e-12);
        EXPECT_NEAR(results[0].block_importance_p_values[j], 0.05, 1e-12);
        EXPECT_TRUE(results[0].block_importance_significant[j]);
    }
    EXPECT_TRUE(std::isnan(results[0].block_importance[4]));
    EXPECT_TRUE(std::isnan(results[0].block_importance_p_values[4]));
    EXPECT_FALSE(results[0].block_importance_significant[4]);
}

TEST(rgcca, component_callback_can_release_bootstrap_results) {
    const std::string data_path = "../data/models/rgcca/";
    constexpr int n_comp_local = 2;

    RGCCA<IndependentSampling>::Options options;
    options.mode = Mode::CovMax;
    options.block_deactivation = true;
    options.max_iter = 2;

    RGCCA<IndependentSampling> rgcca(n_obs, options, n_comp_local);
    add_multivariate_blocks(rgcca, data_path);
    connect_reference_design(rgcca);

    RGCCA<IndependentSampling>::BootstrapConfig bootstrap_config;
    bootstrap_config.max_threads = 1;
    bootstrap_config.B_min = 2;
    bootstrap_config.B_max = 2;
    bootstrap_config.adaptive = false;
    bootstrap_config.fit_max_iter = 2;
    rgcca.set_bootstrap_config(bootstrap_config);

    int callback_count = 0;
    const auto results = rgcca.fit([&](auto& model, const Result& result) {
        EXPECT_EQ(result.h, callback_count);
        EXPECT_FALSE(model.bootstrap_selection_results().empty());
        model.clear_bootstrap_selection_results();
        ++callback_count;
    });

    EXPECT_EQ(callback_count, n_comp_local);
    EXPECT_TRUE(rgcca.bootstrap_selection_results().empty());
    ASSERT_EQ(static_cast<int>(results.size()), n_comp_local);
}

TEST(rgcca, F_GCCA_fem_cov_component_significance) {
    check_fem_cov_component_significance();
}

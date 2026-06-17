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

using namespace fdapde;
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

const RGCCA<IndependentSampling>::BootstrapSelectionResult* find_bootstrap_result(
  const std::vector<RGCCA<IndependentSampling>::BootstrapSelectionResult>& bootstrap_results, int h) {
    for (const auto& result : bootstrap_results) {
        if (result.h == h) return &result;
    }
    return nullptr;
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
    bootstrap_config.B_per_thread_per_batch = 5*12;
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
    bootstrap_config.B_per_thread_per_batch = 5*12;
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

TEST(rgcca, GCCA_NN_cov) {

    std::ofstream("ipopt.opt")
    << "print_level 0\n"
    << "sb yes\n"
    << "print_user_options no\n"
    << "print_timing_statistics no\n";

    check_rgcca_against_first_run("nn_cov",  FunctionalDiscretization::MV, WeightSignConstraint::NonNegative);
}

TEST(rgcca, F_GCCA_fem_cov) {
    check_rgcca_against_first_run("fem_cov", FunctionalDiscretization::FEM);
}

TEST(rgcca, F_GCCA_splines_cov) {
    check_rgcca_against_first_run("splines_cov", FunctionalDiscretization::Splines);
}

TEST(rgcca, F_GCCA_NN_fem_cov) {

    std::ofstream("ipopt.opt")
    << "print_level 0\n"
    << "sb yes\n"
    << "print_user_options no\n"
    << "print_timing_statistics no\n";

    check_rgcca_against_first_run(
      "fem_nn_cov", FunctionalDiscretization::FEM, WeightSignConstraint::NonNegative);
}

TEST(rgcca, F_GCCA_NN_splines_cov) {

    std::ofstream("ipopt.opt")
    << "print_level 0\n"
    << "sb yes\n"
    << "print_user_options no\n"
    << "print_timing_statistics no\n";

    check_rgcca_against_first_run(
      "splines_nn_cov", FunctionalDiscretization::Splines, WeightSignConstraint::NonNegative);
}

TEST(rgcca, F_GCCA_fem_bootstrap_cov) {
    check_fem_bootstrap_rgcca_against_first_run();
}

TEST(rgcca, GCCA_bootstrap_model_selection_without_lambda_grid) {
    check_multivariate_model_selection_without_lambda_grid();
}

TEST(rgcca, F_GCCA_fem_cov_component_significance) {
    check_fem_cov_component_significance();
}

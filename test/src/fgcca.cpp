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

void write_tau_reference(
  const std::vector<Result>& results, const std::string& reference_path, int n_blocks, int n_comp) {
    Eigen::Matrix<double, Dynamic, Dynamic> actual_tau(n_blocks * n_comp, 1);
    for (int h = 0; h < n_comp; ++h) {
        for (int j = 0; j < n_blocks; ++j) {
            actual_tau(h * n_blocks + j, 0) = results[h].tau_values[j];
        }
    }
    write_market_matrix(actual_tau, reference_path + "ref_tau.mtx");
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
        write_tau_reference(results, reference_path, n_blocks, n_comp);
        for (const auto& block : rgcca.blocks()) {
            for (int h = 0; h < n_comp; ++h) {
                write_block_component_reference(*block, reference_path, h);
            }
        }
        return;
    }

    check_tau(results, reference_path + "ref_tau.mtx", n_blocks, n_comp);

    for (const auto& block : rgcca.blocks()) {
        for (int h = 0; h < n_comp; ++h) {
            check_block_component(*block, reference_path, h);
        }
    }
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

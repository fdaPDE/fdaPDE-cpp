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

using namespace fdapde;
using fdapde::test::almost_equal;

namespace {
Eigen::Matrix<double, Dynamic, Dynamic> load_market_matrix(const std::string& path) {
    Eigen::SparseMatrix<double> matrix;
    Eigen::loadMarket(matrix, path);
    return Eigen::Matrix<double, Dynamic, Dynamic>(matrix);
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

void check_rgcca_against_cran(const std::string& reference_case, Mode mode) {
    const std::string data_path = "../data/models/rgcca/";
    const std::string reference_path = data_path + reference_case + "/";

    constexpr int n_blocks = 4;
    constexpr int n_comp = 3;
    constexpr int n_obs = 200;

    RGCCA<IndependentSampling>::Options options;
    options.mode = mode;

    RGCCA<IndependentSampling> rgcca(n_obs, options, n_comp);

    for (int i = 1; i <= n_blocks; ++i) {
        Eigen::Matrix<double, Dynamic, Dynamic> X = read_csv<double>(data_path + "X" + std::to_string(i) + ".csv").as_matrix();
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
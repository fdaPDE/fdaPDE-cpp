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

Eigen::RowVectorXd smooth_mean(
  const Eigen::Matrix<double, Dynamic, Dynamic>& X, fdapde::internals::fe_ls_elliptic& smoother, double lambda) {
    smoother.update_response(X.transpose() * Eigen::VectorXd::Ones(X.rows()) / X.rows());
    smoother.fit(lambda);
    return smoother.fn().transpose();
}

Eigen::RowVectorXd smooth_mean(
  const Eigen::Matrix<double, Dynamic, Dynamic>& X, fdapde::internals::fe_ls_elliptic& smoother,
  const std::vector<double>& lambda_grid, int edf_r, int seed) {
    auto gcv = [&](auto lambda) {
        smoother.update_response(X.transpose() * Eigen::VectorXd::Ones(X.rows()) / X.rows());
        smoother.fit(lambda);
        double dor = X.cols() - smoother.edf(lambda, edf_r, seed);
        return (X.cols() / std::pow(dor, 2)) * (smoother.fn() - smoother.response()).squaredNorm();
    };
    GridSearch<1> optimizer;
    Eigen::Matrix<double, 1, 1> lambda = optimizer.optimize(gcv, lambda_grid);
    smoother.update_response(X.transpose() * Eigen::VectorXd::Ones(X.rows()) / X.rows());
    smoother.fit(lambda);
    return smoother.fn().transpose();
}

void check_fpls_case(const std::string& data_path, double lambda) {
    std::string mesh_path = "../data/mesh/unit_square_60/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);

    Eigen::Matrix<double, Dynamic, Dynamic> X = read_csv<double>(data_path + "X.csv").as_matrix();
    Eigen::Matrix<double, Dynamic, Dynamic> Y = read_csv<double>(data_path + "Y.csv").as_matrix();

    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);

    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);

    fdapde::internals::fe_ls_elliptic x_centering;
    x_centering.discretize(fe_ls_elliptic(a, F).get());
    x_centering.analyze_data(data, Eigen::VectorXd::Ones(data[0].rows()).asDiagonal());
    Eigen::RowVectorXd X_mean = smooth_mean(X, x_centering, lambda);
    Eigen::Matrix<double, Dynamic, Dynamic> X_centered = X.rowwise() - X_mean;
    Eigen::RowVectorXd Y_mean = Y.colwise().mean();
    Eigen::Matrix<double, Dynamic, Dynamic> Y_centered = Y.rowwise() - Y_mean;

    l1.load_blk("X", X_centered.transpose());

    fPLS m("X", Y_centered, data, fe_ls_elliptic(a, F), fe_ls_elliptic(a, F));
    Eigen::Matrix<double, 1, 1> lambda_vec;
    lambda_vec << lambda;
    m.fit(3, lambda_vec, lambda_vec, 20, 1e-2);

    EXPECT_TRUE(almost_equal<double>(m.fitted().rowwise() + Y_mean, data_path + "Y_hat.csv"));
    EXPECT_TRUE(almost_equal<double>(m.reconstructed().rowwise() + X_mean, data_path + "X_hat.csv"));
    EXPECT_TRUE(almost_equal<double>(m.B(), data_path + "B_hat.csv"));
    EXPECT_TRUE(almost_equal<double>(m.fitted(3).rowwise() + Y_mean, data_path + "Y_hat.csv"));
    EXPECT_TRUE(almost_equal<double>(m.reconstructed(3).rowwise() + X_mean, data_path + "X_hat.csv"));
    EXPECT_TRUE(almost_equal<double>(m.B(3), data_path + "B_hat.csv"));
    EXPECT_TRUE(almost_equal<double>(m.Beta(), data_path + "B_hat.csv"));
    EXPECT_TRUE(almost_equal<double>(m.Beta(3), data_path + "B_hat.csv"));
    EXPECT_EQ(m.X_latent_scores().rows(), Y.rows());
    EXPECT_EQ(m.X_latent_scores().cols(), 3);
    EXPECT_EQ(m.Y_latent_scores().rows(), Y.rows());
    EXPECT_EQ(m.Y_latent_scores().cols(), 3);
    EXPECT_EQ(m.fitted(1).cols(), Y.cols());
    EXPECT_EQ(m.reconstructed(1).rows(), Y.rows());
    EXPECT_EQ(m.B(1).cols(), Y.cols());
}

void check_fpls_gcv_case(const std::string& data_path) {
    std::string mesh_path = "../data/mesh/unit_square_60/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);

    Eigen::Matrix<double, Dynamic, Dynamic> X = read_csv<double>(data_path + "X.csv").as_matrix();
    Eigen::Matrix<double, Dynamic, Dynamic> Y = read_csv<double>(data_path + "Y.csv").as_matrix();

    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);

    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);

    std::vector<double> lambda_grid(5);
    for (int i = 0; i < 5; ++i) { lambda_grid[i] = std::pow(10, -4.0 + i); }
    int seed = 476813;

    fdapde::internals::fe_ls_elliptic x_centering;
    x_centering.discretize(fe_ls_elliptic(a, F).get());
    x_centering.analyze_data(data, Eigen::VectorXd::Ones(data[0].rows()).asDiagonal());
    Eigen::RowVectorXd X_mean = smooth_mean(X, x_centering, lambda_grid, 1000, seed);
    Eigen::Matrix<double, Dynamic, Dynamic> X_centered = X.rowwise() - X_mean;
    Eigen::RowVectorXd Y_mean = Y.colwise().mean();
    Eigen::Matrix<double, Dynamic, Dynamic> Y_centered = Y.rowwise() - Y_mean;

    l1.load_blk("X", X_centered.transpose());

    fPLS m("X", Y_centered, data, fe_ls_elliptic(a, F), fe_ls_elliptic(a, F));
    m.fit(3, lambda_grid, lambda_grid, OptimizeGCV, 20, 1e-2, 1000, seed);

    EXPECT_TRUE(almost_equal<double>(m.fitted().rowwise() + Y_mean, data_path + "Y_hat.csv"));
    EXPECT_TRUE(almost_equal<double>(m.reconstructed().rowwise() + X_mean, data_path + "X_hat.csv"));
    EXPECT_TRUE(almost_equal<double>(m.B(), data_path + "B_hat.csv"));
    EXPECT_TRUE(almost_equal<double>(m.fitted(3).rowwise() + Y_mean, data_path + "Y_hat.csv"));
    EXPECT_TRUE(almost_equal<double>(m.reconstructed(3).rowwise() + X_mean, data_path + "X_hat.csv"));
    EXPECT_TRUE(almost_equal<double>(m.B(3), data_path + "B_hat.csv"));
    EXPECT_TRUE(almost_equal<double>(m.Beta(), data_path + "B_hat.csv"));
    EXPECT_TRUE(almost_equal<double>(m.Beta(3), data_path + "B_hat.csv"));
    EXPECT_EQ(m.X_latent_scores().rows(), Y.rows());
    EXPECT_EQ(m.X_latent_scores().cols(), 3);
    EXPECT_EQ(m.Y_latent_scores().rows(), Y.rows());
    EXPECT_EQ(m.Y_latent_scores().cols(), 3);
    EXPECT_EQ(m.fitted(1).cols(), Y.cols());
    EXPECT_EQ(m.reconstructed(1).rows(), Y.rows());
    EXPECT_EQ(m.B(1).cols(), Y.cols());
}

void check_restored_modes_smoke(const std::string& data_path, double lambda) {
    std::string mesh_path = "../data/mesh/unit_square_60/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);

    Eigen::Matrix<double, Dynamic, Dynamic> X = read_csv<double>(data_path + "X.csv").as_matrix();
    Eigen::Matrix<double, Dynamic, Dynamic> Y = read_csv<double>(data_path + "Y.csv").as_matrix();

    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);

    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);

    fdapde::internals::fe_ls_elliptic x_centering;
    x_centering.discretize(fe_ls_elliptic(a, F).get());
    x_centering.analyze_data(data, Eigen::VectorXd::Ones(data[0].rows()).asDiagonal());
    Eigen::RowVectorXd X_mean = smooth_mean(X, x_centering, lambda);
    Eigen::Matrix<double, Dynamic, Dynamic> X_centered = X.rowwise() - X_mean;
    Eigen::RowVectorXd Y_mean = Y.colwise().mean();
    Eigen::Matrix<double, Dynamic, Dynamic> Y_centered = Y.rowwise() - Y_mean;

    l1.load_blk("X", X_centered.transpose());

    Eigen::Matrix<double, 1, 1> lambda_vec;
    lambda_vec << lambda;

    fPLS mode_a("X", Y_centered, data, fe_ls_elliptic(a, F), fe_ls_elliptic(a, F), fPLS_A);
    mode_a.fit(3, lambda_vec, lambda_vec, 20, 1e-2);
    EXPECT_EQ(mode_a.mode(), fPLSMode::ModeA);
    EXPECT_EQ(mode_a.fitted().rows(), Y.rows());
    EXPECT_EQ(mode_a.fitted().cols(), Y.cols());
    EXPECT_EQ(mode_a.reconstructed().rows(), X.rows());
    EXPECT_EQ(mode_a.Y_latent_scores().cols(), 3);
    EXPECT_TRUE(mode_a.fitted().array().isFinite().all());
    EXPECT_TRUE(mode_a.reconstructed().array().isFinite().all());

    fPLS mode_sb("X", Y_centered, data, fe_ls_elliptic(a, F), fe_ls_elliptic(a, F), fPLS_SB);
    mode_sb.fit(3, lambda_vec, lambda_vec, 20, 1e-2);
    EXPECT_EQ(mode_sb.mode(), fPLSMode::SymmetricBlock);
    EXPECT_EQ(mode_sb.fitted().rows(), Y.rows());
    EXPECT_EQ(mode_sb.fitted().cols(), Y.cols());
    EXPECT_EQ(mode_sb.reconstructed().rows(), X.rows());
    EXPECT_EQ(mode_sb.Y_latent_scores().cols(), 3);
    EXPECT_TRUE(mode_sb.fitted().array().isFinite().all());
    EXPECT_TRUE(mode_sb.reconstructed().array().isFinite().all());
}

void check_spline_smoke() {
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 1, 21);
    GeoFrame data(T);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);

    const int n_units = 8;
    const int n_locs = data[0].rows();
    Eigen::Matrix<double, Dynamic, Dynamic> X(n_units, n_locs);
    Eigen::Matrix<double, Dynamic, Dynamic> Y(n_units, 2);
    for (int i = 0; i < n_units; ++i) {
        for (int j = 0; j < n_locs; ++j) {
            double x = T.nodes()(j, 0);
            X(i, j) = std::sin((i + 1) * x) + 0.2 * std::cos((j + 1) * x);
        }
        Y(i, 0) = 0.3 * X.row(i).mean() + i * 0.05;
        Y(i, 1) = X(i, n_locs / 2) - X(i, 0);
    }

    Eigen::RowVectorXd X_mean = X.colwise().mean();
    Eigen::Matrix<double, Dynamic, Dynamic> X_centered = X.rowwise() - X_mean;
    Eigen::RowVectorXd Y_mean = Y.colwise().mean();
    Eigen::Matrix<double, Dynamic, Dynamic> Y_centered = Y.rowwise() - Y_mean;
    l1.load_blk("X", X_centered.transpose());

    BsSpace Bh(T, 3);
    TrialFunction f(Bh);
    TestFunction v(Bh);
    auto a = integral(T)(dxx(f) * dxx(v));
    ZeroField<1> u;
    auto F = integral(T)(u * v);

    fPLS m("X", Y_centered, data, bs_ls_elliptic(a, F), bs_ls_elliptic(a, F));
    Eigen::Matrix<double, 1, 1> lambda_vec;
    lambda_vec << 1e-3;
    m.fit(2, lambda_vec, lambda_vec, 20, 1e-6);

    EXPECT_EQ(m.fitted().rows(), Y.rows());
    EXPECT_EQ(m.fitted().cols(), Y.cols());
    EXPECT_TRUE(m.fitted().array().isFinite().all());
    EXPECT_TRUE(m.reconstructed().array().isFinite().all());
    EXPECT_TRUE(m.Beta().array().isFinite().all());
}

}   // namespace

TEST(fpls, test_01) {
    check_fpls_case("../data/models/fpls/2D_test1/", 10.0);
}

TEST(fpls, test_02) {
    check_fpls_gcv_case("../data/models/fpls/2D_test2/");
}

TEST(fpls, restored_modes_smoke) {
    check_restored_modes_smoke("../data/models/fpls/2D_test1/", 10.0);
}

TEST(fpls, spline_smoke) {
    check_spline_smoke();
}

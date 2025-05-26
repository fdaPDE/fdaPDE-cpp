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

// test 1
//    mesh:         unit_square_60
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
TEST(sr_it, test_01) {
    // geometry
    std::string mesh_path = "../data/mesh/unit_square_60/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);
    l1.load_csv<double>("../data/sr/01/response.csv");
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    SRPDE m("y ~ f", data, fe_ls_elliptic_it(a, F));
    m.fit(1.56206e-08);

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/01/field.mtx"));
}

// test 3
//    mesh:         unit_square_60
//    sampling:     locations = nodes
//    penalization: anisotropic diffusion
//    covariates:   no
//    BC:           no
//    order FE:     1
/*
TEST(sr_it, test_03) {
    // geometry
    std::string mesh_path = "../data/mesh/unit_square_60/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);
    l1.load_csv<double>("../data/sr/03/response.csv");
    // physics: anisotropic diffussion
    Eigen::Matrix<double, 2, 2> K;
    K << 1, 0, 0, 4;
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(K * grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    SRPDE m("y ~ f", data, fe_ls_elliptic_it(a, F));
    m.fit(0.002777777777777778);

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/03/field.mtx"));
}
*/

// test 4
//    mesh:         unit_square_21
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid stochastic
TEST(sr_it, test_04) {
    // geometry
    std::string mesh_path = "../data/mesh/unit_square_21/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);
    l1.load_csv<double>("../data/sr/04/response.csv");
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    SRPDE m("y ~ f", data, fe_ls_elliptic_it(a, F));
    // calibration
    std::vector<double> lambda_grid(13);
    for (int i = 0; i < 13; ++i) { lambda_grid[i] = std::pow(10, -6.0 + 0.25 * i) / data[0].rows(); }
    GridOptimizer<1> optimizer;
    optimizer.optimize(m.gcv(100, 476813), lambda_grid);

    EXPECT_TRUE(almost_equal<double>(optimizer.values(), "../data/sr/04/gcvs.mtx"));
}

// areal test
TEST(sr_it, test_11) {
    // geometry
    std::string mesh_path = "../data/mesh/quasi_circle/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POLYGON>("l1", "../data/sr/11/incidence_mat.csv");
    l1.load_csv<double>("../data/sr/11/response.csv");
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    SRPDE m("y ~ f", data, fe_ls_elliptic_it(a, F));
    m.fit(0.0001428571428571429 * data[0].rows());

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/11/field.mtx"));
}

/*
// areal, non constant coefficient PDE, test
TEST(sr_it, test_12) {
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    // geometry
    std::string mesh_path = "../data/mesh/quasi_circle/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POLYGON>("l1", "../data/sr/12/incidence_mat.csv");
    l1.load_csv<double>("../data/sr/12/response.csv");
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    FeCoeff<2, 2, 2, matrix_t> K(read_csv<double>("../data/sr/12/diffusion.csv").as_matrix());
    FeCoeff<2, 2, 1, matrix_t> b(read_csv<double>("../data/sr/12/transport.csv").as_matrix());
    auto a = integral(D)(dot(K * grad(f), grad(v)) + dot(b, grad(f)) * v);
    FeCoeff<2, 1, 1, vector_t> u(read_csv<double>("../data/sr/12/force.csv").as_matrix());
    auto F = integral(D)(u * v);
    // modeling
    SRPDE m("y ~ f", data, fe_ls_elliptic_it(a, F));
    m.fit(0.0001428571428571429);

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/12/field.mtx"));
}
*/

// test on network (test_4-SR-PDE_no_cov_network.R)
TEST(sr_it, test_15) {
    // geometry
    std::string mesh_path = "../data/mesh/network/";
    Triangulation<1, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);
    l1.load_csv<double>("../data/sr/15/response.csv");
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    SRPDE m("y ~ f", data, fe_ls_elliptic_it(a, F));
    m.fit(1e-4 / data[0].rows());

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/15/field.mtx"));
}
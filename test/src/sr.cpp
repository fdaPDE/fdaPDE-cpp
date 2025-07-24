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
TEST(sr, test_01) {
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
    SRPDE m("y ~ f", data, fe_ls_elliptic(a, F));
    m.fit(1.56206e-08);
    
    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/01/field.mtx"));
}

// test 2
//    mesh:         c_shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
TEST(sr, test_02) {
    // geometry
    std::string mesh_path = "../data/mesh/c_shaped/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/sr/02/locs.csv");
    l1.load_csv<double>("../data/sr/02/response.csv");
    l1.load_csv<double>("../data/sr/02/design_matrix.csv");
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    SRPDE m("y ~ x1 + x2 + f", data, fe_ls_elliptic(a, F));
    m.fit(0.001287161988304094);

    EXPECT_TRUE(almost_equal<double>(m.f()   , "../data/sr/02/field.mtx"));
    EXPECT_TRUE(almost_equal<double>(m.beta(), "../data/sr/02/beta.mtx" ));
}

// test 3
//    mesh:         unit_square_60
//    sampling:     locations = nodes
//    penalization: anisotropic diffusion
//    covariates:   no
//    BC:           no
//    order FE:     1
TEST(sr, test_03) {
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
    TestFunction  v(Vh);
    auto a = integral(D)(dot(K * grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    SRPDE m("y ~ f", data, fe_ls_elliptic(a, F));
    m.fit(0.002777777777777778);
    
    EXPECT_TRUE(almost_equal<double>(m.f() , "../data/sr/03/field.mtx"));
}

// test 4
//    mesh:         unit_square_21
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid stochastic
TEST(sr, test_04) {
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
    TestFunction  v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    SRPDE m("y ~ f", data, fe_ls_elliptic(a, F));
    // calibration
    std::vector<double> lambda_grid(13);
    for (int i = 0; i < 13; ++i) { lambda_grid[i] = std::pow(10, -6.0 + 0.25 * i) / data[0].rows(); }
    GridSearch<1> optimizer;
    optimizer.optimize(m.gcv(100, 476813), lambda_grid);

    EXPECT_TRUE(almost_equal<double>(optimizer.values(), "../data/sr/04/gcvs.mtx"));
}

// test 5
//    mesh:         c_shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    GCV optimization: grid stochastic
TEST(sr, test_05) {
    // geometry
    std::string mesh_path = "../data/mesh/c_shaped/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/sr/05/locs.csv");
    l1.load_csv<double>("../data/sr/05/response.csv");
    l1.load_csv<double>("../data/sr/05/design_matrix.csv");
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);    
    // modeling
    SRPDE m("y ~ x1 + x2 + f", data, fe_ls_elliptic(a, F));
    // calibration
    std::vector<double> lambda_grid(25);
    for (int i = 0; i < 25; ++i) { lambda_grid[i] = std::pow(10, -3.0 + 0.25 * i) / data[0].rows(); }
    GridSearch<1> optimizer;
    optimizer.optimize(m.gcv(100, 66546513), lambda_grid);

    EXPECT_TRUE(almost_equal<double>(optimizer.values(), "../data/sr/05/gcvs.mtx"));
}


TEST(sr, test_06) {
    // geometry
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 2, 11);
    std::string mesh_path = "../data/mesh/unit_square_21/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D, T);
    auto& l1 = data.insert_scalar_layer<POINT, POINT>("l1", std::pair {MESH_NODES, MESH_NODES});
    l1.load_csv<double>("../data/sr/06/response.csv");
    // physics
    FeSpace Vh(D, P1<1>);   // linear finite element in space
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a_D = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u_D;
    auto F_D = integral(D)(u_D * v);

    BsSpace Bh(T, 3);   // cubic B-splines in time
    TrialFunction g(Bh);
    TestFunction  w(Bh);
    auto a_T = integral(T)(dxx(g) * dxx(w));
    ZeroField<1> u_T;
    auto F_T = integral(T)(u_T * w);

    SRPDE m("y ~ f", data, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));
    m.fit(2.06143e-06, 2.06143e-06);

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/06/field.mtx"));
}

TEST(sr, test_07) {
    // geometry
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, std::numbers::pi, 5);
    std::string mesh_path = "../data/mesh/c_shaped/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D, T);
    auto& l1 = data.insert_scalar_layer<POINT, POINT>("l1", std::pair {"../data/sr/07/locs.csv", MESH_NODES});
    l1.load_csv<double>("../data/sr/07/response.csv");
    l1.load_csv<double>("../data/sr/07/design_matrix.csv");
    // physics
    FeSpace Vh(D, P1<1>);   // linear finite element in space
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a_D = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u_D;
    auto F_D = integral(D)(u_D * v);

    BsSpace Bh(T, 3);   // cubic B-splines in time
    TrialFunction g(Bh);
    TestFunction  w(Bh);
    auto a_T = integral(T)(dxx(g) * dxx(w));
    ZeroField<1> u_T;
    auto F_T = integral(T)(u_T * w);
    // modeling
    SRPDE m("y ~ x1 + f", data, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));
    m.fit(1.16959e-05, 1.16959e-05);

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/07/field.mtx"));
    EXPECT_TRUE(almost_equal<double>(m.beta(), "../data/sr/07/beta.mtx"));
}

TEST(sr, test_08) {
    // geometry
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 2, 11);
    std::string mesh_path = "../data/mesh/unit_square_21/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D, T);
    auto& l1 = data.insert_scalar_layer<POINT, POINT>("l1", std::pair {MESH_NODES, "../data/sr/08/locs.csv"});
    l1.load_csv<double>("../data/sr/08/response.csv");
    // physics
    FeSpace Vh(D, P1<1>);   // linear finite element in space
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a_D = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u_D;
    auto F_D = integral(D)(u_D * v);

    BsSpace Bh(T, 3);   // cubic B-splines in time
    TrialFunction g(Bh);
    TestFunction  w(Bh);
    auto a_T = integral(T)(dxx(g) * dxx(w));
    ZeroField<1> u_T;
    auto F_T = integral(T)(u_T * w);    
    // modeling
    SRPDE m("y ~ f", data, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));
    m.fit(0.01 / data[0].rows(), 0.01 / data[0].rows());

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/08/field.mtx"));
}

TEST(sr, test_09) {
    // geometry
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 1, 21);
        std::string mesh_path = "../data/mesh/c_shaped/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D, T);
    auto& l1 = data.insert_scalar_layer<POINT, POINT>(
      "l1", std::pair {"../data/sr/09/space_locs.csv", "../data/sr/09/time_locs.csv"});
    l1.load_csv<double>("../data/sr/09/response.csv");
    // physics
    FeSpace Vh(D, P1<1>);   // linear finite element in space
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a_D = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u_D;
    auto F_D = integral(D)(u_D * v);

    BsSpace Bh(T, 3);   // cubic B-splines in time
    TrialFunction g(Bh);
    TestFunction  w(Bh);
    auto a_T = integral(T)(dxx(g) * dxx(w));
    ZeroField<1> u_T;
    auto F_T = integral(T)(u_T * w);    
    // modeling
    SRPDE m("y ~ f", data, fe_ls_separable_mono(std::pair{a_D, F_D}, std::pair {a_T, F_T}));
    m.fit(4.032258064516129e-07, 4.032258064516129e-07);
    
    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/09/field.mtx"));
}

TEST(sr, test_10) {
    // geometry
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 4, 5);
    std::string mesh_path = "../data/mesh/surface/";
    Triangulation<2, 3> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D, T);
    auto& l1 = data.insert_scalar_layer<POINT, POINT>("l1", std::pair {MESH_NODES, MESH_NODES});
    l1.load_csv<double>("../data/sr/10/response.csv");
    // physics
    FeSpace Vh(D, P1<1>);   // linear finite element in space
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a_D = integral(D)(dot(grad(f), grad(v)));
    ZeroField<3> u_D;
    auto F_D = integral(D)(u_D * v);

    BsSpace Bh(T, 3);   // cubic B-splines in time
    TrialFunction g(Bh);
    TestFunction  w(Bh);
    auto a_T = integral(T)(dxx(g) * dxx(w));
    ZeroField<1> u_T;
    auto F_T = integral(T)(u_T * w);    
    // modeling
    SRPDE m("y ~ f", data, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));
    m.fit(5.882352941176471e-08, 5.882352941176471e-08);
    
    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/10/field.mtx"));
}

// areal test
TEST(sr, test_11) {
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
    TestFunction  v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    SRPDE m("y ~ f", data, fe_ls_elliptic(a, F));
    m.fit(0.0001428571428571429);
    
    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/11/field.mtx"));
}

// areal, non constant coefficient PDE, test
TEST(sr, test_12) {
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
    TestFunction  v(Vh);
    FeCoeff<2, 2, 2, matrix_t> K(read_csv<double>("../data/sr/12/diffusion.csv").as_matrix());
    FeCoeff<2, 2, 1, matrix_t> b(read_csv<double>("../data/sr/12/transport.csv").as_matrix());
    auto a = integral(D)(dot(K * grad(f), grad(v)) + dot(b, grad(f)) * v);
    FeCoeff<2, 1, 1, vector_t> u(read_csv<double>("../data/sr/12/force.csv").as_matrix());
    auto F = integral(D)(u * v);
    // modeling
    SRPDE m("y ~ f", data, fe_ls_elliptic(a, F));
    m.fit(0.0001428571428571429);
    
    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/12/field.mtx"));
}

// // parabolic, areal non-constant coefficient
TEST(sr, test_13) {
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;    
    // geometry
    std::string mesh_path = "../data/mesh/quasi_circle/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 3.6, 10);
    // data
    GeoFrame data(D, T);
    auto& l1 =
      data.insert_scalar_layer<POLYGON, POINT>("l1", std::make_pair("../data/sr/13/incidence_mat.csv", MESH_NODES));
    l1.load_csv<double>("../data/sr/13/response.csv");
    vector_t ic = read_csv<double>("../data/sr/13/ic.csv").as_matrix();
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    FeCoeff<2, 2, 2, matrix_t> K(read_csv<double>("../data/sr/13/diffusion.csv").as_matrix());
    FeCoeff<2, 2, 1, matrix_t> b(read_csv<double>("../data/sr/13/transport.csv").as_matrix());
    auto a = integral(D)(dot(K * grad(f), grad(v)) + dot(b, grad(f)) * v);
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    SRPDE m("y ~ f", data, fe_ls_parabolic_mono(std::pair{a, F}, ic));
    m.fit(1e-6 / data[0].rows(), 1e-6);

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/13/field.mtx"));
}

TEST(sr, test_14) {
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;    
    // geometry
    std::string mesh_path = "../data/mesh/unit_square_21/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 1.8, 10);
    // data
    GeoFrame data(D, T);
    auto& l1 = data.insert_scalar_layer<POINT, POINT>("l1", std::make_pair(MESH_NODES, MESH_NODES));
    l1.load_csv<double>("../data/sr/14/response.csv");
    vector_t ic = read_csv<double>("../data/sr/14/ic.csv").as_matrix();
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    SRPDE m("y ~ f", data, fe_ls_parabolic_ieul(std::pair{a, F}, ic, /* max_iter = */ 50, /* tol = */ 1e-4));
    m.fit(0.0002267573696145125, 1.0);
    
    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/14/field.mtx"));
}

// test on network (test_4-SR-PDE_no_cov_network.R)
TEST(sr, test_15) {
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
    TestFunction  v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    SRPDE m("y ~ f", data, fe_ls_elliptic(a, F));
    m.fit(1e-4 / data[0].rows());

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/15/field.mtx"));
}

// test_5-SR-PDE_with_cov_network.R
TEST(sr, test_16) {
    // geometry
    std::string mesh_path = "../data/mesh/network/";
    Triangulation<1, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);
    l1.load_csv<double>("../data/sr/16/response.csv");
    l1.load_csv<double>("../data/sr/16/design_matrix.csv");
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);    
    // modeling
    SRPDE m("y ~ x1 + f", data, fe_ls_elliptic(a, F));
    m.fit(1e-4 / data[0].rows());

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/16/field.mtx"));
    EXPECT_TRUE(almost_equal<double>(m.beta(), "../data/sr/16/beta.mtx"));    
}

// // test_10-SR-PDE_no_cov_3d.R
TEST(sr, test_17) {
    // geometry
    std::string mesh_path = "../data/mesh/unit_sphere/";
    Triangulation<3, 3> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);
    l1.load_csv<double>("../data/sr/17/response.csv");
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<3> u;
    auto F = integral(D)(u * v);
    // modeling
    SRPDE m("y ~ f", data, fe_ls_elliptic(a, F));
    m.fit(1e-4 / data[0].rows());

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/17/field.mtx"));
}

// test_11-SR-PDE_no_cov_manifold
TEST(sr, test_18) {
    std::string mesh_path = "../data/mesh/surface_horseshoe/";
    Triangulation<2, 3> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);
    l1.load_csv<double>("../data/sr/18/response.csv");
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<3> u;
    auto F = integral(D)(u * v);    
    // modeling
    SRPDE m("y ~ f", data, fe_ls_elliptic(a, F));
    m.fit(1e-2 / data[0].rows());

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/18/field.mtx"));
}

// test_12-tSR-PDE_network.R
TEST(sr, test_19) {
    std::string mesh_path = "../data/mesh/network/";
    Triangulation<1, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    Triangulation<1, 1> T = Triangulation<1, 1>::UnitInterval(6);
    // data
    GeoFrame data(D, T);
    auto& l1 = data.insert_scalar_layer<POINT, POINT>("l1", std::tuple {MESH_NODES, MESH_NODES});
    l1.load_csv<double>("../data/sr/19/response.csv");
    l1.load_csv<double>("../data/sr/19/design_matrix.csv");
    // physics
    FeSpace Vh(D, P1<1>);   // linear finite element in space
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a_D = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u_D;
    auto F_D = integral(D)(u_D * v);

    BsSpace Bh(T, 3);   // cubic B-splines in time
    TrialFunction g(Bh);
    TestFunction  w(Bh);
    auto a_T = integral(T)(dxx(g) * dxx(w));
    ZeroField<1> u_T;
    auto F_T = integral(T)(u_T * w);    
    // modeling
    SRPDE m("y ~ x1 + x2 + f", data, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));
    m.fit(1e-4 / data[0].rows(), 1e-4 / data[0].rows());

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/19/field.mtx"));
    EXPECT_TRUE(almost_equal<double>(m.beta(), "../data/sr/19/beta.mtx"));
}

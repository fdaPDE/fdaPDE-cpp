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
//    mesh:         unit_square_40
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    distribution: poisson
TEST(gsr, test_01) {
    // geometry
    std::string mesh_path = "../data/mesh/unit_square_40/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/gsr/01/locs.csv");
    l1.load_csv<double>("../data/gsr/01/response.csv");
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    GSRPDE m("y ~ f", data, fdapde::poisson_distribution(), fe_ls_elliptic(a, F));
    m.fit(/* lambda = */ 1.25e-06);

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/gsr/01/field.mtx"));
}

// test 2
//    mesh:         unit_square_40
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    distribution: bernulli
TEST(gsr, test_02) {
    // geometry
    std::string mesh_path = "../data/mesh/unit_square_40/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/gsr/02/locs.csv");
    l1.load_csv<double>("../data/gsr/02/response.csv");
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    GSRPDE m("y ~ f", data, fdapde::bernoulli_distribution(), fe_ls_elliptic(a, F));
    m.fit(/* lambda = */ 1.25e-06);

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/gsr/02/field.mtx"));
}

// test 3
//    mesh:         unit_square_40
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    distribution: exponential
TEST(gsr, test_03) {
    // geometry
    std::string mesh_path = "../data/mesh/unit_square_40/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/gsr/03/locs.csv");
    l1.load_csv<double>("../data/gsr/03/response.csv");
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    GSRPDE m("y ~ f", data, fdapde::exponential_distribution(), fe_ls_elliptic(a, F));
    m.fit(/* lambda = */ 1.25e-06);
    
    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/gsr/03/field.mtx"));
}

// test 4
//    mesh:         unit_square_40
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    distribution: gamma
TEST(gsr, test_04) {
    // geometry
    std::string mesh_path = "../data/mesh/unit_square_40/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/gsr/04/locs.csv");
    l1.load_csv<double>("../data/gsr/04/response.csv");
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    GSRPDE m("y ~ f", data, fdapde::gamma_distribution(), fe_ls_elliptic(a, F));
    m.fit(/* lambda = */ 1.25e-06);

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/gsr/04/field.mtx"));
}

TEST(gsr, test_05) {
    // geometry
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 1, 4);
    std::string mesh_path = "../data/mesh/c_shaped/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D, T);
    auto& l1 = data.insert_scalar_layer<POINT, POINT>("l1", std::pair {"../data/gsr/05/locs.csv", MESH_NODES});
    l1.load_csv<double>("../data/gsr/05/response.csv");
    l1.load_csv<double>("../data/gsr/05/design_matrix.csv");
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

    GSRPDE m(
      "y ~ x1 + x2 + f", data, fdapde::gamma_distribution(),
      fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));
    m.fit(/* lambda = */ 1.491640405739802e-06 , 1.491640405739802e-06);
    
    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/gsr/05/field.mtx"));
    EXPECT_TRUE(almost_equal<double>(m.beta(), "../data/gsr/05/beta.mtx"));
}

TEST(gsr, test_06) {
    // geometry
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 2. / 3, 3);
    std::string mesh_path = "../data/mesh/c_shaped/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D, T);
    auto& l1 = data.insert_scalar_layer<POINT, POINT>("l1", std::pair {"../data/gsr/06/locs.csv", MESH_NODES});
    l1.load_csv<double>("../data/gsr/06/response.csv");
    l1.load_csv<double>("../data/gsr/06/design_matrix.csv");
    Eigen::Matrix<double, Dynamic, 1> ic = read_csv<double>("../data/gsr/06/ic.csv").as_matrix();
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    GSRPDE m("y ~ x1 + x2 + f", data, fdapde::gamma_distribution(), fe_ls_parabolic_mono(std::pair {a, F}, ic));
    m.fit(/* lambda = */ std::pow(0.1, 2.5) / data[0].rows(), std::pow(0.1, 2.5));

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/gsr/06/field.mtx"));
    EXPECT_TRUE(almost_equal<double>(m.beta(), "../data/gsr/06/beta.mtx"));
}

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
using fdapde::test::read_mesh;
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
    Triangulation<2, 2> D = read_mesh<2, 2>("../data/mesh/unit_square_40");
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/gsr/01/locs.csv");
    l1.load_csv<double>("../data/gsr/01/response.csv");
    // modeling
    GSRPDE m("y ~ f", data, Poisson, fe_laplace());
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
    Triangulation<2, 2> D = read_mesh<2, 2>("../data/mesh/unit_square_40");
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/gsr/02/locs.csv");
    l1.load_csv<double>("../data/gsr/02/response.csv");
    // modeling
    GSRPDE m("y ~ f", data, Bernoulli, fe_laplace());
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
    Triangulation<2, 2> D = read_mesh<2, 2>("../data/mesh/unit_square_40");
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/gsr/03/locs.csv");
    l1.load_csv<double>("../data/gsr/03/response.csv");
    // modeling
    GSRPDE m("y ~ f", data, Exponential, fe_laplace());
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
    Triangulation<2, 2> D = read_mesh<2, 2>("../data/mesh/unit_square_40");
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/gsr/04/locs.csv");
    l1.load_csv<double>("../data/gsr/04/response.csv");
    // modeling
    GSRPDE m("y ~ f", data, Gamma, fe_laplace());
    m.fit(/* lambda = */ 1.25e-06);

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/gsr/04/field.mtx"));
}

TEST(gsr, test_05) {
    // geometry
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 1, 4);
    Triangulation<2, 2> D = read_mesh<2, 2>("../data/mesh/c_shaped");
    // data
    GeoFrame data(D, T);
    auto& l1 = data.insert_scalar_layer<POINT, POINT>("l1", std::pair {"../data/gsr/05/locs.csv", MESH_NODES});
    l1.load_csv<double>("../data/gsr/05/response.csv");
    l1.load_csv<double>("../data/gsr/05/design_matrix.csv");
    // modeling
    BsSpace Vh(T, 3);
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a = integral(T)(dxx(f) * dxx(v));
    ScalarField<1, decltype([](const Eigen::Matrix<double, 1, 1>& p) { return 0; })> u;
    auto F = integral(T)(u * v);

    GSRPDE m("y ~ x1 + x2 + f", data, Gamma, fe_separable(Direct, fe_laplace(), std::pair {a, F}));
    m.fit(/* lambda = */ 1.491640405739802e-06 , 1.491640405739802e-06);
    
    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/gsr/05/field.mtx"));
    EXPECT_TRUE(almost_equal<double>(m.beta(), "../data/gsr/05/beta.mtx"));
}

TEST(gsr, test_06) {
    // geometry
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 2. / 3, 3);
    Triangulation<2, 2> D = read_mesh<2, 2>("../data/mesh/c_shaped");
    // data
    GeoFrame data(D, T);
    auto& l1 = data.insert_scalar_layer<POINT, POINT>("l1", std::pair {"../data/gsr/06/locs.csv", MESH_NODES});
    l1.load_csv<double>("../data/gsr/06/response.csv");
    l1.load_csv<double>("../data/gsr/06/design_matrix.csv");
    Eigen::Matrix<double, Dynamic, 1> ic = read_csv<double>("../data/gsr/06/ic.csv").as_matrix();
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ScalarField<2, decltype([](const Eigen::Matrix<double, 2, 1>& p) { return 0; })> u;
    auto F = integral(D)(u * v);
    // modeling
    GSRPDE m("y ~ x1 + x2 + f", data, Gamma, fe_parabolic(Direct, std::pair {a, F}, ic));
    m.fit(/* lambda = */ std::pow(0.1, 2.5) / data[0].rows(), std::pow(0.1, 2.5));

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/gsr/06/field.mtx"));
    EXPECT_TRUE(almost_equal<double>(m.beta(), "../data/gsr/06/beta.mtx"));
}

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
//    mesh:         unit_square_60
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
TEST(de, test_01) {
    // geometry
    Triangulation<2, 2> D = read_mesh<2, 2>("../data/mesh/square_density");
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/de/01/points.csv");
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ScalarField<2, decltype([](const Eigen::Matrix<double, 2, 1>&) { return 0; })> u;
    auto F = integral(D)(u * v);
    // modeling
    internals::fe_de_elliptic m(data, std::pair{a, F});
    m.set_llik_tolerance(1e-15);

    Eigen::Matrix<double, Dynamic, 1> g_init = read_csv<double>("../data/de/01/f_init.csv").as_matrix().array().log();
    double lambda = 0.1;
    
    m.fit(lambda, g_init, BFGS<Dynamic> {500, 1e-5, 1e-2});

    EXPECT_TRUE(almost_equal<double>(m.log_density(), "../data/de/01/log_density.mtx"));
}

TEST(de, test_02) {
    // geometry
    Triangulation<2, 2> D = read_mesh<2, 2>("../data/mesh/square_density");
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/de/02/points.csv");
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ScalarField<2, decltype([](const Eigen::Matrix<double, 2, 1>&) { return 0; })> u;
    auto F = integral(D)(u * v);
    // modeling
    internals::fe_de_elliptic m(data, std::pair{a, F});
    Eigen::Matrix<double, Dynamic, 1> g_init = read_csv<double>("../data/de/02/f_init.csv").as_matrix().array().log();
    double lambda = 0.1;
    m.fit(lambda, g_init, GradientDescent<Dynamic, BacktrackingLineSearch> {1000, 1e-5, 1e-2});

    EXPECT_TRUE(almost_equal<double>(m.log_density(), "../data/de/02/log_density.mtx"));
}

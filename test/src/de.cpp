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

TEST(de, test_03) {
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;

    // geometry
    Triangulation<2, 2> D = read_mesh<2, 2>("../data/mesh/unit_square_21");
    Triangulation<1, 1> T = Triangulation<1, 1>::UnitInterval(7);

    matrix_t locs_d = read_csv<double>("../data/de/03/data_space.csv").as_matrix();
    matrix_t locs_t = read_csv<double>("../data/de/03/data_time.csv" ).as_matrix();

    matrix_t locs(locs_d.rows(), 3);
    locs.leftCols(2)  = locs_d;
    locs.rightCols(1) = locs_t;
    
    // data
    GeoFrame data(D, T);
    auto& l1 = data.insert_scalar_layer<POINT, POINT>("l1", locs);
    
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a1 = integral(D)(dot(grad(f), grad(v)));
    ScalarField<2, decltype([](const Eigen::Matrix<double, 2, 1>&) { return 0; })> u1;
    auto F1 = integral(D)(u1 * v);

    BsSpace Qh(T, 3);
    TrialFunction g(Qh);
    TestFunction  h(Qh);
    auto a2 = integral(T)(dxx(g) * dxx(h));
    ScalarField<1, decltype([](const Eigen::Matrix<double, 1, 1>& p) { return 0; })> u2;
    auto F2 = integral(T)(u2 * h);

    // modeling
    internals::fe_de_separable m(data, fe_separable(Direct, std::pair {a1, F1}, std::pair {a2, F2}).get());
    Eigen::Matrix<double, Dynamic, 1> g_init = read_csv<double>("../data/de/03/f_init.csv").as_matrix().array().log();
    double lambda_D = 0.00025, lambda_T = 0.01;
    m.set_tol(1e-15);
    m.fit(lambda_D, lambda_T, g_init, BFGS<Dynamic> {100, 1e-5, 1e-2});

    EXPECT_TRUE(almost_equal<double>(m.log_density(), "../data/de/03/log_density.mtx"));
}

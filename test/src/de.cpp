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
TEST(de, test_01) {
    // geometry
    std::string mesh_path = "../data/mesh/square_density/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    Eigen::Matrix<double, Dynamic, 1> g_init = read_csv<double>("../data/de/01/f_init.csv").as_matrix().array().log();
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/de/01/points.csv");
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    DEPDE m(data, fe_de_elliptic(a, F));
    m.set_llik_tolerance(1e-15);
    m.fit(0.1, g_init, BFGS<Dynamic>(500, 1e-5, 1e-2));

    EXPECT_TRUE(almost_equal<double>(m.log_density(), "../data/de/01/log_density.mtx"));
}

TEST(de, test_02) {
    // geometry
    std::string mesh_path = "../data/mesh/square_density/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    Eigen::Matrix<double, Dynamic, 1> g_init = read_csv<double>("../data/de/02/f_init.csv").as_matrix().array().log();
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/de/02/points.csv");
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    DEPDE m(data, fe_de_elliptic(a, F));
    m.fit(0.1, g_init, GradientDescent<Dynamic>(1000, 1e-5, 1e-2), BacktrackingLineSearch());

    EXPECT_TRUE(almost_equal<double>(m.log_density(), "../data/de/02/log_density.mtx"));
}

TEST(de, test_03) {
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    // geometry
    std::string mesh_path = "../data/mesh/unit_square_21/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    Triangulation<1, 1> T = Triangulation<1, 1>::UnitInterval(7);
    // data
    Eigen::Matrix<double, Dynamic, 1> g_init = read_csv<double>("../data/de/03/f_init.csv").as_matrix().array().log();
    matrix_t locs(500, 3);
    locs.leftCols(2)  = read_csv<double>("../data/de/03/data_space.csv").as_matrix();
    locs.rightCols(1) = read_csv<double>("../data/de/03/data_time.csv" ).as_matrix();
    GeoFrame data(D, T);
    auto& l1 = data.insert_scalar_layer<POINT, POINT>("l1", locs);    
    // physics
    FeSpace Vh(D, P1<1>);   // linear finite element in space
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a_D = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u_D;
    auto F_D = integral(D)(u_D * v);

    BsSpace Qh(T, 3);   // cubic B-spline in time
    TrialFunction g(Qh);
    TestFunction  w(Qh);
    auto a_T = integral(T)(dxx(g) * dxx(w));
    ZeroField<1> u_T;
    auto F_T = integral(T)(u_T * w);
    // modeling
    DEPDE m(data, fe_de_separable(std::pair {a_D, F_D}, std::pair {a_T, F_T}));
    m.set_llik_tolerance(1e-15);
    m.fit(0.00025, 0.01, g_init, BFGS<Dynamic>(100, 1e-5, 1e-2));

    EXPECT_TRUE(almost_equal<double>(m.log_density(), "../data/de/03/log_density.mtx"));
}

// TEST(de, test_04) {
//     using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
//     // geometry
//     Triangulation<1, 2> D = read_mesh<1, 2>("../data/mesh/network");
//     Triangulation<1, 1> T = Triangulation<1, 1>::UnitInterval(3);
//     // data
//     Eigen::Matrix<double, Dynamic, 1> g_init = read_csv<double>("../data/de/04/f_init.csv").as_matrix().array().log();
//     matrix_t locs(500, 3);
//     locs.leftCols(2)  = read_csv<double>("../data/de/04/data_space.csv").as_matrix();
//     locs.rightCols(1) = read_csv<double>("../data/de/04/data_time.csv" ).as_matrix();
//     GeoFrame data(D, T);
//     auto& l1 = data.insert_scalar_layer<POINT, POINT>("l1", locs);

//         std::cout << l1 << std::endl;
    
//     // physics
//     FeSpace Vh(D, P1<1>);   // linear finite element in space
//     TrialFunction f(Vh);
//     TestFunction  v(Vh);
//     auto a_D = integral(D)(dot(grad(f), grad(v)));
//     ZeroField<2> u_D;
//     auto F_D = integral(D)(u_D * v);

//     BsSpace Qh(T, 3);   // cubic B-spline in time
//     TrialFunction g(Qh);
//     TestFunction  w(Qh);
//     auto a_T = integral(T)(dxx(g) * dxx(w));
//     ZeroField<1> u_T;
//     auto F_T = integral(T)(u_T * w);
//     // modeling
//     DEPDE m(data, fe_de_separable(std::pair {a_D, F_D}, std::pair {a_T, F_T}));
//     std::cout << "ci arrivo" << std::endl;
//     m.set_llik_tolerance(1e-15);
//     m.fit(1e-7, 1e-2, g_init, BFGS<Dynamic> {500, 1e-5, 1e-2});

//     std::cout << m.log_density().topRows(20) << std::endl;

//     std::cout << "n: " << m.log_density().rows() << std::endl;
    
    
//     // EXPECT_TRUE(almost_equal<double>(m.log_density(), "../data/de/03/log_density.mtx"));
// }


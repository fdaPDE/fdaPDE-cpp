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
//    mesh:         unit_square_21
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
TEST(qsr, test_01) {
    // geometry
    std::string mesh_path = "../data/mesh/unit_square_21/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);
    l1.load_csv<double>("../data/qsr/01/response.csv");
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    QSRPDE m("y ~ f", data, /* alpha = */ 0.1, fe_ls_elliptic(a, F));
    m.fit(/* lambda = */ 1.778279 * std::pow(0.1, 4));

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/qsr/01/field.mtx"));
}

// test 2
//    mesh:         c_shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
// TEST(qsr, test_02) {
//     // geometry
//     Triangulation<2, 2> D = read_mesh<2, 2>("../data/mesh/c_shaped");
//     // data
//     GeoFrame data(D);
//     auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/qsr/02/locs.csv");
//     l1.load_csv<double>("../data/qsr/02/response.csv");
//     l1.load_csv<double>("../data/qsr/02/design_matrix.csv");
//     // modeling
//     QSRPDE m("y ~ x1 + x2 + f", data, /* alpha = */ 0.9, fe_laplace());
//     m.fit(/* lambda = */ 3.162277 * std::pow(0.1, 4));

//     EXPECT_TRUE(almost_equal<double>(m.f(), "../data/qsr/02/field.mtx"));
//     EXPECT_TRUE(almost_equal<double>(m.beta(), "../data/qsr/02/beta.mtx"));
// }

// // test 3
// //    mesh:         unit_square_21
// //    sampling:     locations = nodes
// //    penalization: anisotropic diffusion
// //    covariates:   no
// //    BC:           no
// //    order FE:     1
// TEST(qsr, test_03) {
//     // geometry
//     Triangulation<2, 2> D = read_mesh<2, 2>("../data/mesh/unit_square_21");
//     // data
//     GeoFrame data(D);
//     auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);
//     l1.load_csv<double>("../data/qsr/03/response.csv");
//     // physics: anisotropic diffussion
//     Eigen::Matrix<double, 2, 2> K;
//     K << 1, 0, 0, 4;
//     FeSpace Vh(D, P1<1>);
//     TrialFunction f(Vh);
//     TestFunction  v(Vh);
//     auto a = integral(D)(dot(K * grad(f), grad(v)));
//     // modeling
//     QSRPDE m("y ~ f", data, /* alpha = */ 0.1, fe_elliptic(a));
//     m.fit(/* lambda = */ 5.623413 * pow(0.1, 4));

//     EXPECT_TRUE(almost_equal<double>(m.f(), "../data/qsr/03/field.mtx"));
// }

// // test 4
// //    mesh:         unit_square_21
// //    sampling:     locations = nodes
// //    penalization: simple laplacian
// //    covariates:   no
// //    BC:           no
// //    order FE:     1
// //    GCV optimization: grid stochastic
// TEST(qsr, test_04) {
//     // geometry
//     Triangulation<2, 2> D = read_mesh<2, 2>("../data/mesh/unit_square_21");
//     // data
//     GeoFrame data(D);
//     auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);
//     l1.load_csv<double>("../data/qsr/04/response.csv");
//     // modeling
//     QSRPDE m("y ~ f", data, /* alpha = */ 0.1, fe_laplace());
//     // calibration
//     std::vector<double> lambda_grid(13);
//     for (int i = 0; i < 13; ++i) { lambda_grid[i] = std::pow(10, -8.0 + 0.25 * i); }
//     GridOptimizer<1> optimizer;
//     optimizer.optimize(m.gcv(1000, 438172), lambda_grid);

//     EXPECT_TRUE(almost_equal<double>(optimizer.values(), "../data/qsr/04/gcvs.mtx"));
// }

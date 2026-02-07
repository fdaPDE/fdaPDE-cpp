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
//    BC:           no
//    order FE:     1
//    solver:       subspace
TEST(fpca, test_01) {
    // geometry
    std::string mesh_path = "../data/mesh/unit_square_40/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);
    // load the data matrix (assume that X is an (n_units,n_locs) data matrix)
    std::string data_path = "../data/fpca/01/";
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> X = read_csv<double>(data_path + "y.csv").as_matrix();
    l1.load_blk("X", X.transpose());
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    fPCA m("X", data, fe_ls_elliptic(a, F));

    // fit
    std::vector<double> lambda_grid(5);
    for (int i = 0; i < 5; ++i) { lambda_grid[i] = std::pow(10, -4.0 +  i);}
    m.fit(
        /* n_comp = */ 3,
        lambda_grid,
        /* options = */ ComputeRandSVD | OptimizeGCV,
        fpca_subspace_solver()
        );
    EXPECT_TRUE(almost_equal<double>(m.F().col(0), data_path + "f1.mtx") || almost_equal<double>(-m.F().col(0), data_path + "f1.mtx"));
    EXPECT_TRUE(almost_equal<double>(m.F().col(1), data_path + "f2.mtx") || almost_equal<double>(-m.F().col(1), data_path + "f2.mtx"));
    EXPECT_TRUE(almost_equal<double>(m.F().col(2), data_path + "f3.mtx") || almost_equal<double>(-m.F().col(2), data_path + "f3.mtx"));
}

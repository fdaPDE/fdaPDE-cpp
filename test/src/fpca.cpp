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
    EXPECT_EQ(m.objective_history().size(), 3);
    EXPECT_EQ(m.iterations().size(), 3);
    EXPECT_EQ(m.monotone().size(), 3);
}

// check vector-valued GCV grids for every fPCA solver policy
TEST(fpca, vector_grid_all_solvers) {
    auto D = Triangulation<2, 2>::Rectangle(0, 1, 0, 1, 4, 4);
    GeoFrame data(D);
    auto& layer = data.insert_scalar_layer<POINT>("l1", MESH_NODES);
    Eigen::MatrixXd X(D.n_nodes(), 8);
    for (int i = 0; i < X.rows(); ++i) {
        for (int j = 0; j < X.cols(); ++j) { X(i, j) = std::sin(0.3 * i + 0.5 * j) + std::cos(0.2 * i - 0.7 * j); }
    }
    X = (X.colwise() - X.rowwise().mean()).eval();
    layer.load_blk("X", X);
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    const std::vector<double> grid {0.01, 0.1};
    auto check_policy = [&](auto policy) {
        fPCA model("X", data, fe_ls_elliptic(a, F));
        model.fit(2, grid, ComputeXactSVD | OptimizeGCV, policy);
        // every loading must be finite after searching the vector grid
        EXPECT_TRUE(model.F().array().isFinite().all());
        // every score must be finite after fitting the selected penalties
        EXPECT_TRUE(model.S().array().isFinite().all());
        for (int i = 0; i < model.lambda().size(); ++i) {
            // the selected penalty for each component must belong to the supplied grid
            EXPECT_TRUE(std::find(grid.begin(), grid.end(), model.lambda().data()[i]) != grid.end());
        }
    };
    check_policy(fpca_power_solver());
    check_policy(fpca_subspace_solver());
    check_policy(fpca_direct_solver());
}

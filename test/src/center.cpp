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

TEST(center, test_01) {
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;

    // data
    matrix_t Y = read_csv<double>("../data/models/centering/2D_test1/X.csv").as_matrix();

    // geometry
    std::string mesh_path = "../data/mesh/unit_square_60/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);
    data[0].data().append_blk("Y", Y.transpose());
    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    FRPDE m("y ~ f", data, fe_ls_elliptic(a, F));

    std::vector<double> lambda_grid;
    for (double x = -6.0; x <= 0.0; x += 0.5) { lambda_grid.push_back(std::pow(10, x) / data[0].rows()); }
    GridOptimizer<1> optimizer;
    optimizer.optimize(m.gcv(100, 66546513), lambda_grid);

    m.fit(optimizer.optimum()[0]);

    // std::cout << m.fitted() << std::endl;

    EXPECT_TRUE(almost_equal<double>(m.residuals(), "../data/models/centering/2D_test1/fitted.mtx"));
}
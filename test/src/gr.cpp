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

void set_random(Eigen::Matrix<double, Dynamic, 1>& v, const double sd = 1.0) {
    std::random_device rd;
    std::mt19937 gen(rd()); // Mersenne Twister engine
    std::normal_distribution<> dist(0.0, sd);

    // Fill the vector with random values
    for (int i = 0; i < v.size(); ++i) v(i) = dist(gen);
}


TEST(gr, test_01) {
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;

    // Topology
    const int n_nodes = 10;
    Graph G = Graph::Path(n_nodes);
    GraphTriangulation GT = GraphTriangulation<1>::FromGraphRegularLayout(G);

    // Data
    GeoFrame data(GT);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);
    matrix_t z_gt = vector_t::Ones(n_nodes);
    for (int i = 0; i < n_nodes; ++i) z_gt(i) = 10 * i / (n_nodes - 1.0);
    matrix_t z = z_gt;
    // z(7) = std::numeric_limits<double>::quiet_NaN();
    vector_t noise(n_nodes); set_random(noise, 0.5);
    z += noise;
    std::cout << "Noise:" << std::endl;
    std::cout << std::setw(12) << noise.transpose() << std::endl;
    std::cout << std::endl;
    std::cout << "GT:" << std::endl;
    std::cout << std::setw(12) << z_gt.transpose() << std::endl;
    std::cout << std::endl;
    std::cout << "Data:" << std::endl;
    std::cout << std::setw(12) << z.transpose() << std::endl;
    std::cout << std::endl;
    l1.load_blk("z", z);

    // Physics
    const GraphSpace Gs(G);
    const TrialFunction f(Gs);
    const TestFunction v(Gs);
    const ZeroField<1> u;
    const auto a = integral(G)( lapG(f) * lapG(v) );
    const auto F = integral(G)( u * v );

    // Solver
    SRPDE model("z ~ f", data, ls_graph(a, F));

    // calibration
    std::vector<double> lambda_grid;
    for (double e = -9.; e <= 1.; e+=0.2) lambda_grid.push_back(std::pow(10, e));
    GridSearch<1> optimizer;
    optimizer.optimize(model.gcv(100, 476813), lambda_grid);

    const double lambda_opt = optimizer.optimum()[0];
    model.fit(lambda_opt);

    std::cout << "Fitted model with lambda = " << lambda_opt << std::endl;
    std::cout << std::setw(12) << model.f().transpose() << std::endl;
    std::cout << std::endl;

    std::cout << "Reconstruction error: " << (z_gt - model.f()).norm()/n_nodes << std::endl;
}
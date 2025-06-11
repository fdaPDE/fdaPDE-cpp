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

TEST(fdti_it, test_01) {
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    // geometry
    std::string mesh_path = "../data/mesh/unit_square_60/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/vsr/tensors/locations.csv");
    vector_t S0 = read_csv<double>("../data/vsr/tensors/S0.csv").as_matrix();
    data[0].data().append_blk("S0", S0);
    matrix_t S = read_csv<double>("../data/vsr/tensors/S.csv").as_matrix();
    data[0].data().append_blk("S", S);
    vector_t b = read_csv<double>("../data/vsr/tensors/b.csv").as_matrix();
    matrix_t g = read_csv<double>("../data/vsr/tensors/g.csv").as_matrix();

    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    // modeling
    FDTI m(b, g, data, fe_it_opt_dti_linearized_gaussian_dirichlet(a, F));
    m.fit(1e-0);

    std::cout << std::endl;
    matrix_t Ln = m.Ln();
    std::cout << Ln.topRows(10) << std::endl;
    std::cout << std::endl;
    matrix_t Ln_true = read_csv<double>("../data/vsr/tensors/L_true_locs.csv").as_matrix();
    std::cout << Ln_true.topRows(10) << std::endl;
    std::cout << std::endl;

    matrix_t exp_data(Ln.rows(), 4);
    for (int i = 0; i < Ln.rows(); ++i) {
        Eigen::Matrix<double, 2, 2> m;
        m(0, 0) = Ln(i, 0);
        m(1, 1) = Ln(i, 1);
        m(0, 1) = Ln(i, 2);   // / std::sqrt(2);
        m(1, 0) = m(0, 1);

        Eigen::Matrix<double, 2, 2> exp_m = expm(m);
        for (int j = 0; j < exp_m.rows(); ++j) {
            for (int h = 0; h < exp_m.cols(); ++h) { exp_data(i, j * exp_m.rows() + h) = exp_m(j, h); }
        }
    }

    write_csv("../data/vsr/tensors/D_est_locs.csv", exp_data);

    // EXPECT_TRUE(almost_equal<double>(Ln, "../data/vsr/tensors/L_true_locs.csv"));
}
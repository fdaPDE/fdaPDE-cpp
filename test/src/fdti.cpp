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

#include "logger.h"
std::ofstream file;

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
    matrix_t Ln_true = read_csv<double>("../data/vsr/tensors/L_true_locs.csv").as_matrix();

    std::string filename = "../data/vsr/tensors/RESULTS/descent.csv";
    {
        std::ofstream clearFile(filename, std::ios::trunc);
        if (!clearFile.is_open()) { std::cerr << "Error clearing file.\n"; }
        // File is cleared here
    }

    file.open(filename, std::ios::app);

    m.fit(1e-18, 1e-10);

    file.close();

    std::cout << std::endl;
    std::cout << std::endl;
    matrix_t Ln = m.Ln();
    std::cout << Ln.topRows(10) << std::endl;
    std::cout << std::endl;
    std::cout << Ln_true.topRows(10) << std::endl;
    std::cout << std::endl;

    matrix_t Dn = m.Dn();

    std::cout << (Ln - Ln_true).norm() / (Ln.size()) << std::endl;

    write_csv("../data/vsr/tensors/D_est_locs.csv", Dn);
    EXPECT_TRUE(almost_equal<double>(Ln, Ln_true));
}
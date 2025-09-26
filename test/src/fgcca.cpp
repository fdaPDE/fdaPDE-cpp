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

TEST(de, test_00) {
    std::string path = "../../../../projects/cca/";

    // data
    Eigen::Matrix<double, Dynamic, Dynamic> X1 = read_csv<double>(path + "X1.csv").as_matrix();
    internals::MultivariateBlock block_1("X1", X1, 0.1);
    std::cout << block_1 << std::endl;

    Eigen::Matrix<double, Dynamic, Dynamic> X2 = read_csv<double>(path + "X2.csv").as_matrix();
    internals::MultivariateBlock block_2("X2", X2, 0.1);
    std::cout << block_2 << std::endl;

    Eigen::Matrix<double, Dynamic, Dynamic> X3 = read_csv<double>(path + "X3.csv").as_matrix();
    internals::MultivariateBlock block_3("X3", X3, 0.1);
    std::cout << block_3 << std::endl;

    Eigen::Matrix<double, Dynamic, Dynamic> X4 = read_csv<double>(path + "X4.csv").as_matrix();
    internals::MultivariateBlock block_4("X4", X4, 0.1);
    std::cout << block_4 << std::endl;

    block_4.l_compute(Eigen::Matrix<double, Eigen::Dynamic, 1>::Ones(201));
    std::cout << block_4.loadings_m().transpose() << "\n" << std::endl;
    std::cout << block_4.components().transpose().leftCols(10) << "\n" << std::endl;

}


TEST(de, test_01) {
    std::string path = "../../../../projects/cca/";

    // geometries
    Triangulation<1, 1> I(0, 1, 21);

    // define physics
    FeSpace Bh(I, P1<1>);
    TrialFunction f(Bh);
    TestFunction  v(Bh);
    auto a = integral(I)(dx(f) * dx(v));
    ZeroField<1> u;
    auto F = integral(I)(u * v);

    // data
    GeoFrame gf_1(I);
    {
        Eigen::Matrix<double, Dynamic, Dynamic> X1 = read_csv<double>(path + "X1.csv").as_matrix();
        auto& level = gf_1.insert_scalar_layer<POINT>("data", path + "locs_1.csv");
        level.load_blk("X1", X1);
    }
    internals::FunctionalBlock block_1("X1", gf_1, fe_ls_elliptic(a, F), 0.1);
    std::cout << block_1 << std::endl;

    GeoFrame gf_2(I);
    {
        Eigen::Matrix<double, Dynamic, Dynamic> X2 = read_csv<double>(path + "X2.csv").as_matrix();
        auto& level = gf_2.insert_scalar_layer<POINT>("data", path + "locs_2.csv");
        level.load_blk("X2", X2);
    }
    internals::FunctionalBlock block_2("X2", gf_2, fe_ls_elliptic(a, F), 0.1);
    std::cout << block_2 << std::endl;

    GeoFrame gf_3(I);
    {
        Eigen::Matrix<double, Dynamic, Dynamic> X3 = read_csv<double>(path + "X3.csv").as_matrix();
        auto& level = gf_3.insert_scalar_layer<POINT>("data", path + "locs_3.csv");
        level.load_blk("X3", X3);
    }
    internals::FunctionalBlock block_3("X3", gf_3, fe_ls_elliptic(a, F), 0.1);
    std::cout << block_3 << std::endl;

    GeoFrame gf_4(I);
    {
        Eigen::Matrix<double, Dynamic, Dynamic> X4 = read_csv<double>(path + "X4.csv").as_matrix();
        auto& level = gf_4.insert_scalar_layer<POINT>("data", path + "locs_4.csv");
        level.load_blk("X4", X4);
    }
    internals::FunctionalBlock block_4("X4", gf_4, fe_ls_elliptic(a, F), 0.1);
    std::cout << block_4 << std::endl;
    block_4.l_compute(Eigen::Matrix<double, Eigen::Dynamic, 1>::Ones(201), 1e-12);
    std::cout << block_4.loadings_m().transpose() << "\n" << std::endl;
    std::cout << block_4.components().transpose().leftCols(10) << "\n" << std::endl;

}
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
//    mesh:         unit_square_40
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    distribution: poisson
TEST(gsr, test_01) {
    // geometry
    Triangulation<2, 2> D = read_mesh<2, 2>("../data/mesh/unit_square_40");
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/gsr/01/locs.csv");
    l1.load_csv<double>("../data/gsr/01/response.csv");
    // modeling
    GSRPDE m("y ~ f", data, Poisson, fe_laplace());
    m.fit(/* lambda = */ 1.25e-06);

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/gsr/01/field.mtx"));
}

// test 2
//    mesh:         unit_square_40
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    distribution: bernulli
TEST(gsr, test_02) {
    // geometry
    Triangulation<2, 2> D = read_mesh<2, 2>("../data/mesh/unit_square_40");
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/gsr/02/locs.csv");
    l1.load_csv<double>("../data/gsr/02/response.csv");
    // modeling
    GSRPDE m("y ~ f", data, Bernoulli, fe_laplace());
    m.fit(/* lambda = */ 1.25e-06);

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/gsr/02/field.mtx"));
}

// test 3
//    mesh:         unit_square_40
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    distribution: exponential
TEST(gsr, test_03) {
    // geometry
    Triangulation<2, 2> D = read_mesh<2, 2>("../data/mesh/unit_square_40");
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/gsr/03/locs.csv");
    l1.load_csv<double>("../data/gsr/03/response.csv");
    // modeling
    GSRPDE m("y ~ f", data, Exponential, fe_laplace());
    m.fit(/* lambda = */ 1.25e-06);
    
    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/gsr/03/field.mtx"));
}

// test 4
//    mesh:         unit_square_40
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    distribution: gamma
TEST(gsr, test_04) {
    // geometry
    Triangulation<2, 2> D = read_mesh<2, 2>("../data/mesh/unit_square_40");
    // data
    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "../data/gsr/04/locs.csv");
    l1.load_csv<double>("../data/gsr/04/response.csv");
    // modeling
    GSRPDE m("y ~ f", data, Gamma, fe_laplace());
    m.fit(/* lambda = */ 1.25e-06);

    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/gsr/04/field.mtx"));
}

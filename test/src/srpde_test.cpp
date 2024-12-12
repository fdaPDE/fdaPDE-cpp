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

#include <cstddef>
#include <gtest/gtest.h>   // testing framework

#include <fdaPDE/core.h>
using fdapde::core::advection;
using fdapde::core::diffusion;
using fdapde::core::FEM;
using fdapde::core::fem_order;
using fdapde::core::laplacian;
using fdapde::core::DiscretizedMatrixField;
using fdapde::core::PDE;
using fdapde::core::DiscretizedVectorField;
using fdapde::core::Triangulation;

#include "../../fdaPDE/models/regression/srpde.h"
#include "../../fdaPDE/models/sampling_design.h"
using fdapde::models::SRPDE;
using fdapde::models::Sampling;

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_csv;

// test 1
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
TEST(srpde_test, laplacian_nonparametric_samplingatnodes) {
    // define domain 
    MeshLoader<Triangulation<2, 2>> domain("unit_square");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/models/srpde/2D_test1/y.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double lambda = 5.623413 * std::pow(0.1, 5);
    SRPDE model(problem, Sampling::mesh_nodes);
    model.set_lambda_D(lambda);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f()  , "../data/models/srpde/2D_test1/sol.mtx"));
}

// test 2
//    domain:       c-shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
TEST(srpde_test, laplacian_semiparametric_samplingatlocations) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("c_shaped");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/models/srpde/2D_test2/locs.csv");
    DMatrix<double> y    = read_csv<double>("../data/models/srpde/2D_test2/y.csv");
    DMatrix<double> X    = read_csv<double>("../data/models/srpde/2D_test2/X.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define statistical model
    double lambda = 0.2201047;
    SRPDE model(problem, Sampling::pointwise);
    model.set_lambda_D(lambda);
    model.set_spatial_locations(locs);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f()   , "../data/models/srpde/2D_test2/sol.mtx" ));
    EXPECT_TRUE(almost_equal(model.beta(), "../data/models/srpde/2D_test2/beta.mtx"));
}

// test 3
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes
//    penalization: costant coefficients PDE
//    covariates:   no
//    BC:           no
//    order FE:     1
TEST(srpde_test, costantcoefficientspde_nonparametric_samplingatnodes) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("unit_square");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/models/srpde/2D_test3/y.csv");
    // define regularizing PDE
    SMatrix<2> K;
    K << 1, 0, 0, 4;
    auto L = -diffusion<FEM>(K);   // anisotropic diffusion
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define  model
    double lambda = 10;
    SRPDE model(problem, Sampling::mesh_nodes);
    model.set_lambda_D(lambda);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f() , "../data/models/srpde/2D_test3/sol.mtx"));
}

// test 4
//    domain:       quasicircular domain
//    sampling:     areal
//    penalization: non-costant coefficients PDE
//    covariates:   no
//    BC:           yes
//    order FE:     1
TEST(srpde_test, noncostantcoefficientspde_nonparametric_samplingareal) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("quasi_circle");
    // import data from files
    DMatrix<double, Eigen::RowMajor> K_data  = read_csv<double>("../data/models/srpde/2D_test4/K.csv");
    DMatrix<double, Eigen::RowMajor> b_data  = read_csv<double>("../data/models/srpde/2D_test4/b.csv");
    DMatrix<double> subdomains = read_csv<double>("../data/models/srpde/2D_test4/incidence_matrix.csv");
    DMatrix<double> u = read_csv<double>("../data/models/srpde/2D_test4/force.csv");
    DMatrix<double> y = read_csv<double>("../data/models/srpde/2D_test4/y.csv"    );
    // define regularizing PDE
    DiscretizedMatrixField<2, 2, 2> K(K_data);
    DiscretizedVectorField<2, 2> b(b_data);
    auto L = -diffusion<FEM>(K) + advection<FEM>(b);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double lambda = std::pow(0.1, 3);
    SRPDE model(problem, Sampling::areal);
    model.set_lambda_D(lambda);
    model.set_spatial_locations(subdomains);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f() , "../data/models/srpde/2D_test4/sol.mtx"));
}

// test 5
//    domain:       c-shaped surface
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
TEST(srpde_test, laplacian_nonparametric_samplingatnodes_surface) {
    // define domain 
  MeshLoader<Triangulation<2, 3>> domain("c_shaped_surface");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/models/srpde/2D_test5/y.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double lambda = 1e-2;
    SRPDE model(problem, Sampling::mesh_nodes);
    model.set_lambda_D(lambda);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f()  , "../data/models/srpde/2D_test5/sol.mtx"));
}

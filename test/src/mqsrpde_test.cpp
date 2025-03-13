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

#include <fdaPDE/core.h>
#include <gtest/gtest.h>   // testing framework

#include <cstddef>
using fdapde::core::advection;
using fdapde::core::diffusion;
using fdapde::core::FEM;
using fdapde::core::fem_order;

using fdapde::core::laplacian;
using fdapde::core::PDE;
using fdapde::core::Triangulation;

#include "../../fdaPDE/models/regression/mqsrpde.h"
#include "../../fdaPDE/models/sampling_design.h"
using fdapde::models::SpaceOnly;
using fdapde::models::MQSRPDE;

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_csv;

// parametri multiple
// double gamma0_ = 5.;                   // crossing penalty   
// double eps_ = 1e-6;                    // crossing tolerance 
// double C_ = 1.5;                       // crossing penalty factor
// double tolerance_ = 1e-5;              // convergence tolerance 
// double tol_weights_ = 1e-6;            // weights tolerance
// std::size_t max_iter_ = 200;           // max number of inner iterations 
// std::size_t max_iter_global_ = 100;    // max number of outer iterations 

// test 1 
//    domain:       unit square
//    sampling:     locations != nodes
//    penalization: constant coefficients PDE
//    covariates:   no
//    BC:           no
//    order FE:     1
TEST(mqsrpde_test1_PullRequest, laplacian_nonparametric_samplingatlocations) {

    // define domain
    MeshLoader<Triangulation<2, 2>> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/models/mqsrpde/2D_test1/locs.csv");
    DMatrix<double> y = read_csv<double>("../data/models/mqsrpde/2D_test1/y.csv");

    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    std::vector<double> alphas = {0.05, 0.10, 0.50, 0.90, 0.95}; 
    DMatrix<double> lambdas;
    lambdas.resize(alphas.size(), 1); 
    double lambda = 1.778279 * std::pow(0.1, 4);
    for(auto i = 0; i < lambdas.rows(); ++i) lambdas(i, 0) = lambda;   
    MQSRPDE<SpaceOnly> model(problem, Sampling::pointwise, alphas);
    model.set_spatial_locations(locs);
    model.setLambdas_D(lambdas);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f(), "../data/models/mqsrpde/2D_test1/sol.mtx"));

}

// test 2 
//    domain:       unit square
//    sampling:     locations != nodes
//    penalization: constant coefficients PDE
//    covariates:   yes
//    BC:           no
//    order FE:     1
TEST(mqsrpde_test2_PullRequest, laplacian_semiparametric_samplingatlocations) {

    // define domain and regularizing PDE
    MeshLoader<Triangulation<2, 2>> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/models/mqsrpde/2D_test1/locs.csv");
    DMatrix<double> y = read_csv<double>("../data/models/mqsrpde/2D_test2/y.csv");
    DMatrix<double> X = read_csv<double>("../data/models/mqsrpde/2D_test2/X.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define statistical model
    std::vector<double> alphas = {0.05, 0.10, 0.50, 0.90, 0.95}; 
    DMatrix<double> lambdas;
    lambdas.resize(alphas.size(), 1); 
    double lambda = 1.778279 * std::pow(0.1, 4);
    for(auto i = 0; i < lambdas.rows(); ++i) lambdas(i, 0) = lambda;   
    MQSRPDE<SpaceOnly> model(problem, Sampling::pointwise, alphas);
    model.set_spatial_locations(locs);
    model.setLambdas_D(lambdas);
    // set model data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f()   , "../data/models/mqsrpde/2D_test2/sol.mtx" ));
    EXPECT_TRUE(almost_equal(model.beta(), "../data/models/mqsrpde/2D_test2/beta.mtx"));

}
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
using fdapde::core::dt;
using fdapde::core::MatrixDataWrapper;
using fdapde::core::PDE;
using fdapde::core::VectorDataWrapper;

#include "../../fdaPDE/models/regression/distributions.h"
#include "../../fdaPDE/models/regression/gsrpde.h"
#include "../../fdaPDE/models/sampling_design.h"
using fdapde::models::GSRPDE;
using fdapde::models::SpaceTimeSeparable;
using fdapde::models::SpaceTimeParabolic;
using fdapde::models::SpaceOnly;
using fdapde::models::GeoStatMeshNodes;
using fdapde::models::GeoStatLocations;
using fdapde::models::Areal;
using fdapde::models::MonolithicSolver;
using fdapde::models::IterativeSolver;
using fdapde::models::Poisson;
using fdapde::models::Bernulli;
using fdapde::models::Exponential;
using fdapde::models::Gamma;

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_mtx;
using fdapde::testing::read_csv;

// test 1
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    distribution: poisson
TEST(gsrpde_test, laplacian_nonparametric_samplingatnodes_poisson) {
    // define domain
    MeshLoader<Mesh2D> domain("unit_square_medium");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/models/gsrpde/2D_test1/locs.csv");
    DMatrix<double> y    = read_csv<double>("../data/models/gsrpde/2D_test1/y.csv"   );
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double lambda_D = 1e-3;
    GSRPDE<decltype(problem), SpaceOnly, GeoStatLocations, MonolithicSolver, Poisson> model(problem);
    model.set_lambda_D(lambda_D);
    model.set_spatial_locations(locs);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f(), "../data/models/gsrpde/2D_test1/sol.mtx"));
}

// test 2
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    distribution: bernulli
TEST(gsrpde_test, laplacian_nonparametric_samplingatlocations_bernulli) {
    // define domain
    MeshLoader<Mesh2D> domain("unit_square_medium");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/models/gsrpde/2D_test2/locs.csv");
    DMatrix<double> y    = read_csv<double>("../data/models/gsrpde/2D_test2/y.csv"   );
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double lambda_D = 1e-3;
    GSRPDE<decltype(problem), SpaceOnly, GeoStatLocations, MonolithicSolver, Bernulli> model(problem);
    model.set_lambda_D(lambda_D);
    model.set_spatial_locations(locs);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f(), "../data/models/gsrpde/2D_test2/sol.mtx"));
}

// test 3
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    distribution: exponential
TEST(gsrpde_test, laplacian_nonparametric_samplingatlocations_exponential) {
    // define domain
    MeshLoader<Mesh2D> domain("unit_square_medium");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/models/gsrpde/2D_test3/locs.csv");
    DMatrix<double> y    = read_csv<double>("../data/models/gsrpde/2D_test3/y.csv"   );
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double lambda_D = 1e-3;
    GSRPDE<decltype(problem), SpaceOnly, GeoStatLocations, MonolithicSolver, Exponential> model(problem);
    model.set_lambda_D(lambda_D);
    model.set_spatial_locations(locs);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f(), "../data/models/gsrpde/2D_test3/sol.mtx"));
}

// test 4
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    distribution: gamma
TEST(gsrpde_test, laplacian_nonparametric_samplingatlocations_gamma) {
    // define domain
    MeshLoader<Mesh2D> domain("unit_square_medium");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/models/gsrpde/2D_test4/locs.csv");
    DMatrix<double> y    = read_csv<double>("../data/models/gsrpde/2D_test4/y.csv"   );
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double lambda_D = 1e-3;
    GSRPDE<decltype(problem), SpaceOnly, GeoStatLocations, MonolithicSolver, Gamma> model(problem);
    model.set_lambda_D(lambda_D);
    model.set_spatial_locations(locs);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f(), "../data/models/gsrpde/2D_test4/sol.mtx"));
}

// test 5
//    domain:       c-shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    time penalization: separable (mass penalization)
//    distribution: gamma
TEST(gsrpde_test, laplacian_semiparametric_samplingatlocations_separable_monolithic_gamma) {
    // define temporal domain
    DVector<double> time_mesh;
    time_mesh.resize(4);
    for (std::size_t i = 0; i < 4; ++i) time_mesh[i] = (1. / 3) * i;
    // define spatial domain
    MeshLoader<Mesh2D> domain("c_shaped");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/models/gsrpde/2D_test5/locs.csv");
    DMatrix<double> y    = read_csv<double>("../data/models/gsrpde/2D_test5/y.csv"   );
    DMatrix<double> X    = read_csv<double>("../data/models/gsrpde/2D_test5/X.csv"   );
    // define regularizing PDE    
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double lambda_D = std::pow(0.1, 2.5);
    double lambda_T = std::pow(0.1, 2.5);
    GSRPDE<decltype(problem), SpaceTimeSeparable, GeoStatLocations, MonolithicSolver, Gamma> model(problem, time_mesh);
    model.set_lambda_D(lambda_D);
    model.set_lambda_T(lambda_T);
    model.set_spatial_locations(locs);
    // set model's data
    BlockFrame<double, int> df;
    df.stack(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f()   , "../data/models/gsrpde/2D_test5/sol.mtx" ));
    EXPECT_TRUE(almost_equal(model.beta(), "../data/models/gsrpde/2D_test5/beta.mtx"));
}

// test 6
//    domain:       c-shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    time penalization: parabolic (monolithic solution)
//    distribution: gamma
TEST(gsrpde_test, laplacian_semiparametric_samplingatlocations_parabolic_monolithic_gamma) {
    // define temporal domain
    DVector<double> time_mesh;
    time_mesh.resize(4);
    for (std::size_t i = 0; i < 4; ++i) time_mesh[i] = (1. / 3) * i;
    // define spatial domain
    MeshLoader<Mesh2D> domain("c_shaped");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/models/gsrpde/2D_test6/locs.csv");
    DMatrix<double> y    = read_csv<double>("../data/models/gsrpde/2D_test6/y.csv"   );
    DMatrix<double> X    = read_csv<double>("../data/models/gsrpde/2D_test6/X.csv"   );
    DMatrix<double> IC   = read_mtx<double>("../data/models/gsrpde/2D_test6/IC.mtx"  );
    // define regularizing PDE    
    auto L = dt<FEM>() - laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, time_mesh.rows());
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double lambda_D = std::pow(0.1, 2.5);
    double lambda_T = std::pow(0.1, 2.5);
    GSRPDE<decltype(problem), SpaceTimeParabolic, GeoStatLocations, MonolithicSolver, Gamma> model(problem, time_mesh);
    model.set_lambda_D(lambda_D);
    model.set_lambda_T(lambda_T);
    model.set_spatial_locations(locs);
    // set model's data
    BlockFrame<double, int> df;
    df.stack(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    model.set_initial_condition(IC);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f()   , "../data/models/gsrpde/2D_test6/sol.mtx" ));
    EXPECT_TRUE(almost_equal(model.beta(), "../data/models/gsrpde/2D_test6/beta.mtx"));    
}

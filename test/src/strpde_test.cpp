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
using fdapde::core::dt;
using fdapde::core::FEM;
using fdapde::core::fem_order;
using fdapde::core::SPLINE;
using fdapde::core::bilaplacian;
using fdapde::core::laplacian;
using fdapde::core::PDE;
using fdapde::core::Triangulation;
using fdapde::core::spline_order;
using fdapde::core::DiscretizedMatrixField; 
using fdapde::core::DiscretizedVectorField; 

#include "../../fdaPDE/models/regression/strpde.h"
#include "../../fdaPDE/models/sampling_design.h"
using fdapde::models::STRPDE;
using fdapde::models::SpaceTimeSeparable;
using fdapde::models::SpaceTimeParabolic;
using fdapde::models::Sampling;

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_mtx;
using fdapde::testing::read_csv;

// test 1
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    time penalization: separable (mass penalization)
TEST(strpde_test, laplacian_nonparametric_samplingatnodes_separable_monolithic) {
    // define temporal and spatial domain
    Triangulation<1, 1> time_mesh(0, 2, 10);
    MeshLoader<Triangulation<2, 2>> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/models/strpde/2D_test1/y.csv");
    // define regularizing PDE in space   
    auto Ld = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3 * time_mesh.n_nodes(), 1);
    PDE<Triangulation<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);
    // define regularizing PDE in time
    auto Lt = -bilaplacian<SPLINE>();
    PDE<Triangulation<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);
    // define model
    double lambda_D = 0.01, lambda_T = 0.01;
    STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::mesh_nodes);
    model.set_lambda_D(lambda_D);
    model.set_lambda_T(lambda_T);
    // set model's data
    BlockFrame<double, int> df;
    df.stack(OBSERVATIONS_BLK, y);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f()  , "../data/models/strpde/2D_test1/sol.mtx"));
}

// test 2
//    domain:       c-shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    time penalization: separable (mass penalization)
TEST(strpde_test, laplacian_semiparametric_samplingatlocations_separable_monolithic) {
    // define temporal and spatial domain
    Triangulation<1, 1> time_mesh(0, fdapde::testing::pi, 4);
    MeshLoader<Triangulation<2, 2>> domain("c_shaped");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/models/strpde/2D_test2/locs.csv");
    DMatrix<double> y    = read_csv<double>("../data/models/strpde/2D_test2/y.csv");
    DMatrix<double> X    = read_csv<double>("../data/models/strpde/2D_test2/X.csv");
    // define regularizing PDE in space
    auto Ld = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<Triangulation<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);
    // define regularizing PDE in time
    auto Lt = -bilaplacian<SPLINE>();
    PDE<Triangulation<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);
    // define model
    double lambda_D = 0.01;
    double lambda_T = 0.01;
    STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::pointwise);
    model.set_lambda_D(lambda_D);
    model.set_lambda_T(lambda_T);
    model.set_spatial_locations(locs);
    // set model's data
    BlockFrame<double, int> df;
    df.stack(OBSERVATIONS_BLK, y);
    df.stack(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f()   , "../data/models/strpde/2D_test2/sol.mtx" ));
    EXPECT_TRUE(almost_equal(model.beta(), "../data/models/strpde/2D_test2/beta.mtx"));
}

// test 3
//    domain:       quasicircular domain
//    sampling:     areal
//    penalization: non-costant coefficients PDE
//    covariates:   no
//    BC:           no
//    order FE:     1
//    time penalization: parabolic (monolithic solution)
TEST(strpde_test, noncostantcoefficientspde_nonparametric_samplingareal_parabolic_monolithic) {
    // define temporal domain
    DVector<double> time_mesh;
    time_mesh.resize(10);
    for (int i = 0; i < time_mesh.size(); ++i) time_mesh[i] = 0.4 * i;
    // define spatial domain
    MeshLoader<Triangulation<2, 2>> domain("quasi_circle");
    // import data from files
    DMatrix<double, Eigen::RowMajor> K_data  = read_csv<double>("../data/models/strpde/2D_test3/K.csv");
    DMatrix<double, Eigen::RowMajor> b_data  = read_csv<double>("../data/models/strpde/2D_test3/b.csv");
    DMatrix<double> subdomains = read_csv<double>("../data/models/strpde/2D_test3/incidence_matrix.csv");
    DMatrix<double> y  = read_csv<double>("../data/models/strpde/2D_test3/y.csv" );
    DMatrix<double> IC = read_csv<double>("../data/models/strpde/2D_test3/IC.csv");   // initial condition
    // define regularizing PDE
    DiscretizedMatrixField<2, 2, 2> K(K_data);
    DiscretizedVectorField<2, 2> b(b_data);
    auto L = dt<FEM>() - diffusion<FEM>(K) + advection<FEM>(b);
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, time_mesh.rows());
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> pde(domain.mesh, time_mesh, L, u);
    pde.set_initial_condition(IC);
    // define model
    double lambda_D = std::pow(0.1, 6);
    double lambda_T = std::pow(0.1, 6);
    STRPDE<SpaceTimeParabolic, fdapde::monolithic> model(pde, Sampling::areal);
    model.set_lambda_D(lambda_D);
    model.set_lambda_T(lambda_T);
    model.set_spatial_locations(subdomains);
    // set model's data
    BlockFrame<double, int> df;
    df.stack(OBSERVATIONS_BLK, y);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f(), "../data/models/strpde/2D_test3/sol.mtx"));
}

// test 4
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    time penalization: parabolic (iterative solver)
TEST(strpde_test, laplacian_nonparametric_samplingatnodes_parabolic_iterative) {
    // define temporal domain
    DVector<double> time_mesh;
    time_mesh.resize(10);
    double x = 0;
    for (int i = 0; i < time_mesh.size(); x += 0.2, ++i) time_mesh[i] = x;
    // define spatial domain
    MeshLoader<Triangulation<2, 2>> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y  = read_mtx<double>("../data/models/strpde/2D_test4/y.mtx" );    
    DMatrix<double> IC = read_mtx<double>("../data/models/strpde/2D_test4/IC.mtx");
    // define regularizing PDE
    auto L = dt<FEM>() - laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, time_mesh.rows());
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> pde(domain.mesh, time_mesh, L, u);
    pde.set_initial_condition(IC);
    // define model
    double lambda_D = 1;
    double lambda_T = 1;
    STRPDE<SpaceTimeParabolic, fdapde::iterative> model(pde, Sampling::mesh_nodes);
    model.set_lambda_D(lambda_D);
    model.set_lambda_T(lambda_T);
    // set model's data
    BlockFrame<double, int> df;
    df.stack(OBSERVATIONS_BLK, y);
    model.set_data(df);
    // set parameters for iterative method
    model.set_tolerance(1e-4);
    model.set_max_iter(50);
    // solve smoothing problem
    model.init();
    model.solve();
    // test corretness
    EXPECT_TRUE(almost_equal(model.f(), "../data/models/strpde/2D_test4/sol.mtx"));
}

// test 5
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes, time locations != time nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    time penalization: separable (mass penalization)
TEST(strpde_test, laplacian_nonparametric_samplingatnodes_timelocations_separable_monolithic) {
    // define temporal and spatial domain
    Triangulation<1, 1> time_mesh(0, 2, 10);
    MeshLoader<Triangulation<2, 2>> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> time_locs = read_csv<double>("../data/models/strpde/2D_test5/time_locations.csv");
    DMatrix<double> y         = read_csv<double>("../data/models/strpde/2D_test5/y.csv");
    // define regularizing PDE in space
    auto Ld = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3 * time_mesh.n_nodes(), 1);
    PDE<Triangulation<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);
    // define regularizing PDE in time
    auto Lt = -bilaplacian<SPLINE>();
    PDE<Triangulation<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

    // define model
    double lambda_D = 0.01;
    double lambda_T = 0.01;
    STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::mesh_nodes);
    model.set_lambda_D(lambda_D);
    model.set_lambda_T(lambda_T);
    model.set_temporal_locations(time_locs);
    // set model's data
    BlockFrame<double, int> df;
    df.stack(OBSERVATIONS_BLK, y);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f(), "../data/models/strpde/2D_test5/sol.mtx"));
}

// test 6
//    domain:         c-shaped
//    space sampling: locations != nodes
//    time sampling:  locations != nodes
//    missing data:   yes
//    penalization:   simple laplacian
//    covariates:     no
//    BC:             no
//    order FE:       1
//    time penalization: separable (mass penalization)
TEST(strpde_test, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {
    // define temporal and spatial domain
    Triangulation<1, 1> time_mesh(0, 1, 20);
    MeshLoader<Triangulation<2, 2>> domain("c_shaped");
    // import data from files
    DMatrix<double> time_locs  = read_csv<double>("../data/models/strpde/2D_test6/time_locations.csv");
    DMatrix<double> space_locs = read_csv<double>("../data/models/strpde/2D_test6/locs.csv");
    DMatrix<double> y          = read_csv<double>("../data/models/strpde/2D_test6/y.csv"   );
    // define regularizing PDE in space
    auto Ld = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3 * time_mesh.n_nodes(), 1);
    PDE<Triangulation<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);
    // define regularizing PDE in time
    auto Lt = -bilaplacian<SPLINE>();
    PDE<Triangulation<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);
    // define model
    STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::pointwise);
    model.set_lambda_D(1e-3);
    model.set_lambda_T(1e-3);
    model.set_spatial_locations(space_locs);
    model.set_temporal_locations(time_locs);
    // set model's data
    BlockFrame<double, int> df;
    df.stack(OBSERVATIONS_BLK, y);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();

    // Save C++ solution 
    DMatrix<double> computedF = model.f();
    const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
    std::ofstream filef("../data/models/strpde/2D_test6/f_STRPDE.csv");
    if (filef.is_open()){
          filef << computedF.format(CSVFormatf);
          filef.close();
    }

    // test correctness
    EXPECT_TRUE(almost_equal(model.f(), "../data/models/strpde/2D_test6/sol.mtx"));
}

// test 7
//    domain:         surface_hub
//    space sampling: locations == nodes
//    time sampling:  locations == nodes
//    missing data:   yes
//    penalization:   simple laplacian
//    covariates:     no
//    BC:             no
//    order FE:       1
//    time penalization: separable (mass penalization)
TEST(strpde_test, laplacian_nonparametric_samplingatnodes_separable_monolithic_surface) {
    // define temporal and spatial domain
    Triangulation<1, 1> time_mesh(0, 4, 4);   // points {0, 1, \ldots, 4}
    MeshLoader<Triangulation<2, 3>> domain("surface");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/models/strpde/2D_test7/y.csv");
    // define regularizing PDE in space
    auto Ld = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3 * time_mesh.n_nodes(), 1);
    PDE<Triangulation<2, 3>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);
    // define regularizing PDE in time
    auto Lt = -bilaplacian<SPLINE>();
    PDE<Triangulation<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);
    // define model
    STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::mesh_nodes);
    model.set_lambda_D(1e-9);
    model.set_lambda_T(1e-6);
    // set model's data
    BlockFrame<double, int> df;
    df.stack(OBSERVATIONS_BLK, y);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f(), "../data/models/strpde/2D_test7/sol.mtx"));
}

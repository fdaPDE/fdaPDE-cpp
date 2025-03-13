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
using fdapde::core::FEM;
using fdapde::core::fem_order;
using fdapde::core::laplacian;
using fdapde::core::diffusion;
using fdapde::core::PDE;
using fdapde::core::Triangulation;
using fdapde::core::bilaplacian;
using fdapde::core::SPLINE;
using fdapde::core::spline_order;
using fdapde::core::Grid;

#include "../../fdaPDE/models/sampling_design.h"
#include "../../fdaPDE/models/regression/qsrpde.h"
#include "../../fdaPDE/models/regression/gcv.h"
using fdapde::models::QSRPDE;
using fdapde::models::SpaceOnly;
using fdapde::models::ExactEDF;
using fdapde::models::StochasticEDF;
using fdapde::models::Sampling;
using fdapde::models::SpaceTime;
using fdapde::models::SpaceTimeSeparable;
#include "../../fdaPDE/calibration/gcv.h"

#include "../../fdaPDE/calibration/kfold_cv.h"
#include "../../fdaPDE/calibration/rmse.h"
using fdapde::calibration::KCV;
using fdapde::calibration::RMSE;

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_mtx;
using fdapde::testing::read_csv;

double RMSE_metric(DVector<double> v1, DVector<double> v2){
    double res = 0.; 
    if(v1.size() != v2.size())
        std::cout << std::endl << "----------ERROR IN RMSE COMPUTATION---------" << std::endl; 
    for(auto i = 0; i < v1.size(); ++i){
        res += (v1[i]-v2[i])*(v1[i]-v2[i]); 
    }
    return std::sqrt(1./(v1.size())*res); 
}

// test 1
//    domain:       unit square [1,1] x [1,1] (coarse)
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
TEST(gcv_qsrpde_test, laplacian_nonparametric_samplingatnodes_spaceonly_gridexact) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test1/y.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double alpha = 0.1;
    QSRPDE<SpaceOnly> model(problem, Sampling::mesh_nodes, alpha);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    auto GCV = model.gcv<ExactEDF>();
    DMatrix<double> lambdas(13, 1);
    for (int i = 0; i < 13; ++i) { lambdas(i, 0) = std::pow(10, -8.0 + 0.25 * i); }
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas);
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/gcv/qsrpde/2D_test1/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/gcv/qsrpde/2D_test1/gcvs.mtx"));
}

// test 2
//    domain:       unit square [1,1] x [1,1] (coarse)
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid stochastic
TEST(gcv_qsrpde_test, laplacian_nonparametric_samplingatnodes_spaceonly_gridstochastic) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test2/y.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double alpha = 0.1;
    QSRPDE<SpaceOnly> model(problem, Sampling::mesh_nodes, alpha);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    std::size_t seed = 438172;
    auto GCV = model.gcv<StochasticEDF>(1000, seed);
    DMatrix<double> lambdas(13, 1);
    for (int i = 0; i < 13; ++i) { lambdas(i, 0) = std::pow(10, -8.0 + 0.25 * i); }
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas);
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/gcv/qsrpde/2D_test2/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/gcv/qsrpde/2D_test2/gcvs.mtx"));
}

// test 3
//    domain:       c-shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
TEST(gcv_qsrpde_test, laplacian_semiparametric_samplingatlocations_gridexact) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("c_shaped");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/gcv/qsrpde/2D_test3/locs.csv");
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test3/y.csv");
    DMatrix<double> X = read_csv<double>("../data/gcv/qsrpde/2D_test3/X.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double alpha = 0.9;
    QSRPDE<SpaceOnly> model(problem, Sampling::pointwise, alpha);
    model.set_spatial_locations(locs);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    auto GCV = model.gcv<ExactEDF>();
    DMatrix<double> lambdas(9, 1);
    for (int i = 0; i < 9; ++i) { lambdas(i, 0) = std::pow(10, -5.0 + 0.25 * i); }
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas);
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/gcv/qsrpde/2D_test3/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/gcv/qsrpde/2D_test3/gcvs.mtx"));
}

// test 4
//    domain:       c-shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    GCV optimization: grid stochastic
TEST(gcv_qsrpde_test, laplacian_semiparametric_samplingatlocations_gridstochastic) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("c_shaped");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/gcv/qsrpde/2D_test4/locs.csv");
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test4/y.csv");
    DMatrix<double> X = read_csv<double>("../data/gcv/qsrpde/2D_test4/X.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double alpha = 0.9;
    QSRPDE<SpaceOnly> model(problem, Sampling::pointwise, alpha);
    model.set_spatial_locations(locs);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D value
    std::size_t seed = 66546513;
    auto GCV = model.gcv<StochasticEDF>(1000, seed);
    DMatrix<double> lambdas(9, 1);
    for (int i = 0; i < 9; ++i) { lambdas(i, 0) = std::pow(10, -5.0 + 0.25 * i); }
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas);
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/gcv/qsrpde/2D_test4/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/gcv/qsrpde/2D_test4/gcvs.mtx"));
}

// test 5
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes
//    penalization: costant coefficients PDE
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
TEST(gcv_qsrpde_test, costantcoefficientspde_nonparametric_samplingatnodes_gridexact) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test5/y.csv");
    // define regularizing PDE
    SMatrix<2> K;
    K << 1, 0, 0, 4;
    auto L = -diffusion<FEM>(K);   // anisotropic diffusion
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double alpha = 0.1;
    QSRPDE<SpaceOnly> model(problem, Sampling::mesh_nodes, alpha);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    auto GCV = model.gcv<ExactEDF>();
    DMatrix<double> lambdas(9, 1);
    for (int i = 0; i < 9; ++i) { lambdas(i, 0) = std::pow(10, -7.0 + 0.25 * i); }
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas);   // optimize gcv field
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/gcv/qsrpde/2D_test5/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/gcv/qsrpde/2D_test5/gcvs.mtx"));
}

// test 6
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes
//    penalization: costant coefficients PDE
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid stochastic
TEST(gcv_qsrpde_test, costantcoefficientspde_nonparametric_samplingatnodes_gridstochastic) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test6/y.csv");
    // define regularizing PDE
    SMatrix<2> K;
    K << 1, 0, 0, 4;
    auto L = -diffusion<FEM>(K);   // anisotropic diffusion
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double alpha = 0.1;
    QSRPDE<SpaceOnly> model(problem, Sampling::mesh_nodes, alpha);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    std::size_t seed = 438172;
    auto GCV = model.gcv<StochasticEDF>(1000, seed);
    DMatrix<double> lambdas(9, 1);
    for (int i = 0; i < 9; ++i) { lambdas(i, 0) = std::pow(10, -7.0 + 0.25 * i); }
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas);   // optimize gcv field
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/gcv/qsrpde/2D_test6/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/gcv/qsrpde/2D_test6/gcvs.mtx"));
}

// test 7
//    domain:       c-shaped
//    sampling:     areal
//    penalization: simple laplacian
//    covariates:   no
//    BC:           yes
//    order FE:     1
//    GCV optimization: grid exact
TEST(gcv_qsrpde_test, laplacian_semiparametric_samplingareal_gridexact) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("c_shaped_areal");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test7/y.csv");
    DMatrix<double> X = read_csv<double>("../data/gcv/qsrpde/2D_test7/X.csv");
    DMatrix<double> subdomains = read_csv<double>("../data/gcv/qsrpde/2D_test7/incidence_matrix.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double alpha = 0.5;
    QSRPDE<SpaceOnly> model(problem, Sampling::areal, alpha);
    model.set_spatial_locations(subdomains);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    auto GCV = model.gcv<ExactEDF>();
    DMatrix<double> lambdas(13, 1);
    for (int i = 0; i < 13; ++i) { lambdas(i, 0) = std::pow(10, -4.0 + 0.25 * i); }
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas);   // optimize gcv field
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/gcv/qsrpde/2D_test7/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/gcv/qsrpde/2D_test7/gcvs.mtx"));
}

// test 8
//    domain:       c-shaped
//    sampling:     areal
//    penalization: simple laplacian
//    covariates:   no
//    BC:           yes
//    order FE:     1
//    GCV optimization: grid stochastic
TEST(gcv_qsrpde_test, laplacian_semiparametric_samplingareal_gridstochastic) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("c_shaped_areal");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test8/y.csv");
    DMatrix<double> X = read_csv<double>("../data/gcv/qsrpde/2D_test8/X.csv");
    DMatrix<double> subdomains = read_csv<double>("../data/gcv/qsrpde/2D_test8/incidence_matrix.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double alpha = 0.5;
    QSRPDE<SpaceOnly> model(problem, Sampling::areal, alpha);
    model.set_spatial_locations(subdomains);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    std::size_t seed = 438172;
    auto GCV = model.gcv<StochasticEDF>(100, seed);
    DMatrix<double> lambdas(13, 1);
    for (int i = 0; i < 13; ++i) { lambdas(i, 0) = std::pow(10, -4.0 + 0.25 * i); }
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas);   // optimize gcv field
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/gcv/qsrpde/2D_test8/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/gcv/qsrpde/2D_test8/gcvs.mtx"));
}

// test 9
//    domain:       c-shaped
//    space sampling: locations != nodes
//    time sampling:  locations != nodes
//    penalization: simple laplacian
//    missing_data: yes
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
//    time penalization: separable (mass penalization)
TEST(gcv_qsrpde_test, laplacian_nonparametric_samplingatlocations_timelocations_separable_gridexact) {
    // define temporal and spatial domain
    Triangulation<1, 1> time_mesh(0, fdapde::testing::pi, 2);   // interval [0, \pi] with 3 knots
    MeshLoader<Triangulation<2, 2>> domain("c_shaped_adj");
    // import data from files
    DMatrix<double> space_locs = read_csv<double>("../data/gcv/qsrpde/2D_test9/locs.csv");
    DMatrix<double> time_locs  = read_csv<double>("../data/gcv/qsrpde/2D_test9/time_locations.csv");
    DMatrix<double> y          = read_csv<double>("../data/gcv/qsrpde/2D_test9/y.csv");
    // define regularizing PDE in space
    auto Ld = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3 * time_mesh.n_nodes(), 1);
    PDE<Triangulation<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);
    // define regularizing PDE in time
    auto Lt = -bilaplacian<SPLINE>();
    PDE<Triangulation<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);
    // define model
    double alpha = 0.5;
    QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);
    model.set_spatial_locations(space_locs);
    model.set_temporal_locations(time_locs);
    // set model's data
    BlockFrame<double, int> df;
    df.stack(OBSERVATIONS_BLK, y);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    auto GCV = model.gcv<ExactEDF>();
    DMatrix<double> lambdas(9, 2);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lambdas(i * 3 + j, 0) = std::pow(10, -4.0 + 1.0 * i);
            lambdas(i * 3 + j, 1) = std::pow(10, -7.0 + 1.0 * j);
        }
    }
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas);
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/gcv/qsrpde/2D_test9/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/gcv/qsrpde/2D_test9/gcvs.mtx"));
}

// test 10
//    domain:       c-shaped
//    space sampling: locations != nodes
//    time sampling:  locations != nodes
//    penalization: simple laplacian
//    missing_data: yes
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid stochastic
//    time penalization: separable (mass penalization)
TEST(gcv_qsrpde_test, laplacian_nonparametric_samplingatlocations_timelocations_separable_gridstochastic) {
    // define temporal and spatial domain
    Triangulation<1, 1> time_mesh(0, fdapde::testing::pi, 2);   // interval [0, \pi] with 3 knots
    MeshLoader<Triangulation<2, 2>> domain("c_shaped_adj");
    // import data from files
    DMatrix<double> space_locs = read_csv<double>("../data/gcv/qsrpde/2D_test10/locs.csv");
    DMatrix<double> time_locs  = read_csv<double>("../data/gcv/qsrpde/2D_test10/time_locations.csv");
    DMatrix<double> y          = read_csv<double>("../data/gcv/qsrpde/2D_test10/y.csv");
    // define regularizing PDE in space 
    auto Ld = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3 * time_mesh.n_nodes(), 1);
    PDE<Triangulation<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);
    // define regularizing PDE in time
    auto Lt = -bilaplacian<SPLINE>();
    PDE<Triangulation<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);
    // define model
    double alpha = 0.5;
    QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);
    model.set_spatial_locations(space_locs);
    model.set_temporal_locations(time_locs);
    // set model's data
    BlockFrame<double, int> df;
    df.stack(OBSERVATIONS_BLK, y);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    std::size_t seed = 66546513;
    auto GCV = model.gcv<StochasticEDF>(100, seed);
    DMatrix<double> lambdas(9, 2);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lambdas(i * 3 + j, 0) = std::pow(10, -4.0 + 1.0 * i);
            lambdas(i * 3 + j, 1) = std::pow(10, -7.0 + 1.0 * j);
        }
    }
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas);
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/gcv/qsrpde/2D_test10/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/gcv/qsrpde/2D_test10/gcvs.mtx"));
    // check consistency with GCV calibrator
    auto GCV_ = fdapde::calibration::GCV<SpaceTime> {Grid<fdapde::Dynamic> {}, StochasticEDF(100, seed)}(lambdas);
    EXPECT_TRUE(GCV_.fit(model) == opt.optimum());
}



/////

// CONFRONTO METODI CV

// test cv
//    domain:       unit square
//    sampling:     locations != nodes
//    penalization: laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
TEST(gcv_sqrpde_test_cv, pde_nonparametric_samplingatlocations_spaceonly_gridexact) {

    std::string test_number = "2";    // "1" "2"

    // path test  
    std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/PhD/Codice/models/QSRPDE/Tests_cv/Test_" + test_number; 

    // define statistical model
    std::vector<double> alphas = {0.5, 0.95}; // {0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99};   

    // define grid of lambda values
    std::vector<std::string> lambda_selection_types = {"gcv_smooth_eps1e-1"}; // {"gcv_smooth_eps1e-3", "gcv_smooth_eps1e-1"};   
    const int eps_power = -1.0;  
    
    bool compute_rmse = true;
    bool compute_gcv = false;     // if you want to compute either gcv, gacv or gacv*
  
    // methods 
    std::vector<std::string> score_types = {"gcv", "k-fold", "10-fold"};   // "gcv" "gacv" "gacv_star" "k-fold" "10-fold"
    std::size_t n_folds; 
    
    const unsigned int n_sim = 10; 

    // define domain
    std::string domain_str; 
    if(test_number == "1"){
        domain_str = "unit_square_15"; 
    }
    if(test_number == "2"){
        domain_str = "unit_square_25"; 
    }
    MeshLoader<Triangulation<2, 2>> domain(domain_str);

    // rhs 
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);

    // define regularizing PDE
    auto L = -laplacian<FEM>();   
    PDE<Triangulation<2, 2>, decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);


    double phi; // gacv parameter
    
    std::vector<double> lambdas_1;
    std::vector<double> lambdas_5; 
    std::vector<double> lambdas_10;
    std::vector<double> lambdas_25; 
    std::vector<double> lambdas_50;
    std::vector<double> lambdas_75; 
    std::vector<double> lambdas_90;
    std::vector<double> lambdas_95; 
    std::vector<double> lambdas_99; 

    std::vector<double> lambdas_longer; 

    if(test_number=="1"){
        for(double x = -9.5; x <= -4.0; x += 0.05) lambdas_1.push_back(std::pow(10, x)); 
        for(double x = -9.5; x <= -4.0; x += 0.05) lambdas_5.push_back(std::pow(10, x));
        for(double x = -9.5; x <= -4.0; x += 0.05) lambdas_10.push_back(std::pow(10, x)); 
        for(double x = -9.5; x <= -4.0; x += 0.05) lambdas_25.push_back(std::pow(10, x));
        for(double x = -9.5; x <= -4.0; x += 0.05) lambdas_50.push_back(std::pow(10, x)); 
        for(double x = -9.5; x <= -4.0; x += 0.05) lambdas_75.push_back(std::pow(10, x)); 
        for(double x = -9.5; x <= -4.0; x += 0.05) lambdas_90.push_back(std::pow(10, x)); 
        for(double x = -9.5; x <= -4.0; x += 0.05) lambdas_95.push_back(std::pow(10, x));
        for(double x = -9.5; x <= -4.0; x += 0.05) lambdas_99.push_back(std::pow(10, x)); 

        for(double x = -4.0; x <= -2.0; x += 0.05) lambdas_longer.push_back(std::pow(10, x)); 
    }
    if(test_number=="2"){
        for(double x = -7.5; x <= -2.0; x += 0.05) lambdas_1.push_back(std::pow(10, x)); 
        for(double x = -7.5; x <= -2.0; x += 0.05) lambdas_5.push_back(std::pow(10, x));
        for(double x = -7.5; x <= -2.0; x += 0.05) lambdas_10.push_back(std::pow(10, x)); 
        for(double x = -7.5; x <= -2.0; x += 0.05) lambdas_25.push_back(std::pow(10, x));
        for(double x = -7.5; x <= -2.0; x += 0.05) lambdas_50.push_back(std::pow(10, x)); 
        for(double x = -7.5; x <= -2.0; x += 0.05) lambdas_75.push_back(std::pow(10, x)); 
        for(double x = -7.5; x <= -2.0; x += 0.05) lambdas_90.push_back(std::pow(10, x)); 
        for(double x = -7.5; x <= -2.0; x += 0.05) lambdas_95.push_back(std::pow(10, x));
        for(double x = -7.5; x <= -2.0; x += 0.05) lambdas_99.push_back(std::pow(10, x)); 

        for(double x = -4.0; x <= -2.0; x += 0.05) lambdas_longer.push_back(std::pow(10, x)); 
    }

    bool force_lambdas_longer = false; 
    double best_lambda; 

    // Read covariates and locations
    DMatrix<double> loc = read_csv<double>(R_path + "/locs.csv"); 

    // Simulations 
    for(auto sim = 1; sim <= n_sim; ++sim){
        std::cout << "--------------------Simulation #" << std::to_string(sim) << "-------------" << std::endl; 
        
        // load data from .csv files
        DMatrix<double> y = read_csv<double>(R_path + "/simulations/sim_" + std::to_string(sim) + "/y.csv");
        BlockFrame<double, int> df;
        df.insert(OBSERVATIONS_BLK, y);

        for(auto alpha : alphas){

            unsigned int alpha_int = alpha*100; 
            std::string alpha_string = std::to_string(alpha_int); 

            std::cout << "------------------alpha=" << alpha_string << "-----------------" << std::endl; 

            std::string solutions_path_rmse = R_path + "/simulations/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + "/RMSE"; 

            std::vector<double> lambdas;
            if(almost_equal(alpha, 0.01)){
                lambdas = lambdas_1; 
            } 
            if(almost_equal(alpha, 0.05)){
                lambdas = lambdas_5; 
            } 
            if(almost_equal(alpha, 0.10)){
                lambdas = lambdas_10; 
            } 
            if(almost_equal(alpha, 0.25)){
                lambdas = lambdas_25; 
            } 
            if(almost_equal(alpha, 0.50)){
                lambdas = lambdas_50;
            }    
            if(almost_equal(alpha, 0.75)){
                lambdas = lambdas_75; 
            }    
            if(almost_equal(alpha, 0.90)){
                lambdas = lambdas_90; 
            }    
            if(almost_equal(alpha, 0.95)){
                lambdas = lambdas_95; 
            } 
            if(almost_equal(alpha, 0.99)){
                lambdas = lambdas_99; 
            } 

            if(force_lambdas_longer){
                lambdas = lambdas_longer; 
            }

            // define lambda sequence as matrix 
            DMatrix<double> lambdas_mat;
            lambdas_mat.resize(lambdas.size(), 1); 
            for(auto i = 0; i < lambdas_mat.rows(); ++i){
                lambdas_mat(i,0) = lambdas[i]; 
            }


            // GCV (o gacv o gacv* o k-fold):
            if(compute_gcv){
                for(auto lambda_selection_type : lambda_selection_types){
                    for(auto score_type : score_types){

                        std::cout << "------------------score=" << score_type << "-----------------" << std::endl;

                        if(score_type == "k-fold"){
                            n_folds = 5; 
                        }
                        if(score_type == "10-fold"){
                            n_folds = 10; 
                        }

                        if(score_type == "k-fold" && n_folds != 5){
                            std::cout << "ATTENTION: k-fold is running with wrong number of folds!" << std::endl; 
                        }

                        std::string solutions_path_gcv = R_path + "/simulations/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + "/" + lambda_selection_type + "/" + score_type; 
                            
                        QSRPDE<SpaceOnly> model_cv(problem, Sampling::pointwise, alpha);
                        model_cv.set_spatial_locations(loc);
                        model_cv.set_eps_power(eps_power); 
                                       
                        model_cv.set_data(df);
                        model_cv.init();

                        if(score_type.length() >= 4 && score_type.substr(score_type.length() - 4) == "fold"){
                            // define KCV engine and search for best lambda which minimizes the model's RMSE
                            KCV kcv(n_folds);
                            auto fit_kcv = kcv.fit(model_cv, lambdas_mat, RMSE(model_cv, alpha, std::pow(10, eps_power)));
                            best_lambda = fit_kcv[0];  // extract the space optimal lambda 

                            // Save k-fold CV score
                            if(force_lambdas_longer){
                                std::ofstream fileGCV_scores(solutions_path_gcv + "/score_ext.csv");
                                for(std::size_t i = 0; i < kcv.avg_scores().size(); ++i) 
                                    fileGCV_scores << std::setprecision(16) << kcv.avg_scores()[i] << "\n"; 
                                fileGCV_scores.close();
                            } else{
                                std::ofstream fileGCV_scores(solutions_path_gcv + "/score.csv");
                                for(std::size_t i = 0; i < kcv.avg_scores().size(); ++i) 
                                    fileGCV_scores << std::setprecision(16) << kcv.avg_scores()[i] << "\n"; 
                                fileGCV_scores.close();
                            }
                            
                        } else{
                            // define GCV function and grid of \lambda_D values
                            auto GCV = model_cv.gcv<ExactEDF>();  
                            phi = std::min(alpha, 1-alpha);
                            if(score_type == "gacv" || score_type == "gacv_star")
                                GCV.set_correction(score_type == "gacv", score_type == "gacv_star", phi); 
            

                            // optimize GCV
                            Grid<fdapde::Dynamic> opt;
                            opt.optimize(GCV, lambdas_mat);
                            
                            best_lambda = opt.optimum()(0,0);

                            // Save GCV score
                            std::ofstream fileGCV_scores(solutions_path_gcv + "/score.csv");
                            for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
                                fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.gcvs()[i]) << "\n"; 
                            fileGCV_scores.close();

                            // Save traces 
                            std::ofstream fileGCV_traces(solutions_path_gcv + "/traces.csv");
                            std::vector<double> vector_traces = GCV.get_trace(); 
                            for(std::size_t i = 0; i < vector_traces.size(); ++i) 
                                fileGCV_traces << std::setprecision(16) << vector_traces[i] << "\n"; 
                            fileGCV_traces.close();
                                
                            // Save edfs
                            std::ofstream fileEDF(solutions_path_gcv + "/edfs.csv");
                            for(auto i = 0; i < GCV.edfs().size(); ++i)
                                fileEDF << std::setprecision(16) << GCV.edfs()[i] << "\n"; 
                            fileEDF.close();

                        }
                    
                        std::cout << "Best lambda is: " << std::setprecision(16) << best_lambda << std::endl; 

                        // Save lambda sequence 
                        if(force_lambdas_longer){
                            std::ofstream fileLambdaS(solutions_path_gcv + "/lambdas_seq_ext.csv");
                            for(std::size_t i = 0; i < lambdas_mat.rows(); ++i) 
                                fileLambdaS << std::setprecision(16) << lambdas_mat(i,0) << "\n"; 
                            fileLambdaS.close();

                            // Save lambda GCVopt for all alphas
                            std::ofstream fileLambdaoptS(solutions_path_gcv + "/lambda_s_opt_ext.csv");
                            if(fileLambdaoptS.is_open()){
                                fileLambdaoptS << std::setprecision(16) << best_lambda;
                                fileLambdaoptS.close();
                            }
                        } else{
                            std::ofstream fileLambdaS(solutions_path_gcv + "/lambdas_seq.csv");
                            for(std::size_t i = 0; i < lambdas_mat.rows(); ++i) 
                                fileLambdaS << std::setprecision(16) << lambdas_mat(i,0) << "\n"; 
                            fileLambdaS.close();

                            // Save lambda GCVopt for all alphas
                            std::ofstream fileLambdaoptS(solutions_path_gcv + "/lambda_s_opt.csv");
                            if(fileLambdaoptS.is_open()){
                                fileLambdaoptS << std::setprecision(16) << best_lambda;
                                fileLambdaoptS.close();
                            }
                        }

                    }
                       
                }
                    
            }

            if(compute_rmse){
                std::cout << "-----RMSE computation-----" << std::endl; 
                // RMSE
                DMatrix<double> f_true = read_csv<double>(R_path + "/true/f_true_" + alpha_string + ".csv");

                std::vector<double> rmse_score; 
                rmse_score.resize(lambdas.size()); 
                double count_l = 0; 
                for(auto lambda : lambdas){
                    QSRPDE<SpaceOnly> model_rmse(problem, Sampling::pointwise, alpha);
                    // set model's data
                    model_rmse.set_spatial_locations(loc);
                    model_rmse.set_lambda_D(lambda);           
                    
                    model_rmse.set_data(df);
                    model_rmse.init();
                    model_rmse.solve();

                    rmse_score[count_l] = RMSE_metric(model_rmse.f(), f_true); 

                    count_l = count_l+1; 
                }

                auto min_idx = std::distance(std::begin(rmse_score), std::min_element(std::begin(rmse_score), std::end(rmse_score))); 
                
                // Save lambda sequence 
                std::ofstream fileLambdaS_rmse(solutions_path_rmse + "/lambdas_seq.csv");
                for(std::size_t i = 0; i < lambdas.size(); ++i) 
                    fileLambdaS_rmse << std::setprecision(16) << lambdas[i] << "\n"; 
                fileLambdaS_rmse.close();

                // Save lambda RMSEopt for all alphas
                std::ofstream fileLambdaoptS_rmse(solutions_path_rmse + "/lambda_s_opt.csv");
                if(fileLambdaoptS_rmse.is_open()){
                    fileLambdaoptS_rmse << std::setprecision(16) << lambdas[min_idx]; ;
                    fileLambdaoptS_rmse.close();
                }

                // Save score 
                std::ofstream fileRMSE_scores(solutions_path_rmse + "/score.csv");
                for(std::size_t i = 0; i < rmse_score.size(); ++i) 
                    fileRMSE_scores << std::setprecision(16) << rmse_score[i] << "\n"; 
                fileRMSE_scores.close();
            
            }

        }


    }
}

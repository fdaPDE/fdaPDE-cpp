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
using fdapde::core::Mesh;
using fdapde::core::bilaplacian;
using fdapde::core::SPLINE;
using fdapde::core::spline_order;
using fdapde::core::Grid;

#include "../../fdaPDE/models/regression/qsrpde.h"
#include "../../fdaPDE/models/regression/gcv.h"
#include "../../fdaPDE/models/sampling_design.h"
using fdapde::models::QSRPDE;
using fdapde::models::SpaceOnly;
using fdapde::models::ExactEDF;
using fdapde::models::StochasticEDF;
using fdapde::models::Sampling;
#include "../../fdaPDE/calibration/gcv.h"

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_mtx;
using fdapde::testing::read_csv;

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
    MeshLoader<Mesh2D> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test1/y.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
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
    std::vector<DVector<double>> lambdas;
    for (double x = -8.0; x <= -5.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
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
    MeshLoader<Mesh2D> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test2/y.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
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
    std::vector<DVector<double>> lambdas;
    for (double x = -8.0; x <= -5.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
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
    MeshLoader<Mesh2D> domain("c_shaped");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/gcv/qsrpde/2D_test3/locs.csv");
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test3/y.csv");
    DMatrix<double> X = read_csv<double>("../data/gcv/qsrpde/2D_test3/X.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
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
    std::vector<DVector<double>> lambdas;
    for (double x = -5.0; x <= -3.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
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
    MeshLoader<Mesh2D> domain("c_shaped");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/gcv/qsrpde/2D_test4/locs.csv");
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test4/y.csv");
    DMatrix<double> X = read_csv<double>("../data/gcv/qsrpde/2D_test4/X.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
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
    std::vector<DVector<double>> lambdas;
    for (double x = -5.0; x <= -3.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
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
    MeshLoader<Mesh2D> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test5/y.csv");
    // define regularizing PDE
    SMatrix<2> K;
    K << 1, 0, 0, 4;
    auto L = -diffusion<FEM>(K);   // anisotropic diffusion
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
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
    std::vector<DVector<double>> lambdas;
    for (double x = -7.0; x <= -5.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
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
    MeshLoader<Mesh2D> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test6/y.csv");
    // define regularizing PDE
    SMatrix<2> K;
    K << 1, 0, 0, 4;
    auto L = -diffusion<FEM>(K);   // anisotropic diffusion
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
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
    std::vector<DVector<double>> lambdas;
    for (double x = -7.0; x <= -5.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
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
    MeshLoader<Mesh2D> domain("c_shaped_areal");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test7/y.csv");
    DMatrix<double> X = read_csv<double>("../data/gcv/qsrpde/2D_test7/X.csv");
    DMatrix<double> subdomains = read_csv<double>("../data/gcv/qsrpde/2D_test7/incidence_matrix.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
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
    std::vector<DVector<double>> lambdas;
    for (double x = -4.0; x <= -1.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
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
    MeshLoader<Mesh2D> domain("c_shaped_areal");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test8/y.csv");
    DMatrix<double> X = read_csv<double>("../data/gcv/qsrpde/2D_test8/X.csv");
    DMatrix<double> subdomains = read_csv<double>("../data/gcv/qsrpde/2D_test8/incidence_matrix.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
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
    std::vector<DVector<double>> lambdas;
    for (double x = -4.0; x <= -1.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
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
    Mesh<1, 1> time_mesh(0, fdapde::testing::pi, 2);   // interval [0, \pi] with 3 knots
    MeshLoader<Mesh2D> domain("c_shaped_adj");
    // import data from files
    DMatrix<double> space_locs = read_csv<double>("../data/gcv/qsrpde/2D_test9/locs.csv");
    DMatrix<double> time_locs  = read_csv<double>("../data/gcv/qsrpde/2D_test9/time_locations.csv");
    DMatrix<double> y          = read_csv<double>("../data/gcv/qsrpde/2D_test9/y.csv");
    // define regularizing PDE in space
    auto Ld = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
    PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);
    // define regularizing PDE in time
    auto Lt = -bilaplacian<SPLINE>();
    PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);
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
    std::vector<DVector<double>> lambda_grid;
    for (double lambda_s = -4.0; lambda_s <= -2.0; lambda_s += 1.0) {
        for (double lambda_t = -7.0; lambda_t <= -5.0; lambda_t += 1.0)
            lambda_grid.push_back(SVector<2>(std::pow(10, lambda_s), std::pow(10, lambda_t)));
    }
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambda_grid);
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
    Mesh<1, 1> time_mesh(0, fdapde::testing::pi, 2);   // interval [0, \pi] with 3 knots
    MeshLoader<Mesh2D> domain("c_shaped_adj");
    // import data from files
    DMatrix<double> space_locs = read_csv<double>("../data/gcv/qsrpde/2D_test10/locs.csv");
    DMatrix<double> time_locs  = read_csv<double>("../data/gcv/qsrpde/2D_test10/time_locations.csv");
    DMatrix<double> y          = read_csv<double>("../data/gcv/qsrpde/2D_test10/y.csv");
    // define regularizing PDE in space 
    auto Ld = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
    PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);
    // define regularizing PDE in time
    auto Lt = -bilaplacian<SPLINE>();
    PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);
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
    std::vector<DVector<double>> lambda_grid;
    for (double lambda_s = -4.0; lambda_s <= -2.0; lambda_s += 1.0) {
        for (double lambda_t = -7.0; lambda_t <= -5.0; lambda_t += 1.0)
            lambda_grid.push_back(SVector<2>(std::pow(10, lambda_s), std::pow(10, lambda_t)));
    }
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambda_grid);
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/gcv/qsrpde/2D_test10/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/gcv/qsrpde/2D_test10/gcvs.mtx"));
    // check consistency with GCV calibrator
    auto GCV_ = fdapde::calibration::GCV {Grid<fdapde::Dynamic> {}, StochasticEDF(100, seed)}(lambda_grid);
    EXPECT_TRUE(GCV_.fit(model) == opt.optimum());
}

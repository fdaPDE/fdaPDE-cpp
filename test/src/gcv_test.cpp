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
using fdapde::core::Grid;
using fdapde::core::laplacian;
using fdapde::core::MatrixDataWrapper;
using fdapde::core::PDE;
using fdapde::core::VectorDataWrapper;

#include "../../fdaPDE/models/regression/srpde.h"
#include "../../fdaPDE/models/sampling_design.h"
using fdapde::models::Areal;
using fdapde::models::GeoStatLocations;
using fdapde::models::GeoStatMeshNodes;
using fdapde::models::SRPDE;

#include "../../fdaPDE/calibration/gcv.h"
using fdapde::calibration::ExactEDF;
using fdapde::calibration::ExactGCV;
using fdapde::calibration::GCV;
using fdapde::calibration::StochasticEDF;

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_mtx;

// test 1
//    domain:       unit square [1,1] x [1,1] (coarse)
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
TEST(gcv_test, laplacian_nonparametric_samplingatnodes_spaceonly_gridexact) {
    // define domain
    MeshLoader<Mesh2D<>> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/models/gcv/2D_test1/y.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM> problem(domain.mesh, L, u);
    // define model
    SRPDE<decltype(problem), GeoStatMeshNodes> model(problem);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
    ScalarField<1, decltype(GCV)> obj(GCV);
    std::vector<SVector<1>> lambdas;
    for (double x = -6.0; x <= -3.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
    // optimize GCV
    Grid<1> opt;
    opt.optimize(obj, lambdas);
    SVector<1> best_lambda = opt.optimum();
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(),   "../data/models/gcv/2D_test1/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.values(), "../data/models/gcv/2D_test1/gcvs.mtx"));
}

// test 2
//    domain:       unit square [1,1] x [1,1] (coarse)
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid stochastic
TEST(gcv_test, laplacian_nonparametric_samplingatnodes_spaceonly_gridstochastic) {
    // define domain
    MeshLoader<Mesh2D<>> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/models/gcv/2D_test2/y.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM> problem(domain.mesh, L, u);
    // define model
    SRPDE<decltype(problem), GeoStatMeshNodes> model(problem);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    std::size_t seed = 476813;
    GCV<decltype(model), StochasticEDF<decltype(model)>> GCV(model, 100, seed);
    ScalarField<1, decltype(GCV)> obj(GCV);
    std::vector<SVector<1>> lambdas;
    for (double x = -6.0; x <= -3.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
    // optimize GCV
    Grid<1> opt;
    opt.optimize(obj, lambdas);
    SVector<1> best_lambda = opt.optimum();
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(),   "../data/models/gcv/2D_test2/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.values(), "../data/models/gcv/2D_test2/gcvs.mtx"));
}

// test 3
//    domain:       c-shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
TEST(gcv_test, laplacian_semiparametric_samplingatlocations_gridexact) {
    // define domain
    MeshLoader<Mesh2D<>> domain("c_shaped");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/models/gcv/2D_test3/locs.csv");
    DMatrix<double> y    = read_csv<double>("../data/models/gcv/2D_test3/y.csv"   );
    DMatrix<double> X    = read_csv<double>("../data/models/gcv/2D_test3/X.csv"   );
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM> problem(domain.mesh, L, u);
    // define model
    SRPDE<decltype(problem), GeoStatLocations> model(problem);
    model.set_spatial_locations(locs);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
    ScalarField<1, decltype(GCV)> obj(GCV);
    std::vector<SVector<1>> lambdas;
    for (double x = -3.0; x <= 3.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
    // optimize GCV
    Grid<1> opt;
    opt.optimize(obj, lambdas);
    SVector<1> best_lambda = opt.optimum();
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(),   "../data/models/gcv/2D_test3/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.values(), "../data/models/gcv/2D_test3/gcvs.mtx"));
}

// test 4
//    domain:       c-shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    GCV optimization: grid stochastic
TEST(gcv_test, laplacian_semiparametric_samplingatlocations_gridstochastic) {
    // define domain
    MeshLoader<Mesh2D<>> domain("c_shaped");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/models/gcv/2D_test4/locs.csv");
    DMatrix<double> y    = read_csv<double>("../data/models/gcv/2D_test4/y.csv"   );
    DMatrix<double> X    = read_csv<double>("../data/models/gcv/2D_test4/X.csv"   );
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM> problem(domain.mesh, L, u);
    // define model
    SRPDE<decltype(problem), GeoStatLocations> model(problem);
    model.set_spatial_locations(locs);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    std::size_t seed = 66546513;
    GCV<decltype(model), StochasticEDF<decltype(model)>> GCV(model, 100, seed);
    ScalarField<1, decltype(GCV)> obj(GCV);
    std::vector<SVector<1>> lambdas;
    for (double x = -3.0; x <= 3.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
    // optimize GCV
    Grid<1> opt;
    opt.optimize(obj, lambdas);
    SVector<1> best_lambda = opt.optimum();
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(),   "../data/models/gcv/2D_test4/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.values(), "../data/models/gcv/2D_test4/gcvs.mtx"));
}

// test 5
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes
//    penalization: costant coefficients PDE
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
TEST(gcv_test, costantcoefficientspde_nonparametric_samplingatnodes_gridexact) {
    // define domain
    MeshLoader<Mesh2D<>> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/models/gcv/2D_test5/y.csv");
    // define regularizing PDE
    SMatrix<2> K;
    K << 1, 0, 0, 4;
    auto L = -diffusion<FEM>(K);   // anisotropic diffusion
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM> problem(domain.mesh, L, u);
    // define model
    SRPDE<decltype(problem), GeoStatMeshNodes> model(problem);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
    ScalarField<1, decltype(GCV)> obj(GCV);
    std::vector<SVector<1>> lambdas;
    for (double x = -6.0; x <= -3.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
    // optimize GCV
    Grid<1> opt;
    opt.optimize(obj, lambdas);   // optimize gcv field
    SVector<1> best_lambda = opt.optimum();
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(),   "../data/models/gcv/2D_test5/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.values(), "../data/models/gcv/2D_test5/gcvs.mtx"));
}

// test 6
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes
//    penalization: costant coefficients PDE
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid stochastic
TEST(gcv_test, costantcoefficientspde_nonparametric_samplingatnodes_gridstochastic) {
    // define domain
    MeshLoader<Mesh2D<>> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/models/gcv/2D_test5/y.csv");
    // define regularizing PDE
    SMatrix<2> K;
    K << 1, 0, 0, 4;
    auto L = -diffusion<FEM>(K);   // anisotropic diffusion
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM> problem(domain.mesh, L, u);
    // define model
    SRPDE<decltype(problem), GeoStatMeshNodes> model(problem);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    std::size_t seed = 4564168;
    GCV<decltype(model), StochasticEDF<decltype(model)>> GCV(model, 100, seed);
    ScalarField<1, decltype(GCV)> obj(GCV);
    std::vector<SVector<1>> lambdas;
    for (double x = -6.0; x <= -3.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
    // optimize GCV
    Grid<1> opt;
    opt.optimize(obj, lambdas);   // optimize gcv field
    SVector<1> best_lambda = opt.optimum();
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(),   "../data/models/gcv/2D_test6/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.values(), "../data/models/gcv/2D_test6/gcvs.mtx"));
}

// test 7
//    domain:       quasicircular domain
//    sampling:     areal
//    penalization: non-costant coefficients PDE
//    covariates:   no
//    BC:           yes
//    order FE:     1
//    GCV optimization: grid exact
TEST(gcv_test, noncostantcoefficientspde_nonparametric_samplingareal_gridexact) {
    // define domain and regularizing PDE
    MeshLoader<Mesh2D<>> domain("quasi_circle");
    // import data from files
    DMatrix<double, Eigen::RowMajor> K_data = read_csv<double>("../data/models/gcv/2D_test7/K.csv");
    DMatrix<double, Eigen::RowMajor> b_data = read_csv<double>("../data/models/gcv/2D_test7/b.csv");
    DMatrix<int> subdomains = read_csv<int>   ("../data/models/gcv/2D_test7/incidence_matrix.csv" );
    DMatrix<double> u       = read_csv<double>("../data/models/gcv/2D_test7/force.csv");
    DMatrix<double> y       = read_csv<double>("../data/models/gcv/2D_test7/y.csv"    );
    // define regularizing PDE
    MatrixDataWrapper<2, 2, 2> K(K_data);
    VectorDataWrapper<2, 2> b(b_data);
    auto L = -diffusion<FEM>(K) + advection<FEM>(b);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM> problem(domain.mesh, L, u);
    // define model
    double lambda = std::pow(0.1, 3);
    SRPDE<decltype(problem), Areal> model(problem);
    model.set_spatial_locations(subdomains);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
    ScalarField<1, decltype(GCV)> obj(GCV);
    std::vector<SVector<1>> lambdas;
    for (double x = -6.0; x <= -3.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
    // optimize GCV
    Grid<1> opt;
    opt.optimize(obj, lambdas);   // optimize gcv field
    SVector<1> best_lambda = opt.optimum();
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(),   "../data/models/gcv/2D_test7/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.values(), "../data/models/gcv/2D_test7/gcvs.mtx"));
}

// test 8
//    domain:       quasicircular domain
//    sampling:     areal
//    penalization: non-costant coefficients PDE
//    covariates:   no
//    BC:           yes
//    order FE:     1
//    GCV optimization: grid stochastic
TEST(gcv_test, noncostantcoefficientspde_nonparametric_samplingareal_gridstochastic) {
    // define domain and regularizing PDE
    MeshLoader<Mesh2D<>> domain("quasi_circle");
    // import data from files
    DMatrix<double, Eigen::RowMajor> K_data = read_csv<double>("../data/models/gcv/2D_test8/K.csv");
    DMatrix<double, Eigen::RowMajor> b_data = read_csv<double>("../data/models/gcv/2D_test8/b.csv");
    DMatrix<int> subdomains = read_csv<int>("../data/models/gcv/2D_test8/incidence_matrix.csv"    );
    DMatrix<double> u       = read_csv<double>("../data/models/gcv/2D_test8/force.csv");
    DMatrix<double> y       = read_csv<double>("../data/models/gcv/2D_test8/y.csv"    );
    // define regularizing PDE
    MatrixDataWrapper<2, 2, 2> K(K_data);
    VectorDataWrapper<2, 2> b(b_data);
    auto L = -diffusion<FEM>(K) + advection<FEM>(b);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM> problem(domain.mesh, L, u);
    // define model
    double lambda = std::pow(0.1, 3);
    SRPDE<decltype(problem), Areal> model(problem);
    model.set_spatial_locations(subdomains);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    std::size_t seed = 438172;
    GCV<decltype(model), StochasticEDF<decltype(model)>> GCV(model, 100, seed);
    ScalarField<1, decltype(GCV)> obj(GCV);
    std::vector<SVector<1>> lambdas;
    for (double x = -6.0; x <= -3.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
    // optimize GCV
    Grid<1> opt;
    opt.optimize(obj, lambdas);
    SVector<1> best_lambda = opt.optimum();
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(),   "../data/models/gcv/2D_test8/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.values(), "../data/models/gcv/2D_test8/gcvs.mtx"));
}

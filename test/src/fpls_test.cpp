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

#include <fdaPDE/core.h>
using fdapde::core::FEM;
using fdapde::core::fem_order;
using fdapde::core::laplacian;
using fdapde::core::PDE;
using fdapde::core::Grid;
using fdapde::core::Triangulation;

#include "../../fdaPDE/models/regression/srpde.h"
using fdapde::models::SRPDE;
#include "../../fdaPDE/models/functional/fpls.h"
#include "../../fdaPDE/models/functional/center.h"
#include "../../fdaPDE/models/sampling_design.h"
using fdapde::models::FPLS;
using fdapde::models::RegularizedSVD;
using fdapde::models::Sampling;
using fdapde::models::center;
using fdapde::models::StochasticEDF;
#include "../../fdaPDE/calibration/symbols.h"
using fdapde::calibration::Calibration;
#include "../../fdaPDE/calibration/off.h"
#include "../../fdaPDE/calibration/gcv.h"

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_csv;
using fdapde::testing::read_mtx;

// test 1
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    BC:           no
//    order FE:     1
//    missing data: no
//    solver: sequential (power iteration) without calibration
TEST(fpls_test, laplacian_samplingatnodes_sequential_off) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("unit_square");
    // import data from files
    DMatrix<double> X = read_csv<double>("../data/models/fpls/2D_test1/X.csv");
    DMatrix<double> Y = read_csv<double>("../data/models/fpls/2D_test1/Y.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> pde(domain.mesh, L, u);
    // define model
    double lambda_D = 10.0;
    RegularizedSVD<fdapde::sequential> rsvd {Calibration::off};
    rsvd.set_tolerance(1e-2);
    rsvd.set_max_iter(20);
    FPLS<SpaceOnly> model(pde, Sampling::mesh_nodes, rsvd);   // functional partial least squares model
    model.set_lambda_D(lambda_D);
    model.set_regression_step_calibrator(fdapde::calibration::Off {}(SVector<1>(lambda_D)));
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, DMatrix<double>(Y.rowwise() - Y.colwise().mean()));   // pointwise centred responses
    // smooth centred functional covariates
    auto centered_covs =
      center(X, SRPDE {pde, Sampling::mesh_nodes}, fdapde::calibration::Off {}(SVector<1>(lambda_D)));
    df.insert(DESIGN_MATRIX_BLK, centered_covs.fitted);
    model.set_data(df);
    // solve FPLS problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.fitted().rowwise() + Y.colwise().mean(), "../data/models/fpls/2D_test1/Y_hat.csv"));
    EXPECT_TRUE(almost_equal(
      model.reconstructed().rowwise() + centered_covs.mean.col(0).transpose(),
      "../data/models/fpls/2D_test1/X_hat.csv"));
    EXPECT_TRUE(almost_equal(model.B(), "../data/models/fpls/2D_test1/B_hat.csv"));
}

// test 2
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    BC:           no
//    order FE:     1
//    missing data: no
//    solver: sequential (power iteration) with GCV calibration
TEST(fpls_test, laplacian_samplingatnodes_sequential_gcv) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("unit_square");
    // import data from files
    DMatrix<double> X = read_csv<double>("../data/models/fpls/2D_test2/X.csv");
    DMatrix<double> Y = read_csv<double>("../data/models/fpls/2D_test2/Y.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> pde(domain.mesh, L, u);
    // define model
    std::size_t seed = 476813;
    // grid for smoothing parameter selection
    DMatrix<double> lambda_grid(5, 1);
    for (int i = 0; i < 5; ++i) lambda_grid(i, 0) = std::pow(10, -4 + 1.0 * i);
    RegularizedSVD<fdapde::sequential> rsvd {Calibration::gcv};
    rsvd.set_tolerance(1e-2);
    rsvd.set_max_iter(20);
    rsvd.set_lambda(lambda_grid);
    rsvd.set_seed(seed);   // for reproducibility purposes in testing
    FPLS<SpaceOnly> model(pde, Sampling::mesh_nodes, rsvd);   // functional partial least square models
    model.set_regression_step_calibrator(
      fdapde::calibration::GCV<SpaceOnly> {Grid<fdapde::Dynamic> {}, StochasticEDF(1000, seed)}(lambda_grid));
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, DMatrix<double>(Y.rowwise() - Y.colwise().mean()));   // pointwise centred responses
    // smooth centred functional covariates (select optimal smoothing)
    auto centered_covs = center(
      X, SRPDE {pde, Sampling::mesh_nodes},
      fdapde::calibration::GCV<SpaceOnly> {Grid<fdapde::Dynamic> {}, StochasticEDF(1000, seed)}(lambda_grid));
    df.insert(DESIGN_MATRIX_BLK, centered_covs.fitted);
    model.set_data(df);
    // solve FPLS problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.fitted().rowwise() + Y.colwise().mean(), "../data/models/fpls/2D_test2/Y_hat.csv"));
    EXPECT_TRUE(almost_equal(
      model.reconstructed().rowwise() + centered_covs.mean.col(0).transpose(),
      "../data/models/fpls/2D_test2/X_hat.csv"));
    EXPECT_TRUE(almost_equal(model.B(), "../data/models/fpls/2D_test2/B_hat.csv"));
}

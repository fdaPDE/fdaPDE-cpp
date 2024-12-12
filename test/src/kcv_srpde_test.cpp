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
using fdapde::core::fem_order;
using fdapde::core::FEM;
using fdapde::core::laplacian;
using fdapde::core::PDE;
using fdapde::core::Triangulation;

#include "../../fdaPDE/models/regression/srpde.h"
#include "../../fdaPDE/models/regression/qsrpde.h"
#include "../../fdaPDE/models/sampling_design.h"
using fdapde::models::SRPDE;
using fdapde::models::QSRPDE;
using fdapde::models::SpaceOnly;
using fdapde::models::Sampling;
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

// test 1
//    domain:       unit square [1,1] x [1,1] (coarse)
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
TEST(kcv_srpde_test, laplacian_nonparametric_samplingatnodes_spaceonly_rmse) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/models/gcv/2D_test1/y.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    SRPDE model(problem, Sampling::mesh_nodes);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    model.init();    
    // define KCV engine and search for best lambda which minimizes the model's RMSE
    std::size_t n_folds = 5;
    KCV kcv(n_folds);
    std::vector<DVector<double>> lambdas;
    for (double x = -6.0; x <= -3.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
    kcv.fit(model, lambdas, RMSE(model));

    auto KCV_ = fdapde::calibration::KCV {n_folds}(lambdas, RMSE());
    KCV_.fit(model);
    
    // test correctness
    // TODO
}


TEST(kcv_srpde_test, qsrpde_laplacian_nonparametric_samplingatnodes_spaceonly_rmse) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/models/qsrpde/2D_test1/y.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double lambda = 1.778279 * std::pow(0.1, 4);
    double alpha = 0.1;
    QSRPDE<SpaceOnly> model(problem, Sampling::mesh_nodes, alpha);
    model.set_lambda_D(lambda);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    // define KCV engine and search for best lambda which minimizes the model's RMSE
    std::size_t n_folds = 5;
    KCV kcv(n_folds);
    std::vector<DVector<double>> lambdas;
    for (double x = -6.0; x <= -3.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
    kcv.fit(model, lambdas, RMSE(model));

    std::cout << kcv.avg_scores() << std::endl;
    
    // calibrator approach
    auto KCV_ = fdapde::calibration::KCV {n_folds}(lambdas, RMSE());
    KCV_.fit(model);
    
    // test correctness
    // TODO
}

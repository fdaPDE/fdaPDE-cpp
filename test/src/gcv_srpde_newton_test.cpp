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
using fdapde::core::Newton;
using fdapde::core::laplacian;
using fdapde::core::PDE;
using fdapde::core::Triangulation;

#include "../../fdaPDE/models/regression/srpde.h"
#include "../../fdaPDE/models/regression/gcv.h"
#include "../../fdaPDE/models/sampling_design.h"
#include "../../fdaPDE/models/regression/regression_type_erasure.h"
using fdapde::models::SRPDE;
using fdapde::models::ExactEDF;
using fdapde::models::GCV;
using fdapde::models::StochasticEDF;
using fdapde::models::Sampling;
using fdapde::models::RegressionView;
#include "../../fdaPDE/calibration/gcv.h"

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_mtx;
using fdapde::testing::read_csv;

/* test 1
   domain:       unit square [1,1] x [1,1] (coarse)
   sampling:     locations = nodes
   penalization: simple laplacian
   covariates:   no
   BC:           no
   order FE:     1
   GCV optimization: newton exact
 */
/*TEST(GCV_SRPDE, Test1_Laplacian_NonParametric_GeostatisticalAtNodes_NewtonExact) {
  // define domain and regularizing PDE
  MeshLoader<Triangulation<2, 2><>> domain("unit_square_coarse");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE

  // define statistical model
  SRPDE<decltype(problem), fdaPDE::models::GeoStatMeshNodes> model(problem);
  
  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/SRPDE/2D_test6/y.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert(OBSERVATIONS_BLK, y);
  model.setData(df);
  model.init(); // init model

  // define GCV function and optimize
  ExactGCV<decltype(model), fdaPDE::models::SpaceOnly> GCV(model);
  NewtonOptimizer<1> opt(10, 0.05, 1);

  Logger<decltype(opt)> opt_logger;
  
  opt.findMinimum(GCV, SVector<1>(6.25e-06), opt_logger); // optimize gcv field
  SVector<1> best_lambda = opt.getSolution();

  // expected values of \lambda explored during optimization
  std::vector<double> expected_lambdas = {
    0.00000625000000000, 0.00000753370873002, 0.00000771088319372, 0.00000771375137800
  };
  for(std::size_t i = 0; i < expected_lambdas.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_lambdas[i], opt_logger.x_vect()[i][0]) );

  // expected value of gcv(\lambda)
  std::vector<double> expected_gcvs = {
    0.03907018620179924, 0.03903868737852587, 0.03903824095062697, 0.03903824083843713
  };
  for(std::size_t i = 0; i < expected_gcvs.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_gcvs[i], opt_logger.y_vect()[i]) );

  // check optimal lambda
  EXPECT_TRUE( almost_equal(best_lambda[0], expected_lambdas[3]) );
  }*/

// test 2
//   domain:       unit square [1,1] x [1,1] (coarse)
//   sampling:     locations = nodes
//   penalization: simple laplacian
//   covariates:   no
//   BC:           no
//   order FE:     1
//   GCV optimization: newton finite differences, exact evaluation of Tr[S]
TEST(gcv_srpde_newton_test, laplacian_nonparametric_samplingatnodes_newton_fd_exact) {
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
  // define GCV function
  auto GCV = model.gcv<ExactEDF>();
  GCV.set_step(4e-08);
  // optimize GCV
  Newton<fdapde::Dynamic> opt(10, 0.05, 1);
  DVector<double> pt = SVector<1>(6.25e-06);
  opt.optimize(GCV, pt);
  auto best_lambda = opt.optimum();
  DVector<double> expected_lambda = SVector<1>(0.00000771375137800);
  // test correctness
  // EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/models/gcv/2D_test1/edfs.mtx"));
  // EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/models/gcv/2D_test1/gcvs.mtx"));

  // check optimal lambda
  EXPECT_TRUE( almost_equal(best_lambda[0], expected_lambda[0]) );
}

// test 3
//   domain:       unit square [1,1] x [1,1] (coarse)
//   sampling:     locations = nodes
//   penalization: simple laplacian
//   covariates:   no
//   BC:           no
//   order FE:     1
//   GCV optimization: newton finite differences, stochastic evaluation of Tr[S]
TEST(gcv_srpde_newton_test, laplacian_nonparametric_samplingatnodes_newton_fd_stochastic) {
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
  // define GCV function
  std::size_t seed = 7896453;
  auto GCV = model.gcv<StochasticEDF>(100, seed);
  GCV.set_step(4e-08);
  // optimize GCV
  Newton<fdapde::Dynamic> opt(10, 0.05, 1);
  DMatrix<double> pt = SVector<1>(6.25e-06);
  opt.optimize(GCV, pt);
  auto best_lambda = opt.optimum();
  DVector<double> expected_lambda = SVector<1>(0.0000075627208132);
  // check optimal lambda
  EXPECT_TRUE(almost_equal(best_lambda[0], expected_lambda[0]));

  // check consistency with GCV calibrator
  auto GCV_ = fdapde::calibration::GCV<SpaceOnly> {Newton<fdapde::Dynamic>(10, 0.05, 1), StochasticEDF(100, seed)};
  GCV_.set_step(4e-8);
  EXPECT_TRUE(GCV_(pt).fit(RegressionView<void>(model)) == opt.optimum());
}

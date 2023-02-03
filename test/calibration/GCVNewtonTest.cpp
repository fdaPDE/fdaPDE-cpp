#include <cstddef>
#include <gtest/gtest.h> // testing framework
#include <unsupported/Eigen/SparseExtra>

#include "../fdaPDE/core/utils/Symbols.h"
#include "../fdaPDE/core/utils/IO/CSVReader.h"
#include "../fdaPDE/core/FEM/PDE.h"
using fdaPDE::core::FEM::PDE;
#include "../fdaPDE/core/FEM/operators/SpaceVaryingFunctors.h"
using fdaPDE::core::FEM::SpaceVaryingDiffusion;
using fdaPDE::core::FEM::SpaceVaryingAdvection;
#include "core/MESH/Mesh.h"
#include "../fdaPDE/models/regression/SRPDE.h"
using fdaPDE::models::SRPDE;
#include "../fdaPDE/models/SamplingDesign.h"
using fdaPDE::models::Sampling;
#include "../fdaPDE/calibration/GCV.h"
using fdaPDE::calibration::FiniteDifferenceGCV;
using fdaPDE::calibration::ExactGCV;
using fdaPDE::calibration::ExactEDF;
using fdaPDE::calibration::StochasticEDF;
#include "../fdaPDE/core/OPT/optimizers/GridOptimizer.h"
using fdaPDE::core::OPT::GridOptimizer;
#include "../fdaPDE/core/OPT/optimizers/Newton.h"
using fdaPDE::core::OPT::NewtonOptimizer;
#include "../fdaPDE/core/OPT/extensions/Logger.h"
using fdaPDE::core::OPT::Logger;

#include "../utils/MeshLoader.h"
using fdaPDE::testing::MeshLoader;
#include "../utils/Constants.h"
using fdaPDE::testing::DOUBLE_TOLERANCE;
#include "../utils/Utils.h"
using fdaPDE::testing::almost_equal;

/* test 1
   domain:       unit square [1,1] x [1,1] (coarse)
   sampling:     locations = nodes
   penalization: simple laplacian
   covariates:   no
   BC:           no
   order FE:     1
   GCV optimization: newton exact
 */
TEST(GCV_SRPDE, Test1_Laplacian_NonParametric_GeostatisticalAtNodes_NewtonExact) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square_coarse");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE

  // define statistical model
  SRPDE<decltype(problem), fdaPDE::models::Sampling::GeoStatMeshNodes> model(problem);
  
  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/SRPDE/2D_test6/y.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert("y", y);
  model.setData(df);
  model.init(); // init model

  // define GCV function and optimize
  ExactGCV<decltype(model), fdaPDE::models::SpaceOnlyTag> GCV(model);
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
}

/* test 2
   domain:       unit square [1,1] x [1,1] (coarse)
   sampling:     locations = nodes
   penalization: simple laplacian
   covariates:   no
   BC:           no
   order FE:     1
   GCV optimization: newton finite differences, exact evaluation of Tr[S]
 */
TEST(GCV_SRPDE, Test2_Laplacian_NonParametric_GeostatisticalAtNodes_NewtonFiniteDifferences_ExactEDF) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square_coarse");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE

  // define statistical model
  SRPDE<decltype(problem), fdaPDE::models::Sampling::GeoStatMeshNodes> model(problem);
  
  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/SRPDE/2D_test6/y.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert("y", y);
  model.setData(df);
  model.init(); // init model

  // define GCV function and optimize
  FiniteDifferenceGCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
  NewtonOptimizer<1> opt(10, 0.05, 1);

  Logger<decltype(opt)> opt_logger;

  ScalarField<1> obj(GCV);
  obj.setStep(4e-08);
  opt.findMinimum(obj, SVector<1>(6.25e-06), opt_logger); // optimize gcv field
  SVector<1> best_lambda = opt.getSolution();
  
  // expected values of \lambda explored during optimization
  std::vector<double> expected_lambdas = {

    0.0000062500000000, 0.0000075337583053, 0.0000077109325588, 0.0000077137989647
  };
  for(std::size_t i = 0; i < expected_lambdas.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_lambdas[i], opt_logger.x_vect()[i][0]) );
    
  // expected value of gcv(\lambda)
  std::vector<double> expected_gcvs = { 
    0.0390701862017992, 0.0390386871313230, 0.0390382409468828, 0.0390382408384679
  };
  for(std::size_t i = 0; i < expected_gcvs.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_gcvs[i], opt_logger.y_vect()[i]) );
    
  // check optimal lambda
  EXPECT_TRUE( almost_equal(best_lambda[0], expected_lambdas[3]) );
}

/* test 3
   domain:       unit square [1,1] x [1,1] (coarse)
   sampling:     locations = nodes
   penalization: simple laplacian
   covariates:   no
   BC:           no
   order FE:     1
   GCV optimization: newton finite differences, stochastic evaluation of Tr[S]
 */
TEST(GCV_SRPDE, Test1_Laplacian_NonParametric_GeostatisticalAtNodes_NewtonFiniteDifferences_StochasticEDF) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square_coarse");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE

  // define statistical model
  SRPDE<decltype(problem), fdaPDE::models::Sampling::GeoStatMeshNodes> model(problem);
  
  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/SRPDE/2D_test6/y.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert("y", y);
  model.setData(df);
  model.init(); // init model

  // define GCV function and optimize
  std::size_t seed = 7896453;
  FiniteDifferenceGCV<decltype(model), StochasticEDF<decltype(model)>> GCV(model, 100, seed);
  NewtonOptimizer<1> opt(10, 0.05, 1);

  Logger<decltype(opt)> opt_logger;

  ScalarField<1> obj(GCV);
  obj.setStep(4e-08);
  opt.findMinimum(obj, SVector<1>(6.25e-06), opt_logger); // optimize gcv field
  SVector<1> best_lambda = opt.getSolution();
  
  // expected values of \lambda explored during optimization
  std::vector<double> expected_lambdas = {
    0.0000062500000000, 0.0000070509504054, 0.0000071167457322
  };
  for(std::size_t i = 0; i < expected_lambdas.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_lambdas[i], opt_logger.x_vect()[i][0]) );
  
  // expected value of gcv(\lambda)
  std::vector<double> expected_gcvs = {
    0.0383343354839899, 0.0383230704663529, 0.0383230077514622
  };
  for(std::size_t i = 0; i < expected_gcvs.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_gcvs[i], opt_logger.y_vect()[i]) );
  
  // check optimal lambda
  EXPECT_TRUE( almost_equal(best_lambda[0], expected_lambdas[2]) );
}

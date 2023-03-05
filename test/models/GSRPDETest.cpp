#include <gtest/gtest.h> // testing framework
#include <unsupported/Eigen/SparseExtra>

#include "../../fdaPDE/core/utils/Symbols.h"
#include "../../fdaPDE/core/utils/IO/CSVReader.h"
#include "../../fdaPDE/core/FEM/PDE.h"
#include "models/ModelBase.h"
using fdaPDE::core::FEM::PDE;
#include "../../fdaPDE/core/FEM/operators/SpaceVaryingFunctors.h"
using fdaPDE::core::FEM::SpaceVaryingDiffusion;
using fdaPDE::core::FEM::SpaceVaryingAdvection;
#include "core/MESH/Mesh.h"
#include "../../fdaPDE/models/regression/GSRPDE.h"
using fdaPDE::models::GSRPDE;
#include "../../fdaPDE/models/ModelTraits.h"
using fdaPDE::models::SolverType;
using fdaPDE::models::SpaceOnlyTag;
#include "../../fdaPDE/models/SamplingDesign.h"
using fdaPDE::models::Sampling;
#include "../../fdaPDE/models/regression/Distributions.h"
#include "../../fdaPDE/preprocess/InitialConditionEstimator.h"
using fdaPDE::preprocess::InitialConditionEstimator;

#include "../utils/MeshLoader.h"
using fdaPDE::testing::MeshLoader;
#include "../utils/Constants.h"
using fdaPDE::testing::DOUBLE_TOLERANCE;
#include "../utils/Utils.h"
using fdaPDE::testing::almost_equal;

// PS: in the following tests the tolerance is set to 10^-10 since FPIRLS makes use of vectorized operations
// the numerical precision of vectorized instructions depend on the specific architecture and compiler used, 10^-10
// is a safety treshold which should work for any architecture and compiler setting.

/* test 1
   domain:       unit square [1,1] x [1,1]
   sampling:     locations != nodes
   penalization: simple laplacian
   covariates:   no
   BC:           no
   order FE:     1
   distribution: poisson
 */
TEST(GSRPDE, Test1_Laplacian_NonParametric_GeostatisticalAtNodes_Poisson) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square_medium");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE

  CSVReader<double> reader{};
  // load locations where data are sampled
  CSVFile<double> locFile;
  locFile = reader.parseFile("data/models/GSRPDE/2D_test1/locs.csv");
  DMatrix<double> loc = locFile.toEigen();
  
  // define statistical model
  double lambda = 1e-3;
  GSRPDE<decltype(problem), SpaceOnlyTag, fdaPDE::models::GeoStatLocations,
	 SolverType::Monolithic, fdaPDE::models::Poisson> model(problem);
  model.setLambdaS(lambda);
  
  // load data from .csv files
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/GSRPDE/2D_test1/y.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert(OBSERVATIONS_BLK, y);
  df.insert(SPACE_LOCATIONS_BLK, loc);
  model.setData(df);

  // solve smoothing problem
  model.init();
  model.solve();

  //   **  test correctness of computed results  **

  // estimate of spatial field \hat f
  SpMatrix<double> expectedSolution;
  Eigen::loadMarket(expectedSolution,   "data/models/GSRPDE/2D_test1/sol.mtx");
  DMatrix<double> computedF = model.f();
  std::size_t N = computedF.rows();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution).topRows(N), computedF, std::pow(0.1, 10)) );
}

/* test 2
   domain:       unit square [1,1] x [1,1]
   sampling:     locations != nodes
   penalization: simple laplacian
   covariates:   no
   BC:           no
   order FE:     1
   distribution: bernulli
 */
TEST(GSRPDE, Test2_Laplacian_NonParametric_GeostatisticalLocations_Bernulli) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square_medium");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE

  CSVReader<double> reader{};
  // load locations where data are sampled
  CSVFile<double> locFile;
  locFile = reader.parseFile("data/models/GSRPDE/2D_test2/locs.csv");
  DMatrix<double> loc = locFile.toEigen();
  
  // define statistical model
  double lambda = 1e-3;
  GSRPDE<decltype(problem), SpaceOnlyTag, fdaPDE::models::GeoStatLocations,
	 SolverType::Monolithic, fdaPDE::models::Bernulli> model(problem);
  model.setLambdaS(lambda);
  
  // load data from .csv files
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/GSRPDE/2D_test2/y.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert(OBSERVATIONS_BLK, y);
  df.insert(SPACE_LOCATIONS_BLK, loc);
  model.setData(df);

  // solve smoothing problem
  model.init();
  model.solve();

  //   **  test correctness of computed results  **
  
  // estimate of spatial field \hat f
  SpMatrix<double> expectedSolution;
  Eigen::loadMarket(expectedSolution,   "data/models/GSRPDE/2D_test2/sol.mtx");
  DMatrix<double> computedF = model.f();
  std::size_t N = computedF.rows();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution).topRows(N), computedF, std::pow(0.1, 10)) );
}

/* test 3
   domain:       unit square [1,1] x [1,1]
   sampling:     locations != nodes
   penalization: simple laplacian
   covariates:   no
   BC:           no
   order FE:     1
   distribution: exponential
 */
TEST(GSRPDE, Test3_Laplacian_NonParametric_GeostatisticalLocations_Exponential) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square_medium");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE

  CSVReader<double> reader{};
  // load locations where data are sampled
  CSVFile<double> locFile;
  locFile = reader.parseFile("data/models/GSRPDE/2D_test3/locs.csv");
  DMatrix<double> loc = locFile.toEigen();
  
  // define statistical model
  double lambda = 1e-3;
  GSRPDE<decltype(problem), SpaceOnlyTag, fdaPDE::models::GeoStatLocations,
	 SolverType::Monolithic, fdaPDE::models::Exponential> model(problem);
  model.setLambdaS(lambda);
  
  // load data from .csv files
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/GSRPDE/2D_test3/y.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert(OBSERVATIONS_BLK, y);
  df.insert(SPACE_LOCATIONS_BLK, loc);
  model.setData(df);

  // solve smoothing problem
  model.init();
  model.solve();

  //   **  test correctness of computed results  **
  
  // estimate of spatial field \hat f
  SpMatrix<double> expectedSolution;
  Eigen::loadMarket(expectedSolution,   "data/models/GSRPDE/2D_test3/sol.mtx");
  DMatrix<double> computedF = model.f();
  std::size_t N = computedF.rows();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution).topRows(N), computedF, std::pow(0.1, 10)) );
}

/* test 4
   domain:       unit square [1,1] x [1,1]
   sampling:     locations != nodes
   penalization: simple laplacian
   covariates:   no
   BC:           no
   order FE:     1
   distribution: gamma
 */
TEST(GSRPDE, Test4_Laplacian_NonParametric_GeostatisticalLocations_Gamma) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square_medium");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE

  CSVReader<double> reader{};
  // load locations where data are sampled
  CSVFile<double> locFile;
  locFile = reader.parseFile("data/models/GSRPDE/2D_test4/locs.csv");
  DMatrix<double> loc = locFile.toEigen();
  
  // define statistical model
  double lambda = 1e-3;
  GSRPDE<decltype(problem), SpaceOnlyTag, fdaPDE::models::GeoStatLocations,
	 SolverType::Monolithic, fdaPDE::models::Gamma> model(problem);
  model.setLambdaS(lambda);
  
  // load data from .csv files
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/GSRPDE/2D_test4/y.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert(OBSERVATIONS_BLK, y);
  df.insert(SPACE_LOCATIONS_BLK, loc);
  model.setData(df);

  // solve smoothing problem
  model.init();
  model.solve();

  //   **  test correctness of computed results  **
  
  // estimate of spatial field \hat f
  SpMatrix<double> expectedSolution;
  Eigen::loadMarket(expectedSolution,   "data/models/GSRPDE/2D_test4/sol.mtx");
  DMatrix<double> computedF = model.f();
  std::size_t N = computedF.rows();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution).topRows(N), computedF, std::pow(0.1, 10)) );
}

/* test 5
   domain:       c-shaped
   sampling:     locations != nodes
   penalization: simple laplacian
   covariates:   yes
   BC:           no
   order FE:     1
   time penalization: separable (mass penalization)
   distribution: gamma
 */
TEST(GSRPDE, Test5_Laplacian_SemiParametric_GeostatisticalAtLocations_Separable_Monolithic_Gamma) {
  // define time domain: unit inteval [0,1] partioned in 4 subintervals
  DVector<double> time_mesh;
  time_mesh.resize(4);
  for(std::size_t i = 0; i < 4; ++i)
    time_mesh[i] = (1./3)*i;

  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("c_shaped");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE
  // define statistical model
  double lambdaS = std::pow(0.1, 2.5); // smoothing in space
  double lambdaT = std::pow(0.1, 2.5); // smoothing in time
  // load sample position
  CSVReader<double> reader{};
  CSVFile<double> locFile; // locations file
  locFile = reader.parseFile("data/models/GSRPDE/2D_test5/locs.csv");
  DMatrix<double> loc = locFile.toEigen();

  GSRPDE<decltype(problem), fdaPDE::models::SpaceTimeSeparableTag,
	 fdaPDE::models::GeoStatLocations, SolverType::Monolithic, fdaPDE::models::Gamma>
    model(problem, time_mesh);
  model.setLambdaS(lambdaS);
  model.setLambdaT(lambdaT);
  
  // load data from .csv files
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile  ("data/models/GSRPDE/2D_test5/y.csv");
  DMatrix<double> y = yFile.toEigen();
  CSVFile<double> XFile; // design matrix
  XFile = reader.parseFile  ("data/models/GSRPDE/2D_test5/X.csv");
  DMatrix<double> X = XFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.stack (OBSERVATIONS_BLK,  y);
  df.insert(DESIGN_MATRIX_BLK, X);
  df.insert(SPACE_LOCATIONS_BLK, loc);
  model.setData(df);
  
  // solve smoothing problem
  model.init();
  model.solve();

  //   **  test correctness of computed results  **   
  
  // estimate of spatial field \hat f
  SpMatrix<double> expectedSolution;
  Eigen::loadMarket(expectedSolution, "data/models/GSRPDE/2D_test5/sol_separable.mtx");
  DMatrix<double> computedF = model.f();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution), computedF, std::pow(0.1, 10)) );  
  // estimate of coefficient vector \hat \beta
  SpMatrix<double> expectedBeta;
  Eigen::loadMarket(expectedBeta, "data/models/GSRPDE/2D_test5/beta_separable.mtx");
  DVector<double> computedBeta = model.beta();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedBeta), computedBeta, std::pow(0.1, 10)) );
}

/* test 6
   domain:       c-shaped
   sampling:     locations != nodes
   penalization: simple laplacian
   covariates:   yes
   BC:           no
   order FE:     1
   time penalization: parabolic (monolithic solution)
   distribution: gamma
 */
TEST(GSRPDE, Test6_Laplacian_SemiParametric_GeostatisticalAtLocations_Parabolic_Monolithic_EstimatedIC_Gamma) {
  // define time domain: unit inteval [0,1] partioned in 4 subintervals
  DVector<double> time_mesh;
  time_mesh.resize(4);
  for(std::size_t i = 0; i < 4; ++i)
    time_mesh[i] = (1./3)*i;

  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("c_shaped");
  auto L = dT() + Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, time_mesh.rows());
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE
  // define statistical model
  double lambdaS = std::pow(0.1, 2.5); // smoothing in space
  double lambdaT = std::pow(0.1, 2.5); // smoothing in time
  // load sample position
  CSVReader<double> reader{};
  CSVFile<double> locFile; // locations file
  locFile = reader.parseFile("data/models/GSRPDE/2D_test5/locs.csv");
  DMatrix<double> loc = locFile.toEigen();
  
  GSRPDE<decltype(problem), fdaPDE::models::SpaceTimeParabolicTag,
	 Sampling::GeoStatLocations, SolverType::Monolithic, fdaPDE::models::Gamma>
    model(problem, time_mesh);
  model.setLambdaS(lambdaS);
  model.setLambdaT(lambdaT);

  // load data from .csv files
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile  ("data/models/GSRPDE/2D_test5/y.csv");
  DMatrix<double> y = yFile.toEigen();
  CSVFile<double> XFile; // design matrix
  XFile = reader.parseFile  ("data/models/GSRPDE/2D_test5/X.csv");
  DMatrix<double> X = XFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.stack (OBSERVATIONS_BLK,  y);
  df.insert(DESIGN_MATRIX_BLK, X);
  df.insert(SPACE_LOCATIONS_BLK, loc);
  model.setData(df);
  
  // define initial condition estimator over grid of lambdas
  InitialConditionEstimator ICestimator(model);
  std::vector<SVector<1>> lambdas;
  lambdas.push_back(SVector<1>(lambdaS)); 
  // compute estimate
  ICestimator.apply(lambdas);
  DMatrix<double> ICestimate = ICestimator.get();
  // set estimated initial condition
  model.setInitialCondition(ICestimate);

  // shift data one time instant forward
  std::size_t n = y.rows();
  model.setData(df.tail(n).extract());
  model.shift_time(1);
  
  model.init();
  model.solve();
  
  //   **  test correctness of computed results  **   

  DMatrix<double> computedF;
  computedF.resize((model.n_time()+1)*model.n_basis(), 1);
  computedF << model.s(), model.f();
  // estimate of spatial field \hat f
  SpMatrix<double> expectedSolution;
  Eigen::loadMarket(expectedSolution, "data/models/GSRPDE/2D_test5/sol_parabolic.mtx");
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution), computedF, std::pow(0.1, 10)) );

  // estimate of coefficient vector \hat \beta
  SpMatrix<double> expectedBeta;
  Eigen::loadMarket(expectedBeta, "data/models/GSRPDE/2D_test5/beta_parabolic.mtx");
  DVector<double> computedBeta = model.beta();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedBeta), computedBeta, std::pow(0.1, 10)) );
}

#include <gtest/gtest.h> // testing framework
#include <iomanip>
#include <unsupported/Eigen/SparseExtra>
#include <unsupported/Eigen/src/SparseExtra/MarketIO.h>

#include "../../fdaPDE/core/utils/Symbols.h"
#include "../../fdaPDE/core/utils/IO/CSVReader.h"
#include "../../fdaPDE/core/FEM/PDE.h"
#include "models/ModelBase.h"
using fdaPDE::core::FEM::PDE;
#include "../../fdaPDE/core/FEM/operators/SpaceVaryingFunctors.h"
using fdaPDE::core::FEM::SpaceVaryingDiffusion;
using fdaPDE::core::FEM::SpaceVaryingAdvection;
#include "core/MESH/Mesh.h"
#include "../../fdaPDE/models/regression/STRPDE.h"
using fdaPDE::models::STRPDE;
#include "../../fdaPDE/models/ModelTraits.h"
using fdaPDE::models::SolverType;
#include "../../fdaPDE/models/SamplingDesign.h"
using fdaPDE::models::Sampling;

#include "../utils/MeshLoader.h"
using fdaPDE::testing::MeshLoader;
#include "../utils/Constants.h"
using fdaPDE::testing::DOUBLE_TOLERANCE;
#include "../utils/Utils.h"
using fdaPDE::testing::almost_equal;

#include "../../fdaPDE/preprocess/InitialConditionEstimator.h"
using fdaPDE::preprocess::InitialConditionEstimator;

/* test 1
   domain:       unit square [1,1] x [1,1]
   sampling:     locations = nodes
   penalization: simple laplacian
   covariates:   no
   BC:           no
   order FE:     1
   time penalization: separable (mass penalization)
 */
TEST(STRPDE, Test1_Laplacian_NonParametric_GeostatisticalAtNodes_Separable_Monolithic) {
  // define time domain
  DVector<double> time_mesh;
  time_mesh.resize(11);
  std::size_t i = 0;
  for(double x = 0; x <= 2; x+=0.2, ++i) time_mesh[i] = x;
  
  // define spatial domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square_coarse");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3*time_mesh.rows(), 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE
  // define statistical model
  double lambdaS = 0.01; // smoothing in space
  double lambdaT = 0.01; // smoothing in time
  STRPDE<decltype(problem), fdaPDE::models::SpaceTimeSeparableTag,
	 Sampling::GeoStatMeshNodes, SolverType::Monolithic> model(problem, time_mesh);
  model.setLambdaS(lambdaS);
  model.setLambdaT(lambdaT);

  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/STRPDE/2D_test1/y.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.stack(OBSERVATIONS_BLK, y);
  model.setData(df);
  
  // // solve smoothing problem
  model.init();
  model.solve();

  //    **  test correctness of computed results  **   
  
  // // \Psi matrix
  SpMatrix<double> expectedPsi;
  Eigen::loadMarket(expectedPsi, "data/models/STRPDE/2D_test1/Psi.mtx");
  SpMatrix<double> computedPsi = model.Psi();
  EXPECT_TRUE( almost_equal(expectedPsi, computedPsi) );

  // R0 matrix (discretization of identity operator)
  SpMatrix<double> expectedR0;
  Eigen::loadMarket(expectedR0,  "data/models/STRPDE/2D_test1/R0.mtx");
  SpMatrix<double> computedR0 = model.R0();
  EXPECT_TRUE( almost_equal(expectedR0, computedR0) );
  
  // R1 matrix (discretization of differential operator)
  SpMatrix<double> expectedR1;
  Eigen::loadMarket(expectedR1,  "data/models/STRPDE/2D_test1/R1.mtx");
  SpMatrix<double> computedR1 = model.R1();
  EXPECT_TRUE( almost_equal(expectedR1, computedR1) );
  
  // estimate of spatial field \hat f
  SpMatrix<double> expectedSolution;
  Eigen::loadMarket(expectedSolution,   "data/models/STRPDE/2D_test1/sol.mtx");
  DMatrix<double> computedF = model.f();
  std::size_t N = computedF.rows();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution).topRows(N), computedF) );
}


/* test 2
   domain:       c-shaped
   sampling:     locations != nodes
   penalization: simple laplacian
   covariates:   yes
   BC:           no
   order FE:     1
   time penalization: separable (mass penalization)
 */
TEST(STRPDE, Test2_Laplacian_SemiParametric_GeostatisticalAtLocations_Separable_Monolithic) {
  // define time domain
  DVector<double> time_mesh;
  time_mesh.resize(5);
  for(std::size_t i = 0; i < 5; ++i)
    time_mesh[i] = (fdaPDE::testing::pi/4)*i;

  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("c_shaped");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE
  // define statistical model
  double lambdaS = 0.01; // smoothing in space
  double lambdaT = 0.01; // smoothing in time
  // load sample position
  CSVReader<double> reader{};
  CSVFile<double> locFile; // locations file
  locFile = reader.parseFile("data/models/STRPDE/2D_test2/locs.csv");
  DMatrix<double> loc = locFile.toEigen();

  STRPDE<decltype(problem), fdaPDE::models::SpaceTimeSeparableTag,
	 Sampling::GeoStatLocations, SolverType::Monolithic> model(problem, time_mesh);
  model.setLambdaS(lambdaS);
  model.setLambdaT(lambdaT);
  
  // load data from .csv files
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile  ("data/models/STRPDE/2D_test2/y.csv");
  DMatrix<double> y = yFile.toEigen();
  CSVFile<double> XFile; // design matrix
  XFile = reader.parseFile  ("data/models/STRPDE/2D_test2/X.csv");
  DMatrix<double> X = XFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.stack (OBSERVATIONS_BLK,  y);
  df.stack (DESIGN_MATRIX_BLK, X);
  df.insert(SPACE_LOCATIONS_BLK, loc);
  model.setData(df);
  
  // solve smoothing problem
  model.init();
  model.solve();

  //   **  test correctness of computed results  **   
  
  // \Psi matrix (sensible to locations != nodes)
  SpMatrix<double> expectedPsi;
  Eigen::loadMarket(expectedPsi, "data/models/STRPDE/2D_test2/Psi.mtx");
  SpMatrix<double> computedPsi = model.Psi();
  EXPECT_TRUE( almost_equal(expectedPsi, computedPsi) );

  // R0 matrix (discretization of identity operator)
  SpMatrix<double> expectedR0;
  Eigen::loadMarket(expectedR0,  "data/models/STRPDE/2D_test2/R0.mtx");
  SpMatrix<double> computedR0 = model.R0();
  EXPECT_TRUE( almost_equal(expectedR0, computedR0) );
  
  // R1 matrix (discretization of differential operator)
  SpMatrix<double> expectedR1;
  Eigen::loadMarket(expectedR1,  "data/models/STRPDE/2D_test2/R1.mtx");
  SpMatrix<double> computedR1 = model.R1();
  EXPECT_TRUE( almost_equal(expectedR1, computedR1) );

  // estimate of spatial field \hat f
  SpMatrix<double> expectedSolution;
  Eigen::loadMarket(expectedSolution, "data/models/STRPDE/2D_test2/sol.mtx");
  DMatrix<double> computedF = model.f();
  std::size_t N = computedF.rows();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution), computedF) );

  // estimate of coefficient vector \hat \beta
  SpMatrix<double> expectedBeta;
  Eigen::loadMarket(expectedBeta, "data/models/STRPDE/2D_test2/beta.mtx");
  DVector<double> computedBeta = model.beta();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedBeta), computedBeta) );
}


/* test 3
   domain:       quasicircular domain
   sampling:     areal
   penalization: non-costant coefficients PDE
   covariates:   no
   BC:           no
   order FE:     1
   time penalization: parabolic (monolithic solution)
 */
TEST(STRPDE, Test3_NonCostantCoefficientsPDE_NonParametric_Areal_Parabolic_Monolithic_EstimatedIC) {
  // define time domain, we skip the first time instant because we are going to use the first block of data
  // for the estimation of the initial condition
  DVector<double> time_mesh;
  time_mesh.resize(11);
  for(std::size_t i = 0; i < 10; ++i) time_mesh[i] = 0.4*i;

  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("quasi_circle");
  // load PDE coefficients data
  CSVReader<double> reader{};
  CSVFile<double> diffFile; // diffusion tensor
  diffFile = reader.parseFile("data/models/STRPDE/2D_test3/K.csv");
  DMatrix<double> diffData = diffFile.toEigen();
  CSVFile<double> adveFile; // transport vector
  adveFile = reader.parseFile("data/models/STRPDE/2D_test3/b.csv");
  DMatrix<double> adveData = adveFile.toEigen();

  // define non-constant coefficients
  SpaceVaryingDiffusion<2> diffCoeff;
  diffCoeff.setData(diffData);
  SpaceVaryingAdvection<2> adveCoeff;
  adveCoeff.setData(adveData);
  // parabolic PDE
  auto L = dT() + Laplacian(diffCoeff.asParameter()) + Gradient(adveCoeff.asParameter());
  
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, time_mesh.rows());
  
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE
  // define statistical model
  double lambdaS = std::pow(0.1, 6); // smoothing in space
  double lambdaT = std::pow(0.1, 6); // smoothing in time

  CSVReader<int> int_reader{};
  CSVFile<int> arealFile; // incidence matrix for specification of subdomains
  arealFile = int_reader.parseFile("data/models/STRPDE/2D_test3/incidence_matrix.csv");
  DMatrix<int> areal = arealFile.toEigen();
  
  STRPDE<decltype(problem), fdaPDE::models::SpaceTimeParabolicTag, Sampling::Areal,
	 SolverType::Monolithic> model(problem, time_mesh);
  model.setLambdaS(lambdaS);
  model.setLambdaT(lambdaT);

  // load data from .csv files
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile  ("data/models/STRPDE/2D_test3/y.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.stack (OBSERVATIONS_BLK, y);
  df.insert(SPACE_AREAL_BLK, areal);
  model.setData(df);

  // define initial condition estimator over grid of lambdas
  InitialConditionEstimator ICestimator(model);
  std::vector<SVector<1>> lambdas;
  for(double x = -9; x <= 3; x += 0.1) lambdas.push_back(SVector<1>(std::pow(10,x))); 
  // compute estimate
  ICestimator.apply(lambdas);
  DMatrix<double> ICestimate = ICestimator.get();
  // test computation initial condition
  CSVFile<double> ICfile;
  ICfile = reader.parseFile("data/models/STRPDE/2D_test3/IC.csv");  
  DMatrix<double> expectedIC = ICfile.toEigen();
  EXPECT_TRUE( almost_equal(expectedIC, ICestimate) );

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
  
  // estimate of spatial field \hat f (with estimatate of initial condition)
  SpMatrix<double> expectedSolution;
  Eigen::loadMarket(expectedSolution, "data/models/STRPDE/2D_test3/sol.mtx");
  std::size_t N = computedF.rows();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution).topRows(N), computedF) );
}

/* test 4
   domain:       unit square [1,1] x [1,1]
   sampling:     locations = nodes
   penalization: simple laplacian
   covariates:   no
   BC:           no
   order FE:     1
   time penalization: parabolic (iterative solver)
 */
TEST(STRPDE, Test4_Laplacian_NonParametric_GeostatisticalAtNodes_Parabolic_Iterative_EstimatedIC) {
  // define time domain
  DVector<double> time_mesh;
  time_mesh.resize(11);
  std::size_t i = 0;
  for(double x = 0; x <= 2; x+=0.2, ++i) time_mesh[i] = x;
  
  // define spatial domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square_coarse");
  auto L = dT() + Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, time_mesh.rows());
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE
  // define statistical model
  double lambdaS = 1;//0.01; // smoothing in space
  double lambdaT = 1;//0.01; // smoothing in time
  STRPDE<decltype(problem), fdaPDE::models::SpaceTimeParabolicTag,
	 Sampling::GeoStatMeshNodes, SolverType::Iterative> model(problem, time_mesh);
  model.setLambdaS(lambdaS);
  model.setLambdaT(lambdaT);

  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/STRPDE/2D_test4/y.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.stack(OBSERVATIONS_BLK, y);
  model.setData(df);

  // define initial condition estimator over grid of lambdas
  InitialConditionEstimator ICestimator(model);
  std::vector<SVector<1>> lambdas;
  for(double x = -9; x <= 3; x += 0.1) lambdas.push_back(SVector<1>(std::pow(10,x))); 
  // compute estimate
  ICestimator.apply(lambdas);

  DMatrix<double> ICestimate = ICestimator.get();

  // set estimated initial condition
  model.setInitialCondition(ICestimate);
  
  // shift data one time instant forward
  std::size_t n = y.rows();
  model.setData(df.tail(n).extract());
  model.shift_time(1);
  // set parameters for iterative method
  model.setTolerance(1e-4);
  model.setMaxIterations(50);

  // solve smoothing problem
  model.init();
  model.solve();

  //   **  test correctness of computed results  **   

  // estimate of spatial field \hat f
  SpMatrix<double> expectedSolution;
  Eigen::loadMarket(expectedSolution, "data/models/STRPDE/2D_test4/sol.mtx");
  DMatrix<double> computedF = model.f();
  std::size_t N = computedF.rows();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution).topRows(N), computedF) );
}


/* test 5
   domain:       unit square [1,1] x [1,1]
   sampling:     locations = nodes, time locations != time nodes
   penalization: simple laplacian
   covariates:   no
   BC:           no
   order FE:     1
   time penalization: separable (mass penalization)
 */
TEST(STRPDE, Test5_Laplacian_NonParametric_GeostatisticalAtNodes_TimeLocations_Separable_Monolithic) {
  // define time domain
  DVector<double> time_mesh;
  time_mesh.resize(11);
  std::size_t i = 0;
  for(double x = 0; x <= 2; x+=0.2, ++i) time_mesh[i] = x;
  
  // define spatial domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square_coarse");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3*time_mesh.rows(), 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE
  // define statistical model
  double lambdaS = 0.01; // smoothing in space
  double lambdaT = 0.01; // smoothing in time
  STRPDE<decltype(problem), fdaPDE::models::SpaceTimeSeparableTag,
	 Sampling::GeoStatMeshNodes, SolverType::Monolithic> model(problem, time_mesh);
  model.setLambdaS(lambdaS);
  model.setLambdaT(lambdaT);

  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/STRPDE/2D_test5/y.csv");
  DMatrix<double> y = yFile.toEigen();
  CSVFile<double> timeLocationsFile;
  timeLocationsFile = reader.parseFile("data/models/STRPDE/2D_test5/time_locations.csv");
  DMatrix<double> timeLocations = timeLocationsFile.toEigen();
  
  // set model data
  BlockFrame<double, int> df;
  df.stack(OBSERVATIONS_BLK, y);
  df.insert(TIME_LOCATIONS_BLK, timeLocations);
  model.setData(df);
  
  // // solve smoothing problem
  model.init();
  model.solve();

  //    **  test correctness of computed results  **   
  
  // estimate of spatial field \hat f
  SpMatrix<double> expectedSolution;
  Eigen::loadMarket(expectedSolution,   "data/models/STRPDE/2D_test5/sol.mtx");
  DMatrix<double> computedF = model.f();
  std::size_t N = computedF.rows();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution).topRows(N), computedF) );
}

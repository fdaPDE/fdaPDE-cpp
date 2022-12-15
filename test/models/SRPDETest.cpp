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

#include "../utils/MeshLoader.h"
using fdaPDE::testing::MeshLoader;
#include "../utils/Constants.h"
using fdaPDE::testing::DOUBLE_TOLERANCE;
#include "../utils/Utils.h"
using fdaPDE::testing::almost_equal;

/* test 1
   domain:       unit square [1,1] x [1,1]
   sampling:     locations = nodes
   penalization: simple laplacian
   covariates:   no
   BC:           no
   order FE:     1
 */
TEST(SRPDE, Test1_Laplacian_NonParametric_GeostatisticalAtNodes) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE
  problem.init();

  // define statistical model
  // use optimal lambda to avoid possible numerical issues
  double lambda = 5.623413 * std::pow(0.1, 5);
  SRPDE<decltype(problem), fdaPDE::models::Sampling::GeoStatMeshNodes>  model(problem);
  model.setLambda(lambda);
  
  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/SRPDE/2D_test1/z.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert("y", y);
  model.setData(df);

  // solve smoothing problem
  model.solve();

  /*   **  test correctness of computed results  **   */
  
  // \Psi matrix
  SpMatrix<double> expectedPsi;
  Eigen::loadMarket(expectedPsi, "data/models/SRPDE/2D_test1/Psi.mtx");
  SpMatrix<double> computedPsi = model.Psi();
  EXPECT_TRUE( almost_equal(expectedPsi, computedPsi) );

  // R0 matrix (discretization of identity operator)
  SpMatrix<double> expectedR0;
  Eigen::loadMarket(expectedR0,  "data/models/SRPDE/2D_test1/R0.mtx");
  SpMatrix<double> computedR0 = model.R0();
  EXPECT_TRUE( almost_equal(expectedR0, computedR0) );
  
  // R1 matrix (discretization of differential operator)
  SpMatrix<double> expectedR1;
  Eigen::loadMarket(expectedR1,  "data/models/SRPDE/2D_test1/R1.mtx");
  SpMatrix<double> computedR1 = model.R1();
  EXPECT_TRUE( almost_equal(expectedR1, computedR1) );
    
  // estimate of spatial field \hat f
  SpMatrix<double> expectedSolution;
  Eigen::loadMarket(expectedSolution,   "data/models/SRPDE/2D_test1/sol.mtx");
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
 */
TEST(SRPDE, Test2_Laplacian_SemiParametric_GeostatisticalAtLocations) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("c_shaped");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE
  problem.init();

  // define statistical model
  CSVReader<double> reader{};
  // load locations where data are sampled
  CSVFile<double> locFile;
  locFile = reader.parseFile("data/models/SRPDE/2D_test2/locs.csv");
  DMatrix<double> loc = locFile.toEigen();

  // use optimal lambda to avoid possible numerical issues
  double lambda = 0.2201047;
  SRPDE<decltype(problem), Sampling::GeoStatLocations> model(problem, loc);
  model.setLambda(lambda);
  
  // load data from .csv files
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile  ("data/models/SRPDE/2D_test2/z.csv");
  DMatrix<double> y = yFile.toEigen();
  CSVFile<double> XFile; // design matrix
  XFile = reader.parseFile  ("data/models/SRPDE/2D_test2/X.csv");
  DMatrix<double> X = XFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert("y", y);
  df.insert("X", X);
  model.setData(df);
  
  // solve smoothing problem
  model.solve();

  /*   **  test correctness of computed results  **   */
  
  // \Psi matrix (sensible to locations != nodes)
  SpMatrix<double> expectedPsi;
  Eigen::loadMarket(expectedPsi, "data/models/SRPDE/2D_test2/Psi.mtx");
  SpMatrix<double> computedPsi = model.Psi();
  EXPECT_TRUE( almost_equal(expectedPsi, computedPsi) );

  // R0 matrix (discretization of identity operator)
  SpMatrix<double> expectedR0;
  Eigen::loadMarket(expectedR0,  "data/models/SRPDE/2D_test2/R0.mtx");
  SpMatrix<double> computedR0 = model.R0();
  EXPECT_TRUE( almost_equal(expectedR0, computedR0) );

  // R1 matrix (discretization of differential operator)
  SpMatrix<double> expectedR1;
  Eigen::loadMarket(expectedR1,  "data/models/SRPDE/2D_test2/R1.mtx");
  SpMatrix<double> computedR1 = model.R1();
  EXPECT_TRUE( almost_equal(expectedR1, computedR1) );

  // estimate of spatial field \hat f
  SpMatrix<double> expectedSolution;
  Eigen::loadMarket(expectedSolution, "data/models/SRPDE/2D_test2/sol.mtx");
  DMatrix<double> computedF = model.f();
  std::size_t N = computedF.rows();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution).topRows(N), computedF) );

  // estimate of coefficient vector \hat \beta
  SpMatrix<double> expectedBeta;
  Eigen::loadMarket(expectedBeta, "data/models/SRPDE/2D_test2/beta.mtx");
  DVector<double> computedBeta = model.beta();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedBeta), computedBeta) );
}

/* test 3
   domain:       unit square [1,1] x [1,1]
   sampling:     locations = nodes
   penalization: costant coefficients PDE
   covariates:   no
   BC:           no
   order FE:     1
 */
TEST(SRPDE, Test3_CostantCoefficientsPDE_NonParametric_GeostatisticalAtNodes) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square");

  // non unitary diffusion tensor
  SMatrix<2> K;
  K << 1,0,0,4;
  auto L = Laplacian(K); // anisotropic diffusion
  
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE
  problem.init();

  // define statistical model
  double lambda = 10;
  SRPDE<decltype(problem), Sampling::GeoStatMeshNodes> model(problem);
  model.setLambda(lambda);
  
  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/SRPDE/2D_test3/z.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert("y", y);
  model.setData(df);

  // solve smoothing problem
  model.solve();

  /*   **  test correctness of computed results  **   */
  
  // \Psi matrix
  SpMatrix<double> expectedPsi;
  Eigen::loadMarket(expectedPsi, "data/models/SRPDE/2D_test3/Psi.mtx");
  SpMatrix<double> computedPsi = model.Psi();
  EXPECT_TRUE( almost_equal(expectedPsi, computedPsi) );

  // R0 matrix (discretization of identity operator)
  SpMatrix<double> expectedR0;
  Eigen::loadMarket(expectedR0,  "data/models/SRPDE/2D_test3/R0.mtx");
  SpMatrix<double> computedR0 = model.R0();
  EXPECT_TRUE( almost_equal(expectedR0, computedR0) );
  
  // R1 matrix (discretization of differential operator)
  SpMatrix<double> expectedR1;
  Eigen::loadMarket(expectedR1,  "data/models/SRPDE/2D_test3/R1.mtx");
  SpMatrix<double> computedR1 = model.R1();
  EXPECT_TRUE( almost_equal(expectedR1, computedR1) );
    
  // estimate of spatial field \hat f
  SpMatrix<double> expectedSolution;
  Eigen::loadMarket(expectedSolution, "data/models/SRPDE/2D_test3/sol.mtx");
  DMatrix<double> computedF = model.f();
  std::size_t N = computedF.rows();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution).topRows(N), computedF) );
}

/* test 4
   domain:       quasicircular domain
   sampling:     areal
   penalization: non-costant coefficients PDE
   covariates:   no
   BC:           yes
   order FE:     1
 */
TEST(SRPDE, Test4_NonCostantCoefficientsPDE_NonParametric_Areal) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("quasi_circle");

  // load PDE coefficients data
  CSVReader<double> reader{};
  CSVFile<double> diffFile; // diffusion tensor
  diffFile = reader.parseFile("data/models/SRPDE/2D_test4/K.csv");
  DMatrix<double> diffData = diffFile.toEigen();
  CSVFile<double> adveFile; // transport vector
  adveFile = reader.parseFile("data/models/SRPDE/2D_test4/b.csv");
  DMatrix<double> adveData = adveFile.toEigen();
  
  // define non-constant coefficients
  SpaceVaryingDiffusion<2> diffCoeff;
  diffCoeff.setData(diffData);
  SpaceVaryingAdvection<2> adveCoeff;
  adveCoeff.setData(adveData);

  auto L = Laplacian(diffCoeff.asParameter()) + Gradient(adveCoeff.asParameter());
  
  // load non-zero forcing term
  CSVFile<double> forceFile; // transport vector
  forceFile = reader.parseFile("data/models/SRPDE/2D_test4/force.csv");
  DMatrix<double> u = forceFile.toEigen();
  
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE
  problem.init();
  
  // define statistical model
  CSVReader<int> int_reader{};
  CSVFile<int> arealFile; // incidence matrix for specification of subdomains
  arealFile = int_reader.parseFile("data/models/SRPDE/2D_test4/incidence_matrix.csv");
  DMatrix<int> areal = arealFile.toEigen();

  double lambda = std::pow(0.1, 3);
  SRPDE<decltype(problem), Sampling::Areal> model(problem, areal);
  model.setLambda(lambda);
  
  // load data from .csv files
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/SRPDE/2D_test4/z.csv");
  DMatrix<double> y = yFile.toEigen();
  
  // set model data
  BlockFrame<double, int> df;
  df.insert("y", y);
  model.setData(df);
  
  // solve smoothing problem
  model.solve();

  /*   **  test correctness of computed results  **   */
 
  // \Psi matrix (sensible to areal sampling)
  SpMatrix<double> expectedPsi;
  Eigen::loadMarket(expectedPsi, "data/models/SRPDE/2D_test4/Psi.mtx");
  SpMatrix<double> computedPsi = model.Psi();
  EXPECT_TRUE( almost_equal(expectedPsi, computedPsi) );
  
  // R0 matrix (discretization of identity operator)
  SpMatrix<double> expectedR0;
  Eigen::loadMarket(expectedR0,  "data/models/SRPDE/2D_test4/R0.mtx");
  SpMatrix<double> computedR0 = model.R0();
  EXPECT_TRUE( almost_equal(expectedR0, computedR0) );
  
  // R1 matrix (discretization of differential operator)
  SpMatrix<double> expectedR1;
  Eigen::loadMarket(expectedR1,  "data/models/SRPDE/2D_test4/R1.mtx");
  SpMatrix<double> computedR1 = model.R1();
  EXPECT_TRUE( almost_equal(expectedR1, computedR1) );
  
  // u vector  (discretization of forcing term)
  SpMatrix<double> expectedU;
  Eigen::loadMarket(expectedU,   "data/models/SRPDE/2D_test4/u.mtx");
  DMatrix<double> computedU = model.u();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedU), computedU) );
 
  // estimate of spatial field \hat f
  SpMatrix<double> expectedSolution;
  Eigen::loadMarket(expectedSolution, "data/models/SRPDE/2D_test4/sol.mtx");
  DMatrix<double> computedF = model.f();
  std::size_t N = computedF.rows();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution).topRows(N), computedF) );
}

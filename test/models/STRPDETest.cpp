#include <Eigen/src/Core/util/Constants.h>
#include <cstddef>
#include <gtest/gtest.h> // testing framework
#include <unsupported/Eigen/SparseExtra>

#include "../../fdaPDE/core/utils/Symbols.h"
#include "../../fdaPDE/core/utils/IO/CSVReader.h"
#include "../../fdaPDE/core/FEM/PDE.h"
#include "models/iStatModel.h"
using fdaPDE::core::FEM::PDE;
#include "../../fdaPDE/core/FEM/operators/SpaceVaryingFunctors.h"
using fdaPDE::core::FEM::SpaceVaryingDiffusion;
using fdaPDE::core::FEM::SpaceVaryingAdvection;
#include "core/MESH/Mesh.h"
#include "../../fdaPDE/models/regression/STRPDE.h"
using fdaPDE::models::STRPDE;

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
   time penalization: separable (mass penalization)
 */
TEST(STRPDE, Test1_Laplacian_NonParametric_GeostatisticalAtNodes_Separable) {
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
  problem.init();

  // define statistical model
  double lambdaS = 0.01; // smoothing in space
  double lambdaT = 0.01; // smoothin in time
  STRPDE<decltype(problem), fdaPDE::models::SpaceTimeSeparable> model(problem, time_mesh);
  model.setLambdaS(lambdaS);
  model.setLambdaT(lambdaT);

  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/STRPDE/2D_test1/y.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.stack("y", y);
  model.setData(df);
  
  // // solve smoothing problem
  model.solve();

  // /*   **  test correctness of computed results  **   */
  
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
TEST(STRPDE, Test2_Laplacian_SemiParametric_GeostatisticalAtLocations_Separable) {
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
  problem.init();

  // define statistical model
  double lambdaS = 0.01; // smoothing in space
  double lambdaT = 0.01; // smoothin in time
  STRPDE<decltype(problem), fdaPDE::models::SpaceTimeSeparable> model(problem, time_mesh);
  model.setLambdaS(lambdaS);
  model.setLambdaT(lambdaT);
  
  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile  ("data/models/STRPDE/2D_test2/y.csv");
  DMatrix<double> y = yFile.toEigen();
  CSVFile<double> XFile; // design matrix
  XFile = reader.parseFile  ("data/models/STRPDE/2D_test2/X.csv");
  DMatrix<double> X = XFile.toEigen();
  CSVFile<double> locFile; // locations file
  locFile = reader.parseFile("data/models/STRPDE/2D_test2/locs.csv");
  DMatrix<double> loc = locFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.stack ("y", y);
  df.stack ("X", X);
  df.repeat("P", loc, time_mesh.rows(), 1);
  model.setData(df);
  
  // solve smoothing problem
  model.solve();

  /*   **  test correctness of computed results  **   */
  
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
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution).topRows(N), computedF) );

  // estimate of coefficient vector \hat \beta
  SpMatrix<double> expectedBeta;
  Eigen::loadMarket(expectedBeta, "data/models/STRPDE/2D_test2/beta.mtx");
  DVector<double> computedBeta = model.beta();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedBeta), computedBeta) );
}

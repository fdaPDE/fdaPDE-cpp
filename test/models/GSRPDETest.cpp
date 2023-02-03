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

#include "../utils/MeshLoader.h"
using fdaPDE::testing::MeshLoader;
#include "../utils/Constants.h"
using fdaPDE::testing::DOUBLE_TOLERANCE;
#include "../utils/Utils.h"
using fdaPDE::testing::almost_equal;

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
	 SolverType::Monolithic, fdaPDE::models::Poisson> model(problem, loc);
  model.setLambdaS(lambda);
  
  // load data from .csv files
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/GSRPDE/2D_test1/y.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert("y", y);
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
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution).topRows(N), computedF) );
}

// << (DMatrix<double>(expectedSolution) - computedF).lpNorm<Eigen::Infinity>();

/* test 2
   domain:       unit square [1,1] x [1,1]
   sampling:     locations != nodes
   penalization: simple laplacian
   covariates:   no
   BC:           no
   order FE:     1
   distribution: bernulli
 */
TEST(GSRPDE, Test2_Laplacian_NonParametric_GeostatisticalAtNodes_Bernulli) {
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
	 SolverType::Monolithic, fdaPDE::models::Bernulli> model(problem, loc);
  model.setLambdaS(lambda);
  
  // load data from .csv files
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/GSRPDE/2D_test2/y.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert("y", y);
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
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution).topRows(N), computedF) );
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
TEST(GSRPDE, Test3_Laplacian_NonParametric_GeostatisticalAtNodes_Exponential) {
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
	 SolverType::Monolithic, fdaPDE::models::Exponential> model(problem, loc);
  model.setLambdaS(lambda);
  
  // load data from .csv files
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/GSRPDE/2D_test3/y.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert("y", y);
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
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution).topRows(N), computedF) );
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
TEST(GSRPDE, Test4_Laplacian_NonParametric_GeostatisticalAtNodes_Gamma) {
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
	 SolverType::Monolithic, fdaPDE::models::Gamma> model(problem, loc);
  model.setLambdaS(lambda);
  
  // load data from .csv files
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/GSRPDE/2D_test4/y.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert("y", y);
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
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution).topRows(N), computedF) );
}

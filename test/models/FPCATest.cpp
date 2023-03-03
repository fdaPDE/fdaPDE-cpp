#include <cstddef>
#include <gtest/gtest.h> // testing framework
#include <unsupported/Eigen/SparseExtra>

#include "../fdaPDE/core/utils/Symbols.h"
#include "../fdaPDE/core/utils/IO/CSVReader.h"
#include "../fdaPDE/core/FEM/PDE.h"
using fdaPDE::core::FEM::PDE;
#include "../fdaPDE/core/FEM/EigenValueProblem.h"
using fdaPDE::core::FEM::EigenValueProblem;
#include "core/MESH/Mesh.h"
#include "../fdaPDE/models/functional/fPCA.h"
using fdaPDE::models::FPCA;
#include "../fdaPDE/models/SamplingDesign.h"
using fdaPDE::models::Sampling;
#include "../../fdaPDE/models/ModelTraits.h"
using fdaPDE::models::SpaceOnlyTag;
using fdaPDE::models::SolverType;

#include "../utils/MeshLoader.h"
using fdaPDE::testing::MeshLoader;
#include "../utils/Constants.h"
using fdaPDE::testing::DOUBLE_TOLERANCE;
#include "../utils/Utils.h"
using fdaPDE::testing::almost_equal;

#include "../fdaPDE/core/OPT/optimizers/Grid.h"

/* test 1
   domain:       unit square [1,1] x [1,1]
   sampling:     locations = nodes
   penalization: simple laplacian
   BC:           no
   order FE:     1
 */
/*TEST(FPCA, Test1_Laplacian_GeostatisticalAtNodes) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE

  // define statistical model
  // use optimal lambda to avoid possible numerical issues
  double lambda = 1e-2;
  FPCA<decltype(problem), SpaceOnlyTag, fdaPDE::models::Sampling::GeoStatMeshNodes,
       fdaPDE::models::gcv_lambda_selection> model(problem);
  model.setLambdaS(lambda);
  
  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/FPCA/2D_test1/y.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert("y", DMatrix<double>(y.transpose()));
  model.setData(df);

  std::vector<SVector<1>> lambdas;
  for(double x = -6.0; x <= -2.0; x++) lambdas.push_back(SVector<1>(std::pow(10,x)));
  model.setLambda(lambdas);
  
  // solve smoothing problem
  model.init();
  model.solve();
  
  //   **  test correctness of computed results  **
  
  // SpMatrix<double> expectedLoadings;
  // Eigen::loadMarket(expectedLoadings, "data/models/FPCA/2D_test1/loadings.mtx");
  // DMatrix<double> computedLoadings = model.loadings();
  // EXPECT_TRUE( almost_equal(DMatrix<double>(expectedLoadings), computedLoadings) );

  // SpMatrix<double> expectedScores;
  // Eigen::loadMarket(expectedScores,   "data/models/FPCA/2D_test1/scores.mtx");
  // DMatrix<double> computedScores = model.scores();
  // EXPECT_TRUE( almost_equal(DMatrix<double>(expectedScores), computedScores) );
  
  }*/


/* test 2
   domain:       unit square [1,1] x [1,1]
   sampling:     locations = nodes
   penalization: simple laplacian
   BC:           no
   order FE:     1
 */
/*TEST(FPCA, Test2_Laplacian_GeostatisticalAtNodes_Separable_Monolithic) {
  // define time domain
  DVector<double> time_mesh;
  time_mesh.resize(10);
  std::size_t i = 0;
  for(double x = 0.5; x <= 0.95; x+=0.05, ++i) time_mesh[i] = x;

  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square05");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE
  problem.init();

  // define statistical model
  // use optimal lambda to avoid possible numerical issues
  double lambdaS = 1e-2;
  double lambdaT = 1e-2;
  // defaults to monolithic solution
  FPCA<decltype(problem), fdaPDE::models::SpaceTimeSeparableTag,
       fdaPDE::models::Sampling::GeoStatMeshNodes, fdaPDE::models::fixed_lambda> model(problem, time_mesh);
  model.setLambdaS(lambdaS);
  model.setLambdaT(lambdaT);
  
  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/FPCA/2D_test2/y.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert("y", DMatrix<double>(y.transpose()));
  model.setData(df);

  // std::vector<SVector<1>> lambdas;
  // for(double x = -6.0; x <= -2.0; x++) lambdas.push_back(SVector<1>(std::pow(10,x)));
  // model.setLambda(lambdas);
  
  // solve smoothing problem
  model.init();
  model.solve();
  
  //   **  test correctness of computed results  **
  
  // SpMatrix<double> expectedLoadings;
  // Eigen::loadMarket(expectedLoadings, "data/models/FPCA/2D_test1/loadings.mtx");
  // DMatrix<double> computedLoadings = model.loadings();
  // EXPECT_TRUE( almost_equal(DMatrix<double>(expectedLoadings), computedLoadings) );

  // SpMatrix<double> expectedScores;
  // Eigen::loadMarket(expectedScores,   "data/models/FPCA/2D_test1/scores.mtx");
  // DMatrix<double> computedScores = model.scores();
  // EXPECT_TRUE( almost_equal(DMatrix<double>(expectedScores), computedScores) );
  
}*/

/* test 2
   domain:       unit square [1,1] x [1,1]
   sampling:     locations = nodes
   penalization: simple laplacian
   BC:           no
   order FE:     1
 */
TEST(FPCA, Test2_Laplacian_GeostatisticalAtNodes_Separable_Monolithic) {
  // define time domain
  DVector<double> time_mesh;
  time_mesh.resize(11);
  std::size_t i = 0;
  for(double x = 0; x <= 0.5; x+=0.05, ++i) time_mesh[i] = x;

  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square_coarse");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE

  // define statistical model
  // use optimal lambda to avoid possible numerical issues
  double lambdaS = 1e-4;
  double lambdaT = 1e-4;
  // defaults to monolithic solution
  FPCA<decltype(problem), fdaPDE::models::SpaceTimeSeparableTag,
       fdaPDE::models::Sampling::GeoStatMeshNodes, fdaPDE::models::gcv_lambda_selection> model;
  model.setPDE(problem);
  model.setTimeDomain(time_mesh);
  model.setLambdaS(lambdaS);
  model.setLambdaT(lambdaT);
  
  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/FPCA/2D_test3/y.csv");
  DMatrix<double> y = yFile.toEigen().leftCols(11*441);
  
  // set model data
  BlockFrame<double, int> df;
  df.insert("y", DMatrix<double>(y.transpose()));
  model.setData(df);

  std::vector<SVector<2>> lambdas;
  for(double x = -4.0; x <= -2.0; x+=0.5) {
    for(double y = -4.0; y <= -2.0; y+=0.5) {
      lambdas.push_back(SVector<2>(std::pow(10,x), std::pow(10,y)));
    }
  }
  model.setLambda(lambdas);
  
  // solve smoothing problem
  model.init();
  model.solve();

  //std::cout << model.loadings() << std::endl;
  
  /*   **  test correctness of computed results  **   */
  
  // SpMatrix<double> expectedLoadings;
  // Eigen::loadMarket(expectedLoadings, "data/models/FPCA/2D_test1/loadings.mtx");
  // DMatrix<double> computedLoadings = model.loadings();
  // EXPECT_TRUE( almost_equal(DMatrix<double>(expectedLoadings), computedLoadings) );

  // SpMatrix<double> expectedScores;
  // Eigen::loadMarket(expectedScores,   "data/models/FPCA/2D_test1/scores.mtx");
  // DMatrix<double> computedScores = model.scores();
  // EXPECT_TRUE( almost_equal(DMatrix<double>(expectedScores), computedScores) );
  
}

#include <cstddef>
#include <gtest/gtest.h> // testing framework
#include <unsupported/Eigen/SparseExtra>

#include "../fdaPDE/core/utils/Symbols.h"
#include "../fdaPDE/core/utils/IO/CSVReader.h"
#include "../fdaPDE/core/FEM/PDE.h"
using fdaPDE::core::FEM::PDE;
#include "core/MESH/Mesh.h"
#include "../fdaPDE/models/functional/fPCA.h"
using fdaPDE::models::FPCA;
#include "../fdaPDE/models/SamplingDesign.h"
using fdaPDE::models::Sampling;
#include "../../fdaPDE/models/ModelTraits.h"
using fdaPDE::models::SpaceOnlyTag;

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
   BC:           no
   order FE:     1
 */
TEST(FPCA, Test1_Laplacian_GeostatisticalAtNodes) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE

  // define statistical model
  // use optimal lambda to avoid possible numerical issues
  double lambda = 1e-2;
  FPCA<decltype(problem), SpaceOnlyTag, fdaPDE::models::Sampling::GeoStatMeshNodes> model(problem);
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

  // solve smoothing problem
  model.init();
  model.solve();
  
  /*   **  test correctness of computed results  **   */
  
  SpMatrix<double> expectedLoadings;
  Eigen::loadMarket(expectedLoadings, "data/models/FPCA/2D_test1/loadings.mtx");
  DMatrix<double> computedLoadings = model.loadings();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedLoadings), computedLoadings) );

  SpMatrix<double> expectedScores;
  Eigen::loadMarket(expectedScores,   "data/models/FPCA/2D_test1/scores.mtx");
  DMatrix<double> computedScores = model.scores();
  EXPECT_TRUE( almost_equal(DMatrix<double>(expectedScores), computedScores) );
  
}

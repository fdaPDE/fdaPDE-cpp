#include <cstddef>
#include <gtest/gtest.h> // testing framework
#include <limits>

#include "../fdaPDE/core/utils/Symbols.h"
#include "../fdaPDE/core/MESH/Mesh.h"
using fdaPDE::core::MESH::Mesh2D;
#include "../fdaPDE/core/MESH/Element.h"
using fdaPDE::core::MESH::Element;
#include "../fdaPDE/core/utils/IO/CSVReader.h"
#include "../fdaPDE/core/FEM/PDE.h"
using fdaPDE::core::FEM::PDE;
#include "../fdaPDE/core/FEM/operators/Laplacian.h"
using fdaPDE::core::FEM::Laplacian;
#include "../fdaPDE/core/FEM/basis/LagrangianBasis.h"
using fdaPDE::core::FEM::LagrangianBasis;
#include "../fdaPDE/core/FEM/integration/Integrator.h"
using fdaPDE::core::FEM::Integrator;
#include "../fdaPDE/regression/SRPDE.h"

// test suite for testing SRPDE model over 2D planar domain.
// Mesh used is the [1,1] x [1,1] unit square: see square_*.csv data series for inspect raw informations

class SRPDETest : public ::testing::Test {
public:
  Mesh2D<> m; // mesh
  double tolerance = 5*std::pow(0.1, 14);
  
  // load mesh from .csv files
  SRPDETest() {
    // compute file names
    std::string point    = "data/square_points.csv";
    std::string edges    = "data/square_edges.csv";
    std::string elements = "data/square_elements.csv";
    std::string neigh    = "data/square_neigh.csv";
    std::string boundary = "data/square_boundary.csv";
    // initialize test objects
    m = Mesh2D<>(point, edges, elements, neigh, boundary);
  };
};

// check that SRPDE correcty computes the field estimate
TEST_F(SRPDETest, Smooth2DNoCov) {
  // load expected parameters from file. Parameters are relative to the correct expected estimation of the field (written in R syntax)
  //
  // f = function(x, y, z = 1){
  //   coe = function(x,y) 1/2*sin(5*pi*x)*exp(-x^2)+1
  //   sin(2*pi*(coe(y,1)*x*cos(z-2)-y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2)+coe(x,1)*y*sin((z-2)*pi/2)))
  // }
  //
  // using L = Laplacian(), u = 0, \Omega = 2D unit square, \lambda = 0.001 and data generated according to the following R code
  //
  // # Exact solution (pointwise at nodes)
  // sol_exact = f(mesh$nodes[,1], mesh$nodes[,2])
  // # generate data with noise
  // set.seed(7893475)
  // ran = range(sol_exact)
  // nnodes = dim(mesh$nodes)[1]
  // data = sol_exact + rnorm(nnodes, mean=0, sd=0.05*abs(ran[2]-ran[1]))
  //
  // generated data are contained in file "data/SRPDE_2Dnocov_data.csv"
  
  CSVReader reader{}; // csv parser
  CSVFile<double> expected_f; // reference estimated parameters
  expected_f = reader.parseFile<double>("data/SRPDE_2Dnocov.csv");
  DVector<double> f = expected_f.toEigen();
  CSVFile<double> dataFile; // observation file
  dataFile = reader.parseFile<double>("data/SRPDE_2Dnocov_data.csv");
  DVector<double> data = dataFile.toEigen();
  
  // perform smoothing
  auto L = Laplacian(); 
  DMatrix<double> u; // forcing term u
  u.resize(m.elements(), 1);
  u.fill(0);
  // define differential problem
  PDE problem(m, L, u);
  DMatrix<double> bc; // set dirichelt BC
  bc.resize(m.elements(), 1);
  bc.fill(0);
  problem.setDirichletBC(bc);
  // define functional basis for FEM approximation of PDE solution
  LagrangianBasis<2, 2, 1> basis{};
  // define integrator for numerical approximations of integrals
  Integrator<2, 6> integrator{};
  // compute R1_ and R0_ matrices
  problem.init(basis, integrator);

  // define spatial regression model
  double lambda = 0.001;
  SRPDE model(problem, lambda);
  model.setObservations(data);
  model.smooth();
  DVector<double> estimated_f = model.f();

  // check estimated parameters equal to expected ones
  for(std::size_t i = 0; i < estimated_f.size(); ++ i)
    EXPECT_NEAR(f[i], estimated_f[i], tolerance);
}

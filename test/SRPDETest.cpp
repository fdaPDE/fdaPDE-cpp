#include <cstddef>
#include <gtest/gtest.h> // testing framework
#include <unsupported/Eigen/SparseExtra>

#include "../fdaPDE/core/utils/Symbols.h"
#include "../fdaPDE/core/utils/IO/CSVReader.h"
#include "../fdaPDE/core/FEM/PDE.h"
using fdaPDE::core::FEM::PDE;
#include "../fdaPDE/regression/SRPDE.h"

#include "utils/MeshLoader.h"
using fdaPDE::testing::MeshLoader;
#include "utils/Constants.h"
using fdaPDE::testing::DOUBLE_TOLERANCE;
using fdaPDE::testing::MACHINE_EPSILON;

// macro to define an SR-PDE model without covariates over a 2D [1,1] x [1,1] unit square domain ("unit_square" sample mesh)
// with L = Laplacian() (isotropic smoothing) , forcing term u = 0, homogeneous boundary conditions
#define DEF_2D_SRPDE_LAPLACIAN_NOCOV(lambda, domain, f_file, z_file)	   \
  /* load mesh file, vector of observations z and expected solution f*/    \
  MeshLoader<Mesh2D<>> UnitSquare(domain);      			   \
  DMatrix<double> f, z;							   \
  readCSV(f, f_file);			                                   \
  readCSV(z, z_file);							   \
  /* define SR-PDE model with simple laplacian regularization */	   \
  auto L = Laplacian();							   \
  DMatrix<double> u(UnitSquare.mesh.elements(), 1);			   \
  u.fill(0);								   \
  PDE problem(UnitSquare.mesh, L, u); /* definition of regularizing PDE*/  \
  LagrangianBasis<2, 2, 1> basis{};   /* definition of FEM basis */        \
  Integrator<2, 3> integrator{};					   \
  problem.init(basis, integrator);    /* compute R1 and R0 */		   \
  /* define model object, set observation vector and perfrom smoothing */  \
  SRPDE model(problem, lambda);	                                           \
  model.setObservations(z);	                                           \
  model.smooth();			                                   \

// check that SRPDE correcty computes the field estimate, case without covariates
TEST(SRPDE_2DNoCov_Test, PsiMatrix) {
  // define statistical model
  DEF_2D_SRPDE_LAPLACIAN_NOCOV(0.001, "unit_square",
			       "data/models/SRPDE_2Dnocov/f.csv", "data/models/SRPDE_2Dnocov/z.csv");
  // load \Psi matrix from file
  SpMatrix<double> expected;
  Eigen::loadMarket(expected, "data/models/SRPDE_2Dnocov/Psi.mtx");
  // request \Psi matrix from model
  SpMatrix<double> computed = *model.Psi();

  // move SpMatrix to DMatrix (to compute infinity norm between two matrices)
  DMatrix<double> x1 = expected;
  DMatrix<double> x2 = computed;
  EXPECT_TRUE((x1 - x2).lpNorm<Eigen::Infinity>() < MACHINE_EPSILON);
}

TEST(SRPDE_2DNoCov_Test, AMatrix) {
  // define statistical model
  DEF_2D_SRPDE_LAPLACIAN_NOCOV(0.001, "unit_square",
			       "data/models/SRPDE_2Dnocov/f.csv", "data/models/SRPDE_2Dnocov/z.csv");
  // load \Psi matrix from file
  SpMatrix<double> expected;
  Eigen::loadMarket(expected, "data/models/SRPDE_2Dnocov/A.mtx");
  // request \Psi matrix from model
  SpMatrix<double> computed = *model.A();

  // move SpMatrix to DMatrix (to compute infinity norm between two matrices)
  DMatrix<double> x1 = -expected; // current fdaPDE library computes A with flipped sign! (need to change here when this version becomes official)
  DMatrix<double> x2 =  computed;
  EXPECT_TRUE((x1 - x2).lpNorm<Eigen::Infinity>() < MACHINE_EPSILON);
}

TEST(SRPDE_2DNoCov_Test, bVector) {
  // define statistical model
  DEF_2D_SRPDE_LAPLACIAN_NOCOV(0.001, "unit_square",
			       "data/models/SRPDE_2Dnocov/f.csv", "data/models/SRPDE_2Dnocov/z.csv");
  // load \Psi matrix from file
  SpMatrix<double> expected;
  Eigen::loadMarket(expected, "data/models/SRPDE_2Dnocov/b.mtx");
  // request \Psi matrix from model
  DVector<double> computed = *model.b();

  // move SpMatrix to DMatrix (to compute infinity norm between two matrices)
  DVector<double> x1 = -expected; // current fdaPDE library computes A with flipped sign! (need to change here when this version becomes official)
  DVector<double> x2 =  computed;
  EXPECT_TRUE((x1 - x2).lpNorm<Eigen::Infinity>() < MACHINE_EPSILON);
}

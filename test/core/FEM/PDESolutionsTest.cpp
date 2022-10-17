// this test suite is totally devoted to high level testing of PDEs, that is FEM module is able to correctly
// discretize differential operators

#include <cstddef>
#include <gtest/gtest.h> // testing framework
#include <unsupported/Eigen/SparseExtra>

#include "../fdaPDE/core/utils/Symbols.h"
#include "../fdaPDE/core/FEM/PDE.h"
using fdaPDE::core::FEM::PDE;

#include "utils/MeshLoader.h"
using fdaPDE::testing::MeshLoader;
#include "utils/Constants.h"
using fdaPDE::testing::DOUBLE_TOLERANCE;
using fdaPDE::testing::MACHINE_EPSILON;

// test if FEM is able to correctly discretize the identity operator on some domain
TEST(PDESolutionsTest, Identity) {
  MeshLoader<Mesh2D<>> UnitSquare("unit_square");
  
  // define identity operator
  auto L = Identity();
  DMatrix<double> u(UnitSquare.mesh.elements()*3, 1);
  u.fill(0);
  LagrangianBasis<2, 2, 1> b{};   
  PDE problem(UnitSquare.mesh, L, b, u); // definition of PDE
  Integrator<2, 3> integrator{};

  // compute discretization matrix 
  problem.init(integrator);
  // load expected result
  SpMatrix<double> expected_R0;
  Eigen::loadMarket(expected_R0, "data/models/SRPDE_2Dnocov/R0.mtx");
  
  // request R0 matrix from pde
  SpMatrix<double> computed_R0 = *problem.R1();

  // move SpMatrix to DMatrix (to compute infinity norm between two matrices)
  DMatrix<double> R0_1 = expected_R0;
  DMatrix<double> R0_2 = computed_R0;
  EXPECT_TRUE((R0_1 - R0_2).lpNorm<Eigen::Infinity>() < DOUBLE_TOLERANCE);
}

// test if FEM is able to correctly discretize the laplacian operator on some domain
TEST(PDESolutionsTest, Laplacian) {
  MeshLoader<Mesh2D<>> UnitSquare("unit_square");
  
  // define identity operator
  auto L = Laplacian();
  DMatrix<double> u(UnitSquare.mesh.elements()*3, 1);
  u.fill(0);
  LagrangianBasis<2, 2, 1> b{};   
  PDE problem(UnitSquare.mesh, L, b, u); // definition of PDE
  Integrator<2, 3> integrator{};

  // compute discretization matrix 
  problem.init(integrator);
  // load expected result
  SpMatrix<double> expected_R0;
  Eigen::loadMarket(expected_R0, "data/models/SRPDE_2Dnocov/R1.mtx");
  
  // request R0 matrix from pde
  SpMatrix<double> computed_R0 = *problem.R1();

  // move SpMatrix to DMatrix (to compute infinity norm between two matrices)
  DMatrix<double> R0_1 = expected_R0;
  DMatrix<double> R0_2 = computed_R0;
  EXPECT_TRUE((R0_1 - R0_2).lpNorm<Eigen::Infinity>() < DOUBLE_TOLERANCE);
}


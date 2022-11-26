#include <cstddef>
#include <gtest/gtest-typed-test.h>
#include <gtest/gtest.h> // testing framework
#include <unsupported/Eigen/SparseExtra>

#include "../../fdaPDE/core/utils/Symbols.h"
#include "../../fdaPDE/core/MESH/Element.h"
using fdaPDE::core::MESH::Element;

#include "../utils/MeshLoader.h"
using fdaPDE::testing::MeshLoader;
#include "../utils/Constants.h"
using fdaPDE::testing::MODEL_TOLERANCE;

#include "../../fdaPDE/core/FEM/basis/LagrangianBasis.h"
using fdaPDE::core::FEM::LagrangianBasis;
#include "../../fdaPDE/core/FEM/operators/Laplacian.h"
using fdaPDE::core::FEM::Laplacian;
#include "../../fdaPDE/core/FEM/PDE.h"
using fdaPDE::core::FEM::PDE;

// test to globally check the correctness of PDE solver

TEST(PDESolutionsTest, FiniteQuadraticElementsGlobalTest) {
  // load sample mesh, request an order 2 basis support
  MeshLoader<Mesh2D<2>> UnitSquare("unit_square");

  // define differential problem
  auto L = Laplacian();
  PDE<2,2,2,decltype(L),DMatrix<double>> problem(UnitSquare.mesh);
  problem.setBilinearForm(L);

  // define forcing term
  DMatrix<double> u; // forcing term u
  u.resize(UnitSquare.mesh.elements()*problem.integrator().numNodes(), 1);

  constexpr double pi = 3.14159265358979323846;
  DMatrix<double> quad_points = problem.quadratureNodes();
  for(std::size_t i = 0; i < quad_points.rows(); ++i){
    u(i,0) = 8.0*pi*pi*std::sin(2.0*pi*quad_points(i,0))*std::sin(2.0*pi*quad_points(i,1));
  }
  problem.setForcing(u);

  // define boundary conditions
  DMatrix<double> b;
  b.resize(UnitSquare.mesh.dof(), 1);
  b.fill(0);
  problem.setDirichletBC(b);

  // solve PDE
  problem.init();
  problem.solve();

  // check expected solution
  SpMatrix<double> expectedSolution;
  Eigen::loadMarket(expectedSolution,  "data/PDEs/FEMorder2_reference_solution.mtx"); 
  DMatrix<double> computedSolution = problem.solution();
  EXPECT_TRUE( (DMatrix<double>(expectedSolution) - computedSolution).lpNorm<Eigen::Infinity>()
	       < MODEL_TOLERANCE);
}

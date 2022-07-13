#include <array>
#include <cmath>
#include <functional>
#include <gtest/gtest.h> // testing framework
#include <limits>
#include <utility>

#include "../src/core/utils/Symbols.h"
#include "../src/core/utils/fields/ScalarField.h"
using fdaPDE::core::ScalarField;
#include "../src/core/OPT/optimizers/GridOptimizer.h"
using fdaPDE::core::OPT::GridOptimizer;

TEST(GridOptimizerTest, OptimizeOver1DGrid) {
  // define objective function: x^3 + 2x^2
  ScalarField<1> field([](SVector<1> x) -> double { return std::pow(x[0], 3) + 2*std::pow(x[0], 2); });
  // define search domain
  auto interval = std::make_pair(-5, 5);
  // define optimizer
  GridOptimizer<1> g;
  // perform optimization
  g.findMinimum(field, {interval}, {0.01});
  // optim must be near to (-5)^3 + 2*5^2 = -75
  double objValue = g.getObjValue();
  EXPECT_TRUE(std::abs(objValue - (-75)) < std::pow(0.1, 8));
  // check optimum point is found: -5
  SVector<1> optPoint = g.getSolution();
  EXPECT_TRUE((optPoint - SVector<1>(-5)).norm() < std::pow(0.1, 4));
}

TEST(GridOptimizerTest, OptimizeOver2DGridConvex) {
  // define objective function: x^2 + y^2
  ScalarField<2> field([](SVector<2> x) -> double { return std::pow(x[0],2) + std::pow(x[1], 2); });
  // define square [-1,1] x [-1,1] domain
  std::array<std::pair<double,double>, 2> domain;
  domain[0] = std::make_pair(-1, 1);
  domain[1] = std::make_pair(-1, 1);
  // define optimizer
  GridOptimizer<2> g;
  // perform optimization, using uniform grid
  g.findMinimum(field, domain, 0.01);
  // optim should be near to 0, true global minimum of given field
  double objValue = g.getObjValue();
  EXPECT_TRUE(std::abs(objValue - 0) < std::pow(0.1, 8));
  // check optimum point is found
  SVector<2> optPoint = g.getSolution();
  EXPECT_TRUE((optPoint - SVector<2>(0,0)).norm() < std::pow(0.1, 4));
}

TEST(GridOptimizerTest, OptimizeOver2DGridNonConvex) {
  // define objective function: x*e^{-x^2 - y^2} + (x^2 + y^2)/20
  ScalarField<2> field([](SVector<2> x) -> double {
    return x[0]*std::exp(- std::pow(x[0],2) - std::pow(x[1], 2)) + (std::pow(x[0],2) + std::pow(x[1], 2))/20;
  });
  // define square [-1,1] x [-1,1] domain
  std::array<std::pair<double,double>, 2> domain;
  domain[0] = std::make_pair(-1, 1);
  domain[1] = std::make_pair(-1, 1);
  // define optimizer
  GridOptimizer<2> g;
  // perform optimization
  g.findMinimum(field, domain, {0.001, 0.001});
  // optim found by matlab using a quasi-newton unconstrained optimization: -0.405236870266690
  double expectedValue = -0.405236870266690;
  double objValue = g.getObjValue();
  EXPECT_TRUE(std::abs(objValue - expectedValue) < std::pow(0.1, 8)); // able to reach matlab result
  // check optimum point is found
  SVector<2> expectedMinimum(-0.669071831647573, 0.000000004602598);
  SVector<2> optPoint = g.getSolution();
  EXPECT_TRUE((optPoint - expectedMinimum).norm() < std::pow(0.1, 4));
}

// check if optimizer is able to optimize over an user defined set of points
TEST(GridOptimizerTest, OptimizeOverCustomdGrid) {
  // define objective function: x^2 + y^2
  ScalarField<2> field([](SVector<2> x) -> double { return std::pow(x[0],2) + std::pow(x[1], 2); });
  // define optimizer
  GridOptimizer<2> g;

  // define custom grid of points
  std::vector<SVector<2>> pointList{SVector<2>(0,0), SVector<2>(0.25, 0.25), SVector<2>(-0.25, 0.5),
				    SVector<2>(5,5), SVector<2>(-0.4, -0.7), SVector<2>(0,1)};
  // perform optimization
  g.findMinimum(field, pointList);
  // optim should be exactly zero, since the exact optimum is included in the point list
  double objValue = g.getObjValue();
  EXPECT_TRUE(std::abs(objValue - 0) == 0);
  // check optimum point is found
  SVector<2> optPoint = g.getSolution();
  EXPECT_TRUE((optPoint - SVector<2>(0,0)).norm() == 0);
}

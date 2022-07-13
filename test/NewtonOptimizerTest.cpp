#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <gtest/gtest.h> // testing framework
#include <limits>
#include <utility>

#include "../src/core/utils/Symbols.h"
#include "../src/core/utils/fields/ScalarField.h"
using fdaPDE::core::ScalarField;
using fdaPDE::core::TwiceDifferentiableScalarField;
#include "../src/core/OPT/optimizers/Newton.h"
using fdaPDE::core::OPT::NewtonOptimizer;

TEST(NewtonOptimizerTest, Optimize1D) {
  // define some objective to optimize
  ScalarField<1> field([](SVector<1> x) -> double { return std::sin(x[0]); });

  // define optimizer
  double tolerance = 0.0001;
  NewtonOptimizer<1> opt(1000, tolerance, 0.01);
  // perform optimization
  // newtons moves toward the direction reducing the gradient vector, is not guaranteed to converge to a minimum if the
  // function is not convex, select a starting point forcing newton to converge to the minimum of the function
  SVector<1> init(-0.1);
  opt.findMinimum(field, init);
  // optimum value should be near to expectedValue, at least below opt tolerance
  double expectedValue = -1;
  double objValue = opt.getObjValue();
  EXPECT_TRUE(std::abs(objValue - expectedValue) < tolerance) << objValue;
}

// optimization using numerical approximation of derivatives works
TEST(NewtonOptimizerTest, Optimize2DConvexNonExact) {
  // define objective function: x^2 + y^2
  ScalarField<2> field([](SVector<2> x) -> double { return std::pow(x[0],2) + std::pow(x[1], 2); });
  // do not define analytical gradient and hessian, instead resort to numerical approximations
  
  // define optimizer
  double tolerance = 0.0001;
  NewtonOptimizer<2> opt(1000, tolerance, 0.01);
  // perform optimization
  SVector<2> init(1,1);
  opt.findMinimum(field, init);
  // optimum value should be near to zero, at least below opt tolerance
  double objValue = opt.getObjValue();
  EXPECT_TRUE(std::abs(objValue - 0) < tolerance);
}

// check that calling .findMinimum after a previous call to .findMinimum returns the same result, this means to check
// that the optimizer is left to a valid state after execution of .findMinimum
TEST(NewtonOptimizerTest, SubsequentCallsSameResult) {
  // define objective function: x^2 + y^2
  ScalarField<2> field([](SVector<2> x) -> double { return std::pow(x[0],2) + std::pow(x[1], 2); });
  // do not define analytical gradient and hessian, instead resort to numerical approximations
  
  // define optimizer
  double tolerance = 0.0001;
  NewtonOptimizer<2> opt(1000, tolerance, 0.01);
  // perform optimization
  SVector<2> init(1,1);
  opt.findMinimum(field, init);
  // optimum value should be near to zero, at least below opt tolerance
  double objValue = opt.getObjValue();
  EXPECT_TRUE(std::abs(objValue - 0) < tolerance);

  // optimize again using same datum, the solution found must be equal
  opt.findMinimum(field, init);
  double objValue2 = opt.getObjValue();
  EXPECT_DOUBLE_EQ(objValue, objValue2);

  // if the above assertion is not met, most likely the optimizer doesn't clean its state after a call to .findMinimum
  // you should always clean the state of the optimizer after the optimization
}

// supply an objective with exact expression of gradient and hessian triggers the optimizer to use these quantites
TEST(NewtonOptimizerTest, Optimize2DConvexExact) {
  // define objective function: x^2 + y^2
  std::function<double(SVector<2>)> fieldExpr = [](SVector<2> x) -> double {
    return std::pow(x[0],2) + std::pow(x[1], 2);
  };
  // define analytical gradient vector
  std::function<double(SVector<2>)> dx = [](SVector<2> x) -> double { return 2*x[0]; };
  std::function<double(SVector<2>)> dy = [](SVector<2> x) -> double { return 2*x[1]; };
  // define analytical hessian matrix
  std::function<SMatrix<2>(SVector<2>)> hess = [](SVector<2> x) -> SMatrix<2> {
    return SMatrix<2> {{2,0}, {0,2}};
  };
  // field approximated using central differences
  ScalarField<2> field_approx(fieldExpr);
  // define the same field but supplying both analytical gradient and hessian
  TwiceDifferentiableScalarField<2> field(fieldExpr, {dx, dy}, hess);
  
  // define optimizer
  double tolerance = 0.0001;
  NewtonOptimizer<2> opt(1000, tolerance, 0.01);
  // perform optimization using approximation 
  SVector<2> init(1,1);
  opt.findMinimum(field_approx, init);
  double objValue_approx = opt.getObjValue();
  std::size_t it_approx  = opt.iterations();
  // perform optimization using analytical gradient and hessian
  opt.findMinimum(field, init);
  double objValue_exact = opt.getObjValue();
  std::size_t it_exact  = opt.iterations();
  
  // expect objValue_exact below tolerance
  EXPECT_TRUE(std::abs(objValue_exact - 0) < tolerance);
  // expect the two values different (this means that the analytical expressions are actually used by the method)
  EXPECT_FALSE(objValue_exact == objValue_approx);
  // expect solution obtained with exact expressions approximately equal to the approximated one
  EXPECT_TRUE(std::abs(objValue_exact - objValue_approx) < tolerance);
}

// use a non convex objective to stress the optimizer but with a unique minimum
TEST(NewtonOptimizerTest, Optimize2DNonConvexUniqueMin) {
  // define non convex objective function having just one minimum : x*e^{-x^2 - y^2} + (x^2 + y^2)/20
  ScalarField<2> field([](SVector<2> x) -> double {
    return x[0]*std::exp(- std::pow(x[0],2) - std::pow(x[1], 2)) + (std::pow(x[0],2) + std::pow(x[1], 2))/20;
  });
  // do not define analytical gradient and hessian, instead resort to numerical approximations
  
  // define optimizer
  double tolerance = 0.0001;
  NewtonOptimizer<2> opt(1000, tolerance, 0.01);
  // perform optimization
  SVector<2> init(-0.5, 0);
  opt.findMinimum(field, init);
  // optimum value should be near to expectedValue, at least below opt tolerance
  double expectedValue = -0.405236870266690;
  double objValue = opt.getObjValue();
  EXPECT_TRUE(std::abs(objValue - expectedValue) < tolerance);
}

// check that if tolerance is not reached maximum iteration is met
TEST(NewtonOptimizerTest, StopOnMaxIter) {
  // select a convex function with a single minimum point: 2*x^2 + 4*y^2
  ScalarField<2> field([](SVector<2> x) -> double { return 2*std::pow(x[0],2) + 4*std::pow(x[1], 2); });

  // define optimizer, set tolerance too low to be reached
  double tolerance = std::pow(0.1, 10);
  NewtonOptimizer<2> opt(1000, tolerance, 0.01);
  // perform optimization
  SVector<2> init(1,1);
  opt.findMinimum(field, init);
  // optimum value should be near to zero, at least below opt tolerance. If this condition is not met, the maximum
  // number of iterations must be reached
  double objValue = opt.getObjValue();
  EXPECT_TRUE(std::abs(objValue - 0) < tolerance || opt.iterations() == 1000);
}

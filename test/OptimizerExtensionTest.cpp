#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <gtest/gtest-typed-test.h>
#include <gtest/gtest.h> // testing framework

#include "../fdaPDE/core/utils/Symbols.h"
#include "../fdaPDE/core/utils/fields/ScalarField.h"
using fdaPDE::core::ScalarField;
using fdaPDE::core::TwiceDifferentiableScalarField;
#include "../fdaPDE/core/OPT/optimizers/BFGS.h"
using fdaPDE::core::OPT::BFGSOptimizer;
#include "../fdaPDE/core/OPT/optimizers/Newton.h"
using fdaPDE::core::OPT::NewtonOptimizer;
#include "../fdaPDE/core/OPT/optimizers/GradientDescent.h"
using fdaPDE::core::OPT::GradientDescentOptimizer;
#include "../fdaPDE/core/OPT/optimizers/GridOptimizer.h"
using fdaPDE::core::OPT::GridOptimizer;
#include "../fdaPDE/core/OPT/extensions/BacktrackingAdaptiveStep.h"
using fdaPDE::core::OPT::BacktrackingAdaptiveStep;

// test if the extension mechanism works as expected.

// create a test extension
struct TestExtension {
  std::size_t counter_ = 0; // internal counter
  std::size_t forcedStop_ = 0; // flag used to stop the optimizer at specific points
  
  // constructor
  TestExtension() = default;
  TestExtension(std::size_t forcedStop) : forcedStop_(forcedStop) {};
  
  // modifiy the counter at different points of the optimization process.
  // given a fixed number of iterations we expect at the end a precise counter value
  template <typename Optimizer, typename Objective>
  bool initOptimization(Optimizer& opt, Objective& obj){
    // executed when optimization starts
    if(forcedStop_ == 1) return true;
    counter_ += 100;
    return false;
  }
  template <typename Optimizer, typename Objective>
  bool initIteration(Optimizer& opt, Objective& obj){
    // executed at the very beginning of each iteration of the iterative method
    if(forcedStop_ == 2 && opt.iterations() == 50) return true;
    counter_ += 2;
    return false;
  }
  template <typename Optimizer, typename Objective>
  bool endIteration(Optimizer& opt, Objective& obj){
    // executed at the very end of each iteration
    if(forcedStop_ == 3) return true;
    if(opt.iterations() < 50)
      counter_--;
    else
      counter_ -= 2;
    return false;
  }
  template <typename Optimizer, typename Objective>
  bool endOptimization(Optimizer& opt, Objective& obj){
    // executed just before the "finalize optimization" step
    if(forcedStop_ == 4){
      counter_ = 1;
      return true;
    }
    counter_ += 5;
    return false;
  }
};

// create templated test fixture for 2D case
template <typename T>
class OptimizerExtensionTest : public ::testing::Test {};

using optList = ::testing::Types<BFGSOptimizer<2>, NewtonOptimizer<2>, GradientDescentOptimizer<2>>;
TYPED_TEST_SUITE(OptimizerExtensionTest, optList);

// check if optimizer is capable to handle an extension and properly change its internal state to collect some result.
// If an optimizer fails this test most likely is receiving extensions as const references or is not handling them in the
// correct way
TYPED_TEST(OptimizerExtensionTest, OptimizerAcceptExtension) {
  // define a dumb objective
  ScalarField<2> field([](SVector<2> x) -> double { return std::pow(x[0], 2) + std::pow(x[1], 2); });

  // define optimizer
  double tolerance = 0.0001;
  TypeParam opt(100, tolerance, 0.01);
  // call findMinimum passing a TestExtension object
  SVector<2> init(1,1);
  TestExtension test_ext;
  opt.findMinimum(field, init, test_ext);
  // after 100 iterations, the value for the test counter must be equal to: 100 + 50 + 5 = 155
  EXPECT_TRUE(test_ext.counter_ == 155);
}

// check if optimizer stops when extension reaches a stopping condition
TYPED_TEST(OptimizerExtensionTest, ForcedStop) {
  // define a dumb objective
  ScalarField<2> field([](SVector<2> x) -> double { return std::pow(x[0], 2) + std::pow(x[1], 2); });

  // define optimizer
  double tolerance = 0.0001;
  TypeParam opt(100, tolerance, 0.01);
  // call findMinimum passing a TestExtension object
  SVector<2> init(1,1);
  TestExtension test_ext1(1);
  opt.findMinimum(field, init, test_ext1);
  // the optimizer stops at iteration 0,   expected counter: 5
  EXPECT_TRUE(test_ext1.counter_ == 5);

  TestExtension test_ext2(2);
  opt.findMinimum(field, init, test_ext2);
  // the optimizer stops at iteration 50,  expected counter: 100 + (2-1)*50 - 2 + 5 = 153
  EXPECT_TRUE(test_ext2.counter_ == 153);

  TestExtension test_ext3(3);
  opt.findMinimum(field, init, test_ext3);
  // the optimizer stops at iteration 1,   expected counter: 100 + 2 + 5 = 107
  EXPECT_TRUE(test_ext3.counter_ == 107);

  TestExtension test_ext4(4);
  opt.findMinimum(field, init, test_ext4);
  // the optimizer stops at iteration 100, expected counter: 1
  EXPECT_TRUE(test_ext4.counter_ == 1);
}

// following tests are specific to extensions

// check if BacktrackingAdaptiveStep produces same result as not Backtracked method, it should just be faster not
// produce different results. Tested on GradientDescent optimizer.
TEST(ExtensionTest, BacktrackingAdaptiveStep) {
  // define non convex objective function having just one minimum : x*e^{-x^2 - y^2} + (x^2 + y^2)/20
  ScalarField<2> field([](SVector<2> x) -> double {
    return x[0]*std::exp(- std::pow(x[0],2) - std::pow(x[1], 2)) + (std::pow(x[0],2) + std::pow(x[1], 2))/20;
  });
  // do not define analytical gradient and hessian, instead resort to numerical approximations
  
  // define optimizer
  double tolerance = 0.0001;
  GradientDescentOptimizer<2> opt(1000, tolerance, 0.01);
  // perform optimization
  SVector<2> init(-0.5, 0);
  opt.findMinimum(field, init);
  // optimum value should be near to expectedValue, at least below opt tolerance
  double expectedValue = -0.405236870266690;
  double objValue = opt.getObjValue();
  EXPECT_TRUE(std::abs(objValue - expectedValue) < tolerance);

  // repeat optimization but now with backtracking on (parameter selected to reach required tolerance)
  BacktrackingAdaptiveStep backtrackExt(1, 0.6, 0.5);
  opt.findMinimum(field, init, backtrackExt);
  // optimum value should be near to expectedValue, at least below opt tolerance
  double objValue_backtracked = opt.getObjValue();
  EXPECT_TRUE(std::abs(objValue_backtracked - expectedValue) < tolerance);

  // we also expect the optima found to be near
  EXPECT_TRUE(std::abs(objValue_backtracked - objValue) < tolerance);
}

// check if the extensions leave the optimizer in a stable state. If this test is not passed most likely the extension
// doesn't restore the optimizer original status after a call to .findMinimum()

// create templated test fixture for extensions
template <typename T>
class ExtensionTest : public ::testing::Test {};

using extensionList = ::testing::Types<BacktrackingAdaptiveStep>;
TYPED_TEST_SUITE(ExtensionTest, extensionList);

TYPED_TEST(ExtensionTest, OptimizerStatusRestored) {
  // define convex objective
  ScalarField<2> field([](SVector<2> x) -> double { return 2*std::pow(x[0],2) + 4*std::pow(x[1], 2); });

  // select an optimizer, we use a GradientDescentOptimizer just as footprint for a general iterative optimizer
  double tolerance = 0.0001;
  GradientDescentOptimizer<2> opt(1000, tolerance, 0.01);
  // perform optimization
  SVector<2> init(-0.5, 0);
  opt.findMinimum(field, init);
  // record obtained solution
  double objValue_before = opt.getObjValue();
  
  // don't care what the extension does in practice, nor the result produced (is not goal of this test to check the
  // correctness of the extension) (extension must be default constructable)
  TypeParam ext;
  opt.findMinimum(field, init, ext);
  
  // now reoptimize same objective without extension on
  opt.findMinimum(field, init);
  double objValue_after = opt.getObjValue();

  // an optimizer is deterministic, the objValue before and after must be exactly equal
  EXPECT_DOUBLE_EQ(objValue_after, objValue_before);
  // if the above assertion is not met, the extension ext did something which changed the initial optimizer configuration, bad.
}

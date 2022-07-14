#include <gtest/gtest.h>

// fields test suites
#include "ScalarFieldTest.cpp"
#include "VectorFieldTest.cpp"
// OPT test suites
#include "GridOptimizerTest.cpp"
//#include "NewtonOptimizerTest.cpp"
//#include "BFGSOptimizerTest.cpp"
#include "IterativeOptimizerTest.cpp"

int main(int argc, char **argv){
  // start testing
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

#include <gtest/gtest.h>

// fields test suites
#include "ScalarFieldTest.cpp"
#include "VectorFieldTest.cpp"
// OPT test suites
#include "GridOptimizerTest.cpp"
#include "IterativeOptimizerTest.cpp"
#include "OptimizerExtensionTest.cpp"
// MESH test suites
#include "MeshTest.cpp"
#include "ElementTest.cpp"

int main(int argc, char **argv){
  // start testing
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

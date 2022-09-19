#include <gtest/gtest.h>

// fields test suites
#include "ScalarFieldTest.cpp"
#include "VectorFieldTest.cpp"
// OPT test suites
// #include "GridOptimizerTest.cpp"
// #include "IterativeOptimizerTest.cpp"
// #include "OptimizerExtensionTest.cpp"
// MESH test suites
#include "core/MESH/MeshTest.cpp"
#include "core/MESH/ElementTest.cpp"
// #include "SearchEngineTest.cpp"
// // NLA test suites
// #include "FSPAITest.cpp"
// #include "VectorSpaceTest.cpp"
// // FEM test suites
// #include "LagrangianBasisTest.cpp"
// #include "IntegratorTest.cpp"
// regression module test suites
// #include "SRPDETest.cpp"

int main(int argc, char **argv){
  // start testing
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

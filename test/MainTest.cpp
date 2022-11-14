#include <gtest/gtest.h> // testing framework
// include eigen now to avoid possible linking errors
#include <Eigen/Dense>
#include <Eigen/Sparse>

// fields test suites
#include "core/ScalarFieldTest.cpp"
#include "core/VectorFieldTest.cpp"
// OPT test suites
// #include "GridOptimizerTest.cpp"
// #include "IterativeOptimizerTest.cpp"
// #include "OptimizerExtensionTest.cpp"
// MESH test suites
#include "core/MeshTest.cpp"
#include "core/ElementTest.cpp"
#include "core/SearchEngineTest.cpp"
// NLA test suites
#include "core/FSPAITest.cpp"
#include "core/VectorSpaceTest.cpp"
// FEM test suites
#include "core/LagrangianBasisTest.cpp"
#include "core/IntegratorTest.cpp"
// regression module test suites
#include "models/SRPDETest.cpp"

int main(int argc, char **argv){
  // start testing
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

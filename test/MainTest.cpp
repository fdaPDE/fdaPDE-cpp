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
#include "core/MESH/MeshTest.cpp"
#include "core/MESH/ElementTest.cpp"
#include "core/MESH/SearchEngineTest.cpp"
// NLA test suites
#include "core/NLA/FSPAITest.cpp"
#include "core/NLA/VectorSpaceTest.cpp"
// // FEM test suites
#include "core/FEM/LagrangianBasisTest.cpp"
// #include "IntegratorTest.cpp"
// regression module test suites
// #include "SRPDETest.cpp"

int main(int argc, char **argv){
  // start testing
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

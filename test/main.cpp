#include <gtest/gtest.h> // testing framework
// include eigen now to avoid possible linking errors
#include <Eigen/Dense>
#include <Eigen/Sparse>

// regression test suites
#include "src/srpde_test.cpp"
#include "src/strpde_test.cpp"
#include "src/gsrpde_test.cpp"
#include "src/sqrpde_test.cpp"
#include "src/fpca_test.cpp"
#include "src/sqrpde_time_test.cpp"
// GCV test suites
#include "src/gcv_srpde_test.cpp"
#include "src/gcv_sqrpde_test.cpp"
#include "src/gcv_sqrpde_time_test.cpp"
// #include "calibration/GCVNewtonTest.cpp"

int main(int argc, char **argv){
  // start testing
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

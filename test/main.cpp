#include <gtest/gtest.h> // testing framework
// include eigen now to avoid possible linking errors
#include <Eigen/Dense>
#include <Eigen/Sparse>

// // regression test suites
// #include "src/srpde_test.cpp"
// #include "src/strpde_test.cpp"
// #include "src/gsrpde_test.cpp"
// #include "src/qsrpde_test.cpp"
// #include "src/stqrpde_test.cpp"
// #include "src/msqrpde_test.cpp"
// #include "src/fpca_test.cpp"
// // GCV test suites
// #include "src/gcv_srpde_test.cpp"
// #include "src/gcv_qsrpde_test.cpp"
// #include "src/gcv_stqrpde_test.cpp"
// #include "src/gcv_msqrpde_test.cpp"
// #include "src/gcv_srpde_newton_test.cpp"
// // Case study 
#include  "src/case_study.cpp"

int main(int argc, char **argv){
  // start testing
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

#include <gtest/gtest.h> // testing framework
// include eigen now to avoid possible linking errors
#include <Eigen/Dense>
#include <Eigen/Sparse>

// // regression test suite
// #include "src/srpde_test.cpp"
// #include "src/strpde_test.cpp"
// #include "src/gsrpde_test.cpp"
// #include "src/qsrpde_test.cpp"
// #include "src/gcv_srpde_test.cpp"
// #include "src/gcv_qsrpde_test.cpp"
// #include "src/gcv_srpde_newton_test.cpp"
// #include "src/kcv_srpde_test.cpp"
// // functional test suite
// #include "src/fpca_test.cpp"
// // regression test suites
// #include "src/qstrpde_test.cpp"   
// #include "src/mqsrpde_test.cpp"
// GCV test suites
// #include "src/gcv_qstrpde_test.cpp"  
// #include "src/gcv_mqsrpde_test.cpp"
// // Case study 
#include  "src/case_study.cpp"

int main(int argc, char **argv){
  // start testing
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

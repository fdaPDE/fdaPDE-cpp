#ifndef __UTILS_H__
#define __UTILS_H__

#include <fdaPDE/utils.h>
#include <unsupported/Eigen/SparseExtra>
#include <string>
#include "constants.h"

// a set of usefull utilities
namespace fdapde {
namespace testing {

  // this function is an implementation of the test for floating point equality based on relative error. There is
  // an huge literature about floating point comparison, refer to it for details
  template <typename T>
  typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
  almost_equal(T a, T b, T epsilon){
    return std::fabs(a-b) < epsilon ||
      std::fabs(a-b) < ((std::fabs(a) < std::fabs(b) ? std::fabs(b) : std::fabs(a)) * epsilon);
  }

  // set default epsilon to DOUBLE_TOLERANCE
  template <typename T>
  typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
  almost_equal(T a, T b){ return almost_equal(a,b, DOUBLE_TOLERANCE); }

  // test if two matrices are equal testing the relative error of the infinte norm of their difference
  inline bool almost_equal(const DMatrix<double>& op1, const DMatrix<double>& op2, double epsilon){
    return (op1-op2).lpNorm<Eigen::Infinity>() < epsilon ||
      (op1-op2).lpNorm<Eigen::Infinity>() < (std::max(op1.lpNorm<Eigen::Infinity>(), op2.lpNorm<Eigen::Infinity>()) * epsilon);
  }
  inline bool almost_equal(const DMatrix<double>& op1, const DMatrix<double>& op2){
    return almost_equal(op1, op2, DOUBLE_TOLERANCE);
  }
  // sparse operands
  inline bool almost_equal(const SpMatrix<double>& op1, const SpMatrix<double>& op2, double epsilon){
    return almost_equal(DMatrix<double>(op1), DMatrix<double>(op2), epsilon);
  }
  inline bool almost_equal(const SpMatrix<double>& op1, const SpMatrix<double>& op2){
    return almost_equal(DMatrix<double>(op1), DMatrix<double>(op2));
  }

  // load rhs from file
  bool almost_equal(const SpMatrix<double>& op1, std::string op2) {
      SpMatrix<double> mem_buff;
      Eigen::loadMarket(mem_buff, op2);
      return almost_equal(op1, mem_buff);
  }
  bool almost_equal(const DMatrix<double>& op1, std::string op2) {
      SpMatrix<double> mem_buff;
      Eigen::loadMarket(mem_buff, op2);
      return almost_equal(op1, DMatrix<double>(mem_buff));
  }

  bool almost_equal(const std::vector<double>& op1, std::string op2){
    DMatrix<double> m;
    m.resize(op1.size(), 1);
    for(std::size_t i = 0; i < op1.size(); ++i) m(i,0) = op1[i];
    return almost_equal(m, op2);
  }

  // utility to import .mtx files
  template <typename T>
  DMatrix<T> read_mtx(const std::string& file_name) {
    SpMatrix<double> buff;
    Eigen::loadMarket(buff, file_name);
    return buff;
  }
  // utility to import .csv files
  template <typename T>
  DMatrix<T> read_csv(const std::string& file_name) {
    CSVReader<T> reader {};
    return reader.template parse_file<DenseStorage>(file_name);
  }
  
}}

#endif // __UTILS_H__

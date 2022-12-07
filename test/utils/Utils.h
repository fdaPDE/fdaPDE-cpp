#ifndef __UTILS_H__
#define __UTILS_H__

#include <limits>
#include "../fdaPDE/core/utils/Symbols.h"

// a set of usefull utilities
namespace fdaPDE{
namespace testing{

  // test if two matrices are equal 
  double spLInfinityNorm(const SpMatrix<double>& op1, const SpMatrix<double>& op2){
    // convert sparse operands into dense ones
    DMatrix<double> d1 = op1, d2 = op2;
    return (d1 - d2).lpNorm<Eigen::Infinity>();
  } 
  
}}

#endif // __UTILS_H__

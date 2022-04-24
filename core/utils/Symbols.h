#ifndef __SYMBOLS_H__
#define __SYMBOLS_H__

// Common symbols and data types used in the Core library
#include <Eigen/Core>

template <unsigned int N> using SVector = Eigen::Matrix<double, N, 1>;
template <unsigned int N> using SMatrix = Eigen::Matrix<double, N, N>;

#endif // __SYMBOLS_H__

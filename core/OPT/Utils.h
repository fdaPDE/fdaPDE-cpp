#ifndef __UTILS_H__
#define __UTILS_H__

#include <Eigen/Core>

template <unsigned int N> using SVector = Eigen::Matrix<double, N, 1>;
template <unsigned int N> using SMatrix = Eigen::Matrix<double, N, N>;

#endif // __UTILS_H__

#ifndef __FDAPDE_H__
#define __FDAPDE_H__

// Insert principal libraries
#define R_NO_REMAP
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include <stdint.h>
#include <iostream>

#include <cstdlib>
//#include <iomanip>
#include <limits>
#include <vector>
#include <array>
#include <stack>
#include <set>
#include "Global_Utilities/Include/Make_Unique.h"

// For debugging purposes
//#include <Eigen/StdVector>
//#include "Eigen/Eigen/Sparse"
//#include "Eigen/Eigen/Dense"
//#define  EIGEN_MPL2_ONLY

//Take the code from the linked RcppEigen
#include <Eigen/StdVector>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#define  EIGEN_MPL2_ONLY

typedef double Real;
typedef int UInt;


typedef Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> MatrixXr;
typedef Eigen::Matrix<UInt,Eigen::Dynamic,Eigen::Dynamic> MatrixXi;
typedef Eigen::Matrix<Real,Eigen::Dynamic,1> VectorXr;
typedef Eigen::Matrix<UInt,Eigen::Dynamic,1> VectorXi;
typedef Eigen::Matrix<VectorXr,Eigen::Dynamic,Eigen::Dynamic> MatrixXv;
typedef Eigen::SparseMatrix<Real> SpMat;
typedef Eigen::SparseVector<Real> SpVec;
typedef Eigen::Triplet<Real> coeff;

#include "Global_Utilities/Include/RObjects.h"

template <bool ... b>
struct multi_bool_type
{};

typedef multi_bool_type<true> t_type;
typedef multi_bool_type<false> f_type;
typedef multi_bool_type<true, true> tt_type;
typedef multi_bool_type<false, true> ft_type;
typedef multi_bool_type<true, false> tf_type;
typedef multi_bool_type<false, false> ff_type;


#endif /* FDAPDE_H_ */

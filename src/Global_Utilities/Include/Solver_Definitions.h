#ifndef __SOLVER_DEFINITIONS_H__
#define __SOLVER_DEFINITIONS_H__

//Take the code from the linked RcppEigen
#include "../../FdaPDE.h"

//#include <Eigen/SuperLUSupport>
//#include <Eigen/UmfPackSupport>

//! Some linear solvers definitions that may be useful for the future

typedef Eigen::SparseLU<SpMat> Sparse_LU;
typedef Eigen::SimplicialLDLT<SpMat> Sparse_Cholesky;
typedef Eigen::ConjugateGradient<SpMat> Sparse_ConjGrad;
typedef Eigen::BiCGSTAB<SpMat> Sparse_BiCGSTAB;
typedef Eigen::BiCGSTAB<SpMat,Eigen::IncompleteLUT<Real>> Sparse_BiCGSTAB_ILUT;
//typedef Eigen::UmfPackLU<SpMat> Sparse_UmfLU;
//typedef Eigen::SuperLU<SpMat> Sparse_SuperLU;

#endif

#ifndef __SRPDE_H__
#define __SRPDE_H__

#include <memory>
#include <Eigen/SparseQR>
#include "../core/utils/Symbols.h"
#include "../core/FEM/PDE.h"
using fdaPDE::core::FEM::PDE;
#include "../core/NLA/SparseBlockMatrix.h"
using fdaPDE::core::NLA::SparseBlockMatrix;
#include "../core/NLA/SMW.h"
using fdaPDE::core::NLA::SMW;
#include "Internals.h"
using fdaPDE::regression::internal::H;
using fdaPDE::regression::internal::Q;
using fdaPDE::regression::internal::psi;

template <unsigned int M, unsigned int N, unsigned int R, typename E>
class SRPDE{
private:
  // PDE of the problem
  const PDE<M,N,R,E>& pde_;
  double lambda_;
  
public:
  // constructor
  SRPDE() = default;
  SRPDE(const PDE<M,N,R,E>& pde, double lambda) : pde_(pde), lambda_(lambda) {};

  // finds a solution to the regression problem
  DVector<double> smooth(const DVector<double>& data);
  DVector<double> smooth(const DVector<double>& data, const DMatrix<double>& covariates);

};

template <unsigned int M, unsigned int N, unsigned int R, typename E>
SRPDE(const PDE<M,N,R,E>& pde_, double lambda_) -> SRPDE<M,N,R,E>;

template <unsigned int M, unsigned int N, unsigned int R, typename E>
DVector<double> SRPDE<M, N, R, E>::smooth(const DVector<double>& observations) {
  // assemble system matrix and solve linear system
  SpMatrix<double> Psi = *(psi(pde_));
  
  // for an SR-PDE model we have to solve A*x = b with
  //
  //         | -\Psi^T*\Psi  \lambda*R1^T |      | -Psi^T*z  |
  //     A = |                            |  b = |           |
  //         | \lambda*R1    \lambda*R0   |      | \lambda*u |
  // 
  // matrices R1 and R0 are obtained from the FEM discretization of the differential operator

  SparseBlockMatrix<double,2,2>
    A(-Psi.transpose()*Psi, lambda_ * pde_.R1().transpose(),
      lambda_ * pde_.R1(),  lambda_ * pde_.R0()            );
  
  DVector<double> b;
  b.resize(A.rows());
  b << -Psi.transpose()*observations,
       lambda_ * pde_.force();

  // define system solver. Matrix A is sparse!
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver{};
  solver.analyzePattern(A.derived());
  solver.factorize(A.derived());
  // solve linear system A*x = b
  DVector<double> sol = solver.solve(b);
  return (Psi*sol).head(A.rows()/2);
}

// smoothing in case covariates are supplied to the model. In this case the efficient solution of the linear system requires the application
// of the SMW decomposition. See NLA module for more details.
template <unsigned int M, unsigned int N, unsigned int R, typename E>
DVector<double> SRPDE<M, N, R, E>::smooth(const DVector<double>& observations, const DMatrix<double>& W) {
  // for an SR-PDE problem with covariates we need to solve M*x = b where
  // 
  //         | -\Psi^T*Q*\Psi  \lambda*R1^T |      | -Psi^T*Q*z |
  //     M = |                              |  b = |            |
  //         | \lambda*R1    \lambda*R0     |      | \lambda*u  |
  //
  // matrices R1 and R0 are obtained from the FEM discretization of the differential operator.
  // Factorize matrix M as the sum of the two following matrices
  //
  //         | -\Psi^T*\Psi  \lambda*R1^T |      | Psi^T*H*Psi    O_N |
  //     A = |                            |  B = |                    |
  //         | \lambda*R1    \lambda*R0   |      |     O_N        O_N |
  //
  // A is the system matrix of a nonparametric SR-PDE problem while B has just one dense block (Psi^T*H*Psi). We can show that
  // M = A + U*C*V where
  //
  //      | Psi^T*W |
  // U  = |         |  C = (W^T*W)^{-1}  V = | W^T*Psi  O_N |
  //      |   O_N   |
  //
  // then an efficient solution to the system M*x = b is obtain using SMW decomposition
  SpMatrix<double> Psi = *(psi(pde_));

  SparseBlockMatrix<double,2,2>
    A(-Psi.transpose()*Psi, lambda_ * pde_.R1().transpose(),
      lambda_ * pde_.R1(),  lambda_ * pde_.R0()            );
  // Define SMW solver
  SMW<> solver{};
  solver.compute(A.derived());

  // compute Q matrix (Q = I - H = I - W*(W*W^T)^{-}1*W^T)
  DMatrix<double> Q_ = *Q(*H(W));
  // compute system vector
  DVector<double> b;
  b.resize(A.rows());
  b << -Psi.transpose()*Q_*observations,
       lambda_ * pde_.force();
  
  Eigen::Index q = W.cols(); // number of covariates
  DMatrix<double> U = DMatrix<double>::Zero(A.rows(), q);
  U.block(0,0, A.rows()/2, q) = Psi.transpose()*W;
  
  // SMW implementation requires directly the inversion of C, compute here directly as W^T*W
  DMatrix<double> invC = W.transpose()*W;
  
  DMatrix<double> V = DMatrix<double>::Zero(q, A.rows());
  V.block(0,0, q, A.rows()/2) = W.transpose()*Psi;
  
  DVector<double> sol = solver.solve(U, invC, V, b); 
  return (Psi*sol).head(A.rows()/2);  
}

#endif // __SRPDE_H__

#ifndef __SRPDE_H__
#define __SRPDE_H__

#include <memory>
#include <Eigen/SparseQR>
#include "../core/utils/Symbols.h"
#include "../core/FEM/PDE.h"
using fdaPDE::core::FEM::PDE;
#include "../core/NLA/SparseBlockMatrix.h"
using fdaPDE::core::NLA::SparseBlockMatrix;
#include "Internals.h"
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
};

template <unsigned int M, unsigned int N, unsigned int R, typename E>
SRPDE(const PDE<M,N,R,E>& pde_, double lambda_) -> SRPDE<M,N,R,E>;

template <unsigned int M, unsigned int N, unsigned int R, typename E>
DVector<double> SRPDE<M, N, R, E>::smooth(const DVector<double>& data) {
  // assemble system matrix and solve linear system
  SpMatrix<double> Psi = *(psi(pde_));
  
  // for an SR-PDE model we have to solve A*x = b with
  //         | -\Psi^T*\Psi  \lambda*R1^T |      | -Psi^T*z  |
  //     A = |                            |  b = |           |
  //         | \lambda*R1    \lambda*R0   |      | \lambda*u |
  // matrices R1 and R0 are obtained from the FEM discretization of the differential operator

  SparseBlockMatrix<double,2,2>
    A(-Psi.transpose()*Psi, lambda_ * pde_.R1().transpose(),
      lambda_ * pde_.R1(),  lambda_ * pde_.R0()            );
  
  DVector<double> b;
  b.resize(A.rows());
  b << -Psi.transpose()*data,
       lambda_ * pde_.force();

  // define system solver. Matrix A is sparse!
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver{};
  solver.analyzePattern(A.derived());
  solver.factorize(A.derived());
  // solve linear system A*x = b
  DVector<double> sol = solver.solve(b);
  return Psi*sol.head(A.rows()/2);
}

#endif // __SRPDE_H__

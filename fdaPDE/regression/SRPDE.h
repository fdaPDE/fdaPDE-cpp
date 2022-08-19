#ifndef __SRPDE_H__
#define __SRPDE_H__

#include "../core/utils/Symbols.h"
#include "../core/FEM/PDE.h"
#include <Eigen/src/SparseCore/SparseMatrix.h>
using fdaPDE::core::FEM::PDE;
#include "../core/NLA/SparseBlockMatrix.h"
using fdaPDE::core::NLA::SparseBlockMatrix;

template <unsigned int M, unsigned int N, unsigned int R, typename E>
class SRPDE{
private:
  // PDE of the problem
  const PDE<M,N,R,E>& pde;
  double lambda; // the smoothing parameter
  
  SpMatrix<double> Psi;
  void computePsi();
public:
  // constructor
  SRPDE() = default;
  SRPDE(const PDE<M,N,R,E>& pde_, double lambda_) : pde(pde_), lambda(lambda_) {
    computePsi();
  };

  // finds a solution to the regression problem
  void smooth(const DVector<double>& data);
};

template <unsigned int M, unsigned int N, unsigned int R, typename E>
SRPDE(const PDE<M,N,R,E>& pde_, double lambda_) -> SRPDE<M,N,R,E>;

template <unsigned int M, unsigned int N, unsigned int R, typename E>
void SRPDE<M, N, R, E>::smooth(const DVector<double>& data) {
  // assemble system matrix and solve linear system
  // for an SR-PDE model we have to solve A*x = b with
  //         | -\Psi^T*\Psi  \lambda*R1^T |      | -Psi^T*z  |
  //     A = |                            |  b = |           |
  //         | \lambda*R1    \lambda*R0   |      | \lambda*u |
  SparseBlockMatrix<double,2,2> A(-Psi.transpose()*Psi, lambda * pde.R1().transpose(),
				  lambda * pde.R1(),    lambda * pde.R0()            );
  
  DVector<double> b;
  b.resize(A.rows());
  b << -Psi.transpose()*data, lambda * pde.force();

  // define system solver. Matrix A is sparse!
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver{};
  solver.analyzePattern(A.derived());
  solver.factorize(A.derived());
  DVector<double> sol = solver.solve(b);
  std::cout << Psi*sol.head(A.rows()/2) << std::endl;
}

// compute Psi matrix assuming locations equal to mesh's nodes (there is 1 only at mesh nodes and 0 elsewhere due to support of lagrangian basis)
// in general it is not diagonal!
template <unsigned int M, unsigned int N, unsigned int R, typename E>
void SRPDE<M, N, R, E>::computePsi() {
  // preallocate space for Psi matrix
  unsigned int locations = pde.domain().nodes();
  unsigned int nbasis = pde.domain().nodes();
  Psi.resize(locations, nbasis);
  
  // fill psi matrix
  std::list<Eigen::Triplet<double>> tripletList;  
  for(std::size_t i = 0; i < locations; ++i){
    tripletList.push_back(Eigen::Triplet<double>(i, i, 1));
  }
  
  Psi.setFromTriplets(tripletList.begin(), tripletList.end());
  Psi.makeCompressed();
  return;
}

#endif // __SRPDE_H__

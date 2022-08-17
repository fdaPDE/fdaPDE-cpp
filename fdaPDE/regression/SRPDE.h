#ifndef __SRPDE_H__
#define __SRPDE_H__

#include "../core/utils/Symbols.h"
#include "../core/FEM/PDE.h"
using fdaPDE::core::FEM::PDE;

template <unsigned int M, unsigned int N, unsigned int R, typename E>
class SRPDE{
private:
  // PDE of the problem
  const PDE<M,N,R,E>& pde;
  double lambda; // the smoothing parameter
  
  SpMatrix<double> Psi{};
  void computePsi();
  void assemble();
public:
  // constructor
  SRPDE() = default;
  SRPDE(const PDE<M,N,R,E>& pde_, double lambda_) : pde(pde_), lambda(lambda_) {};
};

template <unsigned int M, unsigned int N, unsigned int R, typename E>
void SRPDE<M, N, R, E>::assemble() {
  // assemble system matrix and solve linear system

  SpMatrix<double> A = -Psi.transpose()*Psi.transpose();
  SpMatrix<double> B = lambda * pde.R1().transpose();
  SpMatrix<double> C = lambda * pde.R0();

  SpMatrix<double> m;
  unsigned int nbasis = pde.domain().nodes();
  m.resize(nbasis, nbasis);
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
  for(const auto& e : pde.domain()){
    tripletList.push_back(Eigen::Triplet<double>(e->ID(), e->ID(), 1));
  }
  Psi.setFromTriplets(tripletList.begin(), tripletList.end());
  Psi.makeCompressed();
  return;
}

#endif // __SRPDE_H__

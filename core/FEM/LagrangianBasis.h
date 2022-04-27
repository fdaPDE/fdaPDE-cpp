#ifndef __LAGRANGIAN_BASIS_H__
#define __LAGRANGIAN_BASIS_H__

#include "../utils/CompileTime.h"
#include "../utils/Symbols.h"
#include "MultivariatePolynomial.h"
#include <array>
#include <Eigen/QR>

// A class representing a Lagrangian Basis defined over a given set of nodes.
// It uses the Vandermonde matrix to compute coefficients of lagrange polynomials
template <unsigned int N, unsigned int R> class LagrangianBasis {
private:
  // nodes of the lagrangian basis
  std::array<std::array<double, N>, ct_binomial_coefficient(R + N, R)> nodes_;
  // a Lagrangian basis is just a collection of properly defined polynomials
  std::array<MultivariatePolynomial<N,R>, ct_binomial_coefficient(N+R,R)> basis_;
public:
  // constructor
 LagrangianBasis(const std::array<std::array<double, N>, ct_binomial_coefficient(N+R, R)>& nodes) : nodes_(nodes) {

    // build vandermonde matrix
    constexpr unsigned int M = ct_binomial_coefficient(N+R,R);
    constexpr std::array<std::array<unsigned, N>, M> expTable_ = MultivariatePolynomial<N,R>::expTable_;

    // Vandermonde matrix construction
    SMatrix<M> V = Eigen::Matrix<double, M, M>::Ones();
    for(size_t i = 0; i < M; ++i){
      for(size_t j = 1; j < M; ++j){
	V(i,j) = MonomialProduct<N-1, std::array<double, N>, std::array<unsigned, N>>::unfold(nodes_[i], expTable_[j]);
      }
    }    
    // solve system V*a = b with b vector having 1 at position i and 0 everywhere else.
    // Its solution gives the vector of coefficients of the i-th Lagrange polynomial
    Eigen::ColPivHouseholderQR<SMatrix<M>> QRdecomposition(V);
    for(size_t i = 0; i < M; ++i){
      // build rhs vector
      SVector<M> b = Eigen::Matrix<double, M, 1>::Zero(); b[i] = 1;
      SVector<M> coeff = QRdecomposition.solve(b); // solve system

      // cast to array
      std::array<double, M> coeff_array;
      for(size_t j = 0; j < M; ++j) coeff_array[j] = coeff[j];
      
      basis_[i] = MultivariatePolynomial<N, R>(coeff_array); // store basis
    }    
  };

  // get basis element
  MultivariatePolynomial<N, R> getBasisElement(unsigned int n) { return basis_[n]; };
};

#endif // __LAGRANGIAN_BASIS_H__

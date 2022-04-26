#ifndef __LAGRANGIAN_BASIS_H__
#define __LAGRANGIAN_BASIS_H__

#include "../utils/CompileTime.h"
#include "../utils/Symbols.h"
#include <Eigen/src/Core/Matrix.h>
#include <array>
#include <cstddef>
#include <Eigen/QR>

// Content of this file is intended to provide an operative definition of a space of polynomials
// on an N-dimensional simplex of global degree less than or equal to R

/* A polynomial is given by the sum of ct_binomial_coefficient(R+N, R) monomials x1^i1*x2^i2*...*xN^iN. Fixed N and R we can
   precompute the exponents (i1, i2, ..., iN) at compile time.
   For example for a quadratic polinomial p(x) over a 3D space we would obtain a table (named exp table) as follow (R = 2):

   i1  i2  i3 | R
   0   0   0  | 0
   0   0   1  | 1
   0   0   2  | 2
   0   1   0  | 1
   0   1   1  | 2
   0   2   0  | 2
   1   0   0  | 1
   1   0   1  | 2
   1   1   0  | 2
   2   0   0  | 2

   Each row of the above table represent a single monomial of the polynomial p(x). ct_poly_exp() computes the above table for any
   choice of N and R at compile time.
 */

template <unsigned int N, unsigned int R>
constexpr std::array<std::array<unsigned, N>, ct_binomial_coefficient(R+N, R)> ct_poly_exp(){
  const int monomials = ct_binomial_coefficient(R+N, R); // number of monomials in polynomial p(x)
  // initialize empty array
  std::array<std::array<unsigned, N>, monomials> coeff{};

  // auxiliary array. At each iteration this will be a row of the exp table
  std::array<unsigned, N> tmp{};
  int j = 0;
  while(j < monomials){
    // compute exp vector for monomial j
    int i = 0;
    bool found = false;
    while(i < N && !found){
      if(tmp[i] <= R && ct_array_sum<N>(tmp) <= R){
	found = true;
	// copy exp vector for monomial j into result
	for(int k = 0; k < N; ++k) coeff[j] = tmp;
	tmp[0]++; // increment for next iteration
	j++;      // next monomial
      }else{
	// propagate carry to next element
	tmp[i] = 0;
	tmp[++i]++;
      }
    }
  }
  return coeff;
}

// class representing an N-dimensional multivariate polynomial of degree R
template <unsigned int N, unsigned int R>
class MultivariatePolynomial{
private:
  // vector of coefficients
  std::array<double, ct_binomial_coefficient(R+N, R)> coeffVector_;  
public:
  // compute this at compile time once, let publicly accessible
  static const constexpr std::array<std::array<unsigned, N>, ct_binomial_coefficient(R+N, R)> expTable_ = ct_poly_exp<N,R>();
  
  // constructor
  MultivariatePolynomial() = default;
  MultivariatePolynomial(const std::array<double, ct_binomial_coefficient(R+N, R)>& coeffVector) : coeffVector_(coeffVector) {};
  
  double     eval(const SVector<N>& point);    // evaluate polynomial at point
  SVector<N> grad(const SVector<N>& point);    // evaluate gradient at point

  // getter
  std::array<double, ct_binomial_coefficient(R+N, R)> getCoeff() const { return coeffVector_; }
};


// template based unfold of monomial product x1^i1*x2^i2*...*xN^iN. This allow
// to evaluate monomials at point P without explicitly looping over its individual components
// P an std::array<double, N> representing the point where we have to evaluate the monomial, V an array of exp coefficients [i1 i2 ... iN].
template<unsigned int N, typename P, typename V> struct MonomialProduct{
  static double unfold(const P& p, const V& v){
    return v[N] == 0 ? MonomialProduct<N-1, P, V>::unfold(p, v) : std::pow(p[N], v[N]) * MonomialProduct<N-1, P, V>::unfold(p, v);
  }
};

template <typename P, typename V> struct MonomialProduct<0, P, V> { // base case
  static double unfold(const P& p, const V& v){
    return v[0] == 0 ? 1 : std::pow(p[0], v[0]);
  }
};

template <unsigned int N, unsigned int R>
double MultivariatePolynomial<N, R>::eval(const SVector<N> &point) {
  double value = 0;
  // cycle over all monomials
  for(size_t i = 0; i < ct_binomial_coefficient(R+N, R); ++i){
    // this expands at compile time to the product of each factor in the monomial
    value += coeffVector_[i]*MonomialProduct<N-1, SVector<N>, std::array<unsigned, N>>::unfold(point, expTable_[i]);
  }  
  return value;
}

template <unsigned int N, unsigned int R>
SVector<N> MultivariatePolynomial<N, R>::grad(const SVector<N> &point) {
  SVector<N> grad;
  // cycle over dimensions
  for(size_t i = 0; i < N; ++i){
    double value = 0;
    
    // for a given monomial x1^i1*...*xN^iN its partial derivative with respect to dimension K is equal to
    // d/dxK x1^i1*...*xN^iN = iK*x1^i1*...*xK^(iK-1)*...*xN^iN.

    // cycle over monomials
    for(size_t m = 0; m < ct_binomial_coefficient(R+N,R); ++m){
      if(expTable_[m][i] != 0){ // skip powers of zero, their derivative is zero
	double tmp = 1;

	// actual partial derivative computation along dimension i
	for(size_t j = 0; j < N; ++j){
	  tmp *= std::pow(point[j], j != i ? expTable_[m][j] : expTable_[m][j] - 1);
	}
	value += coeffVector_[m]*expTable_[m][i]*tmp;
      }
    }
    // store value of partial derivative in gradient vector
    grad[i] = value;
  }
    
  return grad;
}

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

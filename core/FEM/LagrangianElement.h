#ifndef __LAGRANGIAN_ELEMENT__
#define __LAGRANGIAN_ELEMENT__

#include "../MESH/Element.h"
#include "../utils/CompileTime.h"
#include "../utils/Symbols.h"
#include <Eigen/src/Core/Matrix.h>
#include <array>
#include <cstddef>
#include <Eigen/QR>

using fdaPDE::core::MESH::Element;

// space of polynomials defined on an N-dimensional simplex of global degree less than or equal to R

/* A polynomial is fully described by its coefficients. The equation of a polynomial on a simplex is given as
   p(x) = \sum_{0 \leq i1, ..., iN; i1 + ... + iN \leq R} \alpha_{i1, ..., iN} x1^i1 * x2^i2 * .... * xd^id

      For example for a second order polynomial on a 2D simplex we have
      p(x) = \alpha_{0,0} + \alpha_{1,0} x1 + \alpha_{0,1} + \alpha_{1,1} x1*x2 + \alpha_{2,0} x1^2 + \alpha_{0,2} x2^2

   We store a polynomial by just storing the vector of its coefficients [\alpha_{0,...0}, ..., \alpha{d,...,d}].
   Given this representation is straightforward to evaluate the polynomial at a given point or compute its gradient vector
 */

/* A polynomial is given by the sum of ct_binomial_coefficient(R+N, R) monomials x1^i1*x2^i2*...*xN^iN. Fixed N and R we can
   precompute the exponents (i1, i2, ..., iN) at compile time. For example for a quadratic polinomial p(x) over a 3D space we would
   obtain a table (named exp table) as follow (R = 2):

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
constexpr std::array<unsigned, ct_binomial_coefficient(R+N, R)*N> ct_poly_exp(){
  const int monomials = ct_binomial_coefficient(R+N, R); // number of monomials in polynomial p(x)
  // initialize empty array
  std::array<unsigned, monomials*N> coeff{};

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
	for(int k = 0; k < N; ++k) coeff[j*N + k] = tmp[k];
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

template <unsigned int N, unsigned int R>
class MultivariatePolynomial{
private:
  // vector of coefficients
  std::array<double, ct_binomial_coefficient(R+N, R)> coeffVector_;  
public:
  // compute this at compile time once, let publicly accessible
  static const constexpr std::array<unsigned, ct_binomial_coefficient(R+N, R)*N> expTable_ = ct_poly_exp<N,R>();
  
  // constructor
  MultivariatePolynomial() = default;
  MultivariatePolynomial(const std::array<double, ct_binomial_coefficient(R+N, R)>& coeffVector) : coeffVector_(coeffVector) {};
  
  double     eval(const SVector<N>& point);    // evaluate polynomial at point
  SVector<N> grad(const SVector<N>& point);    // evaluate gradient at point

  // getter
  std::array<double, ct_binomial_coefficient(R+N, R)> getCoeff() const { return coeffVector_; }
};

template <unsigned int N, unsigned int R>
double MultivariatePolynomial<N, R>::eval(const SVector<N> &point) {
  double value = 0;
  // cycle over all monomials
  for(size_t i = 0; i < ct_binomial_coefficient(R+N, R); ++i){
    double tmp = 1;
    for(size_t j = 0; j < N; ++j){
      if(expTable_[i*N + j] != 0) // skip powers of zero
	tmp *= std::pow(point[j], expTable_[i*N + j]);
    }
    value += coeffVector_[i]*tmp;
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
      if(expTable_[i + m*N] != 0){ // skip powers of zero, their derivative is zero
	double tmp = 1;

	// actual partial derivative computation along dimension i
	for(size_t j = 0; j < N; ++j){
	  tmp *= std::pow(point[j], j != i ? expTable_[m*N + j] : expTable_[m*N + j] - 1);
	}
	value += coeffVector_[m]*expTable_[i + m*N]*tmp;
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
  std::array<SVector<N>, ct_binomial_coefficient(R + N, R)> nodes_;
  // a Lagrangian basis is just a collection of properly defined polynomials
  std::array<MultivariatePolynomial<N,R>, ct_binomial_coefficient(N+R,R)> basis_;
public:
  // constructor
  LagrangianBasis(const std::array<SVector<N>, ct_binomial_coefficient(R + N, R)>& nodes) : nodes_(nodes) {

    // build vandermonde matrix
    constexpr unsigned int M = ct_binomial_coefficient(N+R,R);
    constexpr std::array<unsigned, M*N> expTable_ = MultivariatePolynomial<N,R>::expTable_;

    // Vandermonde matrix construction
    SMatrix<M> vander = Eigen::Matrix<double, M, M>::Ones();
    for(size_t i = 0; i < M; ++i){
      for(size_t j = 1; j < M; ++j){
	for(size_t z = 0; z < N; ++z){
	  vander(i,j) *= std::pow(nodes_[i][z], expTable_[j*N + z]);
	}
      }
    }
    
    // solve system V*a = b with V the vandermonde matrix and b the vector having 1 at position i
    // and 0 everywhere else. The solution gives the vector of coefficients of the i-th lagrange polynomial
    Eigen::ColPivHouseholderQR<SMatrix<M>> QRdecomposition(vander);

    for(size_t i = 0; i < M; ++i){
      // build rhs vector
      SVector<M> b = Eigen::Matrix<double, M, 1>::Zero();
      b[i] = 1;

      // system solution
      SVector<M> coeff = QRdecomposition.solve(b);

      // cast to array
      std::array<double, M> coeff_array;
      for(size_t j = 0; j < M; ++j) coeff_array[j] = coeff[j];
      
      // store basis
      basis_[i] = MultivariatePolynomial<N, R>(coeff_array);
    }    
  };

  // get basis element
  MultivariatePolynomial<N, R> getBasisElement(unsigned int n) { return basis_[n]; };
};

// a class representing a Lagrangian finite element
template <unsigned int N, unsigned int M>
class LagrangianElement{

  Element<N,M> domain;
  //LagrangianBasis<N,ORDER> basis;
};


#endif // __LAGRANGIAN_ELEMENT__

#ifndef __LAGRANGIAN_ELEMENT__
#define __LAGRANGIAN_ELEMENT__

#include "../MESH/Element.h"
#include "../utils/CompileTime.h"
#include "../utils/Symbols.h"
#include <array>
#include <cstddef>

using fdaPDE::core::MESH::Element;

// a class representing a Lagrangian finite element
template <unsigned int N, unsigned int M>
class LagrangianElement{

  Element<N,M> domain;
    
};

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
	++i;
	tmp[i]++;
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

  // compute this at compile time once
  static constexpr std::array<unsigned, ct_binomial_coefficient(R+N, R)*N> expTable_ = ct_poly_exp<N,R>();
  
public:
  // constructor
  MultivariatePolynomial(const std::array<double, ct_binomial_coefficient(R+N, R)>& coeffVector) : coeffVector_(coeffVector) {};
  
  double eval(const SVector<N>& point);        // evaluate polynomial at point
  SVector<N> grad(const SVector<N>& point);    // evaluate gradient at point
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


// the reference element coincides with the unit N-dimensional simplex. Over the reference element we define
// a space of polynomials of order R with a particular choice of coefficients which constitute the basis function
class ReferenceElement {

  
};

// given a physical element move to reference for doing computations via an
// affine transformation

// computations involves mainly integrals over elements of kind nabla_phi_i *
// nabla_phi_j or phi_i*phi_j
// here we can implement different quadrature rules for the gaussian case

// hence we can compute stiff and mass matrices for the discretization of an elliptic operator
// first integrate on single elements then exploit linearity of the integral to compute the overall integral

#endif // __LAGRANGIAN_ELEMENT__

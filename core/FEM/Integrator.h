#ifndef __INTEGRATOR_H__
#define __INTEGRATOR_H__

#include "../MESH/Element.h"
#include "../utils/Symbols.h"
#include "../utils/CompileTime.h"
#include "../utils/fields/VectorField.h"
#include "IntegratorTables.h"

using fdaPDE::core::InnerProduct;
using fdaPDE::core::VectorField;
using fdaPDE::core::MESH::Element;

// modify the callable F if this is an InnerProduct by premultiplying each operand by the transpose of the inverse
// of the Jacobian matrix of the transformation which goes from a generic element to the reference one
template <unsigned int M, unsigned int N, typename F>
typename std::enable_if<std::is_same<F, InnerProduct<N>>::value, F>::type
adjusted(const Element<M, N>& e, F& f){
  Eigen::Matrix<double, N, M> invJ = e.getInvBaryMatrix().transpose();
  // matrix-VectorField multiplication (supported by custom arithmetic)
  f.lhs_ = invJ * (f.lhs_);  
  f.rhs_ = invJ * (f.rhs_);
  return f;
}

template <unsigned int M, unsigned int N, typename F>
typename std::enable_if<!std::is_same<F, InnerProduct<N>>::value, F>::type
adjusted(const Element<M, N>& e, F& f){ return f; } // do nothing

// An integrator class. Just integrates a given integrable function (aka
// ScalarField) on a mesh element. N space dimension, M number of nodes in the quadrature rule
template <int N, int M>
class Integrator {
 private:
  // reference to node and weights of the quadrature rule
  IntegratorTable<N, M> integrationTable_;
  
 public:
  // constructor
  Integrator() : integrationTable_(IntegratorTable<N, M>()) {};

  // F is any callable representing the function to integrate. L is the order of the mesh
  template <unsigned int ORDER, typename F>
  double integrate(const Element<ORDER, N>& e, F& f) const;

  // returns the integration table data structure
  const IntegratorTable<N,M>& getTable() const { return integrationTable_; }
};

// integrate a callable F given an integrator table T and a mesh element e.
// Given a function f() its integral \int_K f(x)dx can be approximated using a quadrature rule. The general
// scheme of a qudrature formula for the approximation of integral \int_K f(x)dx is given by a finite sum
// \sum_{i=1}^N [f(x_i) * w_i] where x_i and w_i are properly choosen quadrature nodes and weights.
template <int N, int M>
template <unsigned int ORDER, typename F> double Integrator<N, M>::integrate(const Element<ORDER, N>& e, F& f) const {
  constexpr unsigned int I = integrationTable_.num_nodes;
  F adjusted_f = adjusted(e, f);

  double result = 0;
  // execute quadrature rule
  for(size_t iq = 0; iq < integrationTable_.num_nodes; ++iq){
    result += adjusted_f(SVector<2>(integrationTable_.nodes[iq].data()))*integrationTable_.weights[iq];
  }
  result /= ct_factorial(N); // divide by correction factor

  return result;
}

#endif // __INTEGRATOR_H__

#ifndef __INTEGRATOR_H__
#define __INTEGRATOR_H__

#include "../MESH/Element.h"
#include "../utils/Symbols.h"
#include "../utils/CompileTime.h"
#include "../utils/fields/VectorField.h"
#include "IntegratorTables.h"
#include "../LagrangianBasis.h"
#include "../FunctionalBasis.h"

using fdaPDE::core::InnerProduct;
using fdaPDE::core::VectorField;
using fdaPDE::core::MESH::Element;

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

  template <unsigned int ORDER, typename F>
  typename std::enable_if< std::is_same<F, InnerProduct<N>>::value, double>::type integral(const Element<ORDER, N>& e, F& f) const;

  template <unsigned int ORDER, typename F>
  typename std::enable_if<!std::is_same<F, InnerProduct<N>>::value, double>::type integral(const Element<ORDER, N>& e, F& f) const;
  
};

// modify the callable F if this is an InnerProduct by premultiplying each operand by the transpose of the inverse
// of the Jacobian matrix of the transformation which goes from a generic element to the reference one
template <int N, int M>
template <unsigned int ORDER, typename F>
typename std::enable_if<std::is_same<F, InnerProduct<N>>::value, double>::type
Integrator<N, M>::integral(const Element<ORDER, N>& e, F& f) const {
  // express gradient of f in terms of gradients of basis functions over reference element. This entails to compute
  // (J^{-1})^T * \Nabla phi_i. 

  Eigen::Matrix<double, N, ORDER> invJ = e.getInvBaryMatrix().transpose();
  // Given \Nabla phi_i premultiply it by (J^{-1})^T = invJ.
  // NOTE: f.lhs_ and f.rhs_ are already basis over the reference element
  f.lhs_ = invJ * f.lhs_; 
  f.rhs_ = invJ * f.rhs_; 

  // compute integral
  double value = 0;
  // execute quadrature rule for \Nabla phi_i dot \Nabla phi_j
  for(size_t iq = 0; iq < integrationTable_.num_nodes; ++iq){
    value += f(SVector<N>(integrationTable_.nodes[iq].data()))*integrationTable_.weights[iq];
  }
  
  return value * std::abs(e.getBaryMatrix().determinant())/ct_factorial(N);
}

template <int N, int M>
template <unsigned int ORDER, typename F>
typename std::enable_if<!std::is_same<F, InnerProduct<N>>::value, double>::type
Integrator<N, M>::integral(const Element<ORDER, N>& e, F& f) const {
  double value = 0;
  // execute quadrature rule
  for(size_t iq = 0; iq < integrationTable_.num_nodes; ++iq){
    // map quadrature point to current element e
    SVector<N> p = e.getBaryMatrix()*SVector<N>(integrationTable_.nodes[iq].data()) + e.getCoords()[0];
    value += f(p)*integrationTable_.weights[iq];
  }
  // correct for measure of domain (element e)
  return value * (std::abs(e.getBaryMatrix().determinant())/ct_factorial(N));
}

// integrate a callable F given an integrator table T and a mesh element e.
// Given a function f() its integral \int_K f(x)dx can be approximated using a quadrature rule. The general
// scheme of a qudrature formula for the approximation of integral \int_K f(x)dx is given by a finite sum
// \sum_{i=1}^N [f(x_i) * w_i] where x_i and w_i are properly choosen quadrature nodes and weights.
template <int N, int M>
template <unsigned int ORDER, typename F>
double Integrator<N, M>::integrate(const Element<ORDER, N>& e, F& f) const {
  return integral(e, f);
}

#endif // __INTEGRATOR_H__

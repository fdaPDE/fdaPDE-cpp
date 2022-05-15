#ifndef __INTEGRATOR_H__
#define __INTEGRATOR_H__

#include "../MESH/Element.h"
#include "../utils/Symbols.h"
#include "../utils/CompileTime.h"
#include "../utils/fields/VectorField.h"
#include "IntegratorTables.h"
#include "../LagrangianBasis.h"
#include "../FunctionalBasis.h"
#include <cstddef>

using fdaPDE::core::VectorField;
using fdaPDE::core::MESH::Element;

// An integrator class. Just integrates a given integrable function (i.e. ScalarField or MultivariatePolynomials) on a mesh element.
// Given a function f() its integral \int_K f(x)dx can be approximated using a quadrature rule. The general
// scheme of a qudrature formula for the approximation of integral \int_K f(x)dx is given by a finite sum
// \sum_{i=1}^N [f(x_i) * w_i] where x_i and w_i are properly choosen quadrature nodes and weights.
//N space dimension, K number of nodes in the quadrature rule
template <unsigned int N, unsigned int K>
class Integrator {
 private:
  // reference to node and weights of the quadrature rule
  IntegratorTable<N, K> integrationTable_;
  
 public:
  // constructor
  Integrator() : integrationTable_(IntegratorTable<N, K>()) {};

  // F is any callable representing the function to integrate. L is the order of the mesh
  template <unsigned int ORDER, typename F>
  double integrate(const Element<ORDER, N>& e, F& f) const;

  // integrate a BilinearFormExpr to produce the (i,j)-th element of its discretization
  template <unsigned int ORDER, typename B, typename F>
  double integrate(const B& basis, const Element<ORDER, N>& e, int i , int j, const F& bilinearForm) const;
  
  // returns the integration table data structure
  const IntegratorTable<N, K>& getTable() const { return integrationTable_; }
};

// integrate a BilinearFormExpr over mesh element e using basis B. The pair (i,j) indicates the element position of the produced
// value in the stiff matrix discretizing the form itself. This method is used as part of the assembly loop in the computation of the
// stiff matrix
template <unsigned int N, unsigned int K>
template <unsigned int ORDER, typename B, typename F>
double Integrator<N, K>::integrate(const B& basis, const Element<ORDER, N>& e, int i , int j, const F& bilinearForm) const{
  // apply quadrature rule
  double value = 0;
  for(size_t iq = 0; iq < integrationTable_.num_nodes; ++iq){
    // for a BilinearFormExpr .integrate() is developed in the sum of the single integrals yelding to the discretization of the bilinear form
    SVector<N> p = SVector<N>(integrationTable_.nodes[iq].data());
    value += bilinearForm.integrate(basis, e, i, j, p) * integrationTable_.weights[iq];
  }
  
  return value * std::abs(e.getBaryMatrix().determinant())/ct_factorial(N);
}

// Integrate a callable F given an integrator table T and a mesh element e.
template <unsigned int N, unsigned int K>
template <unsigned int ORDER, typename F>
double Integrator<N, K>::integrate(const Element<ORDER, N> &e, F &f) const {
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


#endif // __INTEGRATOR_H__

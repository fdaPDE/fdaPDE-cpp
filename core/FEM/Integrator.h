#ifndef __INTEGRATOR_H__
#define __INTEGRATOR_H__

#include "../MESH/Element.h"
#include "../utils/Symbols.h"
#include "VectorField.h"
#include "IntegratorTables.h"

using fdaPDE::core::MESH::Element;

// compile time loop unfolding based on template recursion for quadrature rules.
template <unsigned int N,  // template recursion loop variable
	  unsigned int M,  // space dimension
	  typename T,      // type representing an integrator table (or a quadrature rule)
	  typename F>      // callable to integrate, must be evaluable at SVector<M>
struct IntegratorLoop {
  static constexpr double unfold(const T& t, const F& f){
    return t.weights[N]*f(SVector<M>(t.nodes[N].data()))*IntegratorLoop<N-1, M, T, F>::unfold(t,f);
  }
};

template <unsigned int M, typename T, typename F> struct IntegratorLoop<0, M, T, F> {
  static constexpr double unfold(const T& t, const F& f){
    return t.weights[0]*f(SVector<M>(t.nodes[0].data()));
  }
};

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
  const IntegratorTable<N, M>& integrationTable_;
  
 public:
  // constructor
  Integrator() : integrationTable_(IntegratorTable<N, M>()) {};

  // F is any callable representing the function to integrate. L is the order of the mesh
  template <unsigned int L, typename F>
  double integrate(const Element<L, N>& e, F& f) const;

  // returns the integration table data structure
  const IntegratorTable<N,M>& getTable() const { return integrationTable_; }
};

// integrate a callable F given an integrator table T and a mesh element e.
// Given a function f() its integral \int_K f(x)dx can be approximated using a quadrature rule. The general
// scheme of a qudrature formula for the approximation of integral \int_K f(x)dx is given by a finite sum
// \sum_{i=1}^N [f(x_i) * w_i] where x_i and w_i are properly choosen quadrature nodes and weights.
template <int N, int M>
template <unsigned int L, typename F> double Integrator<N, M>::integrate(const Element<L, N>& e, F& f) const {
  constexpr unsigned int I = integrationTable_.num_nodes;
  F adjusted_f = adjusted(e, f);
  double integralValue = IntegratorLoop<I - 1, N, const IntegratorTable<N, M>&, F>::unfold(integrationTable_, adjusted_f); // unfold loop at compile time

  // correct the obtained value by the jacobian of the transformation from e to the reference element
  return integralValue*std::abs(e.getBaryMatrix().determinant())/ct_factorial(N);
}

// NEDDO PASSARE SVECTOR ALLO SCALAR FIELD E NON ARRAY COME NELLE INTEGRATION TABLES!

#endif // __INTEGRATOR_H__

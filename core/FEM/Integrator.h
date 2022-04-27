#ifndef __INTEGRATOR_H__
#define __INTEGRATOR_H__

#include "IntegratorTables.h"

// compile time loop unfolding based on template recursion for quadrature rules.
template <unsigned int N, typename T, typename F> struct IntegratorLoop {
  static constexpr double unfold(const T& t, const F& f){
    return t.weights[N]*f(t.nodes[N])*IntegratorLoop<N-1, T, F>::unfold(t,f);
  }
};

template <typename T, typename F> struct IntegratorLoop<0, T, F> {
  static constexpr double unfold(const T& t, const F& f){
    return t.weights[0]*f(t.nodes[0]);
  }
};

// An integrator class. T is the IntegratorTable which defines the quadrature rule to apply
template <typename T>
class Integrator {
 private:
  // reference to node and weights of the quadrature rule
  const T& integrationTable_;
  
 public:
  // constructor
  Integrator(const T& integrationTable) : integrationTable_(integrationTable) {};

  // F is any callable representing the function to integrate
  template <typename F>
  double integrate(const F& f) const;
};

// integrate a callable F given an integrator table T.
// Given a function f() its integral \int_K f(x)dx can be approximated using a quadrature rule. The general
// scheme of a qudrature formula for the approximation of integral \int_K f(x)dx is given by a finite sum
// \sum_{i=1}^N [f(x_i) * w_i] where x_i and w_i are properly choosen quadrature nodes and weights.
template <typename T> template <typename F> double Integrator<T>::integrate(const F& f) const {
  constexpr unsigned int L = integrationTable_.num_nodes;
  return IntegratorLoop<L - 1, T, F> (integrationTable_, f);
}

#endif // __INTEGRATOR_H__

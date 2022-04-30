#ifndef __POLYNOMIAL_EXPRESSIONS_H__
#define __POLYNOMIAL_EXPRESSIONS_H__

#include "../utils/Symbols.h"

#include <functional>
#include <type_traits>

// MultivariatePolynomial class support expression template based arithmetic.
// The goal is to allow to write expressions of polynomials which are lazily evaluated only when a point evaluation is requested.

// macro to define an arithmetic operator between polynomials.
#define DEF_POLY_EXPR_OPERATOR(OPERATOR, FUNCTOR)			\
  template <typename E1, typename E2>					\
  PolyBinOp<E1, E2, FUNCTOR >						\
  OPERATOR(const PolyExpr<E1>& op1, const PolyExpr<E2>& op2) {		\
    return PolyBinOp<E1, E2, FUNCTOR >					\
      {op1.get(), op2.get(), FUNCTOR()};				\
  }									\
  									\
  template <typename E>							\
  PolyBinOp<E, PolyScalar, FUNCTOR >					\
  OPERATOR(const PolyExpr<E>& op1, double op2) {			\
  return PolyBinOp<E, PolyScalar, FUNCTOR >				\
      (op1.get(), PolyScalar(op2), FUNCTOR());				\
  }									\
  									\
  template <typename E>							\
  PolyBinOp<PolyScalar, E, FUNCTOR >					\
  OPERATOR(double op1, const PolyExpr<E>& op2) {			\
    return PolyBinOp<PolyScalar, E, FUNCTOR >				\
      {PolyScalar(op1), op2.get(), FUNCTOR()};				\
  }									\

// Base class for polynomial expressions
template <typename E> struct PolyExpr {
  // call operator() on the base type E
  template <int N>
  double operator()(const SVector<N>& p) const {
    return static_cast<const E&>(*this)(p);
  }

  // get underyling type composing the expression node
  const E& get() const { return static_cast<const E&>(*this); }
};

// an expression node representing a scalar value (double, int, ... even single valued variables)
class PolyScalar : public PolyExpr<PolyScalar> {
private:
  double value_;
public:
  PolyScalar(double value) : value_(value) { }
  
  // call operator
  template <int N>
  double operator()(const SVector<N>& p) const { return value_; };
};

// expression template based arithmetic
template <typename OP1, typename OP2, typename BinaryOperation>
class PolyBinOp : public PolyExpr<PolyBinOp<OP1, OP2, BinaryOperation>> {
private:
  typename std::remove_reference<OP1>::type op1_;   // first  operand
  typename std::remove_reference<OP2>::type op2_;   // second operand
  BinaryOperation f_;                               // operation to apply

public:
  // constructor
  PolyBinOp(const OP1& op1, const OP2& op2, BinaryOperation f) : op1_(op1), op2_(op2), f_(f) { };

  // call operator, performs the expression evaluation
  template <int N>
  double operator()(const SVector<N>& p) const{
    return f_(op1_(p), op2_(p));
  }
};

DEF_POLY_EXPR_OPERATOR(operator+, std::plus<>      )
DEF_POLY_EXPR_OPERATOR(operator-, std::minus<>     )
DEF_POLY_EXPR_OPERATOR(operator*, std::multiplies<>)

#endif // __POLYNOMIAL_EXPRESSIONS_H__

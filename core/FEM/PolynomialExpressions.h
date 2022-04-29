#ifndef __POLYNOMIAL_EXPRESSIONS_H__
#define __POLYNOMIAL_EXPRESSIONS_H__

#include "../utils/Symbols.h"

#include <functional>

// MultivariatePolynomial class support expression template based arithmetic.
// Base class for expression templates
template <typename E> struct PolyExpr {
  // call operator() on the base type E
  template <int N>
  double operator()(const SVector<N>& p) const {
    return static_cast<const E&>(*this)(p);
  }
};

// expression template based arithmetic
template <typename OP1, typename OP2, typename BinaryOperation>
class PolyBinOp : public PolyExpr<PolyBinOp<OP1, OP2, BinaryOperation>> {
private:
  const OP1& op1_;          // first  operand
  const OP2& op2_;          // second operand
  const BinaryOperation f_; // operation to apply between operands

public:
  // constructor
  PolyBinOp(const OP1& op1, const OP2& op2, const BinaryOperation& f) : op1_(op1), op2_(op2), f_(f) { };

  // call operator
  template <int N>
  double operator()(const SVector<N>& p) const{
    return f_(op1_(p), op2_(p));
  }
};

// sum of 2 expressions
template <typename E1, typename E2>
PolyBinOp<PolyExpr<E1>, PolyExpr<E2>, std::plus<> >
operator+(const PolyExpr<E1>& op1, const PolyExpr<E2>& op2) {
  return PolyBinOp<PolyExpr<E1>, PolyExpr<E2>, std::plus<> >(op1, op2, std::plus<>());
}

// difference of 2 expressions
template <typename E1, typename E2>
PolyBinOp<PolyExpr<E1>, PolyExpr<E2>, std::minus<> >
operator-(const PolyExpr<E1>& op1, const PolyExpr<E2>& op2) {
  return PolyBinOp<PolyExpr<E1>, PolyExpr<E2>, std::minus<> >(op1, op2, std::minus<>());
}

// product of 2 expressions
template <typename E1, typename E2>
PolyBinOp<PolyExpr<E1>, PolyExpr<E2>, std::multiplies<> >
operator*(const PolyExpr<E1>& op1, const PolyExpr<E2>& op2) {
  return PolyBinOp<PolyExpr<E1>, PolyExpr<E2>, std::multiplies<> >(op1, op2, std::multiplies<>());
}

#endif // __POLYNOMIAL_EXPRESSIONS_H__

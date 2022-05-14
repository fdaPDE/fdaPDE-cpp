#ifndef __VECTOR_FIELD_H__
#define __VECTOR_FIELD_H__

#include "../utils/Symbols.h"
#include "VectorFieldExpressions.h"
#include "ScalarField.h"

namespace fdaPDE{
namespace core{

  // forward declaration for InnerProduct functor
  template <int N> class InnerProduct;

  // template class representing a general vector field. Support expression template arithmetic
  template <int N>
  class VectorField : public VectExpr<N, VectorField<N>>{
  private:
    // each array element is a lambda which computes the i-th component of the vector
    std::array<std::function<double(SVector<N>)>, N> field_;
  public:
    // constructor
    VectorField() = default;
    VectorField(std::array<std::function<double(SVector<N>)>, N> field) : field_(field) {}

    // call operator
    SVector<N> operator()(const SVector<N>& point) const;
    // subscript operator
    const std::function<double(SVector<N>)>& operator[](size_t i) const;
    std::function<double(SVector<N>)>& operator[](size_t i);

    // inner product VectorField.dot(VectorField)
    InnerProduct<N> dot(const VectorField<N>& rhs) const;
    ScalarField<N> dot(const SVector<N>& op) const;
  };

  // A functor to represent an inner product.
  template <int N>
  struct InnerProduct{
    // operands of the inner product operation
    VectorField<N> lhs_;
    VectorField<N> rhs_;

    // constructor
    InnerProduct(const VectorField<N>& lhs, const VectorField<N>& rhs) : lhs_(lhs), rhs_(rhs) {}

    // evaluate the form on two different points (this is like evaluating one scalar field at one point,
    // the second one at another point and then perform the inner product between the two numerical results)
    double operator()(const SVector<N>& x, const SVector<N>& y) const;
    // consider the analytical expression of the inner product as a scalar field, evaluate it at point x
    double operator()(const SVector<N>& x) const;
  };

#include "VectorField.tpp"
}}

#endif // __VECTOR_FIELD_H__

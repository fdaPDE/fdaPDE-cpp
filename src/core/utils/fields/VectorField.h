#ifndef __VECTOR_FIELD_H__
#define __VECTOR_FIELD_H__

#include "../utils/Symbols.h"
#include "VectorFieldExpressions.h"
#include "ScalarField.h"
#include "DotProduct.h"

namespace fdaPDE{
namespace core{

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
    DotProduct<VectorField<N>, VectorField<N>> dot(const VectorField<N>& rhs) const;
    // Inner product VectorField.dot(SVector)
    ScalarField<N>  dot(const SVector<N>& op) const;
  };

  // template argument deduction guide
  template <int N> VectorField(std::array<std::function<double(SVector<N>)>, N>) -> VectorField<N>;
  
#include "VectorField.tpp"
}}

#endif // __VECTOR_FIELD_H__

#ifndef __VECTOR_FIELD_H__
#define __VECTOR_FIELD_H__

#include "../utils/Symbols.h"
#include "VectorFieldExpressions.h"

// forward declaration
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

  // inner product
  InnerProduct<N> dot(const VectorField<N>& rhs) const;
};

template <int N>
SVector<N> VectorField<N>::operator()(const SVector<N> &point) const {
  SVector<N> result;
  for(size_t i = 0; i < N; ++i){
    result[i] = field_[i](point);
  }
  return result;
}

template <int N>
const std::function<double(SVector<N>)>& VectorField<N>::operator[](size_t i) const {
  return field_[i];
}

template <int N>
std::function<double(SVector<N>)>& VectorField<N>::operator[](size_t i) {
  return field_[i];
}

// A functor to represent an inner product.
template <int N>
struct InnerProduct{
  // operands of the inner product operation
  VectorField<N> lhs_;
  VectorField<N> rhs_;

  InnerProduct(const VectorField<N>& lhs, const VectorField<N>& rhs) :
    lhs_(lhs), rhs_(rhs) {}

  double operator()(const SVector<N>& x, const SVector<N>& y){
    // implementation of the scalar product operation
    double result;
    for(size_t i = 0; i < N; ++i){
      result += lhs_[i](x)*rhs_[i](y);
    }
    return result;
  }
};

// return a functor representing the inner product operation.
template <int N>
InnerProduct<N> VectorField<N>::dot(const VectorField<N> &rhs) const {  
  return InnerProduct<N>(this, rhs);
}

// matrix vector-field multiplication
template <int N>
VectorField<N> operator*(const SMatrix<N>& op1, const VectorField<N>& op2) {
  VectorField<N> result;
  for(size_t i = 0; i < N; ++i){
    std::function<double(SVector<N>)> f;

    f = [i, op1, op2](const SVector<N>& p) -> double {
      double product = 0;
      for(size_t j = 0; j < N; ++j){
	product += op1(i,j)*op2[j](p);
      }
      return product;
    };
      
    result[i] = f;
  }
  return result;
}  

#endif // __VECTOR_FIELD_H__

#ifndef __VECTOR_FIELD_H__
#define __VECTOR_FIELD_H__

#include <initializer_list>
#include "../Symbols.h"
#include "ScalarField.h"
#include "VectorFieldExpressions.h"
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
    // allow braced-list initialization
    VectorField(std::initializer_list<std::function<double(SVector<N>)>> field) {
      std::size_t i = 0;
      for(auto it = field.begin(); it != field.end(); ++it){
	field_[i] = *it;
	i++;
      }
    }
    
    // this constructor allows for the construction of a VectorField from a single vectorial lambda. Observe that this method
    // is more inefficient wrt the above constructor since evaluating the field along direction i only still requires the evaluation
    // of all other dimesions. Use with care.
    VectorField(const std::function<SVector<N>(SVector<N>)>& field) {
      for(std::size_t i = 0; i < N; ++i){
	field_[i] = [=](SVector<N> x) -> double { return field(x)[i]; };
      }
    }
    
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

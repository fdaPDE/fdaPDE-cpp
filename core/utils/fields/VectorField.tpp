// call operator
template <int N>
SVector<N> VectorField<N>::operator()(const SVector<N> &point) const {
  SVector<N> result;
  for(size_t i = 0; i < N; ++i){
    // call lambda for each dimension of the vector field
    result[i] = field_[i](point);
  }
  return result;
}

// subscript operator, constant and not constant version
template <int N>
const std::function<double(SVector<N>)>& VectorField<N>::operator[](size_t i) const { return field_[i]; }
template <int N>
std::function<double(SVector<N>)>& VectorField<N>::operator[](size_t i) { return field_[i]; }

// out of class definition for VectorField

// matrix-VectorField product operator. Returns a new vector field where each component of the field represent
// the action of the multiplication between a given matrix and a vector field
template <int N, int M>
VectorField<M> operator*(const Eigen::Matrix<double,N,M>& op1, const VectorField<N>& op2) {
  VectorField<M> result;
  for(size_t i = 0; i < M; ++i){
    std::function<double(SVector<N>)> f;

    // create lambda expression to represent the multiplication between a matrix and the field
    f = [i, op1, op2](const SVector<N>& p) -> double {
      double product = 0;
      // matrix row times vector field
      for(size_t j = 0; j < N; ++j){
	product += op1(i,j)*op2[j](p);
      }
      return product;
    };
      
    result[i] = f;
  }
  return result;
}  

// return a functor representing an inner product.
template <int N>
InnerProduct<N> VectorField<N>::dot(const VectorField<N> &rhs) const {  
  return InnerProduct<N>(*this, rhs);
}

// InnerProduct definitions

// evaluate the form on two different points (this is like evaluating one scalar field at one point,
// the second one at another point and then perform the inner product between the two numerical results)
template <int N>
double InnerProduct<N>::operator()(const SVector<N>& x, const SVector<N>& y) const{
  // implementation of the scalar product operation
  double result;
  for(size_t i = 0; i < N; ++i){
    result += lhs_[i](x)*rhs_[i](y);
  }
  return result;
}

// consider the analytical expression of the inner product as a scalar field, evaluate it at point x
template <int N>
double InnerProduct<N>::operator()(const SVector<N>& x) const{
  // implementation of the scalar product operation
  double result;
  for(size_t i = 0; i < N; ++i){
    result += lhs_[i](x)*rhs_[i](x);
  }
  return result;
}

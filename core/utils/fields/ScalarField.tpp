// approximate computation of gradient via central finite differences
template <unsigned int N>
SVector<N> ScalarField<N>::getGradientApprox(const SVector<N>& x, double step) const {
  SVector<N> gradient;
  for(size_t i = 0; i < N; ++i){
    // variation around point x along direction i
    SVector<N> h = SVector<N>::Zero();
    h[i] = step;

    // approximation of i-th partial derivative
    gradient[i] = (f(x + h) - f(x - h))/(2*step);
  }
  
  return gradient;
}

// approximate computation of hessian matrix via central finite differences
template <unsigned int N>
SMatrix<N> ScalarField<N>::getHessianApprox(const SVector<N>& x, double step) const {
  SMatrix<N> hessian = SMatrix<N>::Zero();

  // hessian matrix is symmetric, compute just the lower triangular part
  for(size_t i = 0; i<N; ++i){
    for(size_t j = 0; j<=i; ++j){
      
      SVector<N> h_i = SVector<N>::Zero();
      h_i[i] = step; // variation along i-th direction

      if(i == j){ 
	// central differences
	hessian(i,i) = (-f(x + 2*h_i) + 16*f(x + h_i) - 30*f(x) + 16*f(x - h_i) - f(x - 2*h_i))/(12*pow(step,2));
      }
      else{ 
	SVector<N> h_j = SVector<N>::Zero();
	h_j[j] = step; // variation along j-th direction

	// central differences
	hessian(i,j) = (f(x + h_i + h_j) - f(x + h_i - h_j) - f(x - h_i + h_j) + f(x - h_i - h_j))/(4*pow(step,2));
	
	// exploit symmetry of hessian matrix
	hessian(j,i) = hessian(i,j);
      }
    }
  }
  
  return hessian;
}

// discretize the scalar field over the given mesh
template <unsigned int N>
template <unsigned int L, unsigned int K>
Eigen::Matrix<double, Eigen::Dynamic, 1> ScalarField<N>::discretize(const Mesh<L, K>& mesh) const{
  Eigen::Matrix<double, Eigen::Dynamic, 1> discretizedField;
  discretizedField.resize(mesh.getNumberOfElements(), 1); // allocate space
  discretizedField.fill(0);                               // init result vector to zero

  // perform discretization
  for(const auto& e : mesh) discretizedField[e->getID()] = f(e->computeMidPoint());

  return discretizedField;
}

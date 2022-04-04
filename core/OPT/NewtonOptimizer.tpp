// newton optimization routine
template <unsigned int N>
template <typename... Args>
std::pair<SVector<N>, double> NewtonOptimizer<N>::findMinimum(const ScalarField<N>& objective, const SVector<N>& x0, const Args&... args) { 
  // can be set true by some extension to cause a foced stop at any point of the execution
  bool customStop = false; 
 
  customStop |= Extension::executeInitOptimization(*this, objective, args...);
  
  // algorithm initialization
  x_old = x0;
  error = objective.getGradientApprox(x0, gradient_step).squaredNorm();

  SVector<N> gradientApprox;  // approximated value of the gradient at each step
  SMatrix<N> hessianApprox;   // approximated value of the hessian at each step
  
  while (numIteration < maxIteration && error > tolerance && !customStop){
    customStop |= Extension::executeInitIteration(*this, objective, args...);

    // compute hessian and gradient approximation
    hessianApprox  = objective.getHessianApprox (x_old, hessian_step);
    gradientApprox = objective.getGradientApprox(x_old, gradient_step);

    // solve linear system by using an Householder QR decomposition with column-pivoting: A*P = Q*R
    Eigen::ColPivHouseholderQR<SMatrix<N>> QRdecomposition(hessianApprox);
    update = QRdecomposition.solve(gradientApprox);

    // update step
    x_new = x_old - step*update;
    error = objective.getGradientApprox(x_new, gradient_step).squaredNorm();

    customStop |= Extension::executeEndIteration(*this, objective, args...);
    // prepare next iteration
    x_old = x_new;    
    numIteration++;
  }

  Extension::executeEndOptimization(*this, objective, args...);
  return std::pair<SVector<N>, double>(x_old, objective(x_old));
}

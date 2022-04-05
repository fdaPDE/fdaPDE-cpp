// newton optimization routine
template <unsigned int N>
template <typename... Args>
std::pair<SVector<N>, double> ExactNewtonOptimizer<N>::findMinimum(const TwiceDifferentiableScalarField<N>& objective, // objective to optimize
								   const SVector<N>& x0,                               // starting point
								   const Args&... args) {                              // extensions
  // can be set true by some extension to cause a foced stop at any point of the execution
  bool customStop = false; 
  customStop |= Extension::executeInitOptimization(*this, objective, args...);
  
  // algorithm initialization
  x_old = x0;
  error = objective.derive()(x0).squaredNorm();

  SVector<N> gradientExact;  // exact value of the gradient at each step
  SMatrix<N> hessianExact;   // exact value of the hessian at each step

  // start loop
  while (numIt < maxIt && error > tolerance && !customStop){
    customStop |= Extension::executeInitIteration(*this, objective, args...);

    // compute hessian and gradient exact by resorting to their analytical expression
    hessianExact  = objective.deriveTwice()(x_old);
    gradientExact = objective.derive()(x_old);

    // solve linear system by using an Householder QR decomposition with column-pivoting: A*P = Q*R
    Eigen::ColPivHouseholderQR<SMatrix<N>> QRdecomposition(hessianExact);
    update = QRdecomposition.solve(gradientExact);

    // update step
    x_new = x_old - step*update;
    error = objective.derive()(x_new).squaredNorm();

    customStop |= Extension::executeEndIteration(*this, objective, args...);
    // prepare next iteration
    x_old = x_new;    
    numIt++;
  }
  
  Extension::executeEndOptimization(*this, objective, args...);
  return std::pair<SVector<N>, double>(x_old, objective(x_old));
}

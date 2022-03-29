// newton optimization routine
template <unsigned int N>
std::pair<SVector<N>, double> NewtonOptimizer<N>::findMinimum() {
  CustomizableOptimizer<N>::init();                 // execute custom action
  
  // algorithm initialization
  x_old = x0;
  unsigned int numIteration = 0;
  error = objective.getGradientApprox(x_new, 0.001).squaredNorm();
  
  while (numIteration < maxIteration && error > tolerance && !CustomizableOptimizer<N>::stopCondition()){
    CustomizableOptimizer<N>::beginIteration();     // execute custom action
    // newton step
    hessianApprox  = objective.getHessianApprox(x_old, 0.001);
    gradientApprox = objective.getGradientApprox(x_old, 0.001);

    // solve linear system by using an Householder QR decomposition with column-pivoting: A*P = Q*R
    Eigen::ColPivHouseholderQR<SMatrix<N>> QRdecomposition(hessianApprox);
    SVector<N> update = QRdecomposition.solve(gradientApprox);
    // solution update
    x_new = x_old - step*update;

    // error update
    error = objective.getGradientApprox(x_new, 0.001).squaredNorm();

    CustomizableOptimizer<N>::endIteration();       // execute custom action
    // prepare next iteration
    x_old = x_new;    
    numIteration++;
  }

  CustomizableOptimizer<N>::finalize();             // execute custom action
  return std::pair<SVector<N>, double>(x_old, objective(x_old));
}

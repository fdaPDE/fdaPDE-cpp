// newton optimization routine
template <unsigned int N>
std::pair<SVector<N>, double> ExactNewtonOptimizer<N>::findMinimum() const {
  this->init();                 // execute custom action
  
  // algorithm initialization
  x_old = x0;
  unsigned int numIteration = 0;
  error = objective.derive()(x_old).norm();
  
  while (numIteration < maxIteration && error > tolerance && !this->stopCondition()){
    this->preStep();            // execute custom action
    
    // newton step
    hessianExact  = objective.deriveTwice()(x_old);
    gradientExact = objective.derive()(x_old);

    // solve linear system by using an Householder QR decomposition with column-pivoting: A*P = Q*R
    Eigen::ColPivHouseholderQR<SMatrix<N>> QRdecomposition(hessianExact);
    
    SVector<N> update = QRdecomposition.solve(gradientExact);
    x_new = x_old - NewtonOptimizer<N>::step*update;

    // error update
    error = objective.derive()(x_new).norm();

    this->postStep();           // execute custom action
    
    // prepare next iteration
    x_old = x_new;    
    numIteration++;
  }

  this->finalize();             // execute custom action
  return std::pair<SVector<N>, double>(x_new, objective(x_new));
}

// gradient descent optimization routine
template <unsigned int N>
std::pair<SVector<N>, double> GradientDescentOptimizer<N>::findMinimum(){
  this->init();             // execute custom action
  
  // algorithm initialization
  x_old = x0;
  unsigned int numIteration = 0;
  
  // standard termination criteria based on l^2 norm of the gradient
  error = objective.derive()(x_old).squaredNorm();
  
  while (numIteration < maxIteration && error > tolerance && !this->stopCondition()){
    this->preStep();        // execute custom action
    
    // compute exact gradient
    gradientExact = objective.derive()(x_old);
    // update step    
    x_new = x_old - step*gradientExact;
    // error update: standard termination criteria based on l^2 norm of the gradient
    error = gradientExact.squaredNorm();

    this->postStep();       // execute custom action
    // prepare next iteration    
    x_old = x_new;    
    numIteration++;
  }
  
  this->finalize();         // execute custom action
  return std::pair<SVector<N>, double>(x_old, objective(x_old));
}

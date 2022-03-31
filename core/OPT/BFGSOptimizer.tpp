// BFGS optimization routine
template <unsigned int N>
std::pair<SVector<N>, double> BFGSOptimizer<N>::findMinimum(){
  CustomizableOptimizer<N>::init();             // execute custom action
  
  // algorithm initialization
  x_old = x0;
  unsigned int numIteration = 0;

  // set the hessian approximation at first iteration equal to the identity matrix
  // this will make the first iteration of the algorithm to behave as a gradient descent
  // step. hessian will be corrected as the algorithm goes on.
  hessian = SMatrix<N>::Identity();

  grad_old = objective.derive()(x_old);
  if (grad_old.isApprox(SVector<N>::Zero())) // gradient is zero, already at stationary point
    return std::pair<SVector<N>, double>(x_old, objective(x_old));

  error = grad_old.squaredNorm();
  
  while(numIteration < maxIteration && error > tolerance && !CustomizableOptimizer<N>::stopCondition()){
    CustomizableOptimizer<N>::beginIteration(); // execute custom action
    // compute algorithm update direction
    SVector<N> updateDirection = hessian*grad_old;

    // backtracking based step search
    double ro = 0.5, alpha = 1/ro, gamma = 0.5;
    do{
      // update step
      alpha *= ro;
    }while( objective(x_old + step*updateDirection) > objective(x_old) + gamma*alpha*(grad_old.dot(updateDirection)));
    
    // update along descent direction
    x_new = x_old - alpha*updateDirection;
    // gradient update
    grad_new = objective.derive()(x_new);
    if (grad_new.isApprox(SVector<N>::Zero())){ // gradient is zero, already at stationary point
      return std::pair<SVector<N>, double>(x_new, objective(x_new));
    }
      
    // update inverse hessian approximation
    SVector<N> deltaX    = x_new - x_old;
    SVector<N> deltaGrad = grad_new - grad_old;
    double xg = deltaX.dot(deltaGrad);     // inner product between deltaX and deltaGrad
    SVector<N> hx = hessian*deltaGrad;     // product between hessian matrix and deltaGrad

    // see references for detailed derivation of the equations
    SMatrix<N> U = (1 + (deltaGrad.dot(hx))/xg)*((deltaX*deltaX.transpose())/xg);
    SMatrix<N> V = ((hx*deltaX.transpose() + deltaX*hx.transpose()))/xg;
    
    hessian += U - V; // hessian approximation update

    CustomizableOptimizer<N>::endIteration();   // execute custom action
    // prepare next iteration
    error = grad_new.squaredNorm();
    x_old = x_new;
    grad_old = grad_new;
    numIteration++;
  }

  CustomizableOptimizer<N>::finalize();         // execute custom action
  return std::pair<SVector<N>, double>(x_old, objective(x_old));
}

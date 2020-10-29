#ifndef __OPTIMIZATION_ALGORITHM_IMP_H__
#define __OPTIMIZATION_ALGORITHM_IMP_H__


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>::
  MinimizationAlgorithm(const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dp,
  const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& fp, const std::string& d):
  dataProblem_(dp), funcProblem_(fp){

    direction_ = DescentDirection_factory<Integrator_noPoly, ORDER,  mydim,  ndim>::createDirectionSolver(dp, fp, d);

};


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>::
  MinimizationAlgorithm(const MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>& rhs):
  dataProblem_(rhs.dataProblem_), funcProblem_(rhs.funcProblem_){

    direction_ = rhs.direction_->clone();
};


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>>
FixedStep<Integrator_noPoly, ORDER, mydim, ndim>::clone() const {

  return make_unique<FixedStep<Integrator_noPoly, ORDER, mydim, ndim>>(*this);

}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
VectorXr
FixedStep<Integrator_noPoly,ORDER,mydim,ndim>::apply_core(const SpMat& Psi, Real lambda, const VectorXr& g) const{

  // termination criteria variables
  const Real toll1 = this->dataProblem_.getTol1(), toll2 = this->dataProblem_.getTol2();
  Real norm_grad, dloss = toll1 +1, dllik = toll1 +1, dpen = toll1 +1;

  // to save the current point
  VectorXr g_curr;

  // variables
  VectorXr grad, d;
  Real loss, loss_old, llik, llik_old, pen, pen_old;
  UInt i;

  for(UInt e = 0; e < this->dataProblem_.getNstepProposals(); e++){
    // start always with the initial point
    g_curr = g;

    std::tie(loss, grad, llik, pen) = this->funcProblem_.computeFunctional_g(g_curr, lambda, Psi);
    norm_grad = std::sqrt(grad.dot(grad));

    if(this->dataProblem_.Print()){
      Rprintf("loss %f, llik %f, pen %f, norm_Lp %f\n", loss, llik, pen, norm_grad);
    }

    for(i = 0; i < this->dataProblem_.getNsimulations() && (dloss > toll1 || dllik > toll1 || dpen > toll1) && norm_grad > toll2; i++){

      // Termination criteria variables
      loss_old = loss;
      llik_old = llik;
      pen_old = pen;

      // Compute a descent direction
      d = this->direction_->computeDirection(g_curr, grad);

      // Update the point
      g_curr = g_curr + this->dataProblem_.getStepProposals(e)*d;

      // Update termination criteria variables
      std::tie(loss, grad, llik, pen) = this->funcProblem_.computeFunctional_g(g_curr, lambda, Psi);
      dloss = std::abs((loss - loss_old)/loss_old);
      dllik = std::abs((llik - llik_old)/llik_old);
      dpen = std::abs((pen - pen_old)/pen_old);
      norm_grad = std::sqrt(grad.dot(grad));

      // Check if the step is ok
      if((loss_old - loss) < 0){
        if(this->dataProblem_.Print()){
          Rprintf("The loss function increases: not good. Try decreasing the optimization parameter.\n", norm_grad);
        }
        break;
      }

      if(this->dataProblem_.Print()){
        Rprintf("Iter %d, loss %f, llik %f, pen %f, norm_Lp %f\n", i+1, loss, llik, pen, norm_grad);
      }
    }

    this->direction_->resetParameters();

    if ((loss_old - loss) < 0){}
    else if(dloss <= toll1 && dllik <= toll1 && dpen <= toll1){
      if(this->dataProblem_.Print()){
        Rprintf("The algorithm reaches the tollerance in terms of the functional. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
      }
      return g_curr;
    }
    else if(norm_grad <= toll2){
      if(this->dataProblem_.Print()){
        Rprintf("The algorithm reaches the tollerance in terms of the slope. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
      }
      return g_curr;
    }
    else if(i == this->dataProblem_.getNsimulations()){
      if(this->dataProblem_.Print()){
        Rprintf("The algorithm reaches the maximum number of iterations. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
      }
      return g_curr;
    }
  }

  // If you arrive here you don't have a good gradient parameter
  Rprintf("ERROR: The loss function increases: not good. Try decreasing the optimization parameter");
  //std::abort();
  return VectorXr::Constant(g.size(),0);

}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
VectorXr
AdaptiveStep<Integrator_noPoly,ORDER,mydim,ndim>::apply_core(const SpMat& Psi, Real lambda, const VectorXr& g) const{

  // termination criteria variables
  const Real toll1 = this->dataProblem_.getTol1(), toll2 = this->dataProblem_.getTol2();
  Real norm_grad, dloss = toll1 +1, dllik = toll1 +1, dpen = toll1 +1;

  // to save the current point
  VectorXr g_curr = g;

  // variables
  VectorXr grad, d;
  Real loss, loss_old, llik, llik_old, pen, pen_old, step;
  UInt i;

  std::tie(loss, grad, llik, pen) = this->funcProblem_.computeFunctional_g(g_curr, lambda, Psi);
  norm_grad = std::sqrt(grad.dot(grad));

  if(this->dataProblem_.Print()){
    Rprintf("loss %f, llik %f, pen %f, norm_Lp %f\n", loss, llik, pen, norm_grad);
  }

  for(i = 0; i < this->dataProblem_.getNsimulations() && (dloss > toll1 || dllik > toll1 || dpen > toll1) && norm_grad > toll2; i++){

    // Termination criteria variables
    loss_old = loss;
    llik_old = llik;
    pen_old = pen;

    // Compute a descent direction
    d = this->direction_->computeDirection(g_curr, grad);

    // Compute a step
    step = computeStep(g_curr, loss, grad, d, lambda, Psi);

    // Update the point
    g_curr = g_curr + step*d;

    // Update termination criteria variables
    std::tie(loss, grad, llik, pen) = this->funcProblem_.computeFunctional_g(g_curr, lambda, Psi);
    dloss = std::abs((loss - loss_old)/loss_old);
    dllik = std::abs((llik - llik_old)/llik_old);
    dpen = std::abs((pen - pen_old)/pen_old);
    norm_grad = std::sqrt(grad.dot(grad));

    if(this->dataProblem_.Print()){
      Rprintf("Iter %d, loss %f, llik %f, pen %f, norm_Lp %f\n", i+1, loss, llik, pen, norm_grad);
    }

  }

  this->direction_->resetParameters();

  if(dloss <= toll1 && dllik <= toll1 && dpen <= toll1){
    if(this->dataProblem_.Print()){
      Rprintf("The algorithm reaches the tollerance in terms of the functional. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
    }
    return g_curr;
  }
  else if(norm_grad <= toll2){
    if(this->dataProblem_.Print()){
      Rprintf("The algorithm reaches the tollerance in terms of the slope. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
    }
    return g_curr;
  }
  else{
    if(this->dataProblem_.Print()){
      Rprintf("The algorithm reaches the maximum number of iterations. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
    }
    return g_curr;
  }
}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
Real
BacktrackingMethod<Integrator_noPoly,ORDER,mydim,ndim>::computeStep(const VectorXr& g, Real loss, const VectorXr& grad, const VectorXr& dir, Real lambda, const SpMat& Psi) const{

  Real ro = 0.5, alpha = 1/ro, c = 0.5;

  Real loss_new, llik_new, pen_new, slope, grad_dir;
  VectorXr grad_new, new_point;

  grad_dir =  grad.dot(dir);

  do{
    // update step
    alpha *= ro;

    slope = c*alpha*(grad_dir);

    // Update the point
    new_point = g + alpha*dir;

    // functional in the new point
    std::tie(loss_new, grad_new, llik_new, pen_new) = this->funcProblem_.computeFunctional_g(new_point, lambda, Psi);

  } while(loss_new > (loss + slope));

  return alpha;
}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>>
BacktrackingMethod<Integrator_noPoly, ORDER, mydim, ndim>::clone() const {

  return make_unique<BacktrackingMethod<Integrator_noPoly, ORDER, mydim, ndim>>(*this);

}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
Real
WolfeMethod<Integrator_noPoly,ORDER,mydim,ndim>::computeStep(const VectorXr& g, Real loss, const VectorXr& grad, const VectorXr& dir, Real lambda, const SpMat& Psi) const{

  Real alpha = 1, alphamax = 0, alphamin = 0, c1 = 1e-4, c2 = 0.9;

  Real loss_new, llik_new, pen_new, slope, grad_dir;
  VectorXr grad_new, new_point;

  grad_dir = grad.dot(dir);
  slope = c1*alpha*grad_dir;


  // Update the point
  new_point = g + alpha*dir;

  // functional in the new point
  std::tie(loss_new, grad_new, llik_new, pen_new) = this->funcProblem_.computeFunctional_g(new_point, lambda, Psi);

	bool again = true;

	while(again){

		again = false;

		while(loss_new > (loss + slope)){
      // update step
			alphamax = alpha;
			alpha = 0.5*(alphamin + alphamax);

      // try with the new point
			new_point = g + alpha*dir;
      std::tie(loss_new, grad_new, llik_new, pen_new) = this->funcProblem_.computeFunctional_g(new_point, lambda, Psi);
      slope = c1*alpha*grad_dir;
		}

		if(grad_new.dot(dir) < c2*grad_dir && std::abs(grad_dir)>1e-2){

			again = true;

      // update step
			alphamin = alpha;
			alpha = alphamax==0 ? 2*alphamin : 0.5*(alphamin+alphamax);

      // try with the new point
      new_point = g + alpha*dir;
      std::tie(loss_new, grad_new, llik_new, pen_new) = this->funcProblem_.computeFunctional_g(new_point, lambda, Psi);
			slope =  alpha*c1*grad_dir;
		}
	}

  return alpha;
}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>>
WolfeMethod<Integrator_noPoly, ORDER, mydim, ndim>::clone() const {

  return make_unique<WolfeMethod<Integrator_noPoly, ORDER, mydim, ndim>>(*this);

}



#endif

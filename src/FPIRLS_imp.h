#ifndef FPIRLS_IMP_H
#define FPIRLS_IMP_H

#include "FPIRLS.h"


/*********** FPIRLS_base Methods ************/

// Constructor
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::FPIRLS_Base(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0, bool scale_parameter_flag, Real scale_param):
  mesh_(mesh), inputData_(inputData), regression_(mesh, inputData_), scale_parameter_flag_(scale_parameter_flag), _scale_param(scale_param)
{
  //initialization of mu, current_J_values and past_J_values.
  for(UInt j=0; j< inputData.getLambdaS().size() ; j++){
    mu_.push_back(mu0);
    current_J_values.push_back(std::array<Real,2>{1,1});
    past_J_values.push_back(std::array<Real,2>{1,1});
  }

};

// FPIRLS_base methods 
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::apply( const ForcingTerm& u){
  // f-PRILS implementation
  
  //Initialize the containers size, as LambdaS_len
  const UInt LambdaS_len = mu_.size();
  
  // initialize the algorithm variables
  G_.resize(LambdaS_len);
  WeightsMatrix_.resize(LambdaS_len);
  pseudoObservations_.resize(LambdaS_len);
  n_iterations = std::vector<UInt>(LambdaS_len,0);

  // Initialize the outputs. The temporal dimension is not implemented, for this reason the 2nd dimension is set to 1.
  if( this->inputData_.getCovariates().rows() > 0 )_beta_hat.resize(LambdaS_len,1);
  _fn_hat.resize(LambdaS_len,1);
  _dof.resize(LambdaS_len,1);
  _solution.resize(LambdaS_len,1);
   
  _GCV.resize(LambdaS_len,-1);//If GCV is not computed the vector stores -1
  

  if(isSpaceVarying)
  {
    FiniteElement<Integrator, ORDER, mydim, ndim> fe;
  	Assembler::forcingTerm(mesh_, fe, u, forcingTerm);
  }

  
  for(UInt i=0 ; i < LambdaS_len ; i++){//for-cycle for each spatial penalization (lambdaS). 

    current_J_values[i][0] = past_J_values[i][0] + 2*inputData_.get_treshold();
    current_J_values[i][1] = past_J_values[i][1] + 2*inputData_.get_treshold();

    this->inputData_.setCurrentLambda(i); // set right lambda for the current iteration.
    
    #ifdef R_VERSION_
    Rprintf("Start FPIRLS for the lambda number %d \n", i+1);
    #else
    std::cout<<"Start FPIRLS for the lambda number "<<i+1<<std::endl;
    #endif
    
    
    // start the iterative method for the lambda index i
    while(stopping_criterion(i)){ 
    
      // STEP (1)
      
      compute_G(i);
      compute_Weights(i);
      compute_pseudoObs(i);

      // STEP (2)
      
      this->inputData_.updatePseudodata(pseudoObservations_[i], WeightsMatrix_[i]);
      update_solution(i);

      // STEP (3)
      compute_mu(i);

      // update J
      past_J_values[i] = current_J_values[i];
      current_J_values[i] = compute_J(i);
          
      n_iterations[i]++;

    } //end while
    
    #ifdef R_VERSION_
    Rprintf("\t n. iterations: %d\n \n", n_iterations[i]);
    #else
    std::cout<< "\t n. iterations: "<<n_iterations[i]<<"\n"<<std::endl;
    #endif
    
    _J_minima.push_back(current_J_values[i][0]+current_J_values[i][1]); // compute the minimum value of the J fuctional 

    if(this->inputData_.getGCV_GAM()){ // compute GCV if it is required
      compute_GCV(i);
    }

  }// end for

  // Variance Estimate
  compute_variance_est();
}

template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::update_solution(UInt& lambda_index){
  // performs step (2) of PIRLS. It requires pseudo data after step(1) and mimic regression skeleton behaviour

  // Here we have to solve a weighted regression problem. 
  regression_.recomputeWTW(); // at each iteration of FPIRLS W is updated, so WTW has to be recomputed as well.
  regression_.apply();
  const SpMat& Psi = regression_.getPsi(); // get Psi matrix. It is used for the computation of fn_hat.

  // get the solutions from the regression object.
  _solution(lambda_index,0) = regression_.getSolution()(0,0);
  _dof(lambda_index,0) = regression_.getDOF()(0,0);

  if(inputData_.getCovariates().rows()>0){ 
    _beta_hat(lambda_index,0) = regression_.getBeta()(0,0);
  }

  _fn_hat(lambda_index,0) = Psi *_solution(lambda_index,0).topRows(Psi.cols());


}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_pseudoObs(UInt& lambda_index){
  // compute pseudodata observations
  
  VectorXr first_addendum; // G_ii( z_i - mu_i)
  VectorXr g_mu; // g( mu_i )
  
  VectorXr z = inputData_.getInitialObservations();

  first_addendum.resize(mu_[lambda_index].size());
  g_mu.resize(mu_[lambda_index].size());

  //compute the vecotr first_addendum and g_mu
  for(auto i=0; i < mu_[lambda_index].size(); i++){
    g_mu(i) = link(mu_[lambda_index](i));
    first_addendum(i) = G_[lambda_index](i)*(z(i)-mu_[lambda_index](i));
  }

  pseudoObservations_[lambda_index] = first_addendum + g_mu;

}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_G(UInt& lambda_index){
  // compute the G matrix as G_ii = diag( g'(mu_i))
  
  G_[lambda_index].resize(mu_[lambda_index].size());

  for(UInt i = 0; i<mu_[lambda_index].size(); i++){
    G_[lambda_index](i) = link_deriv(mu_[lambda_index](i));
  }

}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_Weights(UInt& lambda_index){
  // computed W elementwise (it is a diagonal matrix)

  WeightsMatrix_[lambda_index].resize( mu_[lambda_index].size());

  for(auto i=0; i < mu_[lambda_index].size(); i++){
    WeightsMatrix_[lambda_index](i) = 1/(pow(G_[lambda_index](i),2)*(var_function( mu_[lambda_index](i))));
  }

}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_mu(UInt& lambda_index){
  //compute mu as mu_i = g-1( w_ii*beta + fn_hat)

  VectorXr W_beta = VectorXr::Zero(mu_[lambda_index].size()); // initialize the vector w_ii*beta

  if(inputData_.getCovariates().rows()>0)
    W_beta = inputData_.getCovariates()*_beta_hat(lambda_index,0);


  for(UInt j=0; j < W_beta.size(); j++){
      mu_[lambda_index](j) = inv_link(W_beta[j] + _fn_hat(lambda_index,0)(j));
  }

}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
bool FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::stopping_criterion(UInt& lambda_index){
  // return true if the f-PIRLS has to perform another iteration, false if it has to be stopped
  
  bool do_stop_by_iteration = false;  // Do I need to stop becouse n_it > n_max?
  bool do_stop_by_treshold = false; // Do I need to stop becouse |J(k) - J(k+1)| < treshold?

  if(n_iterations[lambda_index] > inputData_.get_maxiter()){
    do_stop_by_iteration = true;
  }

  if(n_iterations[lambda_index] > 1){
    if(abs(past_J_values[lambda_index][0]+past_J_values[lambda_index][1] - current_J_values[lambda_index][0] - current_J_values[lambda_index][1]) < inputData_.get_treshold()){
        do_stop_by_treshold = true;
   }
  }

  return !(do_stop_by_iteration || do_stop_by_treshold );
}

template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
std::array<Real,2> FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_J(UInt& lambda_index){
  // compute the functional J: it is divided in parametric and non parametric part
  Real parametric_value = 0;
  Real non_parametric_value = 0;
  Real tmp;


  VectorXr Lf;

  const VectorXr z = inputData_.getInitialObservations();

  for(UInt i=0; i < mu_.size(); i++){
    tmp = sqrt( var_function( mu_[lambda_index](i)) ) * (z(i) - mu_[lambda_index](i)) ;
    parametric_value += tmp*tmp;
  }

  Lf.resize(_solution(lambda_index,0).size()/2);
  for(UInt i=0; i< Lf.size(); i++){
    Lf(i) = _solution(lambda_index,0)(Lf.size() + i);
  }

  if(isSpaceVarying)
  {
      Lf = Lf - forcingTerm;
  }

  non_parametric_value = Lf.transpose() * regression_.getR0() * Lf;
  non_parametric_value = inputData_.getGlobalLambda()[lambda_index]*non_parametric_value;

  std::array<Real,2> returnObject{parametric_value, non_parametric_value};

  return returnObject;
}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_GCV(UInt& lambda_index){

  //GCV COMPUTATION
  if ( inputData_.computeDOF() ){ // is DOF_matrix to be computed?
    regression_.computeDegreesOfFreedom(0, 0, this->inputData_.getGlobalLambda()[lambda_index], 0);
  }
  _dof(lambda_index,0) = regression_.getDOF()(0,0);

  VectorXr y = inputData_.getInitialObservations();
  Real GCV_value = 0;

  for(UInt j=0; j < y.size();j++)
    GCV_value += dev_function(mu_[lambda_index][j], y[j]); //norm computation

  GCV_value *= y.size();

  GCV_value /= (y.size()-inputData_.getTuneParam()*_dof(lambda_index,0))*(y.size()-inputData_.getTuneParam()*_dof(lambda_index,0));

  _GCV[lambda_index] = GCV_value;

  //best lambda
  if ( GCV_value < _bestGCV)
  {
    bestLambdaS_ = lambda_index;
    _bestGCV = GCV_value;
  }

}

template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_variance_est(){
  Real phi; 
  if(this->scale_parameter_flag_ && !this->inputData_.getGCV_GAM() ){// if scale param should be 
    _variance_estimates.resize(this->mu_.size(),0);
    const UInt n_obs = this->inputData_.getObservations().size();

    //scale parameter computed as: mean((var.link(mu)*phi)/mu), and phi is computed as in Wood IGAM
    for(UInt i=0; i < this->mu_.size();i++){
      phi = (this->scale_parameter_flag_ )? this->current_J_values[i][0]/(n_obs - this->_dof(i,0) ) : _scale_param;
      for(UInt j=0; j < this->mu_[i].size(); j++){
        _variance_estimates[i] += phi* this->var_function(this->mu_[i][j])/this->mu_[i][j];
      }
      _variance_estimates[i] /= this->mu_[i].size();
    }
  }else{
    _variance_estimates.push_back(-1);
  }
}

/*********** FPIRLS apply template specialization ************/

// Laplace or Elliptic case
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS<InputHandler,Integrator,ORDER, mydim, ndim>::apply(){

  FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::apply(ForcingTerm(std::vector<Real>(1)));

}

// SpaceVarying case
template <typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS<GAMDataEllipticSpaceVarying,Integrator,ORDER, mydim, ndim>::apply(){

  this->isSpaceVarying = true;
  FPIRLS_Base<GAMDataEllipticSpaceVarying,Integrator,ORDER, mydim, ndim>::apply(this->inputData_.getU());

}


#endif

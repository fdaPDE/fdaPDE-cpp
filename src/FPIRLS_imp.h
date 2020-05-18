#ifndef FPIRLS_IMP_H
#define FPIRLS_IMP_H

#include "FPIRLS.h"

template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::FPIRLS_Base(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0):
  mesh_(mesh), inputData_(inputData), regression_(mesh, inputData)
{
  Rprintf("Hello I'm FPIRLS_Base constructor \n");

  //mu_(inputData.getLambda().size(),mu0);
  //current_J_values(inputData.getLambda().size(),std::array<Real,2>{1,1});
  //past_J_values(inputData.getLambda().size(),std::array<Real,2>{1,1});

  std::string saving_filename = "TEST_COSTRUTTORE";
  saving_filename = saving_filename + ".txt";
  printer::saveVectorXr(saving_filename,mu0);


  //MU,J
  for(UInt j=0; j< inputData.getLambdaS().size() ; j++){
    mu_.push_back(mu0);
    current_J_values.push_back(std::array<Real,2>{1,1});
    past_J_values.push_back(std::array<Real,2>{1,1});
  }

}; // Constructor


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::apply( const ForcingTerm& u){
  // f-PRILS implementation
  std::string saving_filename = "FPIRLS_applyInit_obs.txt";
  printer::saveVectorXr(saving_filename, this->inputData_.getObservations());

  saving_filename = "FPIRLS_applyInit_WeightMat.txt";
  printer::saveVectorXr(saving_filename, this->inputData_.getWeightsMatrix());


  this->inputData_.copyInitialObservations();
  this->inputData_.copyLambdaVector();
  this->inputData_.setDOFflag();

  const UInt LambdaS_len = mu_.size();
  // STEP 0: initialization
  //  if(mu_.size()==0) //mu can be initialized or not by the user
      //initialize_mu(inputData_.getObservations());

  G_.resize(LambdaS_len);
  WeightsMatrix_.resize(LambdaS_len);
  if( this->inputData_.getCovariates().rows() > 0 ) _beta_hat.resize(LambdaS_len,1);
  _fn_hat.resize(LambdaS_len,1);
  _dof.resize(LambdaS_len,1);
  _solution.resize(LambdaS_len,1);
  pseudoObservations_.resize(LambdaS_len);
  n_iterations = std::vector<UInt>(LambdaS_len,0);
  
  _GCV.resize(LambdaS_len,-1);

  if(isSpaceVarying)
  {
    FiniteElement<Integrator, ORDER, mydim, ndim> fe;
  	Assembler::forcingTerm(mesh_, fe, u, forcingTerm);
  }

  saving_filename = "TEST_1";
  saving_filename = saving_filename + ".txt";
  printer::saveVectorXr(saving_filename,mu_[0]);

  printer::variableInt("lambdaS_len.txt", LambdaS_len);
  printer::variableInt("global_lambda_len.txt", inputData_.getGlobalLambda().size());
  printer::variableInt("local_lambda_len.txt", inputData_.getLambdaS().size());

for(UInt i=0 ; i < LambdaS_len ; i++){

  current_J_values[i][0] = past_J_values[i][0] + 2*inputData_.get_treshold();
  current_J_values[i][1] = past_J_values[i][1] + 2*inputData_.get_treshold();


    saving_filename = "TEST_2";
    saving_filename = saving_filename + ".txt";
    printer::saveVectorXr(saving_filename,mu_[0]);



  this->inputData_.setCurrentLambda(i); // set right lambda for the iteration

  while(stopping_criterion(i))  // n_iterations < inputData_.get_maxiter() && past_J_value - current_J_value < inputData_.get_treshold()
  {

  Rprintf("FPIRLS.apply: while iteration number: %d \n", n_iterations[i]);
  // STEP (1)

  saving_filename = "TEST_3";
  saving_filename = saving_filename + ".txt";
  printer::saveVectorXr(saving_filename,mu_[0]);


    compute_G(i);

    saving_filename = "TEST_4";
    saving_filename = saving_filename + ".txt";
    printer::saveVectorXr(saving_filename,mu_[0]);


    compute_Weights(i);

    saving_filename = "TEST_5";
    saving_filename = saving_filename + ".txt";
    printer::saveVectorXr(saving_filename,mu_[0]);


    compute_pseudoObs(i);


    saving_filename = "TEST_6";
    saving_filename = saving_filename + ".txt";
    printer::saveVectorXr(saving_filename,mu_[0]);


  // STEP (2)
    this->inputData_.updatePseudodata(pseudoObservations_[i], WeightsMatrix_[i]);
    update_solution(i); // here I'm performing FERegression using appropriate objects. I need pseudo data and mesh and it computes solution and dof



    saving_filename = "TEST_7";
    saving_filename = saving_filename + ".txt";
    printer::saveVectorXr(saving_filename,mu_[0]);

  // STEP (3)

    compute_mu(i);


    saving_filename = "TEST_8";
    saving_filename = saving_filename + ".txt";
    printer::saveVectorXr(saving_filename,mu_[0]);

  // update J

    past_J_values[i] = current_J_values[i];
    current_J_values[i] = compute_J(i);
    
    Rprintf("\t current_J_values[0] : %f \n", current_J_values[i][0]);
    Rprintf("\t current_J_values[1] : %f \n", current_J_values[i][1]);

    n_iterations[i]++;

    saving_filename = "TEST_9";
    saving_filename = saving_filename + ".txt";
    printer::saveVectorXr(saving_filename,mu_[0]);


  } //end while

  _J_minima.push_back(current_J_values[i][0]+current_J_values[i][1]);


  if(this->inputData_.getDOF_GAM()){
    saving_filename = "TEST_GCV";
    saving_filename = saving_filename + ".txt";
    printer::saveVectorXr(saving_filename,mu_[0]);
    compute_GCV(i);
  }

}// end for


saving_filename = "TEST_10";
saving_filename = saving_filename + ".txt";
printer::saveVectorXr(saving_filename,mu_[0]);


}// end apply


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::update_solution(UInt& lambda_index){
  // performs step (2) of PIRLS. It requires pseudo data after step(1) and mimic regression skeleton behaviour

  // Here we have to solve a weighted regression problem with laplacian penality
  // COMPOSE PSEUDO DATA
  // I SetMethods di Regression data sono protected, perciò non posso modificare l'oggetto. Tuttavia, potrei creare un metodo ad-hoc in RegressiondDataGAM che faccia ciò di cui ho bisogno, ovvero:
  //      1. l'oggetto di tipo Regression Data abbia tutti i parametri necessari per compiere la regressione
  //      2. observations substituted by z, WeightsMatrix substituted by W

  //Set up regression
  
  //MixedFERegression<InputHandler, Integrator,ORDER, IntegratorGaussP3, 0, 0, mydim, ndim>  regression(mesh_, inputData_);
  std::string saving_filename;

    printer::milestone("first_regression_apply.txt");
    
    regression_.apply();
    const SpMat& Psi = regression_.getPsi();

    saving_filename = std::to_string(n_iterations[lambda_index]) + "_regression_apply.txt"; 
    printer::milestone(saving_filename);

  saving_filename = "bc_index";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::STDvector(  saving_filename ,inputData_.getDirichletIndices());

  saving_filename = "bc_index_dim";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::SaveDimension(saving_filename, inputData_.getDirichletIndices());

  saving_filename = "bc_values_dim";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::SaveDimension(saving_filename, inputData_.getDirichletValues());

  saving_filename = "bc_values";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::STDvector(  saving_filename ,inputData_.getDirichletValues());


  printer::milestone("after_MixedFEregression.txt");

  printer::milestone("after_MixedFEregression_apply.txt");

  _solution(lambda_index,0) = regression_.getSolution()(0,0);
  printer::milestone("after_solution_lambda_index.txt");
  _dof(lambda_index,0) = regression_.getDOF()(0,0);
  printer::milestone("after_dof_lambda_index.txt");
  if(inputData_.getCovariates().rows()>0){
    _beta_hat(lambda_index,0) = regression_.getBeta()(0,0);
    saving_filename = "_beta_hat";
    saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
    printer::saveVectorXr(saving_filename,_beta_hat(lambda_index,0));
  }
  _fn_hat(lambda_index,0) = Psi *_solution(lambda_index,0).topRows(Psi.cols());
    saving_filename = "_fn_hat";
    saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
    printer::saveVectorXr(saving_filename,_fn_hat(lambda_index,0));

  saving_filename = "solution_entire_";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::saveVectorXr(saving_filename,_solution(lambda_index,0));

  saving_filename = "Psi_.txt";
  printer::SaveMatrixXr(saving_filename, regression_.getPsi());


}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_pseudoObs(UInt& lambda_index){

  VectorXr g_mu;
  VectorXr first_addendum;
  VectorXr z = inputData_.getInitialObservations();

  g_mu.resize(mu_[lambda_index].size());
  first_addendum.resize(mu_[lambda_index].size());

  for(auto i=0; i < mu_[lambda_index].size(); i++){
    g_mu(i) = link(mu_[lambda_index](i));
    first_addendum(i) = G_[lambda_index](i)*(z(i)-mu_[lambda_index](i));
  }

  pseudoObservations_[lambda_index] = first_addendum + g_mu;

  std::string saving_filename = "pseudoObservations_";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::saveVectorXr(saving_filename,pseudoObservations_[lambda_index]);

}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_G(UInt& lambda_index){

  G_[lambda_index].resize(mu_[lambda_index].size());

  for(UInt i = 0; i<mu_[lambda_index].size(); i++){
    G_[lambda_index](i) = link_deriv(mu_[lambda_index](i));
  }

  std::string saving_filename = "G_";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::saveVectorXr(saving_filename,G_[lambda_index]);

}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_Weights(UInt& lambda_index){
// computed W elementwise (it is a diagonal matrix)

  WeightsMatrix_[lambda_index].resize( mu_[lambda_index].size());

  for(auto i=0; i < mu_[lambda_index].size(); i++){
    WeightsMatrix_[lambda_index](i) = 1/(pow(G_[lambda_index](i),2)*(var_function( mu_[lambda_index](i))));
  }

  std::string saving_filename = "WeightsMatrix_";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::saveVectorXr(saving_filename,WeightsMatrix_[lambda_index]);

}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_mu(UInt& lambda_index){


  VectorXr W_beta = VectorXr::Zero(mu_[lambda_index].size());

  if(inputData_.getCovariates().rows()>0)
    W_beta = inputData_.getCovariates()*_beta_hat(lambda_index,0);

  std::string saving_filename = "W_beta";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::saveVectorXr(saving_filename,W_beta);

  for(UInt j=0; j < W_beta.size() ; j++){
      mu_[lambda_index](j) = inv_link(W_beta[j] + _fn_hat(lambda_index,0)(j));
  }

  saving_filename = "mu_";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::saveVectorXr(saving_filename,mu_[lambda_index]);

}//end method


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



// Laplace or Elliptic
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS<InputHandler,Integrator,ORDER, mydim, ndim>::apply(){

  FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::apply(ForcingTerm(std::vector<Real>(1)));

}

// SpaceVarying
template <typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS<GAMDataEllipticSpaceVarying,Integrator,ORDER, mydim, ndim>::apply(){

  this->isSpaceVarying = true;
  FPIRLS_Base<GAMDataEllipticSpaceVarying,Integrator,ORDER, mydim, ndim>::apply(this->inputData_.getU());

}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Scaled<InputHandler, Integrator, ORDER, mydim, ndim>::apply(){

  FPIRLS<InputHandler, Integrator, ORDER, mydim, ndim>::apply();

  // Phi estimate
  if(scale_parameter_flag_ && this->_dof(0,0)>=0 ){ //if scale_flag is true and the dofs have been estimated
    compute_scale_param();
  }

}

template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Scaled<InputHandler, Integrator, ORDER, mydim, ndim>::compute_scale_param(){
  Real phi;

  _scale_parameter_estimates.resize(this->mu_.size(),0);
  const UInt n_obs = this->inputData_.getObservations().size();

  //scale parameter computed as: mean((var.link(mu)*phi)/mu), and phi is computed as in Wood IGAM
  for(UInt i=0; i < this->mu_.size();i++){
    phi = this->current_J_values[i][0]/(n_obs - this->_dof(i,0) );
    for(UInt j=0; j < this->mu_[i].size(); j++){
      _scale_parameter_estimates[i] += phi* this->var_function(this->mu_[i][j])/this->mu_[i][j];
    }
    _scale_parameter_estimates[i] /= this->mu_[i].size();
  }
  
}

template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
std::array<Real,2> FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_J(UInt& lambda_index){

  Real parametric_value = 0;
  Real non_parametric_value = 0;
  Real tmp;


  VectorXr Lf;

  const VectorXr z = inputData_.getInitialObservations();
    std::string saving_filename = "z";
    saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
    printer::saveVectorXr(saving_filename,z);

  for(UInt i=0; i < mu_.size(); i++){
    tmp = sqrt( var_function( mu_[lambda_index](i)) ) * (z(i) - mu_[lambda_index](i)) ;
    parametric_value += tmp*tmp;
  }
  Rprintf("\t \t norm_value: %f \n", parametric_value);
  // not sure if the following work/is correct


  saving_filename = "R0_";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::SaveMatrixXr(saving_filename,regression_.getR0());

  Lf.resize(_solution(lambda_index,0).size()/2);
  for(UInt i=0; i< Lf.size(); i++){
    Lf(i) = _solution(lambda_index,0)(Lf.size() + i);
  }

  saving_filename = "Lf_";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::saveVectorXr(saving_filename,Lf);


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
  regression_.computeDegreesOfFreedom(0, 0, this->inputData_.getGlobalLambda()[lambda_index], 0);
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
    bestLambdaS_ = inputData_.getGlobalLambda()[lambda_index];
    _bestGCV = GCV_value;
  }

}

template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
Real FPIRLS_Probit<InputHandler,Integrator,ORDER, mydim, ndim>::erf_inv(const Real& x)const{
   Real tt1, tt2, lnx, sgn;
   sgn = (x < 0) ? -1.0 : 1.0;
   lnx = log((1 - x)*(1 + x));  // ln(x) = ln(1 - x*x);
   tt1 = 2/(M_PI*0.147) + 0.5 * lnx;
   tt2 = 1/(0.147) * lnx;

   return(sgn*sqrt(-tt1 + sqrt(tt1*tt1 - tt2)));
}


/*
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Bernoulli<InputHandler,Integrator,ORDER, mydim, ndim>::initialize_mu(const VectorXr& y){
      this->mu_.resize(y.size());
      for(UInt i = 0; i < y.size(); i++){
          this->mu_[i] = 0.5 * (y[i] + 0.5); // It is different for binary or non-binary outcomes
      }
}
*/




#endif

#ifndef __FPIRLS_IMP_H__
#define __FPIRLS_IMP_H__

#include "FPIRLS.h"


/*********** FPIRLS_base Methods ************/

// Constructor
template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::FPIRLS_Base(const MeshHandler<ORDER,mydim,ndim> & mesh, InputHandler & inputData, OptimizationData & optimizationData,  VectorXr mu0, bool scale_parameter_flag, Real scale_param):
  mesh_(mesh), inputData_(inputData), optimizationData_(optimizationData), regression_(inputData, optimizationData, mesh.num_nodes()), scale_parameter_flag_(scale_parameter_flag), _scale_param(scale_param), lenS_(optimizationData.get_size_S()), lenT_(optimizationData.get_size_T())
{
  //Pre-allocate memory for all quatities
  mu_.resize(lenS_, std::vector<VectorXr>(lenT_));
  pseudoObservations_.resize(lenS_, std::vector<VectorXr>(lenT_));
  G_.resize(lenS_, std::vector<VectorXr>(lenT_));
  WeightsMatrix_.resize(lenS_, std::vector<VectorXr>(lenT_));
  current_J_values.resize(lenS_, std::vector<std::array<Real, 2>>(lenT_));
  past_J_values.resize(lenS_, std::vector<std::array<Real, 2>>(lenT_));
  n_iterations.resize(lenS_, std::vector<UInt>(lenT_));
  _J_minima.resize(lenS_, std::vector<Real>(lenT_));
  _GCV.resize(lenS_, std::vector<Real>(lenT_, -1));
  
  //initialization of mu, current_J_values and past_J_values.
  for(UInt i=0; i<optimizationData_.get_size_S() ; i++){
   for(UInt j=0; j<optimizationData_.get_size_T() ; j++){
    mu_[i][j] = mu0;
    current_J_values[i][j] = std::array<Real,2>{1,1};
    past_J_values[i][j] = std::array<Real,2>{1,1};
   }
  }
};

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::FPIRLS_Base(const MeshHandler<ORDER,mydim,ndim> & mesh, const std::vector<Real>& mesh_time, InputHandler & inputData, OptimizationData & optimizationData,  VectorXr mu0, bool scale_parameter_flag, Real scale_param):
  mesh_(mesh), mesh_time_(mesh_time), inputData_(inputData), optimizationData_(optimizationData), regression_(mesh_time, inputData, optimizationData, mesh.num_nodes()), scale_parameter_flag_(scale_parameter_flag), _scale_param(scale_param), lenS_(optimizationData.get_size_S()), lenT_(optimizationData.get_size_T())
{
  //Pre-allocate memory for all quatities
  mu_.resize(lenS_, std::vector<VectorXr>(lenT_));
  pseudoObservations_.resize(lenS_, std::vector<VectorXr>(lenT_));
  G_.resize(lenS_, std::vector<VectorXr>(lenT_));
  WeightsMatrix_.resize(lenS_, std::vector<VectorXr>(lenT_));
  current_J_values.resize(lenS_, std::vector<std::array<Real, 2>>(lenT_));
  past_J_values.resize(lenS_, std::vector<std::array<Real, 2>>(lenT_));
  n_iterations.resize(lenS_, std::vector<UInt>(lenT_));
  _J_minima.resize(lenS_, std::vector<Real>(lenT_));
  _GCV.resize(lenS_, std::vector<Real>(lenT_, -1));
  
  //initialization of mu, current_J_values and past_J_values.
  for(UInt i=0; i<optimizationData_.get_size_S() ; i++){
   for(UInt j=0; j<optimizationData_.get_size_T() ; j++){
    mu_[i][j] = mu0;
    current_J_values[i][j] = std::array<Real,2>{1,1};
    past_J_values[i][j] = std::array<Real,2>{1,1};
   }
  }
};

// FPIRLS_base methods
template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::apply( const ForcingTerm& u){
  // f-PRILS implementation

  // Initialize the outputs. The temporal dimension is not implemented, for this reason the 2nd dimension is set to 1.
  if( this->inputData_.getCovariates()->rows() > 0 )_beta_hat.resize(lenS_, lenT_);
  _fn_hat.resize(lenS_, lenT_);
  _dof.resize(lenS_, lenT_);
  _solution.resize(lenS_,lenT_);

  if(isSpaceVarying)
  {
    FiniteElement<ORDER, mydim, ndim> fe;
  	Assembler::forcingTerm(mesh_, fe, u, forcingTerm);
  }

  for(UInt i=0 ; i < lenS_ ; i++){//for-cycle for each spatial penalization (lambdaS).
   for(UInt j=0 ; j < lenT_ ; j++){
    current_J_values[i][j][0] = past_J_values[i][j][0] + 2*inputData_.get_treshold();
    current_J_values[i][j][1] = past_J_values[i][j][1] + 2*inputData_.get_treshold();

    this->optimizationData_.setCurrentLambda(i, j); // set right lambda for the current iteration.

    //Rprintf("Start FPIRLS for the lambda number %d \n", i+1);

    // start the iterative method for the lambda index i
    while(stopping_criterion(i, j)){

      // STEP (1)
      compute_G(i, j);
      compute_Weights(i, j);
      compute_pseudoObs(i, j);

      // STEP (2)
      this->inputData_.updatePseudodata(pseudoObservations_[i][j], WeightsMatrix_[i][j]);
      update_solution(i, j);

      // STEP (3)
      compute_mu(i, j);

      // update J
      past_J_values[i][j] = current_J_values[i][j];
      current_J_values[i][j] = compute_J(i, j);

      if( !regression_.isMatrixNoCov_factorized() ) {

         Rprintf("WARNING: System matrix cannot be factorized for optimization parameters in position %d (Space) and  %d (Time). Try increasing optimization parameter.\n", i+1, j+1) ;
          break;
      }
      n_iterations[i][j]++;

    } //end while

    //Rprintf("\t n. iterations: %d\n \n", n_iterations[i]);

    _J_minima[i][j] = current_J_values[i][j][0]+current_J_values[i][j][1]; // compute the minimum value of the J fuctional

    if(this->optimizationData_.get_loss_function()=="GCV"){ // compute GCV if it is required

        if( !regression_.isMatrixNoCov_factorized() ){

            _GCV[i][j] = std::numeric_limits<double>::quiet_NaN();

        }else{

      	   compute_GCV(i,j);

        }

    }
  
  } // end time for
 }// end space for

  // Variance Estimate
  compute_variance_est();
}

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::update_solution(UInt& lambdaS_index, UInt& lambdaT_index){
  // performs step (2) of PIRLS. It requires pseudo data after step(1) and mimic regression skeleton behaviour

  // Here we have to solve a weighted regression problem.
  regression_.recomputeWTW(); // at each iteration of FPIRLS W is updated, so WTW has to be recomputed as well.
  regression_.preapply(this->mesh_);
  regression_.apply();

  // if the system matrix is correctly factorized
  if( regression_.isMatrixNoCov_factorized() ) {
  	const SpMat * Psi = regression_.getpsi_(); // get Psi matrix. It is used for the computation of fn_hat.

  	// get the solutions from the regression object.
  	_solution(lambdaS_index, lambdaT_index) = regression_.getSolution()(0,0);
	_dof(lambdaS_index, lambdaT_index) = regression_.getDOF()(0,0);

  	if(inputData_.getCovariates()->rows()>0){
    	_beta_hat(lambdaS_index, lambdaT_index) = regression_.getBeta()(0,0);
  	}
	_fn_hat(lambdaS_index, lambdaT_index) = (*Psi) *_solution(lambdaS_index,lambdaT_index).topRows(Psi->cols());
  }

}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::compute_pseudoObs(UInt& lambdaS_index, UInt& lambdaT_index){
  // compute pseudodata observations

  VectorXr first_addendum; // G_ii( z_i - mu_i)
  VectorXr g_mu; // g( mu_i )

  const VectorXr * z = inputData_.getInitialObservations();

  first_addendum.resize(mu_[lambdaS_index][lambdaT_index].size());
  g_mu.resize(mu_[lambdaS_index][lambdaT_index].size());

  //compute the vecotr first_addendum and g_mu
  for(auto i=0; i < mu_[lambdaS_index][lambdaT_index].size(); i++){
    g_mu(i) = link(mu_[lambdaS_index][lambdaT_index](i));
    first_addendum(i) = G_[lambdaS_index][lambdaT_index](i)*((*z)(i)-mu_[lambdaS_index][lambdaT_index](i));
  }

  pseudoObservations_[lambdaS_index][lambdaT_index] = first_addendum + g_mu;

}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::compute_G(UInt& lambdaS_index, UInt& lambdaT_index){
  // compute the G matrix as G_ii = diag( g'(mu_i))

  G_[lambdaS_index][lambdaT_index].resize(mu_[lambdaS_index][lambdaT_index].size());

  for(UInt i = 0; i<mu_[lambdaS_index][lambdaT_index].size(); i++){
    G_[lambdaS_index][lambdaT_index](i) = link_deriv(mu_[lambdaS_index][lambdaT_index](i));
  }

}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::compute_Weights(UInt& lambdaS_index, UInt& lambdaT_index){
  // computed W elementwise (it is a diagonal matrix)

  WeightsMatrix_[lambdaS_index][lambdaT_index].resize( mu_[lambdaS_index][lambdaT_index].size());

  for(auto i=0; i < mu_[lambdaS_index][lambdaT_index].size(); i++){
    WeightsMatrix_[lambdaS_index][lambdaT_index](i) = 1/(pow(G_[lambdaS_index][lambdaT_index](i),2)*(var_function( mu_[lambdaS_index][lambdaT_index](i))));
  }

}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::compute_mu(UInt& lambdaS_index, UInt& lambdaT_index){
  //compute mu as mu_i = g-1( w_ii*beta + fn_hat)

  VectorXr W_beta = VectorXr::Zero(mu_[lambdaS_index][lambdaT_index].size()); // initialize the vector w_ii*beta

  if(inputData_.getCovariates()->rows()>0)
    W_beta = (*(inputData_.getCovariates()))*_beta_hat(lambdaS_index,lambdaT_index);

  for(UInt j=0; j < W_beta.size(); j++){
      mu_[lambdaS_index][lambdaT_index](j) = inv_link(W_beta[j] + _fn_hat(lambdaS_index,lambdaT_index)(j));
  }

}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
bool FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::stopping_criterion(UInt& lambdaS_index, UInt& lambdaT_index){
  // return true if the f-PIRLS has to perform another iteration, false if it has to be stopped

  bool do_stop_by_iteration = false;  // Do I need to stop becouse n_it > n_max?
  bool do_stop_by_treshold = false; // Do I need to stop becouse |J(k) - J(k+1)| < treshold?

  if(n_iterations[lambdaS_index][lambdaT_index] > inputData_.get_maxiter()){
    do_stop_by_iteration = true;
  }

  if(n_iterations[lambdaS_index][lambdaT_index] > 1){
    if(abs(past_J_values[lambdaS_index][lambdaT_index][0]+past_J_values[lambdaS_index][lambdaT_index][1] - current_J_values[lambdaS_index][lambdaT_index][0] - current_J_values[lambdaS_index][lambdaT_index][1]) < inputData_.get_treshold()){
        do_stop_by_treshold = true;
   }
  }

  return !(do_stop_by_iteration || do_stop_by_treshold );
}

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
std::array<Real,2> FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::compute_J(UInt& lambdaS_index, UInt& lambdaT_index){
  // compute the functional J: it is divided in parametric and non parametric part
  Real parametric_value = 0;
  Real non_parametric_value = 0;
  Real tmp;


  VectorXr Lf;

  const VectorXr * z = inputData_.getInitialObservations();

  for(UInt i=0; i < mu_.size(); i++){
    tmp =std::sqrt( var_function( mu_[lambdaS_index][lambdaT_index](i)) ) * ((*z)(i) - mu_[lambdaS_index][lambdaT_index](i)) ;
    parametric_value += tmp*tmp;
  }

  Lf.resize(_solution(lambdaS_index,lambdaT_index).size()/2);
  for(UInt i=0; i< Lf.size(); i++){
    Lf(i) = _solution(lambdaS_index, lambdaT_index)(Lf.size() + i);
  }

  if(isSpaceVarying)
  {
      Lf = Lf - forcingTerm;
  }

  non_parametric_value = Lf.transpose() * (*(regression_.getR0_())) * Lf;
  non_parametric_value = (*optimizationData_.get_LambdaS_vector())[lambdaS_index]*non_parametric_value;

  std::array<Real,2> returnObject{parametric_value, non_parametric_value};

  return returnObject;
}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::compute_GCV(UInt& lambdaS_index, UInt& lambdaT_index){

        if (optimizationData_.get_DOF_evaluation() != "not_required") //in this case surely we have already the dofs
        { // is DOF_matrix to be computed?
        regression_.computeDegreesOfFreedom(0, 0, (*optimizationData_.get_LambdaS_vector())[lambdaS_index],
        					  (*optimizationData_.get_LambdaT_vector())[lambdaT_index]);
        _dof(lambdaS_index, lambdaT_index) = regression_.getDOF()(0,0);
        }
        else _dof(lambdaS_index, lambdaT_index) = regression_.getDOF()(lambdaS_index, lambdaT_index);
 
        const VectorXr * y = inputData_.getInitialObservations();
        Real GCV_value = 0;

        for(UInt j=0; j < y->size();j++)
        GCV_value += dev_function(mu_[lambdaS_index][lambdaT_index][j], (*y)[j]); //norm computation

        GCV_value *= y->size();

        GCV_value /= (y->size()-optimizationData_.get_tuning()*_dof(lambdaS_index,lambdaT_index))*(y->size()-optimizationData_.get_tuning()*_dof(lambdaS_index,lambdaT_index));

        _GCV[lambdaS_index][lambdaT_index] = GCV_value;

        //best lambda
        if(GCV_value < optimizationData_.get_best_value())
        {
        optimizationData_.set_best_lambda_S(lambdaS_index);
        optimizationData_.set_best_lambda_T(lambdaT_index);
        optimizationData_.set_best_value(GCV_value);
        }

}

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::compute_variance_est(){
  Real phi;
  if(this->scale_parameter_flag_ && this->optimizationData_.get_loss_function()!="GCV"){// if scale param should be
    _variance_estimates.resize(lenS_, std::vector<Real>(lenT_, 0.0));
    const UInt n_obs = this->inputData_.getObservations()->size();

    //scale parameter computed as: mean((var.link(mu)*phi)/mu), and phi is computed as in Wood IGAM
    for(UInt i=0; i < lenS_;i++){
    	for(UInt j=0; j< lenT_; j++){ 
		phi = (this->scale_parameter_flag_ )? this->current_J_values[i][j][0]/(n_obs - this->_dof(i,j) ) : _scale_param;
		for(UInt k=0; k < this->mu_[i][j].size(); k++){
        	_variance_estimates[i][j] += phi* this->var_function(this->mu_[i][j](k))/this->mu_[i][j](k);
      }
      _variance_estimates[i][j] /= this->mu_[i][j].size();
    }
   }
  }else{
    _variance_estimates.resize(lenS_, std::vector<Real>(lenT_,-1));
  }
}

/*********** FPIRLS apply template specialization ************/

// Laplace or Elliptic case
template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS<InputHandler,ORDER, mydim, ndim>::apply(){

  FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::apply(ForcingTerm());

}

// SpaceVarying case
template <UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS<GAMDataEllipticSpaceVarying,ORDER, mydim, ndim>::apply(){

  this->isSpaceVarying = true;
  FPIRLS_Base<GAMDataEllipticSpaceVarying,ORDER, mydim, ndim>::apply(this->inputData_.getU());

}


#endif

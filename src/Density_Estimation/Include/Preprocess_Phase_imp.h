#ifndef __PREPROCESS_PHASE_IMP_H__
#define __PREPROCESS_PHASE_IMP_H__


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
Preprocess<Integrator_noPoly, ORDER, mydim, ndim>::Preprocess(const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dp,
  const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& fp):
  dataProblem_(dp), funcProblem_(fp){

    densityInit_ = DensityInitialization_factory<Integrator_noPoly, ORDER, mydim, ndim>::createInitializationSolver(dp, fp);

    fInit_.resize(dp.getNlambda());
    fillFInit();

  };


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
void
Preprocess<Integrator_noPoly, ORDER, mydim, ndim>::fillFInit(){

  for(UInt l = 0; l < dataProblem_.getNlambda(); l++){
    fInit_[l] = densityInit_-> chooseInitialization(dataProblem_.getLambda(l));
  }

}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
void
NoCrossValidation<Integrator_noPoly, ORDER, mydim, ndim>::performPreprocessTask(){

  this->bestLambda_ = this->dataProblem_.getLambda(0);
  this->gInit_ = (*(this->fInit_[0])).array().log();

}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
CrossValidation<Integrator_noPoly, ORDER, mydim, ndim>::CrossValidation(const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dp,
  const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& fp,
  std::shared_ptr<MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>> ma):
  Preprocess<Integrator_noPoly, ORDER, mydim, ndim>(dp, fp), minAlgo_(ma), error_(dp){

    K_folds_.resize(dp.getNumberofData());
    CV_errors_.resize(dp.getNlambda(), 0);
    g_sols_.resize(dp.getNlambda());

}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
std::pair<VectorXr, Real>
CrossValidation<Integrator_noPoly, ORDER, mydim, ndim>::performCV(){

  UInt N = this->dataProblem_.getNumberofData();
  UInt K = this->dataProblem_.getNfolds();

  for (UInt i = 0; i< N; i++){
    UInt length = ((i % K) <= (N % K))? (i % K)*(N/K +1) : (N % K) + (i % K)*(N/K);
    K_folds_[length + i/K] = i;
  }

  // cycle on the folds
  for (UInt i = 0; i < K; i++){

    if(this->dataProblem_.Print()){
      Rprintf("X_valid is the fold number %d\n", i);
    }

    std::vector<UInt> x_valid, x_train;

    if (i < N % K){ // fold grossi
      std::set_union(K_folds_.cbegin(), K_folds_.cbegin()+ i*(N/K +1), K_folds_.cbegin()+ (i + 1)*(N/K +1), K_folds_.cend(), std::back_inserter(x_train));
      std::copy(K_folds_.cbegin()+ i*(N/K +1), K_folds_.cbegin()+ (i + 1)*(N/K +1), std::back_inserter(x_valid));
    }
    else{ //fold piccoli
      std::set_union(K_folds_.cbegin(), K_folds_.cbegin()+ (N % K) + i*(N/K), K_folds_.cbegin()+ (N % K) + (i+1)*(N/K) , K_folds_.cend(), std::back_inserter(x_train));
      std::copy(K_folds_.cbegin()+ (N % K) + i*(N/K), K_folds_.cbegin()+ (N % K) + (i+1)*(N/K), std::back_inserter(x_valid));
    }

    SpMat Psi_train = this->dataProblem_.computePsi(x_train);
    SpMat Psi_valid = this->dataProblem_.computePsi(x_valid);

    performCV_core(i, Psi_train, Psi_valid); // it fills g_sols, CV_errors_

  }

  UInt init_best_lambda = std::distance(CV_errors_.cbegin(), std::min_element(CV_errors_.cbegin(), CV_errors_.cend()));

  return std::pair<VectorXr, Real> (g_sols_[init_best_lambda], this->dataProblem_.getLambda(init_best_lambda));

}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
void
CrossValidation<Integrator_noPoly, ORDER, mydim, ndim>::performPreprocessTask(){

  std::tie(this->gInit_, this->bestLambda_) = performCV();

}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
void
SimplifiedCrossValidation<Integrator_noPoly, ORDER, mydim, ndim>::performCV_core(UInt fold, const SpMat& Psi_train, const SpMat& Psi_valid){

  if(this->dataProblem_.Print()){
    Rprintf("lambda: %f\n", this->dataProblem_.getLambda(fold));
  }

  this->g_sols_[fold] = this->minAlgo_->apply_core(Psi_train, this->dataProblem_.getLambda(fold), (*(this->fInit_[fold])).array().log());

  this->CV_errors_[fold] = this->error_(Psi_valid, this->g_sols_[fold]);

}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
RightCrossValidation<Integrator_noPoly, ORDER, mydim, ndim>::
  RightCrossValidation(const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dp,
   const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& fp,
   std::shared_ptr<MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>> ma):
   CrossValidation<Integrator_noPoly, ORDER, mydim, ndim>(dp, fp, ma){

      best_loss_.resize(this->dataProblem_.getNlambda(), std::numeric_limits<double>::max());

}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
void
RightCrossValidation<Integrator_noPoly, ORDER, mydim, ndim>::performCV_core(UInt fold, const SpMat& Psi_train, const SpMat& Psi_valid){

  for (UInt l=0; l < this->dataProblem_.getNlambda(); l++){

     std::unique_ptr<MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>>
     minimizationAlgo = this->minAlgo_->clone();

     if(this->dataProblem_.Print()){
       Rprintf("lambda: %f\n", this->dataProblem_.getLambda(l));
     }

     VectorXr sols;

     sols = minimizationAlgo->apply_core(Psi_train, this->dataProblem_.getLambda(l), (*(this->fInit_[l])).array().log());

     this->CV_errors_[l] += this->error_(Psi_valid, sols);

     Real loss = std::get<0>(this->funcProblem_.computeFunctional_g(sols, this->dataProblem_.getLambda(l), Psi_train));

     if(loss < best_loss_[l]){
       best_loss_[l] = loss;
       this->g_sols_[l] = sols;
     }
  }

}


#endif

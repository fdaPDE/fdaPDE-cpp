#ifndef __FE_DENSITY_ESTIMATION_IMP_H__
#define __FE_DENSITY_ESTIMATION_IMP_H__


template<typename Integrator, typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
FEDE<Integrator, Integrator_noPoly, ORDER, mydim, ndim>::
  FEDE(const DataProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim>& dp,
    const FunctionalProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim>& fp,
    std::shared_ptr<MinimizationAlgorithm<Integrator, Integrator_noPoly, ORDER, mydim, ndim>> ma, const std::string& p):
      dataProblem_(dp), funcProblem_(fp), minAlgo_(ma){

        preprocess_ = Preprocess_factory<Integrator, Integrator_noPoly, ORDER, mydim, ndim>::createPreprocessSolver(dp, fp, ma, p);

}


template<typename Integrator, typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
void
FEDE<Integrator, Integrator_noPoly, ORDER, mydim, ndim>::apply(){

  // perform the preprocess phase
    Rprintf("##### PREPROCESS PHASE #####\n");
  preprocess_ -> performPreprocessTask();

  // collect preprocess results
  VectorXr gInit;
  std::tie(fInit_, gInit, bestLambda_) = preprocess_ -> getPreprocessParameter();

  CV_errors_ = preprocess_ -> getCvError();

  // final minimization descent
    Rprintf("##### FINAL STEP #####\n");

  gcoeff_ = minAlgo_->apply_core(dataProblem_.getGlobalPsi(), bestLambda_, gInit);

}

#endif

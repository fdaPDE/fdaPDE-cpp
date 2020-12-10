#ifndef __FE_DENSITY_ESTIMATION_IMP_H__
#define __FE_DENSITY_ESTIMATION_IMP_H__


template<UInt ORDER, UInt mydim, UInt ndim>
FEDE<ORDER, mydim, ndim>::
  FEDE(const DataProblem<ORDER, mydim, ndim>& dp,
    const FunctionalProblem<ORDER, mydim, ndim>& fp,
    std::shared_ptr<MinimizationAlgorithm<ORDER, mydim, ndim>> ma, const std::string& p):
      dataProblem_(dp), funcProblem_(fp), minAlgo_(ma){

        preprocess_ = Preprocess_factory<ORDER, mydim, ndim>::createPreprocessSolver(dp, fp, ma, p);

}


template<UInt ORDER, UInt mydim, UInt ndim>
void
FEDE<ORDER, mydim, ndim>::apply(){

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

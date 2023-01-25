#ifndef __STOCHASTIC_EDF_H__
#define __STOCHASTIC_EDF_H__

#include <random>
#include "../core/utils/Symbols.h"
#include "../core/NLA/SMW.h"
using fdaPDE::core::NLA::SMW;
#include "../models/regression/RegressionBase.h"
using fdaPDE::models::is_regression_model;

namespace fdaPDE {
namespace calibration{

  // computes an approximation of the trace of S = \Psi*T^{-1}*\Psi^T*Q using a monte carlo approximation.
  template <typename Model>
  class StochasticEDF {
    static_assert(is_regression_model<Model>::value);
  private:
    const Model& model_;
    std::size_t r_;      // number of monte carlo realizations
    std::size_t seed_;
    DMatrix<double> Us_; // sample from Rademacher distribution
    DMatrix<double> Bs_; // \Psi^T*Q*Us_
    DMatrix<double> Y_;  // Us_^T*\Psi 
  public:
    // constructor
    StochasticEDF(const Model& model, std::size_t r, std::size_t seed)
      : model_(model), r_(r), seed_(seed) {
      // compute sample from Rademacher distribution
      std::default_random_engine rng(seed_);
      std::bernoulli_distribution Be(0.5); // bernulli distribution with parameter p = 0.5
      Us_.resize(model_.n_obs(), r_); // preallocate memory for matrix Us
      // fill matrix
      for(std::size_t i = 0; i < model_.n_obs(); ++i){
	for(std::size_t j = 0; j < r_; ++j){
	  if(Be(rng)) Us_(i,j) =  1.0;
	  else        Us_(i,j) = -1.0;
	}
      }
      // prepare matrix Bs_
      Bs_ = DMatrix<double>::Zero(2*model_.n_basis(), r_);
      if(!model_.hasCovariates()) // non-parametric model
	Bs_.topRows(model_.n_basis()) = -model_.PsiTD()*model_.W()*Us_;
      else // semi-parametric model
	Bs_.topRows(model_.n_basis()) = -model_.PsiTD()*model_.lmbQ(Us_);
      // prepare matrix Y
      Y_ = Us_.transpose()*model_.Psi();
    }
    StochasticEDF(const Model& model, std::size_t r)
      : StochasticEDF(model, r, std::random_device()()) {}
    
    // evaluate trace of S exploiting a monte carlo approximation
    double compute() {
      DMatrix<double> sol; // room for problem solution
      if(!model_.hasCovariates()){ // nonparametric case
        sol = model_.invA().solve(Bs_);
      }else{
	// solve system (A+UCV)*x = Bs via woodbury decomposition using matrices U and V cached by model_
	sol = SMW<>().solve(model_.invA(), model_.U(), model_.XtWX(), model_.V(), Bs_);
      }
      // compute approximated Tr[S] using monte carlo mean
      double MCmean = 0;
      for(std::size_t i = 0; i < r_; ++i)
	MCmean += Y_.row(i).dot(sol.col(i).head(model_.n_basis()));
      
      return MCmean/r_;
    }
  };

}}
  
#endif // __STOCHASTIC_EDF_H__

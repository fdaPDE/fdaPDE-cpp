#ifndef __STOCHASTIC_GCV_ENGINE_H__
#define __STOCHASTIC_GCV_ENGINE_H__

#include <memory>
#include <random>
#include <functional>

#include "../core/utils/Symbols.h"
#include "../core/utils/fields/ScalarField.h"
using fdaPDE::core::ScalarField;
#include "../core/NLA/SMW.h"
using fdaPDE::core::NLA::SMW;

namespace fdaPDE{
namespace calibration{

  // computes an approximation of the trace of S = \Psi*T^{-1}*\Psi^T*Q using a monte carlo approximation.
  class StochasticGCVEngine {
  private:
    std::size_t r_; // number of realizations

    // create a sample from a Rademacher distribution to be used for the approximation of Tr[S] as an n x r_ matrix
    std::shared_ptr<DMatrix<double>> Us(std::size_t n) const {
      // define bernulli distribution with parameter p = 0.5 to simulate Rademacher distribution
      std::random_device rng;
      std::bernoulli_distribution Be(0.5);
      // preallocate memory for matrix Us
      std::shared_ptr<DMatrix<double>> Us = std::make_shared<DMatrix<double>>(n, r_);
      // fill matrix
      for(std::size_t i = 0; i < n; ++i){
	for(std::size_t j = 0; j < r_; ++j){
	  if(Be(rng)) Us->coeffRef(i,j) =  1.0;
	  else        Us->coeffRef(i,j) = -1.0;
	}
      }
      return Us;
    }
  
  public:
    // constructor
    StochasticGCVEngine(std::size_t r) : r_(r) { }

    // evaluate trace of S exploiting a monte carlo approximation
    template<typename M>
    double compute(const M& model) const {
      std::size_t n = model.n(); // number of observations
      // apply SMW decomposition
      std::size_t q_ = model.W()->cols(); // number of covariates
      // definition of matrices U and V
      DMatrix<double> U = DMatrix<double>::Zero(model.A()->rows(), q_);
      U.block(0,0, model.A()->rows()/2, q_) = model.Psi()->transpose()*(*model.W());
  
      DMatrix<double> V = DMatrix<double>::Zero(q_, model.A()->rows());
      V.block(0,0, q_, model.A()->rows()/2) = model.W()->transpose()*(*model.Psi());

      // Define system solver. Use SMW solver from NLA module
      SMW<> solver{};
      solver.compute(*model.A());
    
      // define matrix Bs
      DMatrix<double> Bs = DMatrix<double>::Zero(2*model.obs(), r_);
      std::shared_ptr<DMatrix<double>> L = Us(n);
      Bs.topRows(model.obs()) = model.Psi()->transpose() * model.lmbQ(*L);
    
      // solve system Mx = b
      DMatrix<double> solution;
      solution = solver.solve(U, *(model.WTW()), V, Bs);

      // compute approximated Tr[S] using monte carlo mean
      DMatrix<double> Y = L->transpose() * (*model.Psi());
      double MCmean = 0;
      for(std::size_t i = 0; i < r_; ++i){
	MCmean += Y.row(i).dot(solution.col(i).head(n));
      }
      return MCmean/r_;
    }
  };
}}
#endif // __STOCHASTIC_GCV_ENGINE_H__

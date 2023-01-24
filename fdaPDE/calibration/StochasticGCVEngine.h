#ifndef __STOCHASTIC_GCV_ENGINE_H__
#define __STOCHASTIC_GCV_ENGINE_H__

#include <random>
#include "../core/utils/Symbols.h"
#include "../core/NLA/SMW.h"
using fdaPDE::core::NLA::SMW;

namespace fdaPDE{
namespace calibration{

  // computes an approximation of the trace of S = \Psi*T^{-1}*\Psi^T*Q using a monte carlo approximation.
  class StochasticGCVEngine {
  private:
    std::size_t r_; // number of realizations
    DMatrix<double> Us_; // sample from Rademacher distribution

    std::size_t seed_;
  public:
    // constructor
    StochasticGCVEngine(std::size_t r) : r_(r), seed_(std::random_device()()) {}
    StochasticGCVEngine(std::size_t r, std::size_t seed) : r_(r), seed_(seed) {}
    
    // evaluate trace of S exploiting a monte carlo approximation
    template<typename M>
    double compute(const M& model) {
      std::size_t n = model.n_obs(); // number of observations
      std::size_t q = model.q();     // number of covariates
      if(Us_.size() == 0){
	// compute sample from Rademacher distribution
	std::default_random_engine rng(seed_);
	std::bernoulli_distribution Be(0.5); // bernulli distribution with parameter p = 0.5
	Us_.resize(n, r_); // preallocate memory for matrix Us
	// fill matrix
	for(std::size_t i = 0; i < n; ++i){
	  for(std::size_t j = 0; j < r_; ++j){
	    if(Be(rng)) Us_(i,j) =  1.0;
	    else        Us_(i,j) = -1.0;
	  }
	}
      }
      DMatrix<double> sol; // room for problem solution
      if(!model.hasCovariates()){ // nonparametric case
	// define system solver. Use a sparse solver
	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver{};
	solver.compute(model.A());
	// define matrix Bs
	DMatrix<double> Bs = DMatrix<double>::Zero(2*model.n_basis(), r_);
	Bs.topRows(model.n_basis()) = -model.PsiTD()*model.W()*Us_;
	// solve system
        sol = solver.solve(Bs);
      }else{
	// define matrix Bs
        DMatrix<double> Bs = DMatrix<double>::Zero(2*model.n_basis(), r_);
	Bs.topRows(model.n_basis()) = -model.PsiTD()*model.lmbQ(Us_);
	// definition of matrices U and V  for application of woodbury formula
	DMatrix<double> U = DMatrix<double>::Zero(model.A().rows(), q);
	U.block(0,0, model.A().rows()/2, q) = model.PsiTD()*model.W()*model.X();
	DMatrix<double> V = DMatrix<double>::Zero(q, model.A().rows());
	V.block(0,0, q, model.A().rows()/2) = model.X().transpose()*model.W()*model.Psi();
	// Define system solver. Use SMW solver from NLA module
	SMW<> solver{};
	solver.compute(model.A());
	// solve system (A+UCV)*x = Bs
	sol = solver.solve(U, model.XtWX(), V, Bs);
      }
      // compute approximated Tr[S] using monte carlo mean
      DMatrix<double> Y = Us_.transpose()*model.Psi();
      double MCmean = 0;
      for(std::size_t i = 0; i < r_; ++i){
	MCmean += Y.row(i).dot(sol.col(i).head(model.n_basis()));
      }
      return MCmean/r_;
    }
  };
  
}}
#endif // __STOCHASTIC_GCV_ENGINE_H__

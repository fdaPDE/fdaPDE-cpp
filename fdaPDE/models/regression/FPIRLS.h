#ifndef __FPIRLS_H__
#define __FPIRLS_H__

#include "../../core/utils/Symbols.h"
#include "../../core/NLA/SparseBlockMatrix.h"
using fdaPDE::core::NLA::SparseBlockMatrix;
#include "../../core/NLA/SMW.h"
using fdaPDE::core::NLA::SMW;
#include "Distributions.h"
#include <cstddef>
#include "../ModelTraits.h"
#include "SRPDE.h"
#include "STRPDE.h"

#include <chrono>

namespace fdaPDE{
namespace models{

  // trait to select model type to use in the internal loop of FPIRLS
  template <typename Model>
  struct FPRILS_internal_solver {
    using type = typename std::conditional<
      !is_space_time<Model>::value,
      // space-only problem
      SRPDE <typename model_traits<Model>::PDE, model_traits<Model>::sampling>,
      // space-time problem
      STRPDE<typename model_traits<Model>::PDE, typename model_traits<Model>::RegularizationType,
	     model_traits<Model>::sampling, model_traits<Model>::solver>
      >::type;
  };
  
  // a general implementation of the Functional Penalized Iterative Reweighted Least Square (FPIRLS) algorithm
  template <typename Distribution>
  class FPIRLS {
  private:
    // data characterizing the behaviour of the algorithm
    Distribution distribution_{};
    
    double tolerance_;
    std::size_t max_iter_;
    
    // let g() the link function of the considered distribution and y = (y_1, y_2, ..., y_n) the (1 x n) vector of observations
    DVector<double> mu_{};     // \mu^k = [ \mu^k_1, ..., \mu^k_n ] : mean vector at step k
    DVector<double> theta_{};  // \theta^k = [ g(\mu^k_1), ..., g(\mu^k_n) ]
    DVector<double> G_{};      // G^k = diag(g'(\mu^k_1), ..., g'(\mu^k_n))
    DVector<double> py_{};     // \tilde z^k = G^k(y-u^k) + \theta^k : pseudo-observations vector at step k
    DVector<double> V_{};      // V^k = diag(v(\mu^k_1), ..., v(\mu^k_n)) : variance matrix at step k
    DVector<double> W_{};      // W^k = ((G^k)^{-2})*((V^k)^{-1})
    
    DVector<double> f_{};
    DVector<double> g_{};
    DVector<double> beta_{};
    
  public:
    // constructor
    FPIRLS(double tolerance, std::size_t max_iter)
      : tolerance_(tolerance), max_iter_(max_iter) {};
    
    // in general observations are not known at construction time
    template <typename M>
    void compute(M& m_) {
      static_assert(is_regression_model<M>::value);
      // get number of data and preallocate space
      std::size_t n = m_.n_obs();
      theta_.resize(n); G_.resize(n); py_.resize(n); V_.resize(n); W_.resize(n);

      // algorithm initialization
      mu_ = m_.y();
      distribution_.preprocess(mu_);
      std::size_t k = 0; // FPIRLS iteration index
      // define internal problem solver and initialize it
      typename FPRILS_internal_solver<M>::type solver(m_.pde(), m_.locs());
      solver.setLambda(m_.lambda());
      solver.init();
      
      // prepare data for solver, copy covariates if present
      BlockFrame<double, int> df = m_.data();
      if(m_.hasCovariates()) df.insert<double>(DESIGN_MATRIX_BLK, m_.X());
      
      // algorithm stops when an enought small difference between two consecutive values of the J is recordered
      double J_old = tolerance_+1; double J_new = 0;
      // start loop
      while(k < max_iter_ && std::abs(J_new - J_old) > tolerance_){
        theta_ = distribution_.link(mu_);
	G_ = distribution_.der_link(mu_);
	V_ = distribution_.variance(mu_);
	W_ = ((G_.array().pow(2)*V_.array()).inverse()).matrix();
	// compute pseudo observations
	py_ = G_.asDiagonal()*(m_.z() - mu_) + theta_;

	// solve weighted least square problem
	// \argmin_{\beta, f} [ \norm(W^{1/2}(y - X\beta - f_n))^2 + \lambda \int_D (Lf - u)^2 ]
	df.insert<double>(OBSERVATIONS_BLK, py_); // insert should overwrite existing block, if any is present
	df.insert<double>(WEIGHTS_BLK, W_);
	solver.setData(df);
	solver.solve();
	
	// extract estimates from solver
	f_ = solver.f(); g_ = solver.g();
	if(m_.hasCovariates()) beta_ = solver.beta();
	
	// update value of \mu_
	DVector<double> fitted = solver.fitted(); // compute fitted values
	for(std::size_t i = 0; i < n; ++i)
	  mu_[i] = distribution_.inv_link(fitted[i]);
	
	// compute value of functional J for this pair (\beta, f): \norm{V^{-1/2}(y - \mu)}^2 + \int_D (Lf-u)^2
	DVector<double> V = distribution_.variance(mu_).array().sqrt().inverse().matrix();
	double J = (V.asDiagonal()*(m_.y() - mu_)).squaredNorm() + m_.lambda()*g_.dot(m_.R0()*g_); // \int_D (Lf-u)^2

	// prepare for next iteration
	k++;
	J_old = J_new; J_new = J;
      }
      return;
    }

    // getters
    DVector<double> weights() const { return W_; }
    DVector<double> beta() const { return beta_; }
    DVector<double> f() const { return f_; }
  };
  
}}



#endif // __FPIRLS_H__

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
#include "../ModelMacros.h"
#include "SRPDE.h"
#include "STRPDE.h"

#include <chrono>

namespace fdaPDE{
namespace models{

  // trait to select model type to use in the internal loop of FPIRLS
  template <typename Model>
  class FPIRLS_internal_solver {
  private:
    typedef typename std::decay<Model>::type Model_;
    typedef typename model_traits<Model_>::PDE            PDE;
    typedef typename model_traits<Model_>::sampling       sampling;
    typedef typename model_traits<Model_>::solver         solver;
    typedef typename model_traits<Model_>::regularization regularization;
  public:
    using type = typename std::conditional<
      !is_space_time<Model_>::value,
      SRPDE <PDE, sampling>, // space-only problem
      STRPDE<PDE, regularization, sampling, solver> // space-time problem
      >::type;
  };
  
  // a general implementation of the Functional Penalized Iterative Reweighted Least Square (FPIRLS) algorithm
  template <typename Model, typename Distribution>
  class FPIRLS {
  private:
    typedef typename std::decay<Model>::type Model_;
    // data characterizing the behaviour of the algorithm
    Distribution distribution_{};
    Model& m_;
    // algorithm's parameters 
    double tolerance_; 
    std::size_t max_iter_;
    std::size_t k_ = 0; // FPIRLS iteration index
    
    DVector<double> mu_{};    // \mu^k = [ \mu^k_1, ..., \mu^k_n ] : mean vector at step k
    // parameters at convergece
    DVector<double> f_{};     // estimate of non-parametric spatial field
    DVector<double> g_{};     // PDE misfit
    DVector<double> beta_{};  // estimate of coefficient vector
    DVector<double> W_{};     // weight matrix
  public:
    // constructor
    FPIRLS(const Model& m, double tolerance, std::size_t max_iter)
      : m_(m), tolerance_(tolerance), max_iter_(max_iter) {};
    
    // executes the FPIRLS algorithm
    void compute() {
      static_assert(is_regression_model<Model>::value);      
      // algorithm initialization
      mu_ = m_.y();
      distribution_.preprocess(mu_);
      // define internal problem solver
      typename FPIRLS_internal_solver<Model>::type solver;
      if constexpr(!is_space_time<Model>::value) // space-only
	solver = typename FPIRLS_internal_solver<Model>::type(m_.pde());
      else{ // space-time
	solver = typename FPIRLS_internal_solver<Model>::type(m_.pde(), m_.time_domain());
	// in case of parabolic regularization derive initial condition from input model
	if constexpr(is_space_time_parabolic<Model_>::value)
	  solver.setInitialCondition(m_.s());
	// in case of separable regularization set possible temporal locations
	if constexpr(is_space_time_separable<Model_>::value)
	  solver.set_temporal_locations(m_.time_locs());
      }
      // solver initialization
      solver.data() = m_.data();
      solver.setLambda(m_.lambda());
      solver.set_spatial_locations(m_.locs());
      solver.init_pde();
      solver.init_regularization();
      solver.init_sampling();
      solver.init_nan();

      // algorithm stops when an enought small difference between two consecutive values of the J is recordered
      double J_old = tolerance_+1; double J_new = 0;
      // start loop
      while(k_ < max_iter_ && std::abs(J_new - J_old) > tolerance_){
	// request weight matrix W and pseudo-observation vector \tilde y from model
	auto pair = m_.compute(mu_);
	
	// solve weighted least square problem
	// \argmin_{\beta, f} [ \norm(W^{1/2}(y - X\beta - f_n))^2 + \lambda \int_D (Lf - u)^2 ]
	solver.data().template insert<double>(OBSERVATIONS_BLK, std::get<1>(pair));
	solver.data().template insert<double>(WEIGHTS_BLK, std::get<0>(pair));
	// update solver to change in the weight matrix
	solver.init_data();
	solver.init_model();
	solver.solve();
	
	// extract estimates from solver
	f_ = solver.f(); g_ = solver.g();
	if(m_.hasCovariates()) beta_ = solver.beta();
	
	// update value of \mu_
	DVector<double> fitted = solver.fitted(); // compute fitted values
	mu_ = distribution_.inv_link(fitted);
	
	// compute value of functional J for this pair (\beta, f): \norm{V^{-1/2}(y - \mu)}^2 + \int_D (Lf-u)^2
	DVector<double> V = distribution_.variance(mu_).array().sqrt().inverse().matrix();
	double J = (V.asDiagonal()*(m_.y() - mu_)).squaredNorm() + g_.dot(m_.R0()*g_); // \int_D (Lf-u)^2
	// prepare for next iteration
	k_++; J_old = J_new; J_new = J;
      }
      // store weight matrix at convergence
      W_ = std::get<0>(m_.compute(mu_));
      return;
    }

    // getters
    const DVector<double>& mu() const { return mu_; } // mean vector at convergence
    const DVector<double>& weights() const { return W_; } // weights matrix W at convergence
    const DVector<double>& beta() const { return beta_; } // estimate of coefficient vector 
    const DVector<double>& f() const { return f_; } // estimate of spatial field 
    const DVector<double>& g() const { return g_; } // PDE misfit
    std::size_t n_iter() const { return k_ - 1; } // number of iterations
  };
  
}}



#endif // __FPIRLS_H__

#ifndef __GCV_H__
#define __GCV_H__

#include <functional>
#include <memory>
#include <type_traits>

#include "../core/utils/Symbols.h"
#include "../core/utils/fields/ScalarField.h"
using fdaPDE::core::ScalarField;
using fdaPDE::core::TwiceDifferentiableScalarField;
#include "../core/NLA/SMW.h"
using fdaPDE::core::NLA::SMW;
#include "ExactEDF.h"
#include "StochasticEDF.h"

// interfaces
#include "iGCV.h"
#include "../models/ModelTraits.h"
using fdaPDE::models::is_regression_model;

namespace fdaPDE{
namespace calibration{

  // base functor implementing the expression of GCV index for model M. Use type T for evaluation of the expected degrees of freedoms
  template <typename M, typename trS_evaluation_strategy = StochasticEDF<M>>
  class GCV {
    // guarantees statistical model M is compatible with GCV computation
    static_assert(std::is_base_of<iGCV, M>::value,
		  "you are asking to calibrate something which cannot be handled by GCV");
  protected:
    // cache computed values, use pointers to let shallow-copy (optimizers need to access those members)
    typedef std::shared_ptr<std::vector<double>> vect_ptr;
    vect_ptr edfs_   = std::make_shared<std::vector<double>>();
    vect_ptr values_ = std::make_shared<std::vector<double>>();

    M& model_; // model to calibrate
    trS_evaluation_strategy trS_; // strategy used to evaluate the trace of smoothing matrix S

    // cache pairs (lambda, Tr[S]) for fast access if GCV is queried at an already computed point
    std::map<SVector<model_traits<M>::n_lambda>, double, fdaPDE::s_vector_compare<model_traits<M>::n_lambda>> cache_;
  public:
    // SFINAE selection of constructor depending on trace evaluation strategy
    template <typename U = trS_evaluation_strategy, // fake type to enable substitution
	      typename std::enable_if<!std::is_same<U, StochasticEDF<M>>::value,int>::type = 0> 
    GCV(M& model) : model_(model), trS_(model_) {};

    // constructor overloads for stochastic trace approximation
    template <typename U = trS_evaluation_strategy,
	      typename std::enable_if<std::is_same<U, StochasticEDF<M>>::value,int>::type = 0>
    GCV(M& model, std::size_t r) : model_(model), trS_(model_, r) {};
    template <typename U = trS_evaluation_strategy,
	      typename std::enable_if<std::is_same<U, StochasticEDF<M>>::value,int>::type = 0>
    GCV(M& model, std::size_t r, std::size_t seed) : model_(model), trS_(model_, r, seed) {};

    // evaluates the analytical expression of gcv at \lambda (called by any type of GCV optimization)
    //
    // edf = n - (q + Tr[S])
    // GCV(\lambda) = n/(edf^2)*norm(y - \hat y)^2
    double operator()(const SVector<model_traits<M>::n_lambda>& lambda) {
      // fit the model given current lambda
      model_.setLambda(lambda);
      model_.init_model();
      model_.solve();
      // compute equivalent degrees of freedom given current lambda (if not already cached)
      if(cache_.find(lambda) == cache_.end()){
	cache_[lambda] = trS_.compute();
      }
      double trS = cache_[lambda];
      double q = model_.q();          // number of covariates
      std::size_t n = model_.n_obs(); // number of observations      
      double dor = n - (q + trS);     // residual degrees of freedom
      edfs_->emplace_back(q + trS);   // store equivalent degrees of freedom
      
      // return gcv at point
      double gcv_value = (n/std::pow(dor, 2))*( model_.norm(model_.fitted(), model_.y()) ) ;
      values_->emplace_back(gcv_value);
      return gcv_value;
    }

    // returns GCV index of Model in its current state (assume Model already solved)
    double eval() {
      // compute equivalent degrees of freedom given current lambda (if not already cached)
      if(cache_.find(model_.lambda()) == cache_.end()){
	cache_[model_.lambda()] = trS_.compute();
      }
      double trS = cache_[model_.lambda()];
      // GCV(\lambda) = n/((n - (q + Tr[S]))^2)*norm(y - \hat y)^2
      double dor = model_.n_obs() - (model_.q() + trS); // (n - (q + Tr[S])
      return (model_.n_obs()/std::pow(dor, 2))*( model_.norm(model_.fitted(), model_.y()) ) ;
    }
    
    // getters
    const std::vector<double>& edfs() const { return *edfs_; } // equivalent degrees of freedom q + Tr[S]
    const std::vector<double>& values() const { return *values_; } // computed values of GCV index
  };

  // for the optimization of GCV using finite differences we just need the definition of the GCV functional. 
  template <typename M, typename trS_evaluation_strategy = StochasticEDF<M>>
  using FiniteDifferenceGCV = GCV<M, trS_evaluation_strategy>;

  // an optimization of GCV using its exact expression requires the analitycal expression of gradient and hessian matrix.
  // this depends on the type of regularization used by the statistical model M
  template <typename M, typename RegularizationType> class ExactGCV;
  
  // space only specialization of GCV exact derivatives
  // expression of GCV derivatives:
  //    edf = n - (q + Tr[S])
  //    dGCV(\lambda)  = \frac{2n}{edf^2}[ \sigma^2 * Tr[dS] + a ]
  //    ddGCV(\lambda) = \frac{2n}{edf^2}[ \frac{1}{edf}(3*\sigma^2*Tr[dS] + 4*a)*Tr[dS] + \sigma^2*Tr[ddS] + b ]    
  template <typename M>
  class ExactGCV<M, fdaPDE::models::SpaceOnlyTag> : public GCV<M, ExactEDF<M>> {
  private:
    // import symbols from base
    typedef GCV<M, ExactEDF<M>> Base;
    using Base::trS_;
    using Base::model_;
    
    DMatrix<double> L_{}; // T^{-1}*R
    DMatrix<double> F_{}; // (T^{-1}*R)*(T^{-1}*E)
    DVector<double> h_{}; // h = (\lambda*L - I)*T^{-1}*R1^T*R0^{-1}*u
    DVector<double> p_{}; // p = \Psi*h - dS*y

    DMatrix<double> S_  {}; // S = \Psi*T^{-1}*\Psi^T*Q
    DMatrix<double> dS_ {}; // dS = -\Psi*(T^{-1}*R)*(T^{-1}*E)
    DMatrix<double> ddS_{}; // ddS = 2*\Psi*L*F

    // compute first derivative of matrix S: dS = -\Psi*(T^{-1}*R)*(T^{-1}*E)
    const DMatrix<double>& dS() {
      L_ = (trS_.invT_).solve(model_.R()); // T^{-1}*R
      F_ = L_*(trS_.invT_).solve(trS_.E_); // (T^{-1}*R)*(T^{-1}*E)
      dS_ = model_.Psi()*(-F_); 
      return dS_;
    }  
    // compute second derivative of matrix S: ddS = 2*\Psi*L*F
    const DMatrix<double>& ddS(){
      ddS_ = model_.Psi()*2*L_*F_; 
      return ddS_;
    }
    
    // computes the a term in the dGCV expression, given by
    // a = p.dot(y - \hat y)
    //   p = \Psi*h - t
    //     h = (\lambda*L - I)*T^{-1}*g
    //       g = R1^T*R0^{-1}*u
    //     t = dS*y
    double a(){
      DMatrix<double> g = model_.R1().transpose()*model_.invR0().solve(model_.u());
      // cache h and p since needed for computation of second derivative
      h_ = (model_.lambdaS()*L_ - DMatrix<double>::Identity(model_.n_locs(), model_.n_locs()))*(trS_.invT_).solve(g);
      p_ = model_.Psi()*h_ - dS_*model_.y();
      // return a = p.dot(y - \hat y)
      return (( model_.y() - model_.fitted() ).transpose() * p_).coeff(0,0);
    }
  
    // computes the b term in the ddGCV expression, given by
    // b = p.dot(Q*p) + (-ddS*y - 2*\Psi*L*h).dot(y - \hat y) 
    //   p = \Psi*h - t
    //     h = (\lambda*L - I)*T^{-1}*g
    //       g = R1^T*R0^{-1}*u
    //     t = dS*y
    double b(){
      // NB: ddS_ must already contain valid data
      DMatrix<double> C = 2*L_*h_;
      // perform efficient multiplication by permutation matrix Psi
      DMatrix<double> D(model_.n_locs(), 1); // 2*\Psi*L*h
      for(std::size_t k = 0; k < model_.Psi().outerSize(); ++k){
	for(SpMatrix<double>::InnerIterator it(model_.Psi(),k); it; ++it){
	  D.row(it.row()) = C.row(it.col());
	}
      }
      DVector<double> Qp_;
      if(model_.hasCovariates())
	Qp_ = model_.lmbQ(p_); // efficient computation of Q*p
      else Qp_ = model_.W()*p_;
      // return b = p.dot(Q*p) + (-ddS*y - 2*\Psi*L*h).dot(y - \hat y)
      return (( model_.y() - model_.fitted() ).transpose() * ( -ddS_*model_.y() - D )).coeff(0,0) + p_.dot(Qp_);
    }
    
  public:
    // constructor
    ExactGCV(M& model) : GCV<M, ExactEDF<M>>(model) {};
    
    // analytical expression of gcv first derivative (called only by an exact-based GCV optimization)
    //
    // edf      = n - (q + Tr[S])
    // \sigma^2 = \frac{norm(y - \hat y)^2}{n - (q + Tr[S])}
    // a        = p.dot(y - \hat y)
    // dGCV(\lambda) = \frac{2n}{edf^2}[ \sigma^2 * Tr[dS] + a ]
    std::function<SVector<1>(SVector<1>)> derive() {
      return [*this](SVector<1> lambda) mutable -> SVector<1> {
	// fit the model given current lambda
	model_.setLambda(lambda);
	model_.init_model();
	model_.solve();
	// compute trace of matrix S and its first derivative given current lambda
	double trS  = trS_.compute();
	double trdS = dS().trace();
	
	double q = model_.q();           // number of covariates
	std::size_t n = model_.n_locs(); // number of locations
	double edf = n - (q+trS);        // equivalent degrees of freedom
	// \sigma^2 = \frac{(y - \hat y).squaredNorm()}{n - (q + Tr[S])}
	double sigma = ( model_.y() - model_.fitted() ).squaredNorm()/edf;
	// return gradient of GCV at point
	return SVector<1>( 2*n/std::pow(n - (q+trS), 2)*( sigma*trdS + a() ) );
      };
    }

    // analytical expression of gcv second derivative (called only by an exact-based GCV optimization)
    //
    // edf      = n - (q + Tr[S])
    // \sigma^2 = \frac{norm(y - \hat y)^2}{n - (q + Tr[S])}
    // b        = p.dot(Q*p) + (-ddS*y - 2*\Psi*L*h).dot(y - \hat y)
    // ddGCV(\lambda) = \frac{2n}{edf^2}[ \frac{1}{edf}(3*\sigma^2*Tr[dS] + 4*a)*Tr[dS] + \sigma^2*Tr[ddS] + b ]
    std::function<SMatrix<1>(SVector<1>)> deriveTwice() {
      return [*this](SVector<1> lambda) mutable -> SMatrix<1> {
	// fit the model given current lambda
	model_.setLambda(lambda);
	model_.init_model();
	model_.solve();
	// compute trace of matrix S and its first and second derivative given current lambda
	double trS   = trS_.compute();
	double trdS  = dS().trace();
	double trddS = ddS().trace();
      
	double q = model_.q();           // number of covariates
	std::size_t n = model_.n_locs(); // number of locations
	double edf = n - (q+trS);        // equivalent degrees of freedom
	// \sigma^2 = \frac{norm(y - \hat y)^2}{n - (q + Tr[S])}
	double sigma = ( model_.y() - model_.fitted() ).squaredNorm()/edf;
	// return hessian of GCV at point
	return SMatrix<1>( 2*n/std::pow(edf, 2)*( trdS/edf*( 3*sigma*trdS + 4*a() ) + sigma*trddS + b() ) );      
      };
    }
  };

}}
  
#endif // __GCV_H__

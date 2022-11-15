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
#include "ExactGCVEngine.h"
#include "StochasticGCVEngine.h"

// interfaces
#include "iGCV.h"
#include "../models/regression/iRegressionModel.h"
using fdaPDE::models::is_regression_model;

namespace fdaPDE{
namespace calibration{

  // M: statistical model of which the GCV should be optimized
  // T: the trace evaluation strategy, defaulted to stochastic for efficiency reasons.
  template <typename M, typename T = StochasticGCVEngine>
  class GCV {
    // guarantees statistical model M is compatible with GCV computation
    static_assert(is_regression_model<M>::value && std::is_base_of<iGCV, M>::value,
		  "you are asking to calibrate something which cannot be handled by GCV");
  private:
    // analytical expression of gcv field
    std::function<double(SVector<1>)> gcv;
    // analytical expression of first and second derivative of gcv
    std::function<SVector<1>(SVector<1>)> dgcv;
    std::function<SMatrix<1>(SVector<1>)> ddgcv;
  
    M& model_; // the model to which gcv have to be optimized
    T trace;   // Tr[S] evaluation strategy

    // initialize gcv functors
    void init() {
      // analytical expression of gcv (called by any type of GCV optimization)
      //
      // edf      = n - (q + Tr[S])
      // GCV(\lambda) = n/(edf^2)*norm(y - \hat y)^2
      gcv = [*this](SVector<1> lambda) mutable -> double {
	// fit the model given current lambda
	model_.setLambda(lambda[0]);
	model_.solve();
	// compute trace of matrix S given current lambda
	double trS = trace.compute(model_);

	double q = model_.q();        // number of covariates
	std::size_t n = model_.loc(); // number of locations
	double edf = n - (q+trS);     // equivalent degrees of freedom
	// return gcv at point
	return (n/std::pow(edf, 2))*( model_.norm(model_.y(), model_.fitted()) );
      };

      // analytical expression of gcv first derivative (called only by an exact-based GCV optimization)
      //
      // edf      = n - (q + Tr[S])
      // \sigma^2 = \frac{norm(y - \hat y)^2}{n - (q + Tr[S])}
      // a        = p.dot(y - \hat y) (see ExactGCVEngine.h)
      // dGCV(\lambda) = \frac{2n}{edf^2}[ \sigma^2 * Tr[dS] + a ]
      dgcv = [*this](SVector<1> lambda) mutable -> SVector<1> {
	// fit the model given current lambda
	model_.setLambda(lambda[0]);
	model_.solve();
	// compute trace of matrix S and its firstderivative given current lambda
	double trS  = trace.compute(model_);
	double trdS = trace.derive(model_);
      
	double q = model_.q();        // number of covariates
	std::size_t n = model_.loc(); // number of locations
	double edf = n - (q+trS);     // equivalent degrees of freedom
	// \sigma^2 = \frac{(y - \hat y).squaredNorm()}{n - (q + Tr[S])}
	double sigma = ( model_.y() - model_.fitted() ).squaredNorm()/edf;

	double a = trace.a(model_);   // a = p.dot(y - \hat y)
	// return gradient of GCV at point      
	return SVector<1>( 2*n/std::pow(n - (q+trS), 2)*( sigma*trdS + a ) );
      };

      // analytical expression of gcv second derivative (called only by an exact-based GCV optimization)
      //
      // edf      = n - (q + Tr[S])
      // \sigma^2 = \frac{norm(y - \hat y)^2}{n - (q + Tr[S])}
      // b        = p.dot(Q*p) + (-ddS*y - 2*\Psi*L*h).dot(y - \hat y) (see ExactGCVEngine)
      // ddGCV(\lambda) = \frac{2n}{edf^2}[ \frac{1}{edf}(3*\sigma^2*Tr[dS] + 4*a)*Tr[dS] + \sigma^2*Tr[ddS] + b ]
      ddgcv = [*this](SVector<1> lambda) mutable -> SMatrix<1> {
	// fit the model given current lambda
	model_.setLambda(lambda[0]);
	model_.solve();
	// compute trace of matrix S and its first and second derivative given current lambda
	double trS   = trace.compute(model_);
	double trdS  = trace.derive(model_);
	double trddS = trace.deriveTwice(model_);
      
	double q = model_.q();        // number of covariates
	std::size_t n = model_.loc(); // number of locations
	double edf = n - (q+trS);     // equivalent degrees of freedom
	// \sigma^2 = \frac{norm(y - \hat y)^2}{n - (q + Tr[S])}
	double sigma = ( model_.y() - model_.fitted() ).squaredNorm()/edf;

	double a = trace.a(model_);   // a = p.dot(y - \hat y)
	double b = trace.b(model_);   // b = p.dot(Q*p) + (-ddS*y - 2*\Psi*L*h).dot(y - \hat y)
	// return hessian of GCV at point
	return SMatrix<1>( 2*n/std::pow(edf, 2)*( trdS/edf*(3*sigma*trdS + 4*a) + sigma*trddS + b ) );
      };
    }
  
  public:
    // SFINAE selection of constructor depending on trace evaluation strategy
    template <typename U = T, // fake type to enable substitution
	      typename std::enable_if<
		!std::is_same<U, StochasticGCVEngine>::value,
		int>::type = 0> 
    GCV(M& model)
      : model_(model) { init(); };

    template <typename U = T,
	      typename std::enable_if<
		std::is_same<U, StochasticGCVEngine>::value,
		int>::type = 0>
    GCV(M& model, std::size_t r)
      : model_(model), trace(r) { init(); };
  
    // optimizes GCV in an exact way. Requires analytical expression of first and second derivative
    // Because their stochastic estimates are too unreliable, an exact optimization of the GCV is avoided in case
    // a stochastic approximation of Tr[S] is adopted
    template <typename U = T,
	      typename std::enable_if<
		!std::is_same<U, StochasticGCVEngine>::value, int>::type = 0,
	      typename O, typename... Args>
    SVector<1> exact(O& optimizer, Args&... args) {
      // wrap gcv and its derivatives in a TwiceDifferentiableScalarField object.
      // This forces optimizers to employ the analytical expression of gcv's derivatives when it is required
      TwiceDifferentiableScalarField<1> obj(gcv, dgcv, ddgcv);
      // optimize gcv field
      optimizer.findMinimum(obj, args...);
      SVector<1> x = optimizer.getSolution();
      return x;
    }
  
    // optimizes GCV field using finite differences to approximate first and second derivative
    // (doesn't require to evaluate dgcv and ddgcv)
    template <typename O, typename... Args>
    SVector<1> approx(O& optimizer, double gcvFDStep, Args&... args) {
      // wrap gcv in a ScalarField object.
      // This forces optimizers to employ a finite difference approximation of gcv derivatives when it is required
      ScalarField<1> obj(gcv);
      obj.setStep(gcvFDStep);
      // optimize gcv field
      optimizer.findMinimum(obj, args...);
      SVector<1> x = optimizer.getSolution();
      return x;
    };
  };

  // expose some usefull symbols
  template <typename M> using ExactGCV = GCV<M, ExactGCVEngine>;
  template <typename M> using StochsaticGCV = GCV<M, StochasticGCVEngine>;
}}
  
#endif // __GCV_H__

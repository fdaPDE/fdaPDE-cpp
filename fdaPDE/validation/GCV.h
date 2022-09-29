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
#include "ExactTracer.h"
#include "StochasticTracer.h"

// M: statistical model of which the GCV should be optimized
// T: the trace evaluation strategy, defaulted to stochastic for efficiency reasons.
template <typename M, typename T = StochasticTracer>
class GCV {
  // guarantees statistical model M is compatible with GCV computation
  static_assert(std::is_base_of<iGCV, M>::value,
		"can't apply GCV smoothing parameter selection on a model which doesn't extend iGCV interface");
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
    // analytical expression of gcv
    gcv = [*this](SVector<1> lambda) mutable -> double {
      // fit the model given current lambda
      model_.setLambda(lambda[0]);
      model_.smooth();
      // compute trace of matrix S given current lambda
      double trS = trace.compute(model_);
      // GCV(\lambda) = n/(n-(q+trS(\lambda))^2)*norm(z - \hat_z(\lambda))^2
      double q = model_.q();        // number of covariates
      std::size_t n = model_.loc(); // number of locations
      
      return (n/std::pow(n - (q + trS), 2))*( *model_.z() - model_.fitted() ).squaredNorm();
    };

    // analytical expression of gcv first derivative
    dgcv = [*this](SVector<1> lambda) -> SVector<1> {
      // // fit the model given current lambda
      model_.setLambda(lambda[0]);
      model_.smooth();

      // compute derivative of S
      std::shared_ptr<DMatrix<double>> dS = trace.dS(model_);
      dS*model_.fitted();
      
      lambda*trace.L();
      
      return SVector<1>(0);
    };
    
  }
  
public:
  // SFINAE selection of constructor depending on trace evaluation strategy
  template <typename U = T, // fake type to enable substitution
      typename std::enable_if<
	!std::is_same<U, StochasticTracer>::value,
	int>::type = 0> 
  GCV(M& model)
    : model_(model) { init(); };

  template <typename U = T,
      typename std::enable_if<
	std::is_same<U, StochasticTracer>::value,
	int>::type = 0>
  GCV(M& model, std::size_t r)
    : model_(model), trace(r) { init(); };
  
  // optimizes GCV in an exact way. Requires analytical expression of first and second derivative
  // Because their stochastic estimates are too unreliable, an exact optimization of the GCV is avoided in case
  // a stochastic approximation of Tr[S] is adopted
  template <typename U = T,
	    typename std::enable_if<
	      !std::is_same<U, StochasticTracer>::value, int>::type = 0,
	    typename O, typename... Args>
  SVector<1> exact(O& optimizer, Args&... args) {
    // wrap gcv and its derivatives in a TwiceDifferentiableScalarField object.
    // This forces optimizers to employ the analytical expression of gcv's derivatives when it is required
    TwiceDifferentiableScalarField<1> obj(gcv, dgcv, ddgcv);
    // optimize gcv field
    optimizer.findMinimum(obj, args...);
    // recover found solution
    SVector<1> x = optimizer.getSolution();
    return x;
  }
  
  // optimizes GCV field using finite differences to approximate first and second derivative
  template <typename O, typename... Args>
  SVector<1> FDApprox(O& optimizer, double gcvFDStep, Args&... args) {
    // wrap gcv in a ScalarField object.
    // This forces optimizers to employ a finite difference approximation of gcv derivatives when it is required
    ScalarField<1> obj(gcv);
    obj.setStep(gcvFDStep);
    // optimize gcv field
    optimizer.findMinimum(obj, args...);
    // recover found solution
    SVector<1> x = optimizer.getSolution();
    return x;
  };
};

#endif // __GCV_H__

#ifndef __GCV_H__
#define __GCV_H__

#include <functional>

#include "../core/utils/Symbols.h"
#include "../core/utils/fields/ScalarField.h"
using fdaPDE::core::ScalarField;
#include "../core/NLA/SMW.h"
using fdaPDE::core::NLA::SMW;

// M: statistical model of which the GCV should be optimized
// T: the trace evaluation strategy, defaulted to stochastic for efficiency reasons.
template <typename M, typename T = StochasticTracer>
class GCV {
private:
  // analytical expression of gcv field
  std::function<double(SVector<1>)> gcv;
  
  M& model_; // the model to which gcv have to be optimized
  T trace;   // Tr[S] evaluation strategy

  // initialize gcv functors
  void init() {
    // analytical expression of gcv
    gcv = [*this](SVector<1> lambda) -> double {
      // fit the model with this value of lambda
      model_.setLambda(lambda[0]);
      model_.smooth();
      // compute trace of matrix S given current lambda
      double trS = trace.compute(model_);
      // GCV(\lambda) = n/(n-(q+trS(\lambda))^2)*norm(z - \hat_z(\lambda))^2
      double q = model_.q();      // number of covariates
      std::size_t n = model_.n(); // number of locations

      return (n/std::pow(n - (q + trS), 2))*(model_.z() - model_.fitted()).squaredNorm();
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
  
  // optimizes the GCV in an exact way
  double exact();
  // optimizes the GCV in an approximate way using finite differences to approximate first and second derivative
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

#ifndef __GRADIENT_DESCENT_OPTIMIZER__
#define __GRADIENT_DESCENT_OPTIMIZER__

#include "Function.h"
#include "Optimizer.h"
#include "Customizable.h"

#include <tuple>

template <unsigned int N>
class GradientDescentOptimizer : public CustomizableOptimizer<N>{

private:
  double step;
  
  // internal status of the optimizer 
  SVector<N> x_old;                 // value of the optimization point before the update step
  SVector<N> gradientExact;         // exact value of the gradient before the update step
  SVector<N> x_new;                 // value of the optimization point after the update step
  double error;                     // squared l^2 norm of the gradient after the update step

  // optimization problem data
  const SVector<N>& x0;                            // starting point
  unsigned int maxIteration;                       // maximum number of iterations before forced stop
  double tolerance;                                // tolerance on error
  const DifferentiableScalarField<N>& objective;   // objective function to optimize
  
public:
  // constructor
  GradientDescentOptimizer(double step_, const SVector<N>& x0_, unsigned int maxIteration_,
			   double tolerance_, const DifferentiableScalarField<N>& objective_)
    : step(step_), x0(x0_), maxIteration(maxIteration_), tolerance(tolerance_), objective(objective_) {};
  // optimization routine
  std::pair<SVector<N>, double> findMinimum() override;

  // getters
  SVector<N> getX_old() const         { return x_old; };
  SVector<N> getX_new() const         { return x_new; };
  SVector<N> getGradientExact() const { return gradientExact; };
  double     getError() const         { return error; };
};

#include "GradientDescentOptimizer.tpp"

#endif // __GRADIENT_DESCENT_OPTIMIZER__

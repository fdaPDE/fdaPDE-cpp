#ifndef __NEWTON_OPTIMIZER_H__
#define __NEWTON_OPTIMIZER_H__

#include <Eigen/Core>
#include <Eigen/Dense>

#include "Function.h"
#include "Customizable.h"

using std::pair;

// newton method based on approximate gradient and hessian computation. 
template <unsigned int N>
class NewtonOptimizer : public IterativeOptimizer<N>{

protected:
  double step;
  
  // internal status of the optimizer 
  SVector<N> x_old;                 // value of the optimization point before the update step
  SVector<N> x_new;                 // value of the optimization point after the update step
  double error;                     // squared l^2 norm of the gradient after the update step
  SVector<N> gradientApprox;        // approximated value of the gradient before the update step
  SMatrix<N> hessianApprox;         // approximated value of the hessian before the update step
  
  // optimization problem data
  const SVector<N>& x0;             // starting point
  unsigned int maxIteration;        // maximum number of iterations before forced stop
  double tolerance;                 // tolerance on error
  const ScalarField<N>& objective;  // objective function
  
public:
  NewtonOptimizer(double step_, const SVector<N>& x0_, unsigned int maxIteration_,
		  double tolerance_, const ScalarField<N>& objective_)
    : step(step_), x0(x0_), maxIteration(maxIteration_), tolerance(tolerance_), objective(objective_) {};

  // optimization routine
  std::pair<SVector<N>, double> findMinimum() override;
};

#include "NewtonOptimizer.tpp"

#endif // __NEWTON_OPTIMIZER_H__

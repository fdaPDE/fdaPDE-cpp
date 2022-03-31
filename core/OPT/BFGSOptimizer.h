#ifndef __BFGS_OPTIMIZER_H__
#define __BFGS_OPTIMIZER_H__

#include "Function.h"
#include "Customizable.h"

template <unsigned int N>
class BFGSOptimizer : public CustomizableOptimizer<N>{
  
 private:
  double step;

  // internal status of the optimizer 
  SVector<N> x_old;                 // value of the optimization point before the update step
  SVector<N> x_new;                 // value of the optimization point after the update step
  SVector<N> grad_old;              // value of the field's gradient before the update step
  SVector<N> grad_new;              // value of the field's gradient after the update step
  SMatrix<N> hessian;               // value of the hessian matrix approximation at iteration i.
  double error;                     // squared l^2 norm of the gradient after the update step

  // optimization problem data
  const SVector<N>& x0;                            // starting point
  unsigned int maxIteration;                       // maximum number of iterations before forced stop
  double tolerance;                                // tolerance on error
  const DifferentiableScalarField<N>& objective;   // objective function to optimize

 public:
  // constructor
  BFGSOptimizer(double step_, const SVector<N>& x0_, unsigned int maxIteration_,
		double tolerance_, const DifferentiableScalarField<N>& objective_)
    : step(step_), x0(x0_), maxIteration(maxIteration_), tolerance(tolerance_), objective(objective_) {};
  // optimization routine
  std::pair<SVector<N>, double> findMinimum() override;

  // getters
  SVector<N> getX_old() const    { return x_old; };
  SVector<N> getX_new() const    { return x_new; };
  SVector<N> getGrad_old() const { return grad_old; };
  SVector<N> getGrad_new() const { return grad_new; };
  SMatrix<N> getHessian() const  { return hessian; };
  double     getError() const    { return error; };
};

#include "BFGSOptimizer.tpp"

#endif // __BFGS_OPTIMIZER_H_

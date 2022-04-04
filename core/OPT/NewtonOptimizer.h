#ifndef __NEWTON_OPTIMIZER_H__
#define __NEWTON_OPTIMIZER_H__

#include <Eigen/Dense>
#include <chrono>
#include "ScalarField.h"
#include "extensions/Extension.h"
#include "Utils.h"

// newton method based on approximate gradient and hessian computation. 
template <unsigned int N>
class NewtonOptimizer{

private:
  double step;                      // step employed by the optimization scheme.
  
  // internal status of the optimizer 
  SVector<N> x_old;                 // value of the optimization point before the update step
  SVector<N> x_new;                 // value of the optimization point after the update step
  SVector<N> update;                // update vector computed at each step
  double error;                     // squared l^2 norm of the gradient after the update step

  // optimization problem data
  unsigned int maxIteration;        // maximum number of iterations before forced stop
  double tolerance;                 // tolerance on error
  unsigned int numIteration = 0;    // counter to keep track of the number of iterations executed

  double gradient_step;             // step to use in gradient approximation via central differences
  double hessian_step;              // step to use in hessian approximazion via central differences

  using timeType = std::chrono::duration<long int, std::micro>;
  timeType time;
  
public:
  // constructor
  NewtonOptimizer(unsigned int maxIteration_, double tolerance_, double gradient_step_, double hessian_step_)
    : maxIteration(maxIteration_),
      tolerance(tolerance_),
      gradient_step(gradient_step_),
      hessian_step(hessian_step_)
  {};

  // set step size (use this if you want to employ a fixed step method. For adaptive step, use a proper extension)
  void setStepSize(double step_) { step = step_; }
  
  // getters to internal state
  unsigned int getNumIteration() const { return numIteration; }
  double getError()              const { return error;        }
  SVector<N> getXold()           const { return x_old;        }
  SVector<N> getXnew()           const { return x_new;        }
  SVector<N> getUpdate()         const { return update;       }
  timeType getTime()             const { return time;         }
  
  // optimization routine
  template <typename... Args>
  std::pair<SVector<N>, double> findMinimum(const ScalarField<N>& objective, const SVector<N>& x0, const Args&... args);
};

#include "NewtonOptimizer.tpp"

#endif // __NEWTON_OPTIMIZER_H__

#ifndef __EXACT_NEWTON_H__
#define __EXACT_NEWTON_H__

#include <Eigen/Dense>
#include "../../utils/Symbols.h"
#include "../../utils/fields/ScalarField.h"
#include "../extensions/Extension.h"

namespace fdaPDE{
namespace core{
namespace OPT{  

  // newton method based on exact gradient and hessian computation. It requires second differentiabiliy
  // of the scalar field passed as objective
  template <unsigned int N>
  class ExactNewtonOptimizer{

  private:
    double step;                      // step employed by the optimization scheme.
  
    // internal status of the optimizer 
    SVector<N> x_old;                 // value of the optimization point before the update step
    SVector<N> x_new;                 // value of the optimization point after the update step
    SVector<N> update;                // update vector computed at each step
    double error;                     // squared l^2 norm of the gradient after the update step

    // optimization problem data
    unsigned int maxIt;               // maximum number of iterations before forced stop
    double tolerance;                 // tolerance on error
    unsigned int numIt = 0;           // counter to keep track of the number of iterations executed
  
  public:
    // constructor
    ExactNewtonOptimizer(unsigned int maxIt_, double tolerance_)
      : maxIt(maxIt_), tolerance(tolerance_) {};

    // set step size (use this if you want to employ a fixed step method. For adaptive step, use a proper extension)
    void setStepSize(double step_) { step = step_; }

    // getters to internal state
    unsigned int getNumIteration() const { return numIt;  }
    double getError()              const { return error;  }
    SVector<N> getXold()           const { return x_old;  }
    SVector<N> getXnew()           const { return x_new;  }
    SVector<N> getUpdate()         const { return update; }
  
    // optimization routine
    template <typename... Args>
    std::pair<SVector<N>, double> findMinimum(const TwiceDifferentiableScalarField<N>& objective, const SVector<N>& x0, const Args&... args);

    const std::string description = "Newton method with exact gradient and hessian computation";
  };

#include "ExactNewton.tpp"
}}}

#endif // __EXACT_NEWTON_H__

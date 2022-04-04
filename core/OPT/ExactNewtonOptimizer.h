#ifndef __EXACT_NEWTON_OPTIMIZER__
#define __EXACT_NEWTON_OPTIMIZER__

#include <Eigen/Core>
#include <Eigen/Dense>
#include <utility>

#include "ScalarField.h"
#include "IterativeOptimizer.h"

// newton method based on exact gradient and hessian computation. It requires second differentiabiliy
// of the scalar field passed as objective
template <unsigned int N>
class ExactNewtonOptimizer : public IterativeOptimizer<N> {

protected:
  double step;
  
  // internal status of the optimizer 
  SVector<N> x_old;                 // value of the optimization point before the update step
  SVector<N> x_new;                 // value of the optimization point after the update step
  double error;                     // squared l^2 norm of the gradient after the update step
  SVector<N> gradientExact;         // exact value of the gradient before the update step
  SMatrix<N> hessianExact;          // exact value of the hessian before the update step
  
  // optimization problem data
  const SVector<N>& x0;                                // starting point
  unsigned int maxIteration;                           // maximum number of iterations before forced stop
  double tolerance;                                    // tolerance on error
  const TwiceDifferentiableScalarField<N>& objective;  // objective function

 public:
  
  ExactNewtonOptimizer(double step_, const SVector<N>& x0_, unsigned int maxIteration_,
		       double tolerance_, const ScalarField<N>& objective_)
    : step(step_), x0(x0_), maxIteration(maxIteration_), tolerance(tolerance_), objective(objective_) {};

  // optimization routine
  template <typename... Extension>
  std::pair<SVector<N>, double> findMinimum(Extension&&... extensions) const;
  
};


// newton optimization routine
template <unsigned int N>
template <typename... ExtensionList>
std::pair<SVector<N>, double> ExactNewtonOptimizer<N>::findMinimum(ExtensionList&&... extensions) const {
  Extension::initOptimization(extensions...);                 // execute custom action
  
  // algorithm initialization
  x_old = x0;
  unsigned int numIteration = 0;
  error = objective.derive()(x_old).norm();
  
  while (numIteration < maxIteration && error > tolerance && !this->stopCondition()){
    this->preStep();            // execute custom action
    
    // newton step
    hessianExact  = objective.deriveTwice()(x_old);
    gradientExact = objective.derive()(x_old);

    // solve linear system by using an Householder QR decomposition with column-pivoting: A*P = Q*R
    Eigen::ColPivHouseholderQR<SMatrix<N>> QRdecomposition(hessianExact);
    
    SVector<N> update = QRdecomposition.solve(gradientExact);
    x_new = x_old - ExactNewtonOptimizer<N>::step*update;

    // error update
    error = objective.derive()(x_new).norm();

    this->postStep();           // execute custom action
    
    // prepare next iteration
    x_old = x_new;    
    numIteration++;
  }

  this->finalize();             // execute custom action
  return std::pair<SVector<N>, double>(x_new, objective(x_new));
}

#endif // __EXACT_NEWTON_OPTIMIZER__

#ifndef __EXACT_NEWTON_OPTIMIZER__
#define __EXACT_NEWTON_OPTIMIZER__

#include "Function.h"
#include "NewtonOptimizer.h"

// newton method based on exact gradient and hessian computation. It requires second differentiabiliy
// of the scalar field passed as objective
template <unsigned int N>
class ExactNewtonOptimizer : public NewtonOptimizer<N> {

 public:

 ExactNewtonOptimizer(double step_) : NewtonOptimizer<N>(step_) {};
 std::pair<SVector<N>, double> findMinimumExact(const SVector<N>& x0, unsigned int maxIteration, double tolerance, const TwiceDifferentiableScalarField<N>& objective) const;
  
};

template <unsigned int N>
std::pair<SVector<N>, double> ExactNewtonOptimizer<N>::findMinimumExact(const SVector<N>& x0, unsigned int maxIteration, double tolerance,
									const TwiceDifferentiableScalarField<N>& objective) const {

  // algorithm initialization
  SVector<N> x_old = x0;
  SVector<N> x_new;
  unsigned int numIteration = 0;
  double error = objective.derive()(x_old).norm();
  
  while (numIteration < maxIteration && error > tolerance){

    // newton step
    SMatrix<N> hessianExact  = objective.deriveTwice()(x_old);
    SVector<N> gradientExact = objective.derive()(x_old);

    // solve linear system by using an Householder QR decomposition with column-pivoting: A*P = Q*R
    Eigen::ColPivHouseholderQR<SMatrix<N>> QRdecomposition(hessianExact);
    
    SVector<N> update = QRdecomposition.solve(gradientExact);
    x_new = x_old - NewtonOptimizer<N>::step*update;

    // error update
    error = objective.derive()(x_new).norm();

    // prepare next iteration
    x_old = x_new;    
    numIteration++;
  }
  
  return std::pair<SVector<N>, double>(x_new, objective(x_new));
}

#endif // __EXACT_NEWTON_OPTIMIZER__

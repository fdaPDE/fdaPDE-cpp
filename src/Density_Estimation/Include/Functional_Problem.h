#ifndef __FUNCTIONAL_PROBLEM_H__
#define __FUNCTIONAL_PROBLEM_H__

#include "Data_Problem.h"

// This file implements the functionals of the Density Estimation problem

//! @brief A class to store methods regarding the functional of the problem.
template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class FunctionalProblem{
  private:
    // A member to acess data problem methods
    const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dataProblem_;

    //! A method to compute the integrals of the functional.
    std::pair<Real,VectorXr> computeIntegrals(const VectorXr& g) const;

  public:
    //! A constructor
    FunctionalProblem(const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dp): dataProblem_(dp){};
    //! A method to compute the functional for the g-function. Output: loss, gradient, llik, penterm.
    std::tuple<Real, VectorXr, Real, Real> computeFunctional_g(const VectorXr& g, Real lambda, const SpMat& Psi) const;
    //! A method to compute the log-likelihood and the penalization term for the f-function.
    std::pair<Real,Real> computeLlikPen_f(const VectorXr& f) const;

};

#include "Functional_Problem_imp.h"

#endif

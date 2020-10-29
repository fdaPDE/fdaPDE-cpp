#ifndef __OPTIMIZATION_ALGORITHM_H__
#define __OPTIMIZATION_ALGORITHM_H__

#include <memory>
#include "../../Global_Utilities/Include/Make_Unique.h"
#include "Descent_Direction.h"
#include "Descent_Direction_Factory.h"

// This file contains info of the optimization algorithm of the Density Estimation problem

//! @brief An abtract base class to perform the minimization algorithm.
template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class MinimizationAlgorithm{
  protected:
    // A member to access data problem methods
    const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dataProblem_;
    // A member to access functional  methods
    const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& funcProblem_;
    // A pointer to the object which computes the descent direction
    std::unique_ptr<DirectionBase<Integrator_noPoly, ORDER, mydim, ndim>> direction_;

  public:
    //! A constructor.
    MinimizationAlgorithm(const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dp,
      const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& fp, const std::string& d);
    //! A destructor.
    virtual ~MinimizationAlgorithm(){};
    //! A copy constructor.
    MinimizationAlgorithm(const MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>& rhs);
    //! A pure virtual clone method.
    virtual std::unique_ptr<MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>> clone() const = 0;
    //! A pure virtual method to perform the minimization task.
    virtual VectorXr apply_core(const SpMat& Psi, Real lambda, const VectorXr& g) const = 0;

};


//! @brief A class to perform the minimization algorithm when the step parameter is fixed among all the iterations.
template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class FixedStep : public MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>{
  public:
    //! A delegating constructor.
    FixedStep(const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dp,
      const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& fp,
      const std::string& d):
      MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>(dp, fp, d){};

    //! A copy constructor.
    FixedStep(const FixedStep<Integrator_noPoly, ORDER, mydim, ndim>& rhs):
    MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>(rhs){};
    //! Clone method overridden.
    std::unique_ptr<MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>> clone() const override;
    //! A method to perform the minimization algorithm when the step parameter is fixed among all the iterations.
    VectorXr apply_core(const SpMat& Psi, Real lambda, const VectorXr& g) const override;

};


//! @brief An abstract class to perform the minimization algorithm when the step is computed for each iteration.
template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class AdaptiveStep : public MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>{
  protected:
    //! A copy constructor.
    AdaptiveStep(const AdaptiveStep<Integrator_noPoly, ORDER, mydim, ndim>& rhs):
    MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>(rhs){};
    //! A pure virtual method to compute the step.
    virtual Real computeStep (const VectorXr& g,  Real loss, const VectorXr& grad, const VectorXr& dir, Real lambda, const SpMat& Psi) const = 0;

  public:
    //! A delegating constructor.
    AdaptiveStep(const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dp,
      const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& fp,
      const std::string& d):
      MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>(dp, fp, d){};
    //! A pure virtual clone method.
    virtual std::unique_ptr<MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>> clone() const = 0;
    //! A method to perform the minimization algorithm when the step is computed for each iteration.
    VectorXr apply_core(const SpMat& Psi, Real lambda, const VectorXr& g) const override;

};


//! @brief A class to handle the Backtracking Method.
template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class BacktrackingMethod : public AdaptiveStep<Integrator_noPoly, ORDER, mydim, ndim>{
  private:
    //! A method to compute the step using the Backtracking Method.
    Real computeStep(const VectorXr& g, Real loss, const VectorXr& grad, const VectorXr& dir, Real lambda, const SpMat& Psi) const override;
  public:
    //! A delegating constructor.
    BacktrackingMethod(const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dp,
      const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& fp,
      const std::string& d):
      AdaptiveStep<Integrator_noPoly, ORDER, mydim, ndim>(dp, fp, d){};

    //! A copy constructor.
    BacktrackingMethod(const BacktrackingMethod<Integrator_noPoly, ORDER, mydim, ndim>& rhs):
    AdaptiveStep<Integrator_noPoly, ORDER, mydim, ndim>(rhs){};
    //! Clone method overridden.
    std::unique_ptr<MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>> clone() const override;

};


//! @brief A class to handle the Wolfe Method.
template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class WolfeMethod : public AdaptiveStep<Integrator_noPoly, ORDER, mydim, ndim>{
  private:
    //! A method to compute the step using the Wolfe Method.
    Real computeStep(const VectorXr& g, Real loss, const VectorXr& grad, const VectorXr& dir, Real lambda, const SpMat& Psi) const override;
  public:
    //! A delegating constructor.
    WolfeMethod(const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dp,
      const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& fp,
      const std::string& d):
      AdaptiveStep<Integrator_noPoly, ORDER, mydim, ndim>(dp, fp, d){};

    //! A copy constructor.
    WolfeMethod(const WolfeMethod<Integrator_noPoly, ORDER, mydim, ndim>& rhs):
    AdaptiveStep<Integrator_noPoly, ORDER, mydim, ndim>(rhs){};
    //! Clone method overridden.
    std::unique_ptr<MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>> clone() const override;

};

#include "Optimization_Algorithm_imp.h"

#endif

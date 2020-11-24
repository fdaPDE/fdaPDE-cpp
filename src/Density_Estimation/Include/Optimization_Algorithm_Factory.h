#ifndef __OPTIMIZATION_ALGORITHM_FACTORY_H__
#define __OPTIMIZATION_ALGORITHM_FACTORY_H__

#include <memory>

//!brief@ A Factory class: a class for the choice of the step mehod for the optimization algorithm.
template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class MinimizationAlgorithm_factory
{
	public:
		//! A method that builds a pointer to the right object for the step choice, taking as parameters a string and others objects needed for constructor.
	static std::shared_ptr<MinimizationAlgorithm<Integrator_noPoly, ORDER,  mydim,  ndim>>
  createStepSolver(const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dp,
    const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& fp,
		const std::string& d, const std::string& s)
	{
		if(s == "Fixed_Step") return std::make_shared<FixedStep<Integrator_noPoly, ORDER, mydim, ndim>>(dp, fp, d);

    else if(s == "Backtracking_Method") return std::make_shared<BacktrackingMethod<Integrator_noPoly, ORDER, mydim, ndim>>(dp, fp, d);

    else if(s == "Wolfe_Method") return std::make_shared<WolfeMethod<Integrator_noPoly, ORDER, mydim, ndim>>(dp, fp, d);

		else{

      Rprintf("Unknown step option - using fixed step\n");

			return std::make_shared<FixedStep<Integrator_noPoly, ORDER, mydim, ndim>>(dp, fp,  std::move(d));
		}
  }

};

#endif

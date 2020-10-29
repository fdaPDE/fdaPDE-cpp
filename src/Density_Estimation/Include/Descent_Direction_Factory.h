#ifndef __DESCENT_DIRECTION_FACTORY_H__
#define __DESCENT_DIRECTION_FACTORY_H__

#include <memory>
#include "../../Global_Utilities/Include/Make_Unique.h"

//! brief@ A Factory class: a class for the choice of the step mehod for the optimization algorithm.
template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class DescentDirection_factory
{
	public:
	//! A method that builds a pointer to the right object for the direction choice, taking as parameters a string and others objects needed for constructor.
	static std::unique_ptr<DirectionBase<Integrator_noPoly, ORDER,  mydim,  ndim>>
  createDirectionSolver(const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dp,
    const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& fp, const std::string& d)
	{
		if (d=="Gradient")
			return make_unique<DirectionGradient<Integrator_noPoly,ORDER,mydim,ndim>>(fp);
		else if (d=="BFGS")
			return make_unique<DirectionBFGS<Integrator_noPoly,ORDER,mydim,ndim>>(fp, dp.getNumNodes());
		else{

			Rprintf("Unknown direction option - using gradient direction");

			return make_unique<DirectionGradient<Integrator_noPoly,ORDER,mydim,ndim>>(fp);
		}
	}

};

#endif

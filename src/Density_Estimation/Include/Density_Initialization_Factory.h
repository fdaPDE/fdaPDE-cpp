#ifndef __DENSITY_INITIALIZATION_FACTORY_H__
#define __DENSITY_INITIALIZATION_FACTORY_H__

#include <memory>
#include "../../Global_Utilities/Include/Make_Unique.h"

//! brief@ A Factory class: a class for the choice of the density initialization.
template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class DensityInitialization_factory
{
	public:
	//! A method that builds a pointer to the right object for the initialization choice.
	static std::unique_ptr<DensityInitialization<Integrator_noPoly, ORDER,  mydim,  ndim>>
  createInitializationSolver(const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dp,
    const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& fp){

      if(!dp.isFvecEmpty())
    		return make_unique<UserInitialization<Integrator_noPoly, ORDER, mydim, ndim>>(dp);
    	else
    		return make_unique<HeatProcess<Integrator_noPoly, ORDER, mydim, ndim>>(dp, fp);

  }

};

#endif

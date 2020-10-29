#ifndef __PREPROCESS_FACTORY_H__
#define __PREPROCESS_FACTORY_H__

#include <memory>
#include "../../Global_Utilities/Include/Make_Unique.h"

//! @brief A Factory class: a class for the choice of the cross-validation method.
template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class Preprocess_factory
{
	public:
	//! A method that builds a pointer to the right object for the cross-validation method choice, taking as parameters a string and others objects needed for constructor.
	static std::unique_ptr<Preprocess<Integrator_noPoly, ORDER,  mydim,  ndim>>
  createPreprocessSolver(const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dp,
    const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& fp,
    std::shared_ptr<MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>> ma, const std::string& p){

			if(p=="RightCV")
				return make_unique<RightCrossValidation<Integrator_noPoly, ORDER, mydim, ndim>>(dp, fp, ma);
			else if(p=="SimplifiedCV")
				return make_unique<SimplifiedCrossValidation<Integrator_noPoly, ORDER, mydim, ndim>>(dp, fp, ma);
			else if(p=="NoCrossValidation")
      	return make_unique<NoCrossValidation<Integrator_noPoly, ORDER, mydim, ndim>>(dp, fp);
			else
				return make_unique<RightCrossValidation<Integrator_noPoly, ORDER, mydim, ndim>>(dp, fp, ma);

  }

};

#endif

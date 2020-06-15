#ifndef _DESC_DIR_CONST_HPP
#define _DESC_DIR_CONST_HPP

#include <memory>

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

			#ifdef R_VERSION_
			Rprintf("Unknown direction option - using gradient direction");
			#else
			std::cout<<"Unknown direction option - using gradient direction"<<std::endl;
			#endif

			return make_unique<DirectionGradient<Integrator_noPoly,ORDER,mydim,ndim>>(fp);
		}
	}

};

#endif

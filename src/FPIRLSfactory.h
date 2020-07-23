#ifndef __FPIRLSAFACTORY_HPP__
#define __FPIRLSAFACTORY_HPP__

#include "fdaPDE.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "mesh.h"
#include "param_functors.h"
#include "solver.h"
#include "regressionData.h"
#include "FPIRLS.h"
#include "mixedFEFPCAfactory.h"


#include <memory>

/*
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
*/


//! A Factory class. It is used for the choice of the exponential family distribution for the f-PIRLS.
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class FPIRLSfactory
{
	public:
	//! A method that takes as parameter a string and builds a pointer to the right object for the data distribution
	static std::unique_ptr<FPIRLS<InputHandler, Integrator, ORDER,  mydim,  ndim>> createFPIRLSsolver(const std::string &family, const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr& mu0, Real scale_parameter)
	{
		//initial checks: m0 must be initialized correctly for the different distributions 
		
		if(mu0.size() == 0){
		  VectorXr y = inputData.getObservations();
		  if( family == "binomial" ){ //binary outcomes
		    mu0 = VectorXr::Zero(y.size());
		    for(UInt i = 0; i < y.size(); i++){
		        mu0[i] = 0.5 * (y[i] + 0.5);
		    }//end for-i
		  }else{ // not-binary outcome
		    mu0 = y;
		  }
		}//end if

		if(family=="poisson"){
      		for(UInt i = 0; i < mu0.size(); i++){
        		if(mu0[i]<=0) mu0[i] = 1;
			}
      	}

		// Manage scale_parameter
		bool scale_parameter_flag = false;
		if( (family=="gamma") && scale_parameter<0){
			scale_parameter_flag = true;
		}

		if(family=="binomial"){
		    return make_unique<FPIRLS_Bernoulli<InputHandler, Integrator, ORDER,  mydim, ndim>>(mesh, inputData, mu0);
		}else if(family=="poisson"){
		    return make_unique<FPIRLS_Poisson<InputHandler, Integrator, ORDER,  mydim, ndim>>(mesh, inputData, mu0);
		}else if(family=="exponential"){
		    return make_unique<FPIRLS_Exponential<InputHandler, Integrator, ORDER,  mydim, ndim>>(mesh,inputData,mu0);
		}else if(family=="gamma"){
		    return make_unique<FPIRLS_Gamma<InputHandler, Integrator, ORDER,  mydim, ndim>>(mesh,inputData,mu0, scale_parameter, scale_parameter_flag);
		}

		return std::unique_ptr<FPIRLS<InputHandler, Integrator, ORDER,  mydim,  ndim>>(nullptr);
	}


};


#endif

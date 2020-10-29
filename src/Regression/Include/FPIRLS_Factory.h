#ifndef __FPIRLS_FACTORY_H__
#define __FPIRLS_FACTORY_H__

#include "../../FdaPDE.h"
#include "../../Global_Utilities/Include/Make_Unique.h"
#include "../../FE_Assemblers_Solvers/Include/Finite_Element.h"
#include "../../FE_Assemblers_Solvers/Include/Matrix_Assembler.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../FE_Assemblers_Solvers/Include/Param_Functors.h"
#include "../../FE_Assemblers_Solvers/Include/Solver.h"
#include "Regression_Data.h"
#include "../../Lambda_Optimization/Include/Optimization_Data.h"
#include "FPIRLS.h"


#include <memory>


//! A Factory class. It is used for the choice of the exponential family distribution for the f-PIRLS.
template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
class FPIRLSfactory
{
	public:
	//! A method that takes as parameter a string and builds a pointer to the right object for the data distribution
	static std::unique_ptr<FPIRLS<InputHandler, ORDER,  mydim,  ndim>> createFPIRLSsolver(const std::string &family, const MeshHandler<ORDER,mydim,ndim> & mesh, InputHandler & inputData, OptimizationData & optimizationData, VectorXr& mu0, Real scale_parameter)
	{
		//initial checks: m0 must be initialized correctly for the different distributions

		if(mu0.size() == 0){
		  const VectorXr* y = inputData.getObservations();
		  if( family == "binomial" ){ //binary outcomes
		    mu0 = VectorXr::Zero(y->size());
		    for(UInt i = 0; i < y->size(); i++){
		        mu0[i] = 0.5 * ((*y)[i] + 0.5);
		    }//end for-i
		  }else{ // not-binary outcome
		    mu0 = *y;
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
		    return make_unique<FPIRLS_Bernoulli<InputHandler, ORDER,  mydim, ndim>>(mesh, inputData, optimizationData, mu0);
		}else if(family=="poisson"){
		    return make_unique<FPIRLS_Poisson<InputHandler, ORDER,  mydim, ndim>>(mesh, inputData, optimizationData, mu0);
		}else if(family=="exponential"){
		    return make_unique<FPIRLS_Exponential<InputHandler, ORDER,  mydim, ndim>>(mesh, inputData, optimizationData, mu0);
		}else if(family=="gamma"){
		    return make_unique<FPIRLS_Gamma<InputHandler, ORDER,  mydim, ndim>>(mesh,inputData, optimizationData, mu0, scale_parameter, scale_parameter_flag);
		}

		return std::unique_ptr<FPIRLS<InputHandler, ORDER,  mydim,  ndim>>(nullptr);
	}


};


#endif

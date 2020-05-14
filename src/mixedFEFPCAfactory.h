#ifndef __MIXEDFEFPCAFACTORY_HPP__
#define __MIXEDFEFPCAFACTORY_HPP__

#include "fdaPDE.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "mesh.h"
#include "param_functors.h"
#include "solver.h"
#include "FPCAData.h"
#include "FPCAObject.h"
#include "mixedFEFPCA.h"

#include <memory>


template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

//! A Factory class: A class for the choice of the cross-validation method to use for the selection of the parameter lambda for each PC.
template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFEFPCAfactory
{	
	public:
	//! A method that takes as parameter a string and builds a pointer to the right object for the cross-validation
	static std::unique_ptr<MixedFEFPCABase<Integrator, ORDER,  mydim,  ndim>> createFPCAsolver(const std::string &validation, const MeshHandler<ORDER,mydim,ndim>& mesh, const FPCAData& fpcaData){
	
	if(validation=="GCV") 
	    return make_unique<MixedFEFPCAGCV<Integrator, ORDER,  mydim, ndim>>(mesh,fpcaData);
	    
	else if(validation=="KFold") 
	    return make_unique<MixedFEFPCAKFold<Integrator, ORDER,  mydim, ndim>>(mesh,fpcaData);
	
	else if(validation=="NoValidation") 
	    return make_unique<MixedFEFPCA<Integrator, ORDER,  mydim, ndim>>(mesh,fpcaData);
	
	else{
		#ifdef R_VERSION_
		Rprintf("Unknown validation option - using no validation");
		#else
		std::cout<<"Unknown validation option - using no validation";
		#endif
		
		return make_unique<MixedFEFPCA<Integrator, ORDER,  mydim, ndim>>(mesh,fpcaData);
	}
	
	}

};

#endif

#ifndef __MIXED_FE_FPCA_FACTORY_H__
#define __MIXED_FE_FPCA_FACTORY_H__

#include "../../FdaPDE.h"
#include "../../Global_Utilities/Include/Make_Unique.h"
#include "../../FE_Assemblers_Solvers/Include/Finite_Element.h"
#include "../../FE_Assemblers_Solvers/Include/Matrix_Assembler.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../FE_Assemblers_Solvers/Include/Param_Functors.h"
#include "../../FE_Assemblers_Solvers/Include/Solver.h"
#include "FPCA_Data.h"
#include "FPCA_Object.h"
#include "Mixed_FE_FPCA.h"

//! A Factory class: A class for the choice of the cross-validation method to use for the selection of the parameter lambda for each PC.
class MixedFEFPCAfactory
{
	public:
		//! A method that takes as parameter a string and builds a pointer to the right object for the cross-validation
		static std::unique_ptr<MixedFEFPCABase> createFPCAsolver(const std::string & validation, const FPCAData& fpcaData)
		{
			if(validation=="GCV")
			    return make_unique<MixedFEFPCAGCV>(fpcaData);

			else if(validation=="KFold")
			    return make_unique<MixedFEFPCAKFold>(fpcaData);

			else if(validation=="NoValidation")
			    return make_unique<MixedFEFPCA>(fpcaData);

			else{
				Rprintf("Unknown validation option - using no validation");

				return make_unique<MixedFEFPCA>(fpcaData);
			}
		}
};

#endif

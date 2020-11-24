#ifndef __DE_INITIALIZATION_SKELETON_H__
#define __DE_INITIALIZATION_SKELETON_H__

#include "../../FE_Assemblers_Solvers/Include/Finite_Element.h"
#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh_Objects.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../FE_Assemblers_Solvers/Include/Matrix_Assembler.h"
#include "../../Global_Utilities/Include/Solver_Definitions.h"

//Density Estimation
#include "../../Density_Estimation/Include/Data_Problem.h"
#include "../../Density_Estimation/Include/Functional_Problem.h"
#include "../../Density_Estimation/Include/Optimization_Algorithm.h"
#include "../../Density_Estimation/Include/Optimization_Algorithm_Factory.h"
#include "../../Density_Estimation/Include/FE_Density_Estimation.h"

template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
SEXP DE_init_skeleton(SEXP Rdata, SEXP Rorder, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda, SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals,
	SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP Rmesh, SEXP Rsearch, const std::string & init, UInt init_fold)
{
	// Construct data problem object
	DataProblem<Integrator_noPoly, ORDER, mydim, ndim> dataProblem(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rsearch, Rmesh);

	// Construct functional problem object
	FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim> functionalProblem(dataProblem);

	if(init == "Heat"){

		// Construct densityInit object
		std::unique_ptr<DensityInitialization<Integrator_noPoly, ORDER, mydim, ndim>> densityInit = make_unique<HeatProcess<Integrator_noPoly, ORDER, mydim, ndim>>(dataProblem, functionalProblem);

		// fill fInit
		std::vector<VectorXr> fInit(dataProblem.getNlambda());
		for(UInt l = 0; l < dataProblem.getNlambda(); l++){
			fInit[l] = *(densityInit-> chooseInitialization(dataProblem.getLambda(l)));
		}

		// Copy result in R memory
		SEXP result = NILSXP;
		result = PROTECT(Rf_allocVector(VECSXP, 1));
		SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, ((fInit[0])).size(), fInit.size()));

		Real *rans = REAL(VECTOR_ELT(result, 0));
		for(UInt j = 0; j < fInit.size(); j++)
		{
			for(UInt i = 0; i < (fInit[0]).size(); i++)
				rans[i + (fInit[0]).size()*j] = (fInit[j])[i];
		}

		UNPROTECT(1);

		return(result);
	}

	else if(init=="CV"){

		// Construct densityInit object
		std::unique_ptr<Heat_CV<Integrator_noPoly, ORDER, mydim, ndim>> densityInit = make_unique<Heat_CV<Integrator_noPoly, ORDER, mydim, ndim>>(dataProblem, functionalProblem, init_fold);

		// fill fInit
		VectorXr fInit;
		fInit = *(densityInit->chooseInitialization(0));

		// Copy result in R memory
		SEXP result = NILSXP;
		result = PROTECT(Rf_allocVector(VECSXP, 1));
		SET_VECTOR_ELT(result, 0, Rf_allocVector(REALSXP, fInit.size()));

		Real *rans = REAL(VECTOR_ELT(result, 0));
		for(UInt i = 0; i < fInit.size(); i++)
		{
			rans[i] = fInit[i];
		}

		UNPROTECT(1);

		return(result);
	}
	else{

		#ifdef R_VERSION_
		Rprintf("Invalid initialization");
		#endif

		return NILSXP;
	}

}


#endif

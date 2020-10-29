#ifndef __DE_SKELETON_H__
#define __DE_SKELETON_H__

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
SEXP DE_skeleton(SEXP Rdata, SEXP Rorder, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda, SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals,
	SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP Rmesh, SEXP Rsearch,
	const std::string & step_method, const std::string & direction_method, const std::string & preprocess_method)
{
	// Construct data problem object
	DataProblem<Integrator_noPoly, ORDER, mydim, ndim> dataProblem(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rsearch, Rmesh);

	// Construct functional problem object
	FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim> functionalProblem(dataProblem);

	// Construct minimization algorithm object
	std::shared_ptr<MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>> minimizationAlgo =
		MinimizationAlgorithm_factory<Integrator_noPoly, ORDER, mydim, ndim>::createStepSolver(dataProblem, functionalProblem, direction_method, step_method);

	// Construct FEDE object
	FEDE<Integrator_noPoly, ORDER, mydim, ndim> fede(dataProblem, functionalProblem, minimizationAlgo, preprocess_method);

  	// Perform the whole task
	fede.apply();

	// Collect results
	VectorXr g_sol = fede.getDensity_g();
	std::vector<const VectorXr*> f_init = fede.getInitialDensity();
	Real lambda_sol = fede.getBestLambda();
	std::vector<Real> CV_errors = fede.getCvError();

	std::vector<Point<ndim> > data = dataProblem.getData();

	// Copy result in R memory
	SEXP result = NILSXP;
	result = PROTECT(Rf_allocVector(VECSXP, 5 + 5));
	SET_VECTOR_ELT(result, 0, Rf_allocVector(REALSXP, g_sol.size()));
	SET_VECTOR_ELT(result, 1, Rf_allocMatrix(REALSXP, (*(f_init[0])).size(), f_init.size()));
	SET_VECTOR_ELT(result, 2, Rf_allocVector(REALSXP, 1));
	SET_VECTOR_ELT(result, 3, Rf_allocMatrix(REALSXP, data.size(), ndim));
	SET_VECTOR_ELT(result, 4, Rf_allocVector(REALSXP, CV_errors.size()));


	Real *rans = REAL(VECTOR_ELT(result, 0));
	for(UInt i = 0; i < g_sol.size(); i++)
	{
		rans[i] = g_sol[i];
	}

	Real *rans1 = REAL(VECTOR_ELT(result, 1));
	for(UInt j = 0; j < f_init.size(); j++)
	{
		for(UInt i = 0; i < (*(f_init[0])).size(); i++)
			rans1[i + (*(f_init[0])).size()*j] = (*(f_init[j]))[i];
	}

	Real *rans2 = REAL(VECTOR_ELT(result, 2));
	rans2[0] = lambda_sol;

	Real *rans3 = REAL(VECTOR_ELT(result, 3));
	for(UInt j = 0; j < ndim; j++)
	{
		for(UInt i = 0; i < data.size(); i++)
			rans3[i + data.size()*j] = data[i][j];
	}

	Real *rans4 = REAL(VECTOR_ELT(result, 4));
	for(UInt i = 0; i < CV_errors.size(); i++)
	{
		rans4[i] = CV_errors[i];
	}

	if(dataProblem.getSearch()==2){
		//SEND TREE INFORMATION TO R
		SET_VECTOR_ELT(result, 5, Rf_allocVector(INTSXP, 1)); //tree_header information
		int *rans5 = INTEGER(VECTOR_ELT(result, 5));
		rans5[0] = dataProblem.getMesh().getTree().gettreeheader().gettreelev();

		SET_VECTOR_ELT(result, 6, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain origin
		Real *rans6 = REAL(VECTOR_ELT(result, 6));
		for(UInt i = 0; i < ndim*2; i++)
			rans6[i] = dataProblem.getMesh().getTree().gettreeheader().domainorig(i);

		SET_VECTOR_ELT(result, 7, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain scale
		Real *rans7 = REAL(VECTOR_ELT(result, 7));
		for(UInt i = 0; i < ndim*2; i++)
			rans7[i] = dataProblem.getMesh().getTree().gettreeheader().domainscal(i);


		UInt num_tree_nodes = dataProblem.getMesh().num_elements()+1; //Be careful! This is not equal to number of elements
		SET_VECTOR_ELT(result, 8, Rf_allocMatrix(INTSXP, num_tree_nodes, 3)); //treenode information
		int *rans8 = INTEGER(VECTOR_ELT(result, 8));
		for(UInt i = 0; i < num_tree_nodes; i++)
			rans8[i] = dataProblem.getMesh().getTree().gettreenode(i).getid();

		for(UInt i = 0; i < num_tree_nodes; i++)
			rans8[i + num_tree_nodes*1] = dataProblem.getMesh().getTree().gettreenode(i).getchild(0);

		for(UInt i = 0; i < num_tree_nodes; i++)
			rans8[i + num_tree_nodes*2] = dataProblem.getMesh().getTree().gettreenode(i).getchild(1);

		SET_VECTOR_ELT(result, 9, Rf_allocMatrix(REALSXP, num_tree_nodes, ndim*2)); //treenode box coordinate
		Real *rans9 = REAL(VECTOR_ELT(result, 9));
		for(UInt j = 0; j < ndim*2; j++)
		{
			for(UInt i = 0; i < num_tree_nodes; i++)
				rans9[i + num_tree_nodes*j] = dataProblem.getMesh().getTree().gettreenode(i).getbox().get()[j];
		}
	}

	UNPROTECT(1);

	return(result);
}

#endif

#ifndef __GAM_SKELETON_H__
#define __GAM_SKELETON_H__

#include "../../FdaPDE.h"
#include "../../Regression/Include/FPIRLS.h"
#include "../../Regression/Include/FPIRLS_Factory.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../Lambda_Optimization/Include/Optimization_Data.h"

template<typename InputHandler,UInt ORDER, UInt mydim, UInt ndim>
SEXP GAM_skeleton(InputHandler & GAMData, OptimizationData & optimizationData, SEXP Rmesh, SEXP Rmu0, std::string family, SEXP RscaleParam)
{
  MeshHandler<ORDER, mydim, ndim> mesh(Rmesh, GAMData.getSearch());

	// read Rmu0
	VectorXr mu0;
	UInt n_obs_ = Rf_length(Rmu0);
	mu0.resize(n_obs_);

	UInt count = 0;
	for(UInt i=0;i<n_obs_;++i)
		 mu0[i] = REAL(Rmu0)[i];

 	// read scale param
	Real scale_parameter = REAL(RscaleParam)[0];
	// Factory:
	std::unique_ptr<FPIRLS<InputHandler, ORDER, mydim, ndim>> fpirls = FPIRLSfactory<InputHandler, ORDER, mydim, ndim>::createFPIRLSsolver(family, mesh, GAMData, optimizationData, mu0, scale_parameter);


  	fpirls->apply();


  	const MatrixXv& solution = fpirls->getSolution();
  	const MatrixXr& dof = fpirls->getDOF();
  	const std::vector<Real>& J_value = fpirls->get_J();
  	const MatrixXv& fn_hat = fpirls->getFunctionEst();
  	const std::vector<Real> variance_est = fpirls->getVarianceEst();
  	const std::vector<Real>& GCV = fpirls->getGCV();

  	const UInt bestLambda = optimizationData.get_best_lambda_S();

  	MatrixXv beta;
  	if(GAMData.getCovariates()->rows()==0)
   	{
		beta.resize(1,1);
		beta(0,0).resize(1);
		beta(0,0)(0) = 10e20;
	}
	else
		 beta = fpirls->getBetaEst();

	const MatrixXr & barycenters = fpirls->getBarycenters();
	const VectorXi & elementIds = fpirls->getElementIds();

  	// COMPOSIZIONE SEXP result FOR RETURN

	//Copy result in R memory
	SEXP result = R_NilValue;
 	result = PROTECT(Rf_allocVector(VECSXP, 5+3+5+2));
  	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, solution(0).size(), solution.size()));
  	SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, dof.size()));
  	SET_VECTOR_ELT(result, 2, Rf_allocVector(REALSXP, GCV.size()));
  	SET_VECTOR_ELT(result, 3, Rf_allocVector(INTSXP, 1));
  	SET_VECTOR_ELT(result, 4, Rf_allocMatrix(REALSXP, beta(0).size(), beta.size()));

	//return solution
	Real *rans = REAL(VECTOR_ELT(result, 0));
	for(UInt j = 0; j < solution.size(); j++)
	{
		for(UInt i = 0; i < solution(0).size(); i++)
			rans[i + solution(0).size()*j] = solution(j)(i);
	}

	//return DoF
	Real *rans1 = REAL(VECTOR_ELT(result, 1));
	for(UInt i = 0; i < dof.size(); i++)
	{
		rans1[i] = dof(i);
	}

	//return GCV values
  	Real *rans2 = REAL(VECTOR_ELT(result, 2));
	for(UInt i = 0; i < GCV.size(); i++)
	{
		rans2[i] = GCV[i];
	}

	// Copy best lambda
	UInt *rans3 = INTEGER(VECTOR_ELT(result, 3));
	rans3[0] = bestLambda;

	//return beta hat
	Real *rans4 = REAL(VECTOR_ELT(result, 4));
	for(UInt j = 0; j < beta.size(); j++)
	{
		for(UInt i = 0; i < beta(0).size(); i++)
			rans4[i + beta(0).size()*j] = beta(j)(i);
	}

	if(GAMData.getSearch()==2){

		//SEND TREE INFORMATION TO R
		SET_VECTOR_ELT(result, 5, Rf_allocVector(INTSXP, 1)); //tree_header information
		int *rans5 = INTEGER(VECTOR_ELT(result, 5));
		rans5[0] = mesh.getTree().gettreeheader().gettreelev();

		SET_VECTOR_ELT(result, 6, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain origin
		Real *rans6 = REAL(VECTOR_ELT(result, 6));
		for(UInt i = 0; i < ndim*2; i++)
			rans6[i] = mesh.getTree().gettreeheader().domainorig(i);

		SET_VECTOR_ELT(result, 7, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain scale
		Real *rans7 = REAL(VECTOR_ELT(result, 7));
		for(UInt i = 0; i < ndim*2; i++)
			rans7[i] = mesh.getTree().gettreeheader().domainscal(i);


		UInt num_tree_nodes = mesh.num_elements()+1; //Be careful! This is not equal to number of elements
		SET_VECTOR_ELT(result, 8, Rf_allocMatrix(INTSXP, num_tree_nodes, 3)); //treenode information
		int *rans8 = INTEGER(VECTOR_ELT(result, 8));
		for(UInt i = 0; i < num_tree_nodes; i++)
				rans8[i] = mesh.getTree().gettreenode(i).getid();

		for(UInt i = 0; i < num_tree_nodes; i++)
				rans8[i + num_tree_nodes*1] = mesh.getTree().gettreenode(i).getchild(0);

		for(UInt i = 0; i < num_tree_nodes; i++)
				rans8[i + num_tree_nodes*2] = mesh.getTree().gettreenode(i).getchild(1);

		SET_VECTOR_ELT(result, 9, Rf_allocMatrix(REALSXP, num_tree_nodes, ndim*2)); //treenode box coordinate
		Real *rans9 = REAL(VECTOR_ELT(result, 9));
		for(UInt j = 0; j < ndim*2; j++)
		{
			for(UInt i = 0; i < num_tree_nodes; i++)
				rans9[i + num_tree_nodes*j] = mesh.getTree().gettreenode(i).getbox().get()[j];
		}
	}
	
	//SEND BARYCENTER INFORMATION TO R
	SET_VECTOR_ELT(result, 10, Rf_allocVector(INTSXP, elementIds.rows())); //element id of the locations point (vector)
	int *rans10 = INTEGER(VECTOR_ELT(result, 10));
	for(UInt i = 0; i < elementIds.rows(); i++)
		rans10[i] = elementIds(i);

	SET_VECTOR_ELT(result, 11, Rf_allocMatrix(REALSXP, barycenters.rows(), barycenters.cols())); //barycenter information (matrix)
	Real *rans11 = REAL(VECTOR_ELT(result, 11));
	for(UInt j = 0; j < barycenters.cols(); j++)
	{
		for(UInt i = 0; i < barycenters.rows(); i++)
			rans11[i + barycenters.rows()*j] = barycenters(i,j);
	}


	// GAM PARAMETER ESTIMATIONS
	SET_VECTOR_ELT(result, 12, Rf_allocMatrix(REALSXP, fn_hat(0).size(), fn_hat.size()));
	SET_VECTOR_ELT(result, 13, Rf_allocVector(REALSXP, J_value.size()));
	SET_VECTOR_ELT(result, 14, Rf_allocVector(REALSXP, variance_est.size()));

	//return fn hat
	Real *rans12 = REAL(VECTOR_ELT(result, 12));
	for(UInt j = 0; j < fn_hat.size(); j++)
	{
		for(UInt i = 0; i < fn_hat(0).size(); i++)
			rans12[i + fn_hat(0).size()*j] = fn_hat(j)(i);
	}

	//return J_value
  	Real *rans13 = REAL(VECTOR_ELT(result, 13));
  	for(UInt i = 0; i < J_value.size(); i++)
	{
		rans13[i] = J_value[i];
	}

	//return scale parameter
	Real *rans14 = REAL(VECTOR_ELT(result, 14));
	for(UInt j = 0; j < variance_est.size(); j++){
		rans14[j] = variance_est[j];
	}

	UNPROTECT(1);

	return(result);

}

#endif

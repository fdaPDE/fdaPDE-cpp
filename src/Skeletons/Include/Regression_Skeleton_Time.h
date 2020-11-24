#ifndef __REGRESSION_SKELETON_TIME_H__
#define __REGRESSION_SKELETON_TIME_H__

#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../Regression/Include/Mixed_FE_Regression.h"
#include "../../Lambda_Optimization/Include/Optimization_Data.h"

template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
SEXP regression_skeleton_time(InputHandler & regressionData, OptimizationData & optimizationData, SEXP Rmesh, SEXP Rmesh_time)
{
	constexpr UInt SPLINE_DEGREE=MixedSplineRegression<InputHandler>::SPLINE_DEGREE;

	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh, regressionData.getSearch());//! load the mesh
	UInt n_time = Rf_length(Rmesh_time);
	std::vector<Real> mesh_time(n_time);
	for(UInt i=0; i<n_time; ++i)
	{
		mesh_time[i] = REAL(Rmesh_time)[i];
	}
	MixedFERegression<InputHandler> regression(mesh_time, regressionData, optimizationData, mesh.num_nodes(), SPLINE_DEGREE);//! load data in a C++ object

	regression.preapply(mesh); //! solve the problem (compute the _solution, _dof, _GCV, _beta)
    regression.apply();

	//! copy result in R memory
	MatrixXv const & solution = regression.getSolution();
	MatrixXr const & dof = regression.getDOF();
	MatrixXr const & GCV = regression.getGCV();
	UInt bestLambdaS = optimizationData.get_best_lambda_S();
	UInt bestLambdaT = optimizationData.get_best_lambda_T();
	MatrixXv beta;
	if(regressionData.getCovariates()->rows()==0)
	{
		beta.resize(1,1);
		beta(0,0).resize(1);
		beta(0,0)(0) = 10e20;
	}
	else
		 beta = regression.getBeta();

	const MatrixXr & barycenters = regression.getBarycenters();
	const VectorXi & elementIds = regression.getElementIds();
	//!Copy result in R memory
	SEXP result = NILSXP;
	result = PROTECT(Rf_allocVector(VECSXP, 5+5+2));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, solution(0,0).size(), solution.rows()*solution.cols()));
	SET_VECTOR_ELT(result, 1, Rf_allocMatrix(REALSXP, dof.rows(), dof.cols()));
	SET_VECTOR_ELT(result, 2, Rf_allocMatrix(REALSXP, GCV.rows(), GCV.cols()));
	SET_VECTOR_ELT(result, 3, Rf_allocVector(INTSXP, 2));
	SET_VECTOR_ELT(result, 4, Rf_allocMatrix(REALSXP, beta(0,0).size(), beta.rows()*beta.cols()));

	//! Copy solution
	Real *rans = REAL(VECTOR_ELT(result, 0));
	for(UInt i = 0; i < solution.rows(); i++)
	{
		for(UInt j = 0; j < solution.cols(); j++)
		{
			for(UInt k = 0; k < solution(0,0).size(); k++)
				rans[k + solution(0,0).size()*i + solution(0,0).size()*solution.rows()*j] = solution.coeff(i,j)(k);
		}
	}
	//! Copy dof matrix
	Real *rans1 = REAL(VECTOR_ELT(result, 1));
	for(UInt i = 0; i < dof.rows(); i++)
	{
		for(UInt j = 0; j < dof.cols(); j++)
		{
		rans1[i + dof.rows()*j] = dof.coeff(i,j);
		}
	}
	//! Copy GCV matrix
	Real *rans2 = REAL(VECTOR_ELT(result, 2));
	for(UInt i = 0; i < GCV.rows(); i++)
	{
		for(UInt j = 0; j < GCV.cols(); j++)
		{
		rans2[i + GCV.rows()*j] = GCV.coeff(i,j);
		}
	}
	//! Copy best lambdas
	UInt *rans3 = INTEGER(VECTOR_ELT(result, 3));
	rans3[0] = bestLambdaS;
	rans3[1] = bestLambdaT;
	//! Copy betas
	Real *rans4 = REAL(VECTOR_ELT(result, 4));
	for(UInt i = 0; i < beta.rows(); i++)
	{
		for(UInt j = 0; j < beta.cols(); j++)
		{
			for(UInt k = 0; k < beta(0,0).size(); k++)
				rans4[k + beta(0,0).size()*i + beta(0,0).size()*beta.rows()*j] = beta.coeff(i,j)(k);
		}
	}
	if(regressionData.getSearch()==2){
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

	UNPROTECT(1);
	return(result);
}

#endif


#define R_VERSION_

#include "fdaPDE.h"
#include "regressionData.h"
#include "mesh_objects.h"
#include "mesh.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "FPCAData.h"
#include "FPCAObject.h"
#include "solverdefinitions.h"
//#include <chrono>

#include "mixedFEFPCA.h"
#include "mixedFERegression.h"
#include "mixedFEFPCAfactory.h"

//Density Estimation
#include "DataProblem.h"
#include "FunctionalProblem.h"
#include "OptimizationAlgorithm.h"
#include "OptimizationAlgorithm_factory.h"
#include "FEDensityEstimation.h"

// GAM 
#include "FPIRLS.h"
#include "FPIRLSfactory.h"


template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
SEXP regression_skeleton(InputHandler &regressionData, SEXP Rmesh)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);
	MixedFERegression<InputHandler, Integrator,ORDER, IntegratorGaussP3, 0, 0, mydim, ndim> regression(mesh,regressionData);

	regression.apply();

	const MatrixXv& solution = regression.getSolution();
	const MatrixXr& dof = regression.getDOF();
	const MatrixXr & GCV = regression.getGCV();
	UInt bestLambda = regression.getBestLambdaS();
	MatrixXv beta;
	if(regressionData.getCovariates().rows()==0)
	{
		beta.resize(1,1);
		beta(0,0).resize(1);
		beta(0,0)(0) = 10e20;
	}
	else
		 beta = regression.getBeta();
	
	const MatrixXr & barycenters = regression.getBarycenters();
	const VectorXi & elementIds = regression.getElementIds();

	//Copy result in R memory
	SEXP result = NILSXP;
	result = PROTECT(Rf_allocVector(VECSXP, 5+5+2));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, solution(0).size(), solution.size()));
	SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, solution.size()));
	SET_VECTOR_ELT(result, 2, Rf_allocVector(REALSXP, solution.size()));
	SET_VECTOR_ELT(result, 3, Rf_allocVector(INTSXP, 1));
	SET_VECTOR_ELT(result, 4, Rf_allocMatrix(REALSXP, beta(0).size(), beta.size()));

	Real *rans = REAL(VECTOR_ELT(result, 0));
	for(UInt j = 0; j < solution.size(); j++)
	{
		for(UInt i = 0; i < solution(0).size(); i++)
			rans[i + solution(0).size()*j] = solution(j)(i);
	}

	Real *rans1 = REAL(VECTOR_ELT(result, 1));
	for(UInt i = 0; i < solution.size(); i++)
	{
		rans1[i] = dof(i);
	}

	//! Copy GCV vector
	Real *rans2 = REAL(VECTOR_ELT(result, 2));
	for(UInt i = 0; i < solution.size(); i++)
	{
		rans2[i] = GCV(i);
	}

	//! Copy best lambda
	UInt *rans3 = INTEGER(VECTOR_ELT(result, 3));
	rans3[0] = bestLambda;

	//! Copy betas
	Real *rans4 = REAL(VECTOR_ELT(result, 4));
	for(UInt j = 0; j < beta.size(); j++)
	{
		for(UInt i = 0; i < beta(0).size(); i++)
			rans4[i + beta(0).size()*j] = beta(j)(i);
	}

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


template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
SEXP regression_skeleton_time(InputHandler &regressionData, SEXP Rmesh, SEXP Rmesh_time)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);//! load the mesh
	UInt n_time = Rf_length(Rmesh_time);
	std::vector<Real> mesh_time(n_time);
	for(UInt i=0; i<n_time; ++i)
	{
		mesh_time[i] = REAL(Rmesh_time)[i];
	}
	MixedFERegression<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim> regression(mesh, mesh_time,regressionData);//! load data in a C++ object

	regression.apply(); //! solve the problem (compute the _solution, _dof, _GCV, _beta)

	//! copy result in R memory
	MatrixXv const & solution = regression.getSolution();
	MatrixXr const & dof = regression.getDOF();
	MatrixXr const & GCV = regression.getGCV();
	UInt bestLambdaS = regression.getBestLambdaS();
	UInt bestLambdaT = regression.getBestLambdaT();
	MatrixXv beta;
	if(regressionData.getCovariates().rows()==0)
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


template<typename Integrator,UInt ORDER, UInt mydim, UInt ndim>
SEXP FPCA_skeleton(FPCAData &fPCAData, SEXP Rmesh, std::string validation)
{

	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);

	std::unique_ptr<MixedFEFPCABase<Integrator, ORDER, mydim, ndim>> fpca = MixedFEFPCAfactory<Integrator, ORDER, mydim, ndim>::createFPCAsolver(validation, mesh, fPCAData);

	fpca->apply();

	const std::vector<VectorXr>& loadings = fpca->getLoadingsMat();
	const std::vector<VectorXr>& scores = fpca->getScoresMat();
	const std::vector<Real>& lambdas = fpca->getLambdaPC();
	const std::vector<Real>& variance_explained = fpca->getVarianceExplained();
	const std::vector<Real>& cumsum_percentage = fpca->getCumulativePercentage();
	const std::vector<Real>& var = fpca->getVar();
	const MatrixXr & barycenters = fpca->getBarycenters();
	const VectorXi & elementIds = fpca->getElementIds();

	//Copy result in R memory
	SEXP result = NILSXP;
	//result = PROTECT(Rf_allocVector(VECSXP, 7));
	//### why originally 7? Shouldn't it be 6?
	result = PROTECT(Rf_allocVector(VECSXP, 6+5+2));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, loadings[0].size(), loadings.size()));
	SET_VECTOR_ELT(result, 1, Rf_allocMatrix(REALSXP, scores[0].size(), scores.size()));
	SET_VECTOR_ELT(result, 2, Rf_allocVector(REALSXP, lambdas.size()));
	SET_VECTOR_ELT(result, 3, Rf_allocVector(REALSXP, variance_explained.size()));
	SET_VECTOR_ELT(result, 4, Rf_allocVector(REALSXP, cumsum_percentage.size()));
	SET_VECTOR_ELT(result, 5, Rf_allocVector(REALSXP, var.size()));
	Real *rans = REAL(VECTOR_ELT(result, 0));
	for(UInt j = 0; j < loadings.size(); j++)
	{
		for(UInt i = 0; i < loadings[0].size(); i++)
			rans[i + loadings[0].size()*j] = loadings[j][i];
	}

	Real *rans1 = REAL(VECTOR_ELT(result, 1));
	for(UInt j = 0; j < scores.size(); j++)
	{
		for(UInt i = 0; i < scores[0].size(); i++)
			rans1[i + scores[0].size()*j] = scores[j][i];
	}

	Real *rans2 = REAL(VECTOR_ELT(result, 2));
	for(UInt i = 0; i < lambdas.size(); i++)
	{
		rans2[i] = lambdas[i];
	}

	Real *rans3 = REAL(VECTOR_ELT(result, 3));
	for(UInt i = 0; i < variance_explained.size(); i++)
	{
		rans3[i] = variance_explained[i];
	}

	Real *rans4 = REAL(VECTOR_ELT(result, 4));
	for(UInt i = 0; i < cumsum_percentage.size(); i++)
	{
		rans4[i] = cumsum_percentage[i];
	}
	Real *rans5 = REAL(VECTOR_ELT(result, 5));
	for(UInt i = 0; i < var.size(); i++)
	{
		rans5[i] = var[i];
	}

	//TREE INFORMATION
	SET_VECTOR_ELT(result, 6, Rf_allocVector(INTSXP, 7)); //tree_header information
	int *rans6 = INTEGER(VECTOR_ELT(result, 6));
	rans6[0] = mesh.getTree().gettreeheader().gettreeloc();
	rans6[1] = mesh.getTree().gettreeheader().gettreelev();
	rans6[2] = mesh.getTree().gettreeheader().getndimp();
	rans6[3] = mesh.getTree().gettreeheader().getndimt();
	rans6[4] = mesh.getTree().gettreeheader().getnele();
	rans6[5] = mesh.getTree().gettreeheader().getiava();
	rans6[6] = mesh.getTree().gettreeheader().getiend();

	SET_VECTOR_ELT(result, 7, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain origin
	Real *rans7 = REAL(VECTOR_ELT(result, 7));
	for(UInt i = 0; i < ndim*2; i++)
		rans7[i] = mesh.getTree().gettreeheader().domainorig(i);

	SET_VECTOR_ELT(result, 8, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain scale
	Real *rans8 = REAL(VECTOR_ELT(result, 8));
	for(UInt i = 0; i < ndim*2; i++)
		rans8[i] = mesh.getTree().gettreeheader().domainscal(i);


	UInt num_tree_nodes = mesh.num_elements()+1; //Be careful! This is not equal to number of elements
	SET_VECTOR_ELT(result, 9, Rf_allocMatrix(INTSXP, num_tree_nodes, 3)); //treenode information
	int *rans9 = INTEGER(VECTOR_ELT(result, 9));
	for(UInt i = 0; i < num_tree_nodes; i++)
			rans9[i] = mesh.getTree().gettreenode(i).getid();

	for(UInt i = 0; i < num_tree_nodes; i++)
			rans9[i + num_tree_nodes*1] = mesh.getTree().gettreenode(i).getchild(0);

	for(UInt i = 0; i < num_tree_nodes; i++)
			rans9[i + num_tree_nodes*2] = mesh.getTree().gettreenode(i).getchild(1);

	SET_VECTOR_ELT(result, 10, Rf_allocMatrix(REALSXP, num_tree_nodes, ndim*2)); //treenode box coordinate
	Real *rans10 = REAL(VECTOR_ELT(result, 10));
	for(UInt j = 0; j < ndim*2; j++)
	{
		for(UInt i = 0; i < num_tree_nodes; i++)
			rans10[i + num_tree_nodes*j] = mesh.getTree().gettreenode(i).getbox().get()[j];
	}

	//BARYCENTER INFORMATION
	SET_VECTOR_ELT(result, 11, Rf_allocMatrix(REALSXP, barycenters.rows(), barycenters.cols())); //barycenter information (matrix)
	Real *rans11 = REAL(VECTOR_ELT(result, 11));
	for(UInt j = 0; j < barycenters.cols(); j++)
	{
		for(UInt i = 0; i < barycenters.rows(); i++)
			rans11[i + barycenters.rows()*j] = barycenters(i,j);
	}

	SET_VECTOR_ELT(result, 12, Rf_allocVector(INTSXP, elementIds.rows())); //element id of the locations point (vector)
	int *rans12 = INTEGER(VECTOR_ELT(result, 12));
	for(UInt i = 0; i < elementIds.rows(); i++)
		rans12[i] = elementIds(i);

	UNPROTECT(1);

	return(result);
}


template<typename Integrator, typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
SEXP DE_skeleton(SEXP Rdata, SEXP Rorder, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda, SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals,
	SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP Rmesh, SEXP Rsearch,
	const std::string & step_method, const std::string & direction_method, const std::string & preprocess_method)
{
	// Construct data problem object
	DataProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim> dataProblem(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rsearch, Rmesh);

	// Construct functional problem object
	FunctionalProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim> functionalProblem(dataProblem);

	// Construct minimization algorithm object
	std::shared_ptr<MinimizationAlgorithm<Integrator, Integrator_noPoly, ORDER, mydim, ndim>> minimizationAlgo =
		MinimizationAlgorithm_factory<Integrator, Integrator_noPoly, ORDER, mydim, ndim>::createStepSolver(dataProblem, functionalProblem, direction_method, step_method);

	// Construct FEDE object
	FEDE<Integrator, Integrator_noPoly, ORDER, mydim, ndim> fede(dataProblem, functionalProblem, minimizationAlgo, preprocess_method);

  // Perform the whole task
	fede.apply();

	// Collect results
	VectorXr g_sol = fede.getDensity_g();
	std::vector<const VectorXr*> f_init = fede.getInitialDensity();
	Real lambda_sol = fede.getBestLambda();
	std::vector<Real> CV_errors = fede.getCvError();

	std::vector<Point> data = dataProblem.getData();

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

	UNPROTECT(1);

	return(result);
}


template<typename Integrator, typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
SEXP DE_init_skeleton(SEXP Rdata, SEXP Rorder, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda, SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals,
	SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP Rmesh, SEXP Rsearch, const std::string & init, UInt init_fold)
{
	// Construct data problem object
	DataProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim> dataProblem(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rsearch, Rmesh);

	// Construct functional problem object
	FunctionalProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim> functionalProblem(dataProblem);

	if(init == "Heat"){

		// Construct densityInit object
		std::unique_ptr<DensityInitialization<Integrator, Integrator_noPoly, ORDER, mydim, ndim>> densityInit = make_unique<HeatProcess<Integrator, Integrator_noPoly, ORDER, mydim, ndim>>(dataProblem, functionalProblem);

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
		std::unique_ptr<Heat_CV<Integrator, Integrator_noPoly, ORDER, mydim, ndim>> densityInit = make_unique<Heat_CV<Integrator, Integrator_noPoly, ORDER, mydim, ndim>>(dataProblem, functionalProblem, init_fold);

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


template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
SEXP get_integration_points_skeleton(SEXP Rmesh)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);
	FiniteElement<Integrator,ORDER, mydim, ndim> fe;

	SEXP result;
	PROTECT(result=Rf_allocVector(REALSXP, 2*Integrator::NNODES*mesh.num_elements()));
	for(UInt i=0; i<mesh.num_elements(); i++)
	{
		fe.updateElement(mesh.getElement(i));
		for(UInt l = 0;l < Integrator::NNODES; l++)
		{
			Point p = fe.coorQuadPt(l);
			REAL(result)[i*Integrator::NNODES + l] = p[0];
			REAL(result)[mesh.num_elements()*Integrator::NNODES + i*Integrator::NNODES + l] = p[1];
		}
	}

	UNPROTECT(1);
	return(result);
}

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim, typename A>
SEXP get_FEM_Matrix_skeleton(SEXP Rmesh, EOExpr<A> oper)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);

	FiniteElement<Integrator, ORDER, mydim, ndim> fe;

	SpMat AMat;
	Assembler::operKernel(oper, mesh, fe, AMat);

	//Copy result in R memory
	SEXP result;
	result = PROTECT(Rf_allocVector(VECSXP, 2));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(INTSXP, AMat.nonZeros() , 2));
	SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, AMat.nonZeros()));

	int *rans = INTEGER(VECTOR_ELT(result, 0));
	Real  *rans2 = REAL(VECTOR_ELT(result, 1));
	UInt i = 0;
	for (UInt k=0; k < AMat.outerSize(); ++k)
		{
			for (SpMat::InnerIterator it(AMat,k); it; ++it)
			{
				//std::cout << "(" << it.row() <<","<< it.col() <<","<< it.value() <<")\n";
				rans[i] = 1+it.row();
				rans[i + AMat.nonZeros()] = 1+it.col();
				rans2[i] = it.value();
				i++;
			}
		}
	UNPROTECT(1);
	return(result);
}

template<typename InputHandler,typename Integrator,UInt ORDER, UInt mydim, UInt ndim>
SEXP GAM_skeleton(InputHandler &GAMData, SEXP Rmesh, SEXP Rmu0, std::string family, SEXP RscaleParam)
{
  MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);

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
	std::unique_ptr<FPIRLS<InputHandler, Integrator, ORDER, mydim, ndim>> fpirls = FPIRLSfactory<InputHandler, Integrator, ORDER, mydim, ndim>::createFPIRLSsolver(family, mesh, GAMData, mu0, scale_parameter);


  	fpirls->apply();


  	const MatrixXv& solution = fpirls->getSolution();
  	const MatrixXr& dof = fpirls->getDOF();
  	const std::vector<Real>& J_value = fpirls->get_J();
  	const MatrixXv& fn_hat = fpirls->getFunctionEst();
  	const std::vector<Real> variance_est = fpirls->getVarianceEst();
  	const std::vector<Real>& GCV = fpirls->getGCV();

  	const UInt bestLambda = fpirls->getBestLambdaS();

  	MatrixXv beta;
  	if(GAMData.getCovariates().rows()==0)
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




extern "C" {

//! This function manages the various options for Spatial Regression, Sangalli et al version
/*!
	This function is then called from R code.
	\param Robservations an R-vector containing the values of the observations.
	\param Rdesmat an R-matrix containing the design matrix for the regression.
	\param Rmesh an R-object containg the output mesh from Trilibrary
	\param Rorder an R-integer containing the order of the approximating basis.
	\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param Rcovariates an R-matrix of covariates for the regression model
	\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
	\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
			the other are automatically considered in Neumann Condition.
	\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
	\param GCV an R boolean indicating whether dofs of the model have to be computed or not
	\param RGCVmethod an R-integer indicating the method to use to compute the dofs when GCV is TRUE, can be either 1 (exact) or 2 (stochastic)
	\param Rnrealizations the number of random points used in the stochastic computation of the dofs
	\param Rtune a R-double, Tuning parameter used for the estimation of GCV. called 'GCV.inflation.factor' in R code.
	\param RarealDataAvg an R boolean indicating whether the areal data are averaged or not.
	
	\return R-vector containg the coefficients of the solution
*/

SEXP regression_Laplace(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
					SEXP Rlambda, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues,
					SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch, SEXP Rtune, SEXP RarealDataAvg )
{
    //Set input data
	RegressionData regressionData(Rlocations, RbaryLocations, Robservations, Rorder, Rlambda, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, GCV, RGCVmethod, Rnrealizations, DOF, RDOF_matrix, Rsearch, Rtune, RarealDataAvg );

	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

    if(regressionData.getOrder()==1 && mydim==2 && ndim==2)
    	return(regression_skeleton<RegressionData,IntegratorTriangleP2, 1, 2, 2>(regressionData, Rmesh));
    else if(regressionData.getOrder()==2 && mydim==2 && ndim==2)
		return(regression_skeleton<RegressionData,IntegratorTriangleP4, 2, 2, 2>(regressionData, Rmesh));
    else if(regressionData.getOrder()==1 && mydim==2 && ndim==3)
		return(regression_skeleton<RegressionData,IntegratorTriangleP2, 1, 2, 3>(regressionData, Rmesh));
   else if(regressionData.getOrder()==2 && mydim==2 && ndim==3)
		return(regression_skeleton<RegressionData,IntegratorTriangleP4, 2, 2, 3>(regressionData, Rmesh));
	else if(regressionData.getOrder()==1 && mydim==3 && ndim==3)
		return(regression_skeleton<RegressionData,IntegratorTetrahedronP2, 1, 3, 3>(regressionData, Rmesh));
    return(NILSXP);
}

/*!
	This function is then called from R code.
	\param Robservations an R-vector containing the values of the observations.
	\param Rdesmat an R-matrix containing the design matrix for the regression.
	\param Rmesh an R-object containg the output mesh from Trilibrary
	\param Rorder an R-integer containing the order of the approximating basis.
	\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param RK an R-matrix representing the diffusivity matrix of the model
	\param Rbeta an R-vector representing the advection term of the model
	\param Rc an R-double representing the reaction term of the model
	\param Rcovariates an R-matrix of covariates for the regression model
	\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
	\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
			the other are automatically considered in Neumann Condition.
	\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
	\param GCV an R boolean indicating whether dofs of the model have to be computed or not
	\param RGCVmethod an R-integer indicating the method to use to compute the dofs when GCV is TRUE, can be either 1 (exact) or 2 (stochastic)
	\param Rnrealizations the number of random points used in the stochastic computation of the dofs
	\param Rtune a R-double, Tuning parameter used for the estimation of GCV. called 'GCV.inflation.factor' in R code.
	\param RarealDataAvg an R boolean indicating whether the areal data are averaged or not.
	
	\return R-vector containg the coefficients of the solution
*/

SEXP regression_PDE(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
					SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RincidenceMatrix,
					SEXP RBCIndices, SEXP RBCValues, SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch, SEXP Rtune, SEXP RarealDataAvg)
{
	RegressionDataElliptic regressionData(Rlocations, RbaryLocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, GCV, RGCVmethod, Rnrealizations, DOF, RDOF_matrix, Rsearch, Rtune, RarealDataAvg);

	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	if(regressionData.getOrder() == 1 && ndim==2)
		return(regression_skeleton<RegressionDataElliptic,IntegratorTriangleP2, 1, 2, 2>(regressionData, Rmesh));
	else if(regressionData.getOrder() == 2 && ndim==2)
		return(regression_skeleton<RegressionDataElliptic,IntegratorTriangleP4, 2, 2, 2>(regressionData, Rmesh));
	else if(regressionData.getOrder() == 1 && ndim==3)
		return(regression_skeleton<RegressionDataElliptic,IntegratorTriangleP2, 1, 2, 3>(regressionData, Rmesh));
	else if(regressionData.getOrder() == 2 && ndim==3)
		return(regression_skeleton<RegressionDataElliptic,IntegratorTriangleP4, 2, 2, 3>(regressionData, Rmesh));
	return(NILSXP);
}

/*!
	This function is then called from R code.
	\param Robservations an R-vector containing the values of the observations.
	\param Rdesmat an R-matrix containing the design matrix for the regression.
	\param Rmesh an R-object containg the output mesh from Trilibrary
	\param Rorder an R-integer containing the order of the approximating basis.
	\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param RK an R object representing the diffusivity tensor of the model
	\param Rbeta an R object representing the advection function of the model
	\param Rc an R object representing the reaction function of the model
	\param Ru an R object representing the forcing function of the model
	\param Rcovariates an R-matrix of covariates for the regression model
	\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
	\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
			the other are automatically considered in Neumann Condition.
	\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
	\param GCV an R boolean indicating whether dofs of the model have to be computed or not
	\param RGCVmethod an R-integer indicating the method to use to compute the dofs when GCV is TRUE, can be either 1 (exact) or 2 (stochastic)
	\param Rnrealizations the number of random points used in the stochastic computation of the dofs	
	\param Rtune a R-double, Tuning parameter used for the estimation of GCV. called 'GCV.inflation.factor' in R code.
	\param RarealDataAvg an R boolean indicating whether the areal data are averaged or not.
	
	\return R-vector containg the coefficients of the solution
*/


SEXP regression_PDE_space_varying(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
								SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RincidenceMatrix,
								SEXP RBCIndices, SEXP RBCValues, SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch, SEXP Rtune, SEXP RarealDataAvg)
{
    //Set data
	RegressionDataEllipticSpaceVarying regressionData(Rlocations, RbaryLocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Ru, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, GCV,  RGCVmethod, Rnrealizations, DOF, RDOF_matrix, Rsearch, Rtune, RarealDataAvg);

	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	if(regressionData.getOrder() == 1 && ndim==2)
		return(regression_skeleton<RegressionDataEllipticSpaceVarying,IntegratorTriangleP2, 1, 2, 2>(regressionData, Rmesh));
	else if(regressionData.getOrder() == 2 && ndim==2)
		return(regression_skeleton<RegressionDataEllipticSpaceVarying,IntegratorTriangleP4, 2, 2, 2>(regressionData, Rmesh));
	else if(regressionData.getOrder() == 1 && ndim==3)
		return(regression_skeleton<RegressionDataEllipticSpaceVarying,IntegratorTriangleP2, 1, 2, 3>(regressionData, Rmesh));
	else if(regressionData.getOrder() == 2 && ndim==3)
		return(regression_skeleton<RegressionDataEllipticSpaceVarying,IntegratorTriangleP4, 2, 2, 3>(regressionData, Rmesh));
	return(NILSXP);
}

////////////////////////////////////////////////////////////////////////
//												 		  SPACE TIME													 //
//////////////////////////////////////////////////////////////////////

//! This function manages the various options for Spatial Regression, Sangalli et al version
/*!
	This function is then called from R code.
	\param Rlocations an R-matrix containing the spatial locations of the observations
	\param Rtime_locations an R-vector containing the temporal locations of the observations
	\param Robservations an R-vector containing the values of the observations.
	\param Rmesh an R-object containing the spatial mesh
	\param Rmesh_time an R-vector containing the temporal mesh
	\param Rorder an R-integer containing the order of the approximating basis in space.
	\param Rmydim an R-integer specifying if the mesh nodes lie in R^2 or R^3
	\param Rndim  an R-integer specifying if the "local dimension" is 2 or 3
	\param RlambdaS an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param RlambdaT an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param Rcovariates an R-matrix of covariates for the regression model
	\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
	\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
			the other are automatically considered in Neumann Condition.
	\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
	\param Rflag_mass an R-integer that in case of separable problem specifies whether to use mass discretization or identity discretization
	\param Rflag_parabolic an R-integer specifying if the problem is parabolic or separable
	\param Ric an R-vector containing the initial condition needed in case of parabolic problem
	\param GCV an R-integer indicating if the GCV has to be computed or not
	\param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic)
	\param DOF an R boolean indicating whether dofs of the model have to be computed or not
	\param RDOF_matrix a R-matrix containing the dofs (for every combination of the values in RlambdaS and RlambdaT) if they are already known from precedent computations
	\param Rnrealizations the number of random points used in the stochastic computation of the dofs
	\param Rtune a R-double, Tuning parameter used for the estimation of GCV. called 'GCV.inflation.factor' in R code.
	\param RarealDataAvg an R boolean indicating whether the areal data are averaged or not.
	
	\return R-vector containg the coefficients of the solution
*/

SEXP regression_Laplace_time(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rmesh, SEXP Rmesh_time, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
					SEXP RlambdaS, SEXP RlambdaT, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric,
					SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch, SEXP Rtune, SEXP RarealDataAvg)
{
    //Set input data
	RegressionData regressionData(Rlocations, RbaryLocations, Rtime_locations, Robservations, Rorder, RlambdaS, RlambdaT, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, Rflag_mass, Rflag_parabolic, Ric, GCV, RGCVmethod, Rnrealizations, DOF, RDOF_matrix, Rsearch, Rtune, RarealDataAvg);

	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

    if(regressionData.getOrder()==1 && mydim==2 && ndim==2)
    	return(regression_skeleton_time<RegressionData,IntegratorTriangleP2, 1, IntegratorGaussP5, 3, 2, 2, 2>(regressionData, Rmesh, Rmesh_time));
    else if(regressionData.getOrder()==2 && mydim==2 && ndim==2)
		return(regression_skeleton_time<RegressionData,IntegratorTriangleP4, 2, IntegratorGaussP5, 3, 2, 2, 2>(regressionData, Rmesh, Rmesh_time));
    else if(regressionData.getOrder()==1 && mydim==2 && ndim==3)
		return(regression_skeleton_time<RegressionData,IntegratorTriangleP2, 1, IntegratorGaussP5, 3, 2, 2, 3>(regressionData, Rmesh, Rmesh_time));
   else if(regressionData.getOrder()==2 && mydim==2 && ndim==3)
		return(regression_skeleton_time<RegressionData,IntegratorTriangleP4, 2, IntegratorGaussP5, 3, 2, 2, 3>(regressionData, Rmesh, Rmesh_time));
	else if(regressionData.getOrder()==1 && mydim==3 && ndim==3)
		return(regression_skeleton_time<RegressionData,IntegratorTetrahedronP2, 1, IntegratorGaussP5, 3, 2, 3, 3>(regressionData, Rmesh, Rmesh_time));
    return(NILSXP);
}

/*!
	This function is then called from R code.
	\param Rlocations an R-matrix containing the spatial locations of the observations
	\param Rtime_locations an R-vector containing the temporal locations of the observations
	\param Robservations an R-vector containing the values of the observations.
	\param Rmesh an R-object containing the spatial mesh
	\param Rmesh_time an R-vector containing the temporal mesh
	\param Rorder an R-integer containing the order of the approximating basis in space.
	\param Rmydim an R-integer specifying if the mesh nodes lie in R^2 or R^3
	\param Rndim  an R-integer specifying if the "local dimension" is 2 or 3
	\param RlambdaS an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param RlambdaT an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param RK an R-matrix representing the diffusivity matrix of the model
	\param Rbeta an R-vector representing the advection term of the model
	\param Rc an R-double representing the reaction term of the model
	\param Rcovariates an R-matrix of covariates for the regression model
	\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
	\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
			the other are automatically considered in Neumann Condition.
	\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
	\param Rflag_mass an R-integer that in case of separable problem specifies whether to use mass discretization or identity discretization
	\param Rflag_parabolic an R-integer specifying if the problem is parabolic or separable
	\param Ric an R-vector containing the initial condition needed in case of parabolic problem
	\param GCV an R-integer indicating if the GCV has to be computed or not
	\param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic)
	\param DOF an R boolean indicating whether dofs of the model have to be computed or not
	\param RDOF_matrix a R-matrix containing the dofs (for every combination of the values in RlambdaS and RlambdaT) if they are already known from precedent computations
	\param Rnrealizations the number of random points used in the stochastic computation of the dofs
	\param Rtune a R-double, Tuning parameter used for the estimation of GCV. called 'GCV.inflation.factor' in R code.
	\param RarealDataAvg an R boolean indicating whether the areal data are averaged or not.
	
	\return R-vector containg the coefficients of the solution
*/

SEXP regression_PDE_time(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rmesh, SEXP Rmesh_time, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
					SEXP RlambdaS, SEXP RlambdaT, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RincidenceMatrix,
					SEXP RBCIndices, SEXP RBCValues, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric, SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch, SEXP Rtune, SEXP RarealDataAvg)
{
	RegressionDataElliptic regressionData(Rlocations, RbaryLocations, Rtime_locations, Robservations, Rorder, RlambdaS, RlambdaT, RK, Rbeta, Rc, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, Rflag_mass, Rflag_parabolic, Ric, GCV, RGCVmethod, Rnrealizations, DOF, RDOF_matrix, Rsearch, Rtune, RarealDataAvg);

	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	if(regressionData.getOrder() == 1 && ndim==2)
		return(regression_skeleton_time<RegressionDataElliptic,IntegratorTriangleP2, 1, IntegratorGaussP5, 3, 2, 2, 2>(regressionData, Rmesh, Rmesh_time));
	else if(regressionData.getOrder() == 2 && ndim==2)
		return(regression_skeleton_time<RegressionDataElliptic,IntegratorTriangleP4, 2, IntegratorGaussP5, 3, 2, 2, 2>(regressionData, Rmesh, Rmesh_time));
	else if(regressionData.getOrder() == 1 && ndim==3)
		return(regression_skeleton_time<RegressionDataElliptic,IntegratorTriangleP2, 1, IntegratorGaussP5, 3, 2, 2, 3>(regressionData, Rmesh, Rmesh_time));
	else if(regressionData.getOrder() == 2 && ndim==3)
		return(regression_skeleton_time<RegressionDataElliptic,IntegratorTriangleP4, 2, IntegratorGaussP5, 3, 2, 2, 3>(regressionData, Rmesh, Rmesh_time));
	return(NILSXP);
}


/*!
	This function is then called from R code.
	\param Rlocations an R-matrix containing the spatial locations of the observations
	\param Rtime_locations an R-vector containing the temporal locations of the observations
	\param Robservations an R-vector containing the values of the observations.
	\param Rmesh an R-object containing the spatial mesh
	\param Rmesh_time an R-vector containing the temporal mesh
	\param Rorder an R-integer containing the order of the approximating basis in space.
	\param Rmydim an R-integer specifying if the mesh nodes lie in R^2 or R^3
	\param Rndim  an R-integer specifying if the "local dimension" is 2 or 3
	\param RlambdaS an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param RlambdaT an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param RK an R object representing the diffusivity tensor of the model
	\param Rbeta an R object representing the advection function of the model
	\param Rc an R object representing the reaction function of the model
	\param Ru an R object representing the forcing function of the model
	\param Rcovariates an R-matrix of covariates for the regression model
	\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
	\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
			the other are automatically considered in Neumann Condition.
	\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
	\param Rflag_mass an R-integer that in case of separable problem specifies whether to use mass discretization or identity discretization
	\param Rflag_parabolic an R-integer specifying if the problem is parabolic or separable
	\param Ric an R-vector containing the initial condition needed in case of parabolic problem
	\param GCV an R-integer indicating if the GCV has to be computed or not
	\param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic)
	\param DOF an R boolean indicating whether dofs of the model have to be computed or not
	\param RDOF_matrix a R-matrix containing the dofs (for every combination of the values in RlambdaS and RlambdaT) if they are already known from precedent computations
	\param Rnrealizations the number of random points used in the stochastic computation of the dofs
	\param Rtune a R-double, Tuning parameter used for the estimation of GCV. called 'GCV.inflation.factor' in R code.
	\param RarealDataAvg an R boolean indicating whether the areal data are averaged or not.
	
	\return R-vector containg the coefficients of the solution
*/


SEXP regression_PDE_space_varying_time(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rmesh, SEXP Rmesh_time, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
					SEXP RlambdaS, SEXP RlambdaT, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RincidenceMatrix,
					SEXP RBCIndices, SEXP RBCValues, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric, SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch, SEXP Rtune, SEXP RarealDataAvg)
{
    //Set data
	RegressionDataEllipticSpaceVarying regressionData(Rlocations, RbaryLocations, Rtime_locations, Robservations, Rorder, RlambdaS, RlambdaT, RK, Rbeta, Rc, Ru, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, Rflag_mass, Rflag_parabolic, Ric, GCV, RGCVmethod, Rnrealizations, DOF, RDOF_matrix, Rsearch, Rtune, RarealDataAvg);

	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	if(regressionData.getOrder() == 1 && ndim==2)
		return(regression_skeleton_time<RegressionDataEllipticSpaceVarying,IntegratorTriangleP2, 1, IntegratorGaussP5, 3, 2, 2, 2>(regressionData, Rmesh, Rmesh_time));
	else if(regressionData.getOrder() == 2 && ndim==2)
		return(regression_skeleton_time<RegressionDataEllipticSpaceVarying,IntegratorTriangleP4, 2, IntegratorGaussP5, 3, 2, 2, 2>(regressionData, Rmesh, Rmesh_time));
	else if(regressionData.getOrder() == 1 && ndim==3)
		return(regression_skeleton_time<RegressionDataEllipticSpaceVarying,IntegratorTriangleP2, 1, IntegratorGaussP5, 3, 2, 2, 3>(regressionData, Rmesh, Rmesh_time));
	else if(regressionData.getOrder() == 2 && ndim==3)
		return(regression_skeleton_time<RegressionDataEllipticSpaceVarying,IntegratorTriangleP4, 2, IntegratorGaussP5, 3, 2, 2, 3>(regressionData, Rmesh, Rmesh_time));
	return(NILSXP);
}

//! A function required for anysotropic and nonstationary regression (only 2D)
/*!
    \return points where the PDE space-varying params are evaluated in the R code
*/
SEXP get_integration_points(SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim)
{
	//Declare pointer to access data from C++
	int order = INTEGER(Rorder)[0];

	//Get mydim and ndim
	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];
//Not implemented for ndim==3
    if(order == 1 && ndim ==2)
    	return(get_integration_points_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh));
    else if(order == 2 && ndim==2)
    	return(get_integration_points_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh));
    return(NILSXP);
}

//! A utility, not used for system solution, may be used for debugging

SEXP get_FEM_mass_matrix(SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim)
{
	int order = INTEGER(Rorder)[0];

	//Get mydim and ndim
	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);

    if(order==1 && ndim==2)
    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh, mass));
	if(order==2 && ndim==2)
		return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh, mass));
	return(NILSXP);
}

//! A utility, not used for system solution, may be used for debugging
SEXP get_FEM_stiff_matrix(SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim)
{
	int order = INTEGER(Rorder)[0];

	//Get mydim and ndim
	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	typedef EOExpr<Stiff> ETMass;   Stiff EStiff;   ETMass stiff(EStiff);

    if(order==1 && ndim==2)
    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh, stiff));
	if(order==2 && ndim==2)
		return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh, stiff));
	return(NILSXP);
}

//! A utility, not used for system solution, may be used for debugging
SEXP get_FEM_PDE_matrix(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder,SEXP Rmydim, SEXP Rndim, SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc,
				   SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP GCV,SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch, SEXP Rtune, SEXP RarealDataAvg)
{
	RegressionDataElliptic regressionData(Rlocations, RbaryLocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, GCV, RGCVmethod, Rnrealizations, DOF, RDOF_matrix, Rsearch, Rtune, RarealDataAvg);

	//Get mydim and ndim
	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
	typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

	const Real& c = regressionData.getC();
	const Eigen::Matrix<Real,2,2>& K = regressionData.getK();
	const Eigen::Matrix<Real,2,1>& beta = regressionData.getBeta();

    if(regressionData.getOrder()==1 && ndim==2)
    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
	if(regressionData.getOrder()==2 && ndim==2)
		return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
	return(NILSXP);
}

//! A utility, not used for system solution, may be used for debugging
SEXP get_FEM_PDE_space_varying_matrix(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim, SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru,
		   SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP GCV,SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch, SEXP Rtune, SEXP RarealDataAvg)
{
	RegressionDataEllipticSpaceVarying regressionData(Rlocations, RbaryLocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Ru, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, GCV, RGCVmethod, Rnrealizations, DOF, RDOF_matrix, Rsearch, Rtune, RarealDataAvg);

	//Get mydim and ndim
	//UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
	typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

	const Reaction& c = regressionData.getC();
	const Diffusivity& K = regressionData.getK();
	const Advection& beta = regressionData.getBeta();

    if(regressionData.getOrder()==1 && ndim==2)
    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
	if(regressionData.getOrder()==2 && ndim==2)
		return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
	return(NILSXP);
}



//! This function manages the various options for SF-PCA
/*!
	This function is than called from R code.
	\param Rdatamatrix an R-matrix containing the datamatrix of the problem.
	\param Rlocations an R-matrix containing the location of the observations.
	\param Rmesh an R-object containg the output mesh from Trilibrary
	\param Rorder an R-integer containing the order of the approximating basis.
	\param RincidenceMatrix an R-matrix representing the incidence matrix defining regions in the model with areal data
	\param Rmydim an R-integer containing the dimension of the problem we are considering.
	\param Rndim an R-integer containing the dimension of the space in which the location are.
	\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param RnPC an R-integer specifying the number of principal components to compute.
	\param Rvalidation an R-string containing the method to use for the cross-validation of the penalization term lambda.
	\param RnFolds an R-integer specifying the number of folds to use if K-Fold cross validation method is chosen.
	\param RGCVmethod an R-integer specifying if the GCV computation has to be exact(if = 1) or stochastic (if = 2).
	\param Rnrealizations an R-integer specifying the number of realizations to use when computing the GCV stochastically.
	\param Rsearch an R-integer to decide the search algorithm type (tree or naive or walking search algorithm).
	\return R-vector containg the coefficients of the solution
*/
SEXP Smooth_FPCA(SEXP Rlocations, SEXP RbaryLocations, SEXP Rdatamatrix, SEXP Rmesh, SEXP Rorder, SEXP RincidenceMatrix, SEXP Rmydim, SEXP Rndim, SEXP Rlambda, SEXP RnPC, SEXP Rvalidation, SEXP RnFolds, SEXP RGCVmethod, SEXP Rnrealizations, SEXP Rsearch){
//Set data

	FPCAData fPCAdata(Rlocations, RbaryLocations, Rdatamatrix, Rorder, RincidenceMatrix, Rlambda, RnPC, RnFolds, RGCVmethod, Rnrealizations, Rsearch);

	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	std::string validation=CHAR(STRING_ELT(Rvalidation,0));

	if(fPCAdata.getOrder() == 1 && mydim==2 && ndim==2)
		return(FPCA_skeleton<IntegratorTriangleP2, 1, 2, 2>(fPCAdata, Rmesh, validation));
	else if(fPCAdata.getOrder() == 2 && mydim==2 && ndim==2)
		return(FPCA_skeleton<IntegratorTriangleP4, 2, 2, 2>(fPCAdata, Rmesh, validation));
	else if(fPCAdata.getOrder() == 1 && mydim==2 && ndim==3)
		return(FPCA_skeleton<IntegratorTriangleP2, 1, 2, 3>(fPCAdata, Rmesh, validation));
	else if(fPCAdata.getOrder() == 2 && mydim==2 && ndim==3)
		return(FPCA_skeleton<IntegratorTriangleP4, 2, 2, 3>(fPCAdata, Rmesh, validation));
	else if(fPCAdata.getOrder() == 1 && mydim==3 && ndim==3)
		return(FPCA_skeleton<IntegratorTetrahedronP2, 1, 3, 3>(fPCAdata, Rmesh, validation));
	return(NILSXP);
	 }


	 // Density estimation (interface with R)

//! This function manages the various options for DE-PDE algorithm
/*!
	This function is than called from R code.
	\param Rdata an R-matrix containing the data.
	\param Rmesh an R-object containg the output mesh from Trilibrary
	\param Rorder an R-integer containing the order of the approximating basis.
	\param Rmydim an R-integer containing the dimension of the problem we are considering.
	\param Rndim an R-integer containing the dimension of the space in which the location are.
	\param Rfvec an R-vector containing the the initial solution coefficients given by the user.
	\param RheatStep an R-double containing the step for the heat equation initialization.
	\para, RheatIter an R-integer containing the number of iterations to perfrom the heat equation initialization.
	\param Rlambda an R-vector containing the penalization terms.
	\param Rnfolds an R-integer specifying the number of folds for cross validation.
	\param Rnsim an R-integer specifying the number of iterations to use in the optimization algorithm.
	\param RstepProposals an R-vector containing the step parameters useful for the descent algotihm.
	\param Rtol1 an R-double specifying the tolerance to use for the termination criterion based on the percentage differences.
	\param Rtol2 an R-double specifying the tolerance to use for the termination criterion based on the norm of the gradient.
	\param Rprint and R-integer specifying if print on console.
	\param RstepMethod an R-string containing the method to use to choose the step in the optimization algorithm.
	\param RdirectionMethod an R-string containing the descent direction to use in the optimization algorithm.
	\param RpreprocessMethod an R-string containing the cross-validation method to use.
	\param Rsearch an R-integer to decide the search algorithm type (tree or naive search algorithm).

	\return R-list containg solutions.
*/

SEXP Density_Estimation(SEXP Rdata, SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda,
	 SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals, SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP RstepMethod, SEXP RdirectionMethod, SEXP RpreprocessMethod, SEXP Rsearch)
{
	UInt order= INTEGER(Rorder)[0];
  UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	std::string step_method=CHAR(STRING_ELT(RstepMethod, 0));
	std::string direction_method=CHAR(STRING_ELT(RdirectionMethod, 0));
	std::string preprocess_method=CHAR(STRING_ELT(RpreprocessMethod, 0));

  if(order== 1 && mydim==2 && ndim==2)
		return(DE_skeleton<IntegratorTriangleP2, IntegratorGaussTriangle3, 1, 2, 2>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, step_method, direction_method, preprocess_method));
	else if(order== 2 && mydim==2 && ndim==2)
		return(DE_skeleton<IntegratorTriangleP4, IntegratorGaussTriangle3, 2, 2, 2>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, step_method, direction_method, preprocess_method));
	else if(order== 1 && mydim==2 && ndim==3)
		return(DE_skeleton<IntegratorTriangleP2, IntegratorGaussTriangle3, 1, 2, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, step_method, direction_method, preprocess_method));
	else if(order== 2 && mydim==2 && ndim==3)
		return(DE_skeleton<IntegratorTriangleP4, IntegratorGaussTriangle3, 2, 2, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, step_method, direction_method, preprocess_method));
	else if(order == 1 && mydim==3 && ndim==3)
		return(DE_skeleton<IntegratorTetrahedronP2, IntegratorGaussTetra3, 1, 3, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, step_method, direction_method, preprocess_method));
	// else if(order == 1 && mydim==3 && ndim==3)
	// 	return(DE_skeleton<IntegratorTetrahedronP2, IntegratorGaussTetra3, 1, 3, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, step_method, direction_method, preprocess_method));

	return(NILSXP);
}
  
  
  SEXP Density_Initialization(SEXP Rdata, SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda,
	 SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals, SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP Rsearch, SEXP Rinit, SEXP Rinit_fold)
{
	UInt order= INTEGER(Rorder)[0];
  UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	UInt init_fold=INTEGER(Rinit_fold)[0];

	std::string init=CHAR(STRING_ELT(Rinit, 0));

  if(order== 1 && mydim==2 && ndim==2)
		return(DE_init_skeleton<IntegratorTriangleP2, IntegratorGaussTriangle3, 1, 2, 2>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, init, init_fold));
	else if(order== 2 && mydim==2 && ndim==2)
		return(DE_init_skeleton<IntegratorTriangleP4, IntegratorGaussTriangle3, 2, 2, 2>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, init, init_fold));
	else if(order== 1 && mydim==2 && ndim==3)
		return(DE_init_skeleton<IntegratorTriangleP2, IntegratorGaussTriangle3, 1, 2, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, init, init_fold));
	else if(order== 2 && mydim==2 && ndim==3)
		return(DE_init_skeleton<IntegratorTriangleP4, IntegratorGaussTriangle3, 2, 2, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, init, init_fold));
	else if(order == 1 && mydim==3 && ndim==3)
		return(DE_init_skeleton<IntegratorTetrahedronP2, IntegratorGaussTetra3, 1, 3, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, init, init_fold));
	// else if(order == 1 && mydim==3 && ndim==3)
	// 	return(DE_init_skeleton<IntegratorTetrahedronP2, IntegratorTetrahedronP2, 1, 3, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, init, init_fold));

	return(NILSXP);
}

/*!
	This function is then called from R code.
	\param Rlocations an R-matrix containing the location of the observations.
	\param Robservations an R-vector containing the values of the observations.
	\param Rmesh an R-object containg the output mesh from Trilibrary.
	\param Rorder an R-integer containing the order of the approximating basis.
	\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param Rcovariates an R-matrix of covariates for the regression model.
	\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data.
	\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
			the other are automatically considered in Neumann Condition.
	\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices.
	\param DOF an R boolean indicating whether dofs of the model have to be computed or not.
	\param GCV an R boolean indicating whether GCV of the model have to be computed or not.
	\param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic).
	\param Rnrealizations the number of random points used in the stochastic computation of the dofs.
	\param Rfamily Denotes the distribution of the data, within the exponential family.
	\param Rmax_num_iteration Maximum number of steps run in the PIRLS algorithm, set to 15 by default.
	\param Rtreshold an R-double used for arresting FPIRLS algorithm. Algorithm stops when two successive iterations lead to improvement in penalized log-likelihood smaller than threshold.
	\param Rtune It is usually set to 1, but may be higher. It gives more weight to the equivalent degrees of freedom in the computation of the value of the GCV.
	\param Rmu0 Initial value of the mean (natural parameter). There exists a default value for each familiy
	\param RscaleParam If necessary and not supplied, the scale parameter \phi is estimated. See method.phi for details.
	\param Rsearch an R-integer to decide the search algorithm type (tree or naive search algorithm).
	\param Rtune a R-double, Tuning parameter used for the estimation of GCV. called 'GCV.inflation.factor' in R code.
	\param RarealDataAvg an R boolean indicating whether the areal data are averaged or not.

	\return R-vector containg the outputs.
*/


 SEXP gam_Laplace(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
  					SEXP Rlambda, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues,
  					SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations , SEXP Rfamily, SEXP Rmax_num_iteration, SEXP Rtreshold, SEXP Rtune, SEXP Rmu0, SEXP RscaleParam, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch, SEXP RarealDataAvg )
{
    // set up the GAMdata structure for the laplacian case
	GAMDataLaplace regressionData(Rlocations, RbaryLocations, Robservations, Rorder, Rlambda, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, GCV, RGCVmethod, Rnrealizations, DOF, RDOF_matrix, Rsearch, Rmax_num_iteration, Rtreshold, Rtune, RarealDataAvg);

 	UInt mydim=INTEGER(Rmydim)[0];// set the mesh dimension form R to C++
	UInt ndim=INTEGER(Rndim)[0];// set the mesh space dimension form R to C++


  	std::string family = CHAR(STRING_ELT(Rfamily,0));

    if(regressionData.getOrder()==1 && mydim==2 && ndim==2)
    	return(GAM_skeleton<GAMDataLaplace,IntegratorTriangleP2, 1, 2, 2>(regressionData, Rmesh, Rmu0 , family, RscaleParam));
    else if(regressionData.getOrder()==2 && mydim==2 && ndim==2)
		return(GAM_skeleton<GAMDataLaplace,IntegratorTriangleP4, 2, 2, 2>(regressionData, Rmesh, Rmu0, family, RscaleParam));
    else if(regressionData.getOrder()==1 && mydim==2 && ndim==3)
    	return(GAM_skeleton<GAMDataLaplace,IntegratorTriangleP2, 1, 2, 3>(regressionData, Rmesh, Rmu0, family, RscaleParam));
   	else if(regressionData.getOrder()==2 && mydim==2 && ndim==3)
   		return(GAM_skeleton<GAMDataLaplace,IntegratorTriangleP4, 2, 2, 3>(regressionData, Rmesh, Rmu0, family, RscaleParam));
	else if(regressionData.getOrder()==1 && mydim==3 && ndim==3)
		return(GAM_skeleton<GAMDataLaplace,IntegratorTetrahedronP2, 1, 3, 3>(regressionData, Rmesh, Rmu0, family, RscaleParam));
    return(R_NilValue);
}


/*!
	This function is then called from R code.
	\param Rlocations an R-matrix containing the location of the observations.
	\param Robservations an R-vector containing the values of the observations.
	\param Rmesh an R-object containg the output mesh from Trilibrary.
	\param Rorder an R-integer containing the order of the approximating basis.
	\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param RK an R object representing the diffusivity tensor of the model
	\param Rbeta an R object representing the advection function of the model
	\param Rc an R object representing the reaction function of the model
	\param Rcovariates an R-matrix of covariates for the regression model.
	\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data.
	\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
			the other are automatically considered in Neumann Condition.
	\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices.
	\param DOF an R boolean indicating whether dofs of the model have to be computed or not.
	\param GCV an R boolean indicating whether GCV of the model have to be computed or not.
	\param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic).
	\param Rnrealizations the number of random points used in the stochastic computation of the dofs.
	\param Rfamily Denotes the distribution of the data, within the exponential family.
	\param Rmax_num_iteration Maximum number of steps run in the PIRLS algorithm, set to 15 by default.
	\param Rtreshold an R-double used for arresting FPIRLS algorithm. Algorithm stops when two successive iterations lead to improvement in penalized log-likelihood smaller than threshold.
	\param Rmu0 Initial value of the mean (natural parameter). There exists a default value for each familiy
	\param RscaleParam If necessary and not supplied, the scale parameter \phi is estimated. See method.phi for details.
	\param Rsearch an R-integer to decide the search algorithm type (tree or naive search algorithm).
	\param Rtune a R-double, Tuning parameter used for the estimation of GCV. called 'GCV.inflation.factor' in R code.
	\param RarealDataAvg an R boolean indicating whether the areal data are averaged or not.

	\return R-vector containg the outputs.
*/

  SEXP gam_PDE(SEXP Rlocations, SEXP RbaryLocations ,SEXP Robservations, SEXP Rmesh, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
  					SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues,
  					SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations , SEXP Rfamily, SEXP Rmax_num_iteration, SEXP Rtreshold, SEXP Rtune, SEXP Rmu0, SEXP RscaleParam, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch, SEXP RarealDataAvg)
{
    
	GAMDataElliptic regressionData(Rlocations, RbaryLocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, GCV, RGCVmethod, Rnrealizations, DOF, RDOF_matrix, Rsearch, Rmax_num_iteration, Rtreshold, Rtune, RarealDataAvg);

	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];


  	std::string family = CHAR(STRING_ELT(Rfamily,0));

    if(regressionData.getOrder()==1 && mydim==2 && ndim==2)
    	return(GAM_skeleton<GAMDataElliptic,IntegratorTriangleP2, 1, 2, 2>(regressionData, Rmesh, Rmu0 , family, RscaleParam));
    else if(regressionData.getOrder()==2 && mydim==2 && ndim==2)
		return(GAM_skeleton<GAMDataElliptic,IntegratorTriangleP4, 2, 2, 2>(regressionData, Rmesh, Rmu0, family, RscaleParam));

    return(R_NilValue);
}

/*!
	This function is then called from R code.
	\param Rlocations an R-matrix containing the location of the observations.
	\param Robservations an R-vector containing the values of the observations.
	\param Rmesh an R-object containg the output mesh from Trilibrary.
	\param Rorder an R-integer containing the order of the approximating basis.
	\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param RK an R object representing the diffusivity tensor of the model
	\param Rbeta an R object representing the advection function of the model
	\param Rc an R object representing the reaction function of the model
	\param Ru an R object representing the forcing function of the model
	\param Rcovariates an R-matrix of covariates for the regression model.
	\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data.
	\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
			the other are automatically considered in Neumann Condition.
	\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices.
	\param DOF an R boolean indicating whether dofs of the model have to be computed or not.
	\param GCV an R boolean indicating whether GCV of the model have to be computed or not.
	\param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic).
	\param Rnrealizations the number of random points used in the stochastic computation of the dofs.
	\param Rfamily Denotes the distribution of the data, within the exponential family.
	\param Rmax_num_iteration Maximum number of steps run in the PIRLS algorithm, set to 15 by default.
	\param Rtreshold an R-double used for arresting FPIRLS algorithm. Algorithm stops when two successive iterations lead to improvement in penalized log-likelihood smaller than threshold.
	\param Rmu0 Initial value of the mean (natural parameter). There exists a default value for each familiy
	\param RscaleParam If necessary and not supplied, the scale parameter \phi is estimated. See method.phi for details.
	\param Rsearch an R-integer to decide the search algorithm type (tree or naive search algorithm).
	\param Rtune a R-double, Tuning parameter used for the estimation of GCV. called 'GCV.inflation.factor' in R code.
	\param RarealDataAvg an R boolean indicating whether the areal data are averaged or not.

	\return R-vector containg the outputs.
*/

  SEXP gam_PDE_space_varying(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
  					SEXP Rlambda,SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues,
  					SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations , SEXP Rfamily, SEXP Rmax_num_iteration, SEXP Rtreshold, SEXP Rtune, SEXP Rmu0, SEXP RscaleParam, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch, SEXP RarealDataAvg )
{
    
	GAMDataEllipticSpaceVarying regressionData(Rlocations, RbaryLocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Ru, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, GCV,  RGCVmethod, Rnrealizations, DOF, RDOF_matrix, Rsearch, Rmax_num_iteration, Rtreshold, Rtune, RarealDataAvg);

	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];


  	std::string family = CHAR(STRING_ELT(Rfamily,0));

    if(regressionData.getOrder()==1 && mydim==2 && ndim==2)
    	return(GAM_skeleton<GAMDataEllipticSpaceVarying,IntegratorTriangleP2, 1, 2, 2>(regressionData, Rmesh, Rmu0 , family, RscaleParam));
    else if(regressionData.getOrder()==2 && mydim==2 && ndim==2)
		return(GAM_skeleton<GAMDataEllipticSpaceVarying,IntegratorTriangleP4, 2, 2, 2>(regressionData, Rmesh, Rmu0, family, RscaleParam));
    return(R_NilValue);
}


}

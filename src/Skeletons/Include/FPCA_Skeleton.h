#ifndef __FPCA_SKELETON_H__
#define __FPCA_SKELETON_H__

#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../FPCA/Include/Mixed_FE_FPCA.h"
#include "../../FPCA/Include/Mixed_FE_FPCA_Factory.h"

template<UInt ORDER, UInt mydim, UInt ndim>
SEXP FPCA_skeleton(FPCAData &fPCAData, SEXP Rmesh, std::string validation)
{

	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh, fPCAData.getSearch());

	std::unique_ptr<MixedFEFPCABase> fpca = MixedFEFPCAfactory::createFPCAsolver(validation, fPCAData);

	fpca->template SetAndFixParameters<ORDER,mydim,ndim>(mesh);
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

	if(fPCAData.getSearch()==2){
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

#endif

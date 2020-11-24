/*
 * FEMeval.cpp
 *
 *  Created on: Aug 16, 2015
 *      Author: eardi
 */


#include "../../FdaPDE.h"
//#include "IO_handler.h"
#include "../../FE_Assemblers_Solvers/Include/Spline.h"
#include "../../Mesh/Include/Mesh_Objects.h"
#include "../../Mesh/Include/Mesh.h"
#include "../Include/Evaluator.h"
#include "../Include/Projection.h"

template<UInt ORDER, UInt mydim, UInt ndim>
SEXP tree_mesh_skeleton(SEXP Rmesh) {
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh, 2);

	//Copy result in R memory
	SEXP result = NILSXP;
	result = PROTECT(Rf_allocVector(VECSXP, 5));


	//SEND TREE INFORMATION TO R
	SET_VECTOR_ELT(result, 0, Rf_allocVector(INTSXP, 1)); //tree_header information
	int *rans = INTEGER(VECTOR_ELT(result, 0));
	rans[0] = mesh.getTree().gettreeheader().gettreelev();

	SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain origin
	Real *rans1 = REAL(VECTOR_ELT(result, 1));
	for(UInt i = 0; i < ndim*2; i++)
		rans1[i] = mesh.getTree().gettreeheader().domainorig(i);

	SET_VECTOR_ELT(result, 2, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain scale
	Real *rans2 = REAL(VECTOR_ELT(result, 2));
	for(UInt i = 0; i < ndim*2; i++)
		rans2[i] = mesh.getTree().gettreeheader().domainscal(i);


	UInt num_tree_nodes = mesh.num_elements()+1; //Be careful! This is not equal to number of elements
	SET_VECTOR_ELT(result, 3, Rf_allocMatrix(INTSXP, num_tree_nodes, 3)); //treenode information
	int *rans3 = INTEGER(VECTOR_ELT(result, 3));
	for(UInt i = 0; i < num_tree_nodes; i++)
		rans3[i] = mesh.getTree().gettreenode(i).getid();

	for(UInt i = 0; i < num_tree_nodes; i++)
		rans3[i + num_tree_nodes*1] = mesh.getTree().gettreenode(i).getchild(0);

	for(UInt i = 0; i < num_tree_nodes; i++)
		rans3[i + num_tree_nodes*2] = mesh.getTree().gettreenode(i).getchild(1);

	SET_VECTOR_ELT(result, 4, Rf_allocMatrix(REALSXP, num_tree_nodes, ndim*2)); //treenode box coordinate
	Real *rans4 = REAL(VECTOR_ELT(result, 4));
	for(UInt j = 0; j < ndim*2; j++)
	{
		for(UInt i = 0; i < num_tree_nodes; i++)
			rans4[i + num_tree_nodes*j] = mesh.getTree().gettreenode(i).getbox().get()[j];
	}


	UNPROTECT(1);
	return(result);
}

SEXP CPP_eval_FEM_fd(SEXP Rmesh, double* X,  double* Y,  double* Z, UInt n_X, UInt** incidenceMatrix, UInt nRegions, UInt nElements, double* coef, UInt order, UInt fast, UInt mydim, UInt ndim, int search, SEXP RbaryLocations)
{
	SEXP result;

	std::vector<UInt> element_id;
	Real **barycenters;

	//RECIEVE BARYCENTER INFORMATION FROM R
	if (TYPEOF(RbaryLocations) != 0) { //have location information
		element_id.assign(INTEGER(VECTOR_ELT(RbaryLocations, 1)), INTEGER(VECTOR_ELT(RbaryLocations, 1))+n_X);
		UInt n_ = INTEGER(Rf_getAttrib(VECTOR_ELT(RbaryLocations, 2), R_DimSymbol))[0]; //barycenter rows (number of locations)
		UInt p_ = INTEGER(Rf_getAttrib(VECTOR_ELT(RbaryLocations, 2), R_DimSymbol))[1]; //barycenter columns (number of vertices)

		barycenters = (Real**) malloc(sizeof(Real*)*n_);
		for (int i=0; i<n_; i++)
		{
			barycenters[i] = (Real*) malloc(sizeof(Real)*p_);
			for (int j=0; j<p_; j++)
			{
				barycenters[i][j] = REAL(VECTOR_ELT(RbaryLocations, 2))[i+n_*j];
			}
		}
	}


	if (n_X>0) //pointwise data
	{
		PROTECT(result = Rf_allocVector(REALSXP, n_X));
		std::vector<bool> isinside(n_X);
		//Set the mesh
		if(order==1 && mydim==2 && ndim==2)
		{
			MeshHandler<1,2,2> mesh(Rmesh, search);
			Evaluator<1,2,2> evaluator(mesh);
			if (TYPEOF(RbaryLocations) == 0) { //doesn't have location information
				evaluator.eval(X, Y, n_X, coef, fast, REAL(result), isinside);
			} else { //have location information
				evaluator.evalWithInfo(X, Y, n_X, coef, fast, REAL(result), isinside, element_id, barycenters);
			}
		}
		else if(order==2 && mydim==2 && ndim==2)
		{
			MeshHandler<2,2,2> mesh(Rmesh, search);
			Evaluator<2,2,2> evaluator(mesh);
			if (TYPEOF(RbaryLocations) == 0) { //doesn't have location information
				evaluator.eval(X, Y, n_X, coef, fast, REAL(result), isinside);
			} else { //have location information
				evaluator.evalWithInfo(X, Y, n_X, coef, fast, REAL(result), isinside, element_id, barycenters);
			}
		}
		else if(order==1 && mydim==2 && ndim==3)
		{
			MeshHandler<1,2,3> mesh(Rmesh, search);
			Evaluator<1,2,3> evaluator(mesh);
			if (TYPEOF(RbaryLocations) == 0) { //doesn't have location information
				evaluator.eval(X, Y, Z, n_X, coef, fast, REAL(result), isinside);
			} else { //have location information
				evaluator.evalWithInfo(X, Y, Z, n_X, coef, fast, REAL(result), isinside, element_id, barycenters);
			}

		}
		else if(order==2 && mydim==2 && ndim==3)
		{
			MeshHandler<2,2,3> mesh(Rmesh, search);
			Evaluator<2,2,3> evaluator(mesh);
			if (TYPEOF(RbaryLocations) == 0) { //doesn't have location information
				evaluator.eval(X, Y, Z, n_X, coef, fast, REAL(result), isinside);
			} else { //have location information
				evaluator.evalWithInfo(X, Y, Z, n_X, coef, fast, REAL(result), isinside, element_id, barycenters);
			}
		}
		else if(order==1 && mydim==3 && ndim==3)
		{
			MeshHandler<1,3,3> mesh(Rmesh, search);
			Evaluator<1,3,3> evaluator(mesh);
			if (TYPEOF(RbaryLocations) == 0) { //doesn't have location information
				evaluator.eval(X, Y, Z, n_X, coef, fast, REAL(result), isinside);
			} else { //have location information
				evaluator.evalWithInfo(X, Y, Z, n_X, coef, fast, REAL(result), isinside, element_id, barycenters);
			}
		}
		else if(order==2 && mydim==3 && ndim==3)
		{
			MeshHandler<2,3,3> mesh(Rmesh, search);
			Evaluator<2,3,3> evaluator(mesh);
			if (TYPEOF(RbaryLocations) == 0) { //doesn't have location information
				evaluator.eval(X, Y, Z, n_X, coef, fast, REAL(result), isinside);
			} else { //have location information
				evaluator.evalWithInfo(X, Y, Z, n_X, coef, fast, REAL(result), isinside, element_id, barycenters);
			}
		}


		for (int i=0; i<n_X; ++i)
		{
			if(!(isinside[i]))
			{
				REAL(result)[i]=NA_REAL;
			}
		}
	}
	else //areal data
	{
		PROTECT(result = Rf_allocVector(REALSXP, nRegions));
		if(order==1 && mydim==2 && ndim==2)
		{
			MeshHandler<1,2,2> mesh(Rmesh);
			Evaluator<1,2,2> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==2 && mydim==2 && ndim==2)
		{
			MeshHandler<2,2,2> mesh(Rmesh);
			Evaluator<2,2,2> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==1 && mydim==2 && ndim==3)
		{
			MeshHandler<1,2,3> mesh(Rmesh);
			Evaluator<1,2,3> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==2 && mydim==2 && ndim==3)
		{
			MeshHandler<2,2,3> mesh(Rmesh);
			Evaluator<2,2,3> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==1 && mydim==3 && ndim==3)
		{
			MeshHandler<1,3,3> mesh(Rmesh);
			Evaluator<1,3,3> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));

		}
		else if(order==2 && mydim==3 && ndim==3)
		{
			MeshHandler<2,3,3> mesh(Rmesh);
			Evaluator<2,3,3> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));

		}

	}

	UNPROTECT(1);
    // result list
	return(result);
}

extern "C" {
//! This function manages the various option for the solution evaluation.
/*!
	This function is then called from R code.
	Calls the walking algoritm for efficient point location inside the mesh in 2D.

	\param Rmesh an R-object containg the output mesh from Trilibrary
	\param Rlocations an R-matrix (seen as an array) containing the xyz coordinates of the points where the solution has to be evaluated
	\param RincidenceMatrix an R-matrix for the incidence matrix defining the regions in the case of areal data
	\param Rcoef an R-vector the coefficients of the solution
	\param Rorder an R integer containg the order of the solution
	\param Rfast an R integer 0 for Naive location algorithm, 1 for Walking Algorithm (can miss location for non convex meshes)
*/


SEXP eval_FEM_fd(SEXP Rmesh, SEXP Rlocations, SEXP RincidenceMatrix, SEXP Rcoef, SEXP Rorder, SEXP Rfast, SEXP Rmydim, SEXP Rndim, SEXP Rsearch, SEXP RbaryLocations)
{
	int n_X = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[0];
	int nRegions = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[0];
	int nElements = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[1]; //number of triangles/tetrahedron if areal data

	std::vector<UInt> element_id;
	Real **barycenters;

	//RECIEVE BARYCENTER INFORMATION FROM R
	if (TYPEOF(RbaryLocations) != 0) { //have location information
		element_id.assign(INTEGER(VECTOR_ELT(RbaryLocations, 1)), INTEGER(VECTOR_ELT(RbaryLocations, 1))+n_X);
		UInt n_ = INTEGER(Rf_getAttrib(VECTOR_ELT(RbaryLocations, 2), R_DimSymbol))[0]; //barycenter rows (number of locations)
		UInt p_ = INTEGER(Rf_getAttrib(VECTOR_ELT(RbaryLocations, 2), R_DimSymbol))[1]; //barycenter columns (number of vertices)

		barycenters = (Real**) malloc(sizeof(Real*)*n_);
		for (int i=0; i<n_; i++)
		{
			barycenters[i] = (Real*) malloc(sizeof(Real)*p_);
			for (int j=0; j<p_; j++)
			{
				barycenters[i][j] = REAL(VECTOR_ELT(RbaryLocations, 2))[i+n_*j];
			}
		}
	}

	//Declare pointer to access data from C++
	double *X, *Y, *Z;
	UInt **incidenceMatrix;
	double *coef;
	int order, mydim, ndim, search;
	bool fast;

	coef  = REAL(Rcoef);
	order = INTEGER(Rorder)[0];
	fast  = INTEGER(Rfast)[0];
	mydim = INTEGER(Rmydim)[0];
	ndim  = INTEGER(Rndim)[0];
	search  = INTEGER(Rsearch)[0];

	incidenceMatrix = (UInt**) malloc(sizeof(UInt*)*nRegions);

    // Cast all computation parameters
	if (ndim==3)
	{
		X = REAL(Rlocations);
		Y = REAL(Rlocations)+n_X;
		Z = REAL(Rlocations)+2*n_X;
	}
	else //ndim==2
	{
		X = REAL(Rlocations);
		Y = REAL(Rlocations)+n_X;
	}
	for (int i=0; i<nRegions; i++)
	{
		incidenceMatrix[i] = (UInt*) malloc(sizeof(UInt)*nElements);
		for (int j=0; j<nElements; j++)
		{
			incidenceMatrix[i][j] = INTEGER(RincidenceMatrix)[i+nRegions*j];
		}
	}

	SEXP result;

	if (n_X>0) //pointwise data
	{
		PROTECT(result = Rf_allocVector(REALSXP, n_X));
		std::vector<bool> isinside(n_X);
		//Set the mesh
		if(order==1 && mydim==2 && ndim==2)
		{
			MeshHandler<1,2,2> mesh(Rmesh, search);
			Evaluator<1,2,2> evaluator(mesh);
			if (TYPEOF(RbaryLocations) == 0) { //doesn't have location information
				evaluator.eval(X, Y, n_X, coef, fast, REAL(result), isinside);
			} else { //have location information
				evaluator.evalWithInfo(X, Y, n_X, coef, fast, REAL(result), isinside, element_id, barycenters);
			}
		}
		else if(order==2 && mydim==2 && ndim==2)
		{
			MeshHandler<2,2,2> mesh(Rmesh, search);
			Evaluator<2,2,2> evaluator(mesh);
			if (TYPEOF(RbaryLocations) == 0) { //doesn't have location information
				evaluator.eval(X, Y, n_X, coef, fast, REAL(result), isinside);
			} else { //have location information
				evaluator.evalWithInfo(X, Y, n_X, coef, fast, REAL(result), isinside, element_id, barycenters);
			}
		}
		else if(order==1 && mydim==2 && ndim==3)
		{
			MeshHandler<1,2,3> mesh(Rmesh, search);
			Evaluator<1,2,3> evaluator(mesh);
			if (TYPEOF(RbaryLocations) == 0) { //doesn't have location information
				evaluator.eval(X, Y, Z, n_X, coef, fast, REAL(result), isinside);
			} else { //have location information
				evaluator.evalWithInfo(X, Y, Z, n_X, coef, fast, REAL(result), isinside, element_id, barycenters);
			}
		}
		else if(order==2 && mydim==2 && ndim==3)
		{
			MeshHandler<2,2,3> mesh(Rmesh, search);
			Evaluator<2,2,3> evaluator(mesh);
			if (TYPEOF(RbaryLocations) == 0) { //doesn't have location information
				evaluator.eval(X, Y, Z, n_X, coef, fast, REAL(result), isinside);
			} else { //have location information
				evaluator.evalWithInfo(X, Y, Z, n_X, coef, fast, REAL(result), isinside, element_id, barycenters);
			}
		}
		else if(order==1 && mydim==3 && ndim==3)
		{
			MeshHandler<1,3,3> mesh(Rmesh, search);
			Evaluator<1,3,3> evaluator(mesh);
			if (TYPEOF(RbaryLocations) == 0) { //doesn't have location information
				evaluator.eval(X, Y, Z, n_X, coef, fast, REAL(result), isinside);
			} else { //have location information
				evaluator.evalWithInfo(X, Y, Z, n_X, coef, fast, REAL(result), isinside, element_id, barycenters);
			}
		}
		else if(order==2 && mydim==3 && ndim==3)
		{
			MeshHandler<2,3,3> mesh(Rmesh, search);
			Evaluator<2,3,3> evaluator(mesh);
			if (TYPEOF(RbaryLocations) == 0) { //doesn't have location information
				evaluator.eval(X, Y, Z, n_X, coef, fast, REAL(result), isinside);
			} else { //have location information
				evaluator.evalWithInfo(X, Y, Z, n_X, coef, fast, REAL(result), isinside, element_id, barycenters);
			}
		}


		for (int i=0; i<n_X; ++i)
		{
			if(!(isinside[i]))
			{
				REAL(result)[i]=NA_REAL;
			}
		}
	}
	else //areal data
	{
		PROTECT(result = Rf_allocVector(REALSXP, nRegions));
		if(order==1 && mydim==2 && ndim==2)
		{
			MeshHandler<1,2,2> mesh(Rmesh);
			Evaluator<1,2,2> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==2 && mydim==2 && ndim==2)
		{
			MeshHandler<2,2,2> mesh(Rmesh);
			Evaluator<2,2,2> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==1 && mydim==2 && ndim==3)
		{
			MeshHandler<1,2,3> mesh(Rmesh);
			Evaluator<1,2,3> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==2 && mydim==2 && ndim==3)
		{
			MeshHandler<2,2,3> mesh(Rmesh);
			Evaluator<2,2,3> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==1 && mydim==3 && ndim==3)
		{
			MeshHandler<1,3,3> mesh(Rmesh);
			Evaluator<1,3,3> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));

		}
		else if(order==2 && mydim==3 && ndim==3)
		{
			MeshHandler<2,3,3> mesh(Rmesh);
			Evaluator<2,3,3> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));

		}

	}

	for (int i=0; i<nRegions; i++)
	{
		free(incidenceMatrix[i]);
	}
	free(incidenceMatrix);


	UNPROTECT(1);
    // result list
	return(result);
}


//! This function evaluates the solution on a set of given points by evaluating the tensorial basis expansion.
/*!
	This function is then called from R code.
	Calls the walking algoritm for efficient point location inside the mesh in 2D.

	\param Rmesh an R-object containg the output mesh from Trilibrary
  \param Rmesh_time an R-vector containg the time mesh
  \param Rlocations an R-matrix (seen as an array) containing the xyz coordinates of the points where the solution has to be evaluated
  \param Rtime_locations an R-vector (seen as an array) containing the xyz coordinates of the points where the solution has to be evaluated
	\param RincidenceMatrix an R-matrix for the incidence matrix defining the regions in the case of areal data
	\param Rcoef an R-vector the coefficients of the solution
	\param Rorder an R integer containg the order of the solution
	\param Rfast an R integer 0 for Naive location algorithm, 1 for Walking Algorithm (can miss location for non convex meshes)
  \param Rflag_parabolic an R logical (seen as an integer) 1 if parabolic smoothing, 0 otherwise
  \param Rmydim an R integer containing the mesh space size, 2 for 2D and 2.5D, 3 for 3D
  \param Rmydim an R integer containing the space size, 2 for 2D , 3 for 2.5D and 3D
*/
SEXP eval_FEM_time(SEXP Rmesh, SEXP Rmesh_time, SEXP Rlocations, SEXP Rtime_locations, SEXP RincidenceMatrix, SEXP Rcoef, SEXP Rorder, SEXP Rfast, SEXP Rflag_parabolic, SEXP Rmydim, SEXP Rndim, SEXP Rsearch, SEXP RbaryLocations)
{
  UInt mydim = INTEGER(Rmydim)[0];
  UInt ndim  = INTEGER(Rndim)[0];
  UInt n = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[0];
  
  UInt ns = INTEGER(Rf_getAttrib(VECTOR_ELT(Rmesh, 0), R_DimSymbol))[0];

  UInt nt = Rf_length(Rmesh_time);
	UInt nRegions = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[0];
	UInt nElements = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[1]; //number of triangles/tetrahedron if areal data
	//Declare pointer to access data from C++
	Real *X, *Y, *Z, *mesh_time, *t;
	UInt **incidenceMatrix;
	double *coef;
	int order, search;
	bool fast,flag_par;

	coef  = REAL(Rcoef);
  order = INTEGER(Rorder)[0];
  search  = INTEGER(Rsearch)[0];
  fast  = INTEGER(Rfast)[0];
	flag_par = INTEGER(Rflag_parabolic)[0];
	mesh_time = REAL(Rmesh_time);
	t = REAL(Rtime_locations);

	incidenceMatrix = (UInt**) malloc(sizeof(UInt*)*nRegions);

    // Cast all computation parameters
	if (ndim==3)
	{
		for (int i=0; i<n; i++)
		{
			X = REAL(Rlocations);
			Y = REAL(Rlocations)+n;
			Z = REAL(Rlocations)+2*n;
		}
	}
	else //ndim==2
	{
		for (int i=0; i<n; i++)
		{
			X = REAL(Rlocations);
			Y = REAL(Rlocations)+n;
		}
	}
	for (int i=0; i<nRegions; i++)
	{
		incidenceMatrix[i] = (UInt*) malloc(sizeof(UInt)*nElements);
		for (int j=0; j<nElements; j++)
		{
			incidenceMatrix[i][j] = INTEGER(RincidenceMatrix)[i+nRegions*j];
		}
	}

  // Compute the matrix of temporal basis evaluation in the given points
	UInt DEGREE = flag_par ? 1 : 3;
	UInt M = nt + DEGREE - 1;
	SpMat phi(n,M);
  UInt N = nRegions==0 ? n : nRegions;
	if(flag_par)
	{
		Spline<1,0>spline(mesh_time,nt);
		Real value;
		for (UInt i = 0; i < N; ++i)
		{
			for (UInt j = 0; j < M; ++j)
			{
				value = spline.BasisFunction(j, t[i]);
				if (value!=0)
				{
					phi.coeffRef(i,j) = value;
				}
			}
		}
	}
	else
	{
		Spline<3,2>spline(mesh_time,nt);
		Real value;
		for (UInt i = 0; i < N; ++i)
		{
			for (UInt j = 0; j < M; ++j)
			{
				value = spline.BasisFunction(j, t[i]);
				if (value!=0)
				{
					phi.coeffRef(i,j) = value;
				}
			}
		}
	}
	phi.makeCompressed();

	SEXP result;

  PROTECT(result=Rf_allocVector(REALSXP, N));
	Real* COEFF;
	COEFF = (double*) malloc(sizeof(double)*ns);
	std::vector<Real> XX,YY,ZZ;
  std::vector<UInt> indices;

  //!evaluates the solution on the given points location at the first
  //!node of the time mesh to initialize the array of results and retrieve the points out of mesh (NA)
	for(UInt j=0; j<ns; ++j)
	{
		COEFF[j] = coef[j];
	}

	SEXP temp = CPP_eval_FEM_fd(Rmesh, X, Y, Z, n, incidenceMatrix, nRegions, nElements, COEFF, order, fast, mydim, ndim, search, RbaryLocations);
	for(UInt k=0; k < N; k++)
	{
		REAL(result)[k] = REAL(temp)[k];
		if(!ISNA(REAL(result)[k]))
			REAL(result)[k] = REAL(result)[k]*phi.coeff(k,0);
	}

  //! loop over time b-splines basis and evaluate the solution only on the points that have
  //! the coefficient corresponding to that basis different from 0
	for(UInt i=1; i<M; ++i)
	{
		for(UInt j=0; j<ns; ++j)
		{
			COEFF[j] = coef[i*ns+j];
		}
    if (ndim==3)
    {
  		for(UInt k=0; k<n; k++)
  		{
  			if(phi.coeff(k,i)!=0 && !ISNA(REAL(result)[k]))
  			{
  				XX.push_back(X[k]);
  				YY.push_back(Y[k]);
  				ZZ.push_back(Z[k]);
  				indices.push_back(k);
  			}
      }
    }
  	else //ndim==2
    {
  		for(UInt k=0; k<n; k++)
  		{
  			if(phi.coeff(k,i)!=0 && !ISNA(REAL(result)[k]))
				{
					XX.push_back(X[k]);
					YY.push_back(Y[k]);
					indices.push_back(k);
				}
			}
		}
    UInt count=0;
    UInt **INCIDENCE_MATRIX;
    INCIDENCE_MATRIX = (UInt**)malloc(sizeof(UInt*)*phi.col(i).nonZeros());

  	for (UInt k=0; k<nRegions; k++)
  	{
      if(phi.coeff(k,i)!=0 && !ISNA(REAL(result)[k]))
      {
        INCIDENCE_MATRIX[count] = (UInt*) malloc(sizeof(UInt)*nElements);
    		for (UInt j=0; j<nElements; j++)
    		{
    			INCIDENCE_MATRIX[count][j] = incidenceMatrix[k][j];
    		}
        indices.push_back(k);
        count++;
      }
    }
		temp = CPP_eval_FEM_fd(Rmesh, XX.data(), YY.data(), ZZ.data(), XX.size(), INCIDENCE_MATRIX, phi.col(i).nonZeros(), nElements, COEFF, order, fast, mydim, ndim, search, RbaryLocations);
		for(UInt k=0; k<indices.size(); ++k)
		{
			REAL(result)[indices[k]] = REAL(result)[indices[k]] + REAL(temp)[k]*phi.coeff(indices[k],i);
		}
		XX.clear();YY.clear();ZZ.clear();indices.clear();

    if(nRegions!=0)
    {
      for (int l=0; l<phi.col(i).nonZeros(); l++)
      {
        free(INCIDENCE_MATRIX[l]);
      }
    }
    free(INCIDENCE_MATRIX);
	}

	free(COEFF);
	for (int i=0; i<nRegions; i++)
	{
		free(incidenceMatrix[i]);
	}
	free(incidenceMatrix);

	UNPROTECT(1);
	return(result);
}



//! This function evaluates the solution on the mesh nodes at a given time with the purpose of plotting it.
/*!
	This function is then called from R code.

  \param Rns an R integer containing the number of mesh nodes
	\param Rmesh_time an R-vector containg the time mesh
	\param Rtime an R double containing the time at which the solution has to be evaluated
	\param Rcoef an R-vector the coefficients of the solution
  \param Rflag_parabolic an R logical TRUE for parabolic smoothing, FALSE otherwise
*/
  SEXP eval_FEM_time_nodes(SEXP Rns, SEXP Rmesh_time, SEXP Rtime, SEXP Rcoef, SEXP Rflag_parabolic)
  {
  	UInt ns = INTEGER(Rns)[0];
  	UInt nt = Rf_length(Rmesh_time);
  	UInt n = Rf_length(Rtime);

  	Real *mesh_time = REAL(Rmesh_time);
  	Real *t = REAL(Rtime);
  	bool flag_par = INTEGER(Rflag_parabolic)[0];

  	UInt DEGREE = flag_par ? 1 : 3;
  	UInt M = nt + DEGREE - 1;
  	MatrixXr phi(M,n);

  	if(flag_par)
  	{
  		Spline<1,0>spline(mesh_time,nt);
  		for (UInt i=0; i < n; ++i)
  		{
  			for (UInt j = 0; j < M; ++j)
  			{
  				phi(j,i) = spline.BasisFunction(j, t[i]);
  			}
  		}
  	}
  	else
  	{
  		Spline<3,2>spline(mesh_time,nt);
  		for (UInt i=0; i < n; ++i)
  		{
  			for (UInt j = 0; j < M; ++j)
  			{
  				phi(j,i) = spline.BasisFunction(j, t[i]);
  			}
  		}
  	}

  	SEXP result;

  	PROTECT(result=Rf_allocVector(REALSXP, ns*n));

  	for(UInt j=0; j<n; ++j)
  	{
  		for(UInt k=0; k<ns; ++k)
  		{
  			REAL(result)[k+j*ns] = REAL(Rcoef)[k]*phi(0,j);
  		}
  	}
  	for(UInt i=1; i < M; i++)
  	{
  		for(UInt j=0; j<n; ++j)
  		{
  			if(phi(i,j)!=0)
  			{
  				for(UInt k=0; k<ns; ++k)
  				{
  					REAL(result)[k+j*ns] = REAL(result)[k+j*ns] + REAL(Rcoef)[k+ns*i]*phi(i,j);
  				}
  			}
  		}
  	}

  	UNPROTECT(1);
  	return(result);
  }

  SEXP points_projection(SEXP Rmesh, SEXP Rlocations)
  {
  	int n_X = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[0];
	//Declare pointer to access data from C++
  	double X, Y, Z;

    // Cast all computation parameters
    std::vector<Point<3> > deData_(n_X); // the points to be projected
    std::vector<Point<3> > prjData_(n_X); // the projected points

    //RECIEVE PROJECTION INFORMATION FROM R
    for (int i=0; i<n_X; i++)
    {
    	X = REAL(Rlocations)[i + n_X*0];
    	Y = REAL(Rlocations)[i + n_X*1];
    	Z = REAL(Rlocations)[i + n_X*2];
    	deData_[i]=Point<3>({X,Y,Z});
    }

    SEXP result;

	if (n_X>0) //pointwise data
	{
		PROTECT(result = Rf_allocMatrix(REALSXP, n_X, 3));
		UInt order = INTEGER(VECTOR_ELT(Rmesh,10))[0];

		if (order == 1) {
			MeshHandler<1,2,3> mesh(Rmesh);
			projection<1,2,3> projector(mesh, deData_);
			prjData_ = projector.computeProjection();
		}

		// if (order == 2) {
		// MeshHandler<2,2,3> mesh(Rmesh);
		// projection<2,2,3> projector(mesh, deData_);
		// prjData_ = projector.computeProjection();
		// }
	}

	for (int i=0; i<n_X; ++i)
	{
		REAL(result)[i + n_X*0]=prjData_[i][0];
		REAL(result)[i + n_X*1]=prjData_[i][1];
		REAL(result)[i + n_X*2]=prjData_[i][2];
	}

	UNPROTECT(1);
    // result matrix
	return(result);
}

SEXP tree_mesh_construction(SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim) {
	UInt ORDER=INTEGER(Rorder)[0];
	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	if(ORDER == 1 && mydim==2 && ndim==2)
		return(tree_mesh_skeleton<1, 2, 2>(Rmesh));
	else if(ORDER == 2 && mydim==2 && ndim==2)
		return(tree_mesh_skeleton<2, 2, 2>(Rmesh));
	else if(ORDER == 1 && mydim==2 && ndim==3)
		return(tree_mesh_skeleton<1, 2, 3>(Rmesh));
	else if(ORDER == 2 && mydim==2 && ndim==3)
		return(tree_mesh_skeleton<2, 2, 3>(Rmesh));
	else if(ORDER == 1 && mydim==3 && ndim==3)
		return(tree_mesh_skeleton<1, 3, 3>(Rmesh));
	else if(ORDER == 2 && mydim==3 && ndim==3)
		return(tree_mesh_skeleton<2, 3, 3>(Rmesh));

	return(NILSXP);
}

}

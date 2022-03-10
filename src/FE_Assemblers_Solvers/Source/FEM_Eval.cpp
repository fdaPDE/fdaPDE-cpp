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
#include "../Include/FE_Skeleton.h"

extern "C" {
//! This function manages the various option for the solution evaluation.
/*!
	This function is then called from R code.
	Calls the walking algoritm for efficient point location inside the mesh in 2D.

	\param Rmesh an R-object containg the output mesh from Trilibrary in 2D (in 1.5D, 2.5D and 3D R functions can produce a compatible object)
	\param Rlocations an R-matrix (seen as an array) containing the xyz coordinates of the points where the solution has to be evaluated
	\param RincidenceMatrix an R-matrix for the incidence matrix defining the regions in the case of areal data
	\param Rcoef an R-vector the coefficients of the solution
	\param Rorder an R integer containg the order of the solution
	\param Rfast an R integer to enforce verbose search for Walking Algorithm (can miss location for non convex meshes)
	\param Rmydim an R integer containing the mesh space size, 1 for 1.5D, 2 for 2D and 2.5D, 3 for 3D
    	\param Rndim an R integer containing the space size, 2 for 1.5D and 2D , 3 for 2.5D and 3D
    	\param Rsearch an R integer 0 for Naive location algorithm, 1 for Walking Algorithm (can miss location for non convex meshes), 3 for Tree search
*/
SEXP eval_FEM_fd(SEXP Rmesh, SEXP Rlocations, SEXP RincidenceMatrix, SEXP Rcoef, SEXP Rorder, SEXP Rfast, SEXP Rmydim, SEXP Rndim, SEXP Rsearch, SEXP RbaryLocations){

    UInt order = INTEGER(Rorder)[0];
    UInt mydim = INTEGER(Rmydim)[0];
    UInt ndim  = INTEGER(Rndim)[0];

    if(order==1 && mydim==1 && ndim==2)
        return Eval_FEM_fd_skeleton<1,1,2>(Rmesh, Rlocations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
    else if(order==2 && mydim==1 && ndim==2)
        return Eval_FEM_fd_skeleton<2,1,2>(Rmesh, Rlocations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
    else if(order==1 && mydim==2 && ndim==2)
        return Eval_FEM_fd_skeleton<1,2,2>(Rmesh, Rlocations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
    else if(order==2 && mydim==2 && ndim==2)
        return Eval_FEM_fd_skeleton<2,2,2>(Rmesh, Rlocations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
    else if(order==1 && mydim==2 && ndim==3)
        return Eval_FEM_fd_skeleton<1,2,3>(Rmesh, Rlocations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
    else if(order==2 && mydim==2 && ndim==3)
        return Eval_FEM_fd_skeleton<2,2,3>(Rmesh, Rlocations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
    else if(order==1 && mydim==3 && ndim==3)
        return Eval_FEM_fd_skeleton<1,3,3>(Rmesh, Rlocations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
    else if(order==2 && mydim==3 && ndim==3)
        return Eval_FEM_fd_skeleton<2,3,3>(Rmesh, Rlocations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);

    return NILSXP;
}

//! This function evaluates the solution on a set of given points by evaluating the tensorial basis expansion.
/*!
	This function is then called from R code.

	\param Rmesh an R-object containg the output mesh from Trilibrary in 2D (in 1.5D, 2.5D and 3D R functions can produce a compatible object)
  \param Rmesh_time an R-vector containg the time mesh
  \param Rlocations an R-matrix containing the xyz coordinates of the points where the solution has to be evaluated
  \param Rtime_locations an R-vector containing the time coordinates of the points where the solution has to be evaluated
	\param RincidenceMatrix an R-matrix for the incidence matrix defining the regions in the case of areal data
	\param Rcoef an R-vector the coefficients of the solution
	\param Rorder an R integer containg the order of the solution
	\param Rfast an R integer to enforce verbose search for Walking Algorithm (can miss location for non convex meshes)
	\param Rsearch an R integer 0 for Naive location algorithm, 1 for Walking Algorithm (can miss location for non convex meshes), 3 for Tree search
  \param Rflag_parabolic an R logical (seen as an integer) 1 if parabolic smoothing, 0 otherwise
  \param Rmydim an R integer containing the mesh space size, 1 for 1.5D, 2 for 2D and 2.5D, 3 for 3D
  \param Rndim an R integer containing the space size, 2 for 1.5D and 2D, 3 for 2.5D and 3D
*/
SEXP eval_FEM_time(SEXP Rmesh, SEXP Rmesh_time, SEXP Rlocations, SEXP Rtime_locations, SEXP RincidenceMatrix, SEXP Rcoef, SEXP Rorder, SEXP Rfast, SEXP Rflag_parabolic, SEXP Rmydim, SEXP Rndim, SEXP Rsearch, SEXP RbaryLocations)
    {
        UInt order = INTEGER(Rorder)[0];
        UInt mydim = INTEGER(Rmydim)[0];
        UInt ndim  = INTEGER(Rndim)[0];
        UInt flag_par = INTEGER(Rflag_parabolic)[0]; // UInt DEGREE =  flag_par ? 1 : 3

        if(order==1 && mydim==2 && ndim==2 && flag_par==1)
            return Eval_FEM_time_skeleton<1,2,2,1>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==1 && mydim==2 && ndim==2 && flag_par!=1)
            return Eval_FEM_time_skeleton<1,2,2,3>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==2 && mydim==2 && ndim==2 && flag_par==1)
            return Eval_FEM_time_skeleton<2,2,2,1>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==2 && mydim==2 && ndim==2 && flag_par!=1)
            return Eval_FEM_time_skeleton<2,2,2,3>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==1 && mydim==2 && ndim==3 && flag_par==1)
            return Eval_FEM_time_skeleton<1,2,3,1>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==1 && mydim==2 && ndim==3 && flag_par!=1)
            return Eval_FEM_time_skeleton<1,2,3,3>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==2 && mydim==2 && ndim==3 && flag_par==1)
            return Eval_FEM_time_skeleton<2,2,3,1>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==2 && mydim==2 && ndim==2 && flag_par!=1)
            return Eval_FEM_time_skeleton<1,2,3,3>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==1 && mydim==3 && ndim==3 && flag_par==1)
            return Eval_FEM_time_skeleton<1,3,3,1>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==1 && mydim==2 && ndim==2 && flag_par!=1)
            return Eval_FEM_time_skeleton<1,3,3,3>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==2 && mydim==3 && ndim==3 && flag_par==1)
            return Eval_FEM_time_skeleton<2,3,3,1>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==2 && mydim==3 && ndim==3 && flag_par!=1)
            return Eval_FEM_time_skeleton<2,3,3,3>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==1 && mydim==1 && ndim==2 && flag_par==1)
            return Eval_FEM_time_skeleton<1,1,2,1>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==1 && mydim==1 && ndim==2 && flag_par!=1)
            return Eval_FEM_time_skeleton<1,1,2,3>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==2 && mydim==1 && ndim==2 && flag_par==1)
            return Eval_FEM_time_skeleton<2,1,2,1>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);
        else if(order==2 && mydim==1 && ndim==2 && flag_par!=1)
            return Eval_FEM_time_skeleton<2,1,2,3>(Rmesh, Rmesh_time, Rlocations, Rtime_locations, RincidenceMatrix, Rcoef, Rfast, Rsearch, RbaryLocations);

        return NILSXP;
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

//! This function build the binary-tree object in the tree search algorithm.
/*!
	This function is then called from R code.

	\param Rmesh an R-object containg the output mesh from Trilibrary in 2D (in 1.5D, 2.5D and 3D R functions can produce a compatible object)
	\param Rorder an R integer representing the order of the element of the mesh
  	\param Rmydim an R integer containing the mesh space size, 1 for 1.5D, 2 for 2D and 2.5D, 3 for 3D
    \param Rndim an R integer containing the space size, 2 for 1.5D and 2D , 3 for 2.5D and 3D
*/
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
	else if(ORDER == 1 && mydim==1 && ndim==2)
        	return tree_mesh_skeleton<1,1,2>(Rmesh);
    	else if(ORDER == 2 && mydim==1 && ndim==2)
        	return tree_mesh_skeleton<2,1,2>(Rmesh);

	return(NILSXP);
}

//! This function projects the points on the given mesh.
/*!
	This function is then called from R code.

  \param Rmesh an R-object containg the output mesh from Trilibrary in 2D (in 1.5D, 2.5D and 3D R functions can produce a compatible object)
  \param Rlocation an R matrix containing the coordinates of the points where the solution has to be evaluated
  \param Rmydim an R integer containing the mesh space size, 1 for 1.5D, 2 for 2D and 2.5D, 3 for 3D
  \param Rndim an R integer containing the space size, 2 for 1.5D and 2D , 3 for 2.5D and 3D
*/
SEXP points_projection(SEXP Rmesh, SEXP Rlocations,SEXP Rmydim, SEXP Rndim){
    UInt order = INTEGER(VECTOR_ELT(Rmesh,10))[0];
    UInt mydim = INTEGER(Rmydim)[0];
    UInt ndim = INTEGER(Rndim)[0];

    if(order==1 && mydim==1 && ndim==2)
        return points_projection_skeleton<1,1,2>(Rmesh,Rlocations);
    else if(order==2 && mydim==1 && ndim==2)
        return points_projection_skeleton<2,1,2>(Rmesh,Rlocations);
    else if(order==1 && mydim==2 && ndim==3)
        return points_projection_skeleton<1,2,3>(Rmesh,Rlocations);
    else if(order==2 && mydim==2 && ndim==3)
        return points_projection_skeleton<2,2,3>(Rmesh,Rlocations);

    return(NILSXP);
}

}

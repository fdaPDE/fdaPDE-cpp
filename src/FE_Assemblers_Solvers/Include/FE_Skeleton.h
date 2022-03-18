#ifndef __FE_SKELETON_H__
#define __FE_SKELETON_H__

#include "../../FdaPDE.h"
#include "Projection.h"
#include "../../Mesh/Include/Mesh.h"

//! This function build the binary-tree object used in the tree search algorithm according to template the parameter.
/*!
	This function is then called from <tree_mesh_construction>"()" function.

	\param Rmesh an R-object containg the output mesh from Trilibrary in 2D (in 1.5D, 2.5D and 3D R functions can produce a compatible object)
*/
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

//! This function projects the points on the mesh according to the template the parameter.
/*!
	This function is then called from  code.

	\param Rmesh an R-object containg the output mesh from Trilibrary in 2D (in 2.5D and 3D R functions can produce a compatible object if the triangulation is already available)
  	\param Rlocations R matrix containing the coordinates of the points to project
*/
template<UInt ORDER,UInt mydim,UInt ndim>
SEXP points_projection_skeleton(SEXP Rmesh, SEXP Rlocations){
    //RECIEVE PROJECTION INFORMATION FROM R
    RNumericMatrix locations(Rlocations);
    UInt n_X = locations.nrows();

    // Cast all computation parameters
    std::vector<Point<ndim> > deData_(n_X); // the points to be projected
    std::vector<Point<ndim> > prjData_(n_X); // the projected points

    std::array<Real,ndim> coords;
    for(UInt i=0; i<n_X; ++i) {
        for(UInt n=0; n<ndim; ++n)
            coords[n] = locations(i,n);
        deData_[i] = Point<ndim>(coords);

    }
    SEXP result;

    if (n_X>0) //pointwise data
    {
        PROTECT(result = Rf_allocMatrix(REALSXP, n_X, ndim));
        MeshHandler<ORDER,mydim,ndim> mesh(Rmesh);
        projection<ORDER,mydim,ndim> projector(mesh, deData_);
        prjData_ = projector.computeProjection();

        RNumericMatrix res(result);
        for(UInt i=0; i<n_X; ++i){
            for(UInt n=0; n<ndim; ++n)
                res(i,n) = prjData_[i][n];
        }

        UNPROTECT(1);
        return(result);
    }

    return(NILSXP);
}

//! This function manages the solution evaluation on given locations according to template parameters.
/*!
	This function is called from <eval_FEM_fd>"()" function.

	\param Rmesh an R-object containg the output mesh from Trilibrary in 2D (in 1.5D, 2.5D and 3D R functions can produce a compatible object)
	\param Rlocations an R-matrix containing the coordinates of the points where the solution has to be evaluated
	\param RincidenceMatrix an R-matrix for the incidence matrix defining the regions in the case of areal data
	\param Rcoef an R-vector the coefficients of the solution
	\param Rorder an R integer containg the order of the solution
	\param Rfast an R integer to enforce verbose search for Walking Algorithm (can miss location for non convex meshes)
	\param Rsearch an R integer 0 for Naive location algorithm, 1 for Walking Algorithm (can miss location for non convex meshes), 3 for Tree search
 */
template<UInt ORDER, UInt mydim, UInt ndim>
SEXP Eval_FEM_fd_skeleton(SEXP Rmesh, SEXP Rlocations, SEXP RincidenceMatrix, SEXP Rcoef, SEXP Rfast, SEXP Rsearch, SEXP RbaryLocations){

    RNumericMatrix barycenters( VECTOR_ELT(RbaryLocations,2));
    RIntegerMatrix id_element( VECTOR_ELT(RbaryLocations,1));
    RIntegerMatrix incidenceMatrix( RincidenceMatrix );
    RNumericMatrix locations(Rlocations);

    UInt n_X = locations.nrows();
    UInt nRegions = incidenceMatrix.nrows();
    RNumericMatrix coef(Rcoef);
    UInt search;
    bool fast;

    fast  = INTEGER(Rfast)[0];
    search  = INTEGER(Rsearch)[0];
    MeshHandler<ORDER, mydim, ndim> mesh(Rmesh, search);

    Evaluator<ORDER, mydim, ndim> evaluator(mesh);

    SEXP result;

    if(n_X >0) {
        PROTECT(result=Rf_allocMatrix(REALSXP,n_X,1));
        RNumericMatrix result_(result);

        std::vector<bool> isinside(n_X);
        if (barycenters.nrows() == 0) {
            evaluator.eval(locations, coef, fast, result_, isinside);
        } else {
            evaluator.evalWithInfo(locations, coef, fast, result_, isinside, id_element, barycenters);
        }

        for (int i = 0; i < n_X; ++i) {
            if (!(isinside[i])) {
                result_[i] = NA_REAL;
            }
        }
    }
    else{
        PROTECT(result=Rf_allocMatrix(REALSXP, nRegions,1));
        RNumericMatrix result_(result);
        evaluator.integrate(incidenceMatrix,coef,result_);

    }

    UNPROTECT(1);
    return result;

}

//! This function evaluates the solution on a set of given points by evaluating the tensorial basis expansion.
/*!
	This function is then called from <eval_FEM_time>"()" function.

	\param Rmesh an R-object containg the output mesh from Trilibrary in 2D (in 1.5D, 2.5D and 3D R functions can produce a compatible object)
  \param Rmesh_time an R-vector containg the time mesh
  \param Rlocations an R-matrix containing the xyz coordinates of the points where the solution has to be evaluated
  \param Rtime_locations an R-vector containing the coordinates of the points where the solution has to be evaluated
	\param RincidenceMatrix an R-matrix for the incidence matrix defining the regions in the case of areal data
	\param Rcoef an R-vector the coefficients of the solution
	\param Rorder an R integer containg the order of the solution
	\param Rfast an R integer to enforce verbose search for Walking Algorithm (can miss location for non convex meshes)
	\param Rsearch an R integer 0 for Naive location algorithm, 1 for Walking Algorithm (can miss location for non convex meshes), 3 for Tree search
  \param Rflag_parabolic an R logical (seen as an integer) 1 if parabolic smoothing, 0 otherwise
*/
template<UInt ORDER, UInt mydim, UInt ndim, UInt DEGREE>
SEXP Eval_FEM_time_skeleton (SEXP Rmesh, SEXP Rmesh_time, SEXP Rlocations, SEXP Rtime_locations, SEXP RincidenceMatrix, SEXP Rcoef, SEXP Rfast, SEXP Rsearch, SEXP RbaryLocations)
{
    RNumericMatrix coef(Rcoef);
    RNumericMatrix locations(Rlocations);
    RIntegerMatrix incidenceMatrix(RincidenceMatrix);

    UInt n_X = locations.nrows();
    UInt n_nodes = INTEGER(Rf_getAttrib(VECTOR_ELT(Rmesh, 0), R_DimSymbol))[0];
    UInt n_t = Rf_length(Rmesh_time);
    UInt nRegions = incidenceMatrix.nrows();
    UInt nElements = incidenceMatrix.ncols();

    Real *mesh_time, *t;
    mesh_time = REAL(Rmesh_time);
    t = REAL(Rtime_locations);

    UInt M = n_t + DEGREE - 1;
    UInt N = nRegions==0 ? n_X : nRegions;
    SpMat phi(N,M);
    Spline<DEGREE,DEGREE-1>spline(mesh_time,n_t);
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
    phi.makeCompressed();

    SEXP result;
    PROTECT(result=Rf_allocVector(REALSXP, N));  // 1

    SEXP Rcoef_0;
    PROTECT(Rcoef_0=Rf_allocMatrix(REALSXP, n_nodes,1)); // 2
    RNumericMatrix coef_0(Rcoef_0);

    //!evaluates the solution on the given points location at the first
    //!node of the time mesh to initialize the array of results and retrieve the points out of mesh (NA)
    for(UInt j=0; j<n_nodes; ++j)
    {
        coef_0[j] = coef[j];
    }

    SEXP temp = Eval_FEM_fd_skeleton<ORDER,mydim,ndim>(Rmesh,Rlocations, RincidenceMatrix,Rcoef_0, Rfast,Rsearch, RbaryLocations);
    UNPROTECT(1);

    for(UInt k=0; k < N; k++) {
        REAL(result)[k] = REAL(temp)[k];
        if (!ISNA(REAL(result)[k]))
            REAL(result)[k] = REAL(result)[k] * phi.coeff(k, 0);
    }

    //! loop over time b-splines basis and evaluate the solution only on the points that have
    //! the coefficient corresponding to that basis different from 0
    std::vector<UInt> indices;
    for(UInt i=1; i<M; ++i)
    {

        SEXP Rcoef_i;
        PROTECT(Rcoef_i=Rf_allocMatrix(REALSXP, n_nodes,1)); //1
        RNumericMatrix coef_i(Rcoef_i);
        for(UInt j=0; j<n_nodes; ++j)
        {
            coef_i[j] = coef[i*n_nodes+j];
        }
        // if n_X > 0
        UInt n_i = 0;
        for(UInt k=0; k<n_X; k++) {
            if (phi.coeff(k, i) != 0 && !ISNA(REAL(result)[k])) {
                ++n_i;
                indices.push_back(k);
            }
        }

        SEXP Rlocations_i;
        PROTECT(Rlocations_i = Rf_allocMatrix(REALSXP,n_i,ndim)); //2
        RNumericMatrix locations_i(Rlocations_i);
        n_i = 0;
        for(UInt k=0; k<indices.size(); ++k){
            for(UInt l=0; l<ndim;++l)
                locations_i(n_i,l) = locations(indices[k],l);
            ++n_i;
        }

        SEXP RincidenceMatrix_i;
        PROTECT(RincidenceMatrix_i=Rf_allocMatrix(INTSXP,phi.col(i).nonZeros(), nElements)); //3
        RIntegerMatrix incidenceMatrix_i(RincidenceMatrix_i);
        UInt nRegions_i = 0;
        for (UInt k=0; k<nRegions; ++k)
        {
            if(phi.coeff(k,i)!=0 && !ISNA(REAL(result)[k]))
            {
                for (UInt j=0; j<nElements; j++)
                {
                    incidenceMatrix_i(nRegions_i , j) = incidenceMatrix(k,j);
                }
                indices.push_back(k);
                ++nRegions_i;
            }
        }
        temp = Eval_FEM_fd_skeleton<ORDER,mydim,ndim>(Rmesh,Rlocations_i, RincidenceMatrix_i,Rcoef_i, Rfast,Rsearch, RbaryLocations);
        UNPROTECT(3); // UNPROTECT Rcoef_i, Rlocations_i RincidenceMatrix_i

        for(UInt k=0; k<indices.size(); ++k)
        {
            REAL(result)[indices[k]] = REAL(result)[indices[k]] + REAL(temp)[k]*phi.coeff(indices[k],i);
        }
        indices.clear();

    }

    UNPROTECT(1); //UNPROTECT result;
    return result;
}

#endif

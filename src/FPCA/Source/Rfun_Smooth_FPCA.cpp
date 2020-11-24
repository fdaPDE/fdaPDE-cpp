#include "../../FdaPDE.h"
#include "../../Skeletons/Include/FPCA_Skeleton.h"
#include "../Include/FPCA_Data.h"
#include "../../FE_Assemblers_Solvers/Include/Integration.h"

extern "C" {
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
        		return(FPCA_skeleton<1, 2, 2>(fPCAdata, Rmesh, validation));
        	else if(fPCAdata.getOrder() == 2 && mydim==2 && ndim==2)
        		return(FPCA_skeleton<2, 2, 2>(fPCAdata, Rmesh, validation));
        	else if(fPCAdata.getOrder() == 1 && mydim==2 && ndim==3)
        		return(FPCA_skeleton<1, 2, 3>(fPCAdata, Rmesh, validation));
        	else if(fPCAdata.getOrder() == 2 && mydim==2 && ndim==3)
        		return(FPCA_skeleton<2, 2, 3>(fPCAdata, Rmesh, validation));
        	else if(fPCAdata.getOrder() == 1 && mydim==3 && ndim==3)
        		return(FPCA_skeleton<1, 3, 3>(fPCAdata, Rmesh, validation));
                else if(fPCAdata.getOrder() == 2 && mydim==3 && ndim==3)
                        return(FPCA_skeleton<2, 3, 3>(fPCAdata, Rmesh, validation));

        	return(NILSXP);
        	 }
}

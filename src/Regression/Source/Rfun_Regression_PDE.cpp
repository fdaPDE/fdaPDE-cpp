#include "../../FdaPDE.h"
#include "../../Skeletons/Include/Auxiliary_Skeleton.h"
#include "../../Skeletons/Include/Regression_Skeleton.h"
#include "../../Skeletons/Include/Regression_Skeleton_Time.h"
#include "../../Skeletons/Include/GAM_Skeleton.h"
#include "../Include/Regression_Data.h"
#include "../../FE_Assemblers_Solvers/Include/Integration.h"
#include "../../Lambda_Optimization/Include/Optimization_Data.h"

extern "C"
{
        //! This function manages the various options for Spatial Regression
        /*!
                This function is then called from R code.
                \param Rlocations an R-matrix containing the spatial locations of the observations
                \param RbaryLocations A list with three vectors:
                                location points which are same as the given locations options (to checks whether both locations are the same),
                                a vector of element id of the points from the mesh where they are located,
                                a vector of barycenter of points from the located element.
                \param Robservations an R-vector containing the values of the observations.
                \param Rmesh an R-object containg the output mesh from Trilibrary
                \param Rorder an R-integer containing the order of the approximating basis.
                \param Rmydim an R-integer specifying if the mesh nodes lie in R^2 or R^3
                \param Rndim  an R-integer specifying if the "local dimension" is 2 or 3
                \param RK an R-matrix representing the diffusivity matrix of the model
                \param Rbeta an R-vector representing the advection term of the model
                \param Rc an R-double representing the reaction term of the model
                \param Rcovariates an R-matrix of covariates for the regression model
                \param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
                                the other are automatically considered in Neumann Condition.
                \param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
                \param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
                \param RarealDataAvg an R boolean indicating whether the areal data are averaged or not.
                \param Rsearch an R-integer to decide the search algorithm type (tree or naive search algorithm).
                \param Roptim optimzation type, DOF evaluation and loss function used coded as integer vector
                \param Rlambda a vector containing the penalization term of the empirical evidence respect to the prior one. or initial codition for optimized methods
                \param Rnrealizations integer, the number of random points used in the stochastic computation of the dofs
                \param Rseed integer, user defined seed for stochastic DOF computation methods
                \param RDOF_matrix user provided DOF matrix for GCV computation
                \param Rtune a R-double, Tuning parameter used for the estimation of GCV. called 'GCV.inflation.factor' in R code.
                \param Rsct user defined stopping criterion tolerance for optimized methods (newton or newton with finite differences)
                \return R-vectors containg the coefficients of the solution, prediction of the values, optimization data and much more
        */
        SEXP regression_PDE(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
                SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch,
                SEXP Roptim, SEXP Rlambda, SEXP Rnrealizations, SEXP Rseed, SEXP RDOF_matrix, SEXP Rtune, SEXP Rsct)
        {
        	RegressionDataElliptic regressionData(Rlocations, RbaryLocations, Robservations, Rorder, RK, Rbeta, Rc, Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch);
                OptimizationData optimizationData(Roptim, Rlambda, Rnrealizations, Rseed, RDOF_matrix, Rtune, Rsct);

        	UInt mydim = INTEGER(Rmydim)[0];
        	UInt ndim = INTEGER(Rndim)[0];

        	if(regressionData.getOrder()==1 && ndim==2)
        		return(regression_skeleton<RegressionDataElliptic, 1, 2, 2>(regressionData, optimizationData, Rmesh));
        	else if(regressionData.getOrder()==2 && ndim==2)
        		return(regression_skeleton<RegressionDataElliptic, 2, 2, 2>(regressionData, optimizationData, Rmesh));
        	else if(regressionData.getOrder()==1 && mydim==2 && ndim==3)
        		return(regression_skeleton<RegressionDataElliptic, 1, 2, 3>(regressionData, optimizationData, Rmesh));
        	else if(regressionData.getOrder()==2 && mydim==2 && ndim==3)
        		return(regression_skeleton<RegressionDataElliptic, 2, 2, 3>(regressionData, optimizationData, Rmesh));
            else if(regressionData.getOrder()==1 && mydim==3 && ndim==3)
                return(regression_skeleton<RegressionDataElliptic, 1, 3, 3>(regressionData, optimizationData, Rmesh));
            else if(regressionData.getOrder()==2 && mydim==3 && ndim==3)
                return(regression_skeleton<RegressionDataElliptic, 2, 3, 3>(regressionData, optimizationData, Rmesh));

        	return(NILSXP);
        }

        //! This function manages the various options for Spatio-Temporal Regression
	/*!
		This function is then called from R code.
		\param Rlocations an R-matrix containing the spatial locations of the observations
		\param RbaryLocations A list with three vectors:
				location points which are same as the given locations options (to checks whether both locations are the same),
		 		a vector of element id of the points from the mesh where they are located,
				a vector of barycenter of points from the located element.
		\param Rtime_locations an R-vector containing the temporal locations of the observations
		\param Robservations an R-vector containing the values of the observations.
		\param Rmesh an R-object containing the spatial mesh
		\param Rmesh_time an R-vector containing the temporal mesh
		\param Rorder an R-integer containing the order of the approximating basis in space.
		\param Rmydim an R-integer specifying if the mesh nodes lie in R^2 or R^3
		\param Rndim  an R-integer specifying if the "local dimension" is 2 or 3
                \param RK an R-matrix representing the diffusivity matrix of the model
                \param Rbeta an R-vector representing the advection term of the model
                \param Rc an R-double representing the reaction term of the model
		\param Rcovariates an R-matrix of covariates for the regression model
		\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
				the other are automatically considered in Neumann Condition.
		\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
		\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
		\param RarealDataAvg an R boolean indicating whether the areal data are averaged or not.
		\param Rflag_mass an R-integer that in case of separable problem specifies whether to use mass discretization or identity discretization
		\param Rflag_parabolic an R-integer specifying if the problem is parabolic or separable
		\param Ric an R-vector containing the initial condition needed in case of parabolic problem
		\param Rsearch an R-integer to decide the search algorithm type (tree or naive search algorithm).
		\param Roptim optimzation type, DOF evaluation and loss function used coded as integer vector
		\param Rlambda_S a vector containing the penalization term of the empirical evidence respect to the prior one. or initial codition for optimized methods
		\param Rlambda_T a vector containing the temporal penalization term of the empirical evidence respect to the prior one. or initial codition for optimized methods in eparable context
		\param Rnrealizations integer, the number of random points used in the stochastic computation of the dofs
		\param Rseed integer, user defined seed for stochastic DOF computation methods
		\param RDOF_matrix user provided DOF matrix for GCV computation
		\param Rtune a R-double, Tuning parameter used for the estimation of GCV. called 'GCV.inflation.factor' in R code.
		\param Rsct user defined stopping criterion tolerance for optimized methods (newton or newton with finite differences)
		\return R-vectors containg the coefficients of the solution, prediction of the values, optimization data and much more
	*/
        SEXP regression_PDE_time(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rmesh, SEXP Rmesh_time, SEXP Rorder, SEXP Rmydim, SEXP Rndim,
		SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates,	SEXP RBCIndices, SEXP RBCValues,  SEXP RincidenceMatrix, SEXP RarealDataAvg,
                SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric, SEXP Rsearch, SEXP Roptim, SEXP Rlambda_S, SEXP Rlambda_T, SEXP Rnrealizations, SEXP Rseed, SEXP RDOF_matrix, SEXP Rtune, SEXP Rsct)
        {
        	RegressionDataElliptic regressionData(Rlocations, RbaryLocations, Rtime_locations, Robservations, Rorder, RK, Rbeta, Rc,
                        Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rflag_mass, Rflag_parabolic, Ric, Rsearch);
                OptimizationData optimizationData(Roptim, Rlambda_S, Rlambda_T, Rflag_parabolic, Rnrealizations, Rseed, RDOF_matrix, Rtune, Rsct);

        	UInt mydim = INTEGER(Rmydim)[0];
        	UInt ndim = INTEGER(Rndim)[0];

        	if(regressionData.getOrder()==1 && ndim==2)
        		return(regression_skeleton_time<RegressionDataElliptic, 1, 2, 2>(regressionData, optimizationData, Rmesh, Rmesh_time));
        	else if(regressionData.getOrder()==2 && ndim==2)
        		return(regression_skeleton_time<RegressionDataElliptic, 2, 2, 2>(regressionData, optimizationData, Rmesh, Rmesh_time));
        	else if(regressionData.getOrder()==1 && mydim==2 && ndim==3)
        		return(regression_skeleton_time<RegressionDataElliptic, 1, 2, 3>(regressionData, optimizationData, Rmesh, Rmesh_time));
        	else if(regressionData.getOrder()==2 && mydim==2 && ndim==3)
        		return(regression_skeleton_time<RegressionDataElliptic, 2, 2, 3>(regressionData, optimizationData, Rmesh, Rmesh_time));
            else if(regressionData.getOrder()==1 && mydim==3 && ndim==3)
                return(regression_skeleton_time<RegressionDataElliptic, 1, 3, 3>(regressionData, optimizationData, Rmesh, Rmesh_time));
            else if(regressionData.getOrder()==2 && mydim==3 && ndim==3)
                return(regression_skeleton_time<RegressionDataElliptic, 2, 3, 3>(regressionData, optimizationData, Rmesh, Rmesh_time));

        	return(NILSXP);
        }

        //! This function manages the various options for GAM Spatial Regression
	/*!
		This function is then called from R code.
		\param Rlocations an R-matrix containing the spatial locations of the observations
		\param RbaryLocations A list with three vectors:
				location points which are same as the given locations options (to checks whether both locations are the same),
				a vector of element id of the points from the mesh where they are located,
				a vector of barycenter of points from the located element.
		\param Robservations an R-vector containing the values of the observations.
		\param Rmesh an R-object containg the output mesh from Trilibrary
		\param Rorder an R-integer containing the order of the approximating basis.
		\param Rmydim an R-integer specifying if the mesh nodes lie in R^2 or R^3
		\param Rndim  an R-integer specifying if the "local dimension" is 2 or 3
                \param RK an R-matrix representing the diffusivity matrix of the model
                \param Rbeta an R-vector representing the advection term of the model
                \param Rc an R-double representing the reaction term of the model
		\param Rcovariates an R-matrix of covariates for the regression model
		\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
				the other are automatically considered in Neumann Condition.
		\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
		\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
		\param RarealDataAvg an R boolean indicating whether the areal data are averaged or not.
		\param Rfamily Denotes the distribution of the data, within the exponential family.
		\param Rmax_num_iteration Maximum number of steps run in the PIRLS algorithm, set to 15 by default.
		\param Rtreshold an R-double used for arresting FPIRLS algorithm. Algorithm stops when two successive iterations lead to improvement in penalized log-likelihood smaller than threshold.
		\param Rtune It is usually set to 1, but may be higher. It gives more weight to the equivalent degrees of freedom in the computation of the value of the GCV.
		\param Rmu0 Initial value of the mean (natural parameter). There exists a default value for each familiy
		\param RscaleParam If necessary and not supplied, the scale parameter \phi is estimated. See method.phi for details.
		\param Rsearch an R-integer to decide the search algorithm type (tree or naive search algorithm).
		\param Roptim optimzation type, DOF evaluation and loss function used coded as integer vector
		\param Rlambda a vector containing the penalization term of the empirical evidence respect to the prior one. or initial codition for optimized methods
		\param Rnrealizations integer, the number of random points used in the stochastic computation of the dofs
		\param Rseed integer, user defined seed for stochastic DOF computation methods
		\param RDOF_matrix user provided DOF matrix for GCV computation
		\param Rtune a R-double, Tuning parameter used for the estimation of GCV. called 'GCV.inflation.factor' in R code.
		\param Rsct user defined stopping criterion tolerance for optimized methods (newton or newton with finite differences)
		\return R-vectors containg the coefficients of the solution, prediction of the values, optimization data and much more
	*/
        SEXP gam_PDE(SEXP Rlocations, SEXP RbaryLocations ,SEXP Robservations, SEXP Rmesh, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
                SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg,
          	SEXP Rfamily, SEXP Rmax_num_iteration, SEXP Rtreshold, SEXP Rmu0, SEXP RscaleParam, SEXP Rsearch,
                SEXP Roptim, SEXP Rlambda, SEXP Rnrealizations, SEXP Rseed, SEXP RDOF_matrix, SEXP Rtune, SEXP Rsct)
        {
        	GAMDataElliptic regressionData(Rlocations, RbaryLocations, Robservations, Rorder, RK, Rbeta, Rc, Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch, Rmax_num_iteration, Rtreshold);
                OptimizationData optimizationData(Roptim, Rlambda, Rnrealizations, Rseed, RDOF_matrix, Rtune, Rsct);

        	UInt mydim = INTEGER(Rmydim)[0];
        	UInt ndim = INTEGER(Rndim)[0];

          	std::string family = CHAR(STRING_ELT(Rfamily,0));

                if(regressionData.getOrder()==1 && mydim==2 && ndim==2)
                	return(GAM_skeleton<GAMDataElliptic, 1, 2, 2>(regressionData, optimizationData, Rmesh, Rmu0, family, RscaleParam));
                else if(regressionData.getOrder()==2 && mydim==2 && ndim==2)
                	return(GAM_skeleton<GAMDataElliptic, 2, 2, 2>(regressionData, optimizationData, Rmesh, Rmu0, family, RscaleParam));
                else if(regressionData.getOrder()==1 && mydim==2 && ndim==3)
                    return(GAM_skeleton<GAMDataElliptic, 1, 2, 3>(regressionData, optimizationData, Rmesh, Rmu0, family, RscaleParam));
                else if(regressionData.getOrder()==2 && mydim==2 && ndim==3)
                    return(GAM_skeleton<GAMDataElliptic, 2, 2, 3>(regressionData, optimizationData, Rmesh, Rmu0, family, RscaleParam));
                else if(regressionData.getOrder()==1 && mydim==3 && ndim==3)
                    return(GAM_skeleton<GAMDataElliptic, 1, 3, 3>(regressionData, optimizationData, Rmesh, Rmu0, family, RscaleParam));
                else if(regressionData.getOrder()==2 && mydim==3 && ndim==3)
                    return(GAM_skeleton<GAMDataElliptic, 2, 3, 3>(regressionData, optimizationData, Rmesh, Rmu0, family, RscaleParam));

                return(R_NilValue);
        }

        //! A utility, not used for system solution, may be used for debugging
        /*!
                This function is then called from R code.
                \param Rlocations an R-matrix containing the spatial locations of the observations
                \param RbaryLocations A list with three vectors:
                                location points which are same as the given locations options (to checks whether both locations are the same),
                                a vector of element id of the points from the mesh where they are located,
                                a vector of barycenter of points from the located element.
                \param Robservations an R-vector containing the values of the observations.
                \param Rmesh an R-object containg the output mesh from Trilibrary
                \param Rorder an R-integer containing the order of the approximating basis.
                \param Rmydim an R-integer specifying if the mesh nodes lie in R^2 or R^3
                \param Rndim  an R-integer specifying if the "local dimension" is 2 or 3
                \param RK an R-matrix representing the diffusivity matrix of the model
                \param Rbeta an R-vector representing the advection term of the model
                \param Rc an R-double representing the reaction term of the model
                \param Rcovariates an R-matrix of covariates for the regression model
                \param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
                                the other are automatically considered in Neumann Condition.
                \param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
                \param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
                \param RarealDataAvg an R boolean indicating whether the areal data are averaged or not.
                \param Rsearch an R-integer to decide the search algorithm type (tree or naive search algorithm).
                \return R matrix produced by FEM skeleton
        */
        SEXP get_FEM_PDE_matrix(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder,SEXP Rmydim, SEXP Rndim, SEXP RK, SEXP Rbeta, SEXP Rc,
                SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch)
        {
                RegressionDataElliptic regressionData(Rlocations, RbaryLocations, Robservations, Rorder, RK, Rbeta, Rc, Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch);

                // Get mydim and ndim
                UInt mydim = INTEGER(Rmydim)[0];
                UInt ndim = INTEGER(Rndim)[0];

                typedef EOExpr<Mass>  ETMass;  Mass  EMass;  ETMass  mass(EMass);
                typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
                typedef EOExpr<Grad>  ETGrad;  Grad  EGrad;  ETGrad  grad(EGrad);

                const Real& c = regressionData.getC();
                const Diffusion<PDEParameterOptions::Constant>& K = regressionData.getK();
                const Advection<PDEParameterOptions::Constant>& beta = regressionData.getBeta();

                if(regressionData.getOrder()==1 && ndim==2)
                        return(get_FEM_Matrix_skeleton<1,2,2>(Rmesh, c*mass+stiff[K]+beta.dot(grad)));
                if(regressionData.getOrder()==2 && ndim==2)
                        return(get_FEM_Matrix_skeleton<2,2,2>(Rmesh, c*mass+stiff[K]+beta.dot(grad)));
                return(NILSXP);
        }
}

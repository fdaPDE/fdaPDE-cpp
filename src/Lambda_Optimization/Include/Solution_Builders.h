#ifndef __SOLUTION_BUILDERS_H__
#define __SOLUTION_BUILDERS_H__

// HEADERS
#include <cmath>
#include <limits>
#include <utility>
#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../Regression/Include/Regression_Data.h"
#include "../../Global_Utilities/Include/Lambda.h"

//! Output struct to be used to return values in R
template<UInt num_params>
struct output_Data
{
        std::string             content{"Empty"};          //!< Suggests what the output is containing and how it should be used
        MatrixXr                z_hat;                     //!< Model predicted values in the locations
        std::vector<Real>       rmse;                      //!< Model root mean squared error
        Real                    sigma_hat_sq    = -1.0;    //!< Model estimated variance of errors
        std::vector<Real>       dof             = {};      //!< tr(S) + q, degrees of freedom of the model
	lambda::type<num_params>	lambda_sol = lambda::init<num_params>(0.0); //!< Lambda obtained in the solution
        UInt                    lambda_pos      = 0;       //!< Position of optimal lambda, only for grid evaluation, in R numebring starting from 1 (0 means no grid used)
        UInt                    n_it            = 0;       //!< Number of iterations for the method
        Real                    time_partial    = 0.0;     //!< Time, from beginning to end of the optimization method
        std::vector<Real>       GCV_evals       = {-1};    //!< GCV evaluations vector of explored lambda, with the optimization iterative method or grid
        std::vector<lambda::type<num_params>> lambda_vec = {lambda::init<num_params>(-1)}; //!< Vector of explored lambda with with the optimization iterative method or grid
        Real                    GCV_opt         = -1;      //!< GCV optimal comptued in the vector of lambdas
        int                     termination     = -2;      //!< Reason of termination of the iterative optimization method (reached tolerance or max number of iterations)
        MatrixXv                betas;                     //!< Regression coefficients of the optimal solution
        //the following two are to be able to deduce the single lambdaS and lambdaT indices from the pair index (e.g. lambda_pos) in the temporal case
	UInt                    size_S = 0;		   //!< Size of the lambdaS vector
        UInt                    size_T = 0;		   //!< Size of the lambdaT vector      
};

//! Unique namespace to manage the output
namespace Solution_Builders
{
        //! Function to build the output of regression problems
        /*!
         \tparam InputHandler type of data used
         \tparam ORDER order of the finite elements in mesh
         \tparam mydim dimension of elements in mesh
         \tparam ndim dimension of space in mesh
         \param solution matrix collecting the output of apply function
         \param output output_Data struct coming from a Lambda_Optimizer type class
         \param mesh to be returned to the user
         \param regressionData the original data passed by the user
         \return SEXP containg all the data that will be managed by R code
        */
        template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
        static SEXP build_solution_plain_regression(const MatrixXr & solution, const output_Data<1> & output, const MeshHandler<ORDER, mydim, ndim> & mesh, const InputHandler & regressionData, const MixedFERegression<InputHandler>& regression);
        
        template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
        static SEXP build_solution_temporal_regression(const MatrixXr & solution, const output_Data<2> & output, const MeshHandler<ORDER, mydim, ndim> & mesh, const InputHandler & regressionData, const MixedFERegression<InputHandler>& regression);
};

#include "Solution_Builders_imp.h"

#endif

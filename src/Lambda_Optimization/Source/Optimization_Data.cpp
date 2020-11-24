#include "../Include/Optimization_Data.h"

//! Utility used by the constructor to set the parameters common to all methods or that do not not_require specific treatment
/*!
 \param Roptim optimization method (used to fill criterion, DOF_evaluation and loss_function)
 \param Rnrealizations number of realizations for the stochastic gcv computation
 \param Rseed seed to be stored for reproducibility of stochastic gcv computation
 \param RDOF_MATRIX matrix of dof possibly passed by the user
 \param Rtune tuning parameter for gcv, used in GAM methods
 \param Rstc stopping criterion for optimized methods
*/
void OptimizationData::builder_utility(SEXP Roptim, SEXP Rnrealizations, SEXP Rseed, SEXP RDOF_matrix, SEXP Rtune, SEXP Rsct)
{
        UInt criterion = INTEGER(Roptim)[0]; // Decipher the Roptim sequence of numbers, first criterion
        if(criterion == 2)
        {
                this->set_criterion("newton_fd");
                this->set_stopping_criterion_tol(REAL(Rsct)[0]);
        }
        else if(criterion == 1)
        {
                this->set_criterion("newton");
                this->set_stopping_criterion_tol(REAL(Rsct)[0]);
        }
        else if(criterion == 0)
                this->set_criterion("grid");

        UInt DOF_evaluation = INTEGER(Roptim)[1]; // Decipher the Roptim sequence of numbers, second DOF_evaluaton
        if(DOF_evaluation == 0)
        {
                this->set_DOF_evaluation("not_required");
        }
        else if(DOF_evaluation == 1)
        {
                this->set_DOF_evaluation("stochastic");              // since the method is stochastic
                this->set_nrealizations(INTEGER(Rnrealizations)[0]); // set nrealizations
                this->set_seed(INTEGER(Rseed)[0]);                   // set seed
        }
        else
        {
                this->set_DOF_evaluation("exact");
        }

        UInt loss_function = INTEGER(Roptim)[2]; // Decipher the Roptim sequence of numbers, third loss function
        if(loss_function == 0)
        {
                this->set_loss_function("unused");
        }
        else if(loss_function == 1)
        {
                this->set_loss_function("GCV");
        }

        // Tuning parameter, set from R
        this->set_tuning(REAL(Rtune)[0]);

        // DOF_matrix, set from R
        UInt n_ = INTEGER(Rf_getAttrib(RDOF_matrix, R_DimSymbol))[0]; // #Rows
        UInt p_ = INTEGER(Rf_getAttrib(RDOF_matrix, R_DimSymbol))[1]; // #Columns
        DOF_matrix.resize(n_, p_);
        for(auto i=0; i<n_; ++i) // Fill the matrix if not empty
        {
                for(auto j=0; j<p_ ; ++j)
                {
                        DOF_matrix(i,j) = REAL(RDOF_matrix)[i+ n_*j];
                }
        }

        if(n_==0 || p_==0) // If one of the dimensions is null, put the other also to 0 for consistency
        {
                n_ = 0;
                p_ = 0;
        }
        DOF_matrix.resize(n_,p_);
}

//! Utility used by the constructor to set the vector of lambdas (space or time)
/*!
 \param Rlambda the vector of lambdas to be set passed by the user
 \param vect the vector of lambdas to be filled with Rlambda
 \param size the size of vect, to be stored
*/
void OptimizationData::fill_lambda(SEXP Rlambda, std::vector<Real> & vect, UInt & size)
{
        size = Rf_length(Rlambda); // stoe teh size
        vect.resize(size);

        for(UInt i=0; i<size; ++i) // fill the vector of lambdas
        {
                vect[i] = REAL(Rlambda)[i];
        }
}

//! Utility used by the constructor to set initial value of an opimized method
/*!
 \param Rlambda the initial lambda passed by the user
 \param init the initial value to be set
*/
void OptimizationData::initialize_lambda(SEXP Rlambda, Real & init)
{
        if(Rf_length(Rlambda)>0) // check if an initial value is passed
                init = REAL(Rlambda)[0]; // store the value
}

//! Main constructor of the class for pure spatial problems
/*!
 \param Roptim optimization method (used to fill criterion, DOF_evaluation and loss_function)
 \param Rlambda is the initial value if the Roptim[0] is 1 or 2 else is a vector of spatial lambdas
 \param Rnrealizations number of realizations for the stochastic gcv computation
 \param Rseed seed to be stored for reproducibility of stochastic gcv computation
 \param RDOF_MATRIX matrix of dof possibly passed by the user
 \param Rtune tuning parameter for gcv, used in GAM methods
 \param Rstc stopping criterion for optimized methods
*/
OptimizationData::OptimizationData(SEXP Roptim, SEXP Rlambda, SEXP Rnrealizations, SEXP Rseed, SEXP RDOF_matrix, SEXP Rtune, SEXP Rsct)
{
        builder_utility(Roptim, Rnrealizations, Rseed, RDOF_matrix, Rtune, Rsct); // build common terms

        // Lambda
        if(this->criterion == "grid") // Rlambda is a vector
        {
                fill_lambda(Rlambda, this->lambda_S, this->size_S);
                set_lambdaS_backup(); // save a backup
        }
        else // Rlambda is an initial value (or null)
        {
                initialize_lambda(Rlambda, this->initial_lambda_S);
        }
}

//! Main constructor of the class for pure spatio-temporal problems
/*!
 \param Roptim optimization method (used to fill criterion, DOF_evaluation and loss_function)
 \param Rlambda_S is the initial value if the Roptim[0] is 1 or 2 else is a vector of spatial lambdas
 \param Rlambda_T is the initial value if the Roptim[0] is 1 or 2 else is a vector of temporal lambdas
 \param Rflag_parabolic check if the spatio-temporal problem is parabolic or separable
 \param Rnrealizations number of realizations for the stochastic gcv computation
 \param Rseed seed to be stored for reproducibility of stochastic gcv computation
 \param RDOF_MATRIX matrix of dof possibly passed by the user
 \param Rtune tuning parameter for gcv, used in GAM methods
 \param Rstc stopping criterion for optimized methods
*/
OptimizationData::OptimizationData(SEXP Roptim, SEXP Rlambda_S, SEXP Rlambda_T, SEXP Rflag_parabolic, SEXP Rnrealizations, SEXP Rseed, SEXP RDOF_matrix, SEXP Rtune, SEXP Rsct)
{
        builder_utility(Roptim, Rnrealizations, Rseed, RDOF_matrix, Rtune, Rsct); // build common terms

        // Lambda
        if(this->criterion == "grid") // Rlambda is a vector
        {
                fill_lambda(Rlambda_S, this->lambda_S, this->size_S);
                fill_lambda(Rlambda_T, this->lambda_T, this->size_T);
                set_lambdaS_backup(); // save a backup
        }
        else // Rlambda is an initial value (or null)
        {
                initialize_lambda(Rlambda_S, this->initial_lambda_S);
                initialize_lambda(Rlambda_T, this->initial_lambda_T);
        }
}

//! Simple method to print the characteristics of an optimization method
void OptimizationData::print_opt_data(void) const
{
        Rprintf("\nOptimization data:\n");
        Rprintf("Criterion: %s\n", criterion.c_str());
        Rprintf("DOF valuation: %s\n", DOF_evaluation.c_str());
        Rprintf("Loss Function: %s\n",loss_function.c_str());
}

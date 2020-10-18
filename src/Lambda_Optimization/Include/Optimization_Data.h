#ifndef __OPTIMIZATION_DATA_H__
#define __OPTIMIZATION_DATA_H__

#include "../../FdaPDE.h"
#include <string>

//!  Class to collect data for optimization
/*!
 * This class collects all the data used in the optimization framework, is constructed
 * using the input from R and stores a structure flexible enough for any problem currently
 * implemented in the library
*/
class  OptimizationData
{
        private:
                std::string criterion      = "grid";            //!< grid [default], newton or newton_fd
                std::string DOF_evaluation = "not_required";    //!< not_required [default], stochastic or exact
                std::string loss_function  = "unused";          //!< unused [default] or GCV

                // For grid
                std::vector<Real> lambda_S = {-1.};             //!< Stores the vector of spatial lambdas to be evaluated in case of criterion = grid
                std::vector<Real> lambda_T = {-1.};             //!< Stores the vector of temporal lambdas to be evaluated in case of criterion = grid
                UInt size_S = 1;                                //!< Size of vector lambda_S
                UInt size_T = 1;                                //!< Size of vector lambda_T

                // For grid fixed method
                UInt best_lambda_S = 0;	                        //!< Stores the index of the best lambdaS according to method
                UInt best_lambda_T = 0;	                        //!< Stores the index of the best lambdaT according to method
                Real best_value    = std::numeric_limits<Real>::max();	//!< Stores the value of the best loss function

                // For optimized methods
                Real initial_lambda_S = 0.;                     //!< Initial lambda_S for optimized methods (newton or newton_fd)
                Real initial_lambda_T = 0.;                     //!< Initial lambda_T for optimized methods (newton or newton_fd)
                UInt seed             = 0;                      //!< The seed of random points used in the stochastic computation of the dofs [default 0]
                UInt nrealizations    = 100;                    //!< The number of random points used in the stochastic computation of the dofs [default 100]

                // To keep track of optimization
                Real last_lS_used = std::numeric_limits<Real>::infinity();      //!< last lambda_S used in optimization
                Real last_lT_used = std::numeric_limits<Real>::infinity();      //!< last lambda_T used in optimization
                Real current_lambdaS = -1.;                                     //!< Value of the lambda_S for which we are currently performing the computation

                // If already present
                MatrixXr DOF_matrix;                            //!< Matrix of dof (if passed by the user no need to compute dofs, we can use this)

                // Tuning parameter
                Real tuning = 1.;                               //!< To tune gcv value computation (useful in GAM methods)

                // For GAM
                std::vector<Real> lambdaS_backup;               //!< Backup vector of lambda_S (as passed by the user), used in GAM methods

                // For optimized methods
                Real stopping_criterion_tol = 0.05;             //!< Contains the user defined tolerance for optimized methods


                void builder_utility(SEXP Roptim, SEXP Rnrealizations, SEXP Rseed, SEXP RDOF_matrix, SEXP Rtune, SEXP Rsct);
                void fill_lambda(SEXP Rlambda, std::vector<Real> & vect, UInt & size);
                void initialize_lambda(SEXP Rlambda, Real & init);

        public:
                //! Default constructor of the class
                OptimizationData() = default;

                OptimizationData(SEXP Roptim, SEXP Rlambda, SEXP Rnrealizations, SEXP Rseed, SEXP RDOF_matrix, SEXP Rtune, SEXP Rsct);
                OptimizationData(SEXP Roptim, SEXP Rlambda_S, SEXP Rlambda_T, SEXP Rflag_parabolic, SEXP Rnrealizations, SEXP Rseed, SEXP RDOF_matrix, SEXP Rtune, SEXP Rsct);

                // Setters
                inline void set_criterion(const std::string && criterion_) {criterion = criterion_;}                            //!< Setter of criterion \param criterion_ new criterion
                inline void set_DOF_evaluation(const std::string && DOF_evaluation_) {DOF_evaluation = DOF_evaluation_;}        //!< Setter of DOF_evaluation \param DOF_evaluation_ new DOF_evaluation
                inline void set_loss_function(const std::string && loss_function_) {loss_function = loss_function_;}            //!< Setter of loss_function \param loss_function_ new \loss_function
                inline void set_lambda_S(const std::vector<Real> & lambda_S_) {lambda_S = lambda_S_;}                           //!< Setter of lambda_S vector \param lambda_S_ new lambda_S
                inline void set_lambda_T(const std::vector<Real> & lambda_T_) {lambda_T = lambda_T_;}                           //!< Setter of lambda_T vector \param lambda_T_ new lambda_T
                inline void set_best_lambda_S(const UInt best_lambda_S_) {best_lambda_S = best_lambda_S_;}                      //!< Setter of best_lambda_S \param best_lambda_S_ new best_lambda_S
                inline void set_best_lambda_T(const UInt best_lambda_T_) {best_lambda_T = best_lambda_T_;}                      //!< Setter of best_lambda_T \param best_lambda_T_ new best_lambda_T
                inline void set_best_value(const Real best_value_) {best_value = best_value_;}                                  //!< Setter of best_value \param best_value_ new best_value
                inline void set_initial_lambda_S(const Real initial_lambda_S_) {initial_lambda_S = initial_lambda_S_;}          //!< Setter of initial_lambda_S \param initial_lambda_S_ new initial_lambda_S
                inline void set_initial_lambda_T(const Real initial_lambda_T_) {initial_lambda_T = initial_lambda_T_;}          //!< Setter of initial_lambda_T \param initial_lambda_T_ new initial_lambda_T
                inline void set_seed(const UInt seed_){seed = seed_;}                                                           //!< Setter of seed \param seed_ new seed
                inline void set_nrealizations(const UInt nrealizations_) {nrealizations = nrealizations_;}                      //!< Setter of nrealizations \param nrealizations_ new nrealizations
                inline void set_last_lS_used(const Real last_lS_used_) {last_lS_used = last_lS_used_;}                          //!< Setter of last_lS_used \param last_lS_used_ new last_lS_used
                inline void set_last_lT_used(const Real last_lT_used_) {last_lT_used = last_lT_used_;}                          //!< Setter of last_lT_used \param last_lT_used_ new last_lT_used
                inline void set_DOF_matrix(const MatrixXr & DOF_matrix_) {DOF_matrix = DOF_matrix_;}                            //!< Setter of DOF_matrix \param DOF_matrix_ new DOF_matrix
                inline void set_stopping_criterion_tol(Real stc_) {stopping_criterion_tol = stc_;}                              //!< Setter of stopping_criterion_tol \param stc_ new stopping_criterion_tol
                inline void set_tuning(const Real tuning_) {tuning = tuning_;}                                                  //!< Setter of tuning \param tuning_ new tuning
                inline void set_current_lambdaS(const Real new_lambdaS) {current_lambdaS = new_lambdaS;}                        //!< Utility for GAM problems, that always need a vector \param new_lambdaS, new current_lambdaS
                inline void setCurrentLambda(UInt lambda_index) {lambda_S = std::vector<Real>(1,lambdaS_backup[lambda_index]);} //!< Setter of a backup of lambda_S manpualted in setCurrentLambda
                inline void set_lambdaS_backup(void) {lambdaS_backup = lambda_S;}

                // Getters
                inline std::string get_criterion(void) const {return criterion;}                        //!< Getter of criterion \return criterion
                inline std::string get_DOF_evaluation(void) const {return DOF_evaluation;}              //!< Getter of DOF_evaluation \return DOF_evaluation
                inline std::string get_loss_function(void) const {return loss_function;}                //!< Getter of loss_function \return loss_function
                inline std::vector<Real> get_lambda_S(void) const {return lambda_S;}                    //!< Getter of lambda_S \return lambda_S
                inline std::vector<Real> get_lambda_T(void) const {return lambda_T;}                    //!< Getter of lambda_T \return lambda_T
                inline UInt get_size_S(void) const {return lambda_S.size();}                            //!< Getter of size_S \return size_S
                inline UInt get_size_T(void) const {return lambda_T.size();}                            //!< Getter of size_T \return size_T
                inline UInt get_best_lambda_S(void) const {return best_lambda_S;}                       //!< Getter of best_lambda_S \return best_lambda_S
                inline UInt get_best_lambda_T(void) const {return best_lambda_T;}                       //!< Getter of best_lambda_T \return best_lambda_T
                inline Real get_best_value(void) const {return best_value;}                             //!< Getter of best_value \return best_value
                inline Real get_initial_lambda_S(void) const {return initial_lambda_S;}                 //!< Getter of initial_lambda_S \return initial_lambda_S
                inline Real get_initial_lambda_T(void) const {return initial_lambda_T;}                 //!< Getter of initial_lambda_T \return initial_lambda_T
                inline UInt get_seed(void) const {return seed;}                                         //!< Getter of seed \return seed
                inline UInt get_nrealizations(void) const {return nrealizations;}                       //!< Getter of nrealizations  \return nrealizations
                inline Real get_last_lS_used(void) const {return last_lS_used;}                         //!< Getter of last_lS_used \return last_lS_used
                inline Real get_last_lT_used(void) const {return last_lT_used;}                         //!< Getter of last_lT_used \return last_lT_used
                inline MatrixXr const & get_DOF_matrix(void) const {return DOF_matrix;}                 //!< Getter of DOF_matrix \return DOF_matrix
                inline Real get_tuning(void) const {return tuning;}                                     //!< Getter of tuning \return tuning
                inline Real get_stopping_criterion_tol(void) const {return stopping_criterion_tol;}     //!< Getter of stopping_criterion_tol \return stopping_criterion_tol
                inline Real get_current_lambdaS(void) const {return current_lambdaS;}                   //!< Getter of current_lambdaS \return current_lambdaS
                inline const std::vector<Real> * get_LambdaS_vector() const {return &lambdaS_backup;}   //!< Getter of backup lamnda_S vector for GAM problems \return &lambdaS_backup

                // Debugging
                void print_opt_data(void) const;
};

#endif

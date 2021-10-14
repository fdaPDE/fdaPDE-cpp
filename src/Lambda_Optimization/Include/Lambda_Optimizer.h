#ifndef __LAMBDA_OPTIMIZER_H__
#define __LAMBDA_OPTIMIZER_H__

// HEADERS
#include "../../FdaPDE.h"
#include "Auxiliary_Optimizer.h"
#include "Carrier.h"
#include "Gof_Updater.h"
#include "../../FE_Assemblers_Solvers/Include/Solver.h"
#include <algorithm>
#include "../../Global_Utilities/Include/Lambda.h"

// CLASSES
// **** GENERAL METHODS ***

//! Father class used for multidiemnsional lambda optimization purposes
/*!
 This virtual class stores the model to be used by all its children,
 i. e. the classes actually instantiating a real optimization method and
 performing evaluations.
 \tparam InputCarrier Carrier-type parameter that contains insight about the problem to be solved
 \tparam size specialization parameter used to characterize the size of the lambda to be used
*/
template <typename InputCarrier, UInt size>
class Lambda_optimizer
{
        protected:
                //! Model containing all the information necessary for the computation of the optimal value
                InputCarrier & the_carrier;

                // CONSTRUCTORS
                //! Constructor of the class given the InputCarrier
                /*!
                 \param the_carrier the structure from which to take all the data for the derived classes
                */
                Lambda_optimizer(InputCarrier & the_carrier_):
                        the_carrier(the_carrier_) {}

                // UPDATERS
                //! A pure virtual member used in children classes for updates of internal data
                /*!
                 \param lambda the value of lambda with which to perform the update
                */
        virtual void update_parameters(lambda::type<size> lambda) = 0;

        public:
                //! Virtual Destuctor
                virtual ~Lambda_optimizer(){};
};

//----------------------------------------------------------------------------//
// *** GCV-BASED ***

//! Father class used for multidimensional lambda gcv-based methods
/*!
 This virtual class inherits from the generic multidimensional optimizer Lambda_optimizer
 and contains the methods used to compute the main statistics
 related to the optimization problem, like predicted output, errors and variances
 also stores a container to retrieve the data and pass them to the user.
 \tparam InputCarrier Carrier-type parameter that contains insight about the problem to be solved
 \tparam size specialization parameter used to characterize the size of the lambda to be used
*/
template <typename InputCarrier, UInt size>
class GCV_Family: Lambda_optimizer<InputCarrier, size>
{
        protected:
                //! Model containing all the information necessary for the computation of the optimal value
                using Lambda_optimizer<InputCarrier, size>::the_carrier;

                // Output data
                VectorXr        z_hat;                  //!< Model predicted values in the locations [size s]
                VectorXr        eps_hat;                //!< Model predicted error in the locations (residuals) [size s]
                Real            SS_res = 0.0;           //!< Model predicted sum of squares of the residuals
                Real            rmse = 0.0;             //!< Model root mean squared error
                Real            sigma_hat_sq = 0.0;     //!< Model estimated variance of error
                UInt            s;                      //!< Model number of observations (i.e. #locations)
                output_Data<size>     output;		//!< Output, needed to be user-available, necessarily public

                // Degrees of freedom
                Real            dof = 0.0;              //!< tr(S) + q, degrees of freedom of the model
                Real            dor = 0.0;              //!< s - dof, degrees of freedom of the residuals

                UInt            use_index = -1;         //!< Index of the DOF_matrix to be used, if non empty

                // SETTERS of the output data
        virtual void compute_z_hat(lambda::type<size> lambda) = 0;    //!< Utility to compute the size of predicted value in the locations
                void compute_z_hat_from_f_hat(const VectorXr & f_hat);
                void compute_eps_hat(void);
                void compute_SS_res(void);
                void compute_rmse(void);
                void compute_sigma_hat_sq(void);
                void compute_s(void);

                // UPDATERS
                void update_errors(lambda::type<size> lambda);

                // DOF methods
        virtual void update_dof(lambda::type<size> lambda) = 0;       //!< Utility to compute the degrees of freedom of the model
        virtual void update_dor(lambda::type<size> lambda) = 0;       //!< Utility to compute the degrees of freedom of the residuals

                // CONSTRUCTORS
                //! Constructor of the class given the InputCarrier
                /*!
                 \param the_carrier the structure from which to take all the data for the derived classes
                 \sa compute_s()
                */
                GCV_Family(InputCarrier & the_carrier_):
                        Lambda_optimizer<InputCarrier, size>(the_carrier_)
                        {
                                this->compute_s();      // stores immediately the number of locations
                                output.size_S = the_carrier_.get_opt_data()->get_size_S();
				output.size_T = (size == 1) ? 0 : the_carrier_.get_opt_data()->get_size_T();
                        }

        public:
                // UTILITY FOR DOF MATRIX
        inline  void set_index(UInt index){this->use_index = index;}
                // PUBLIC UPDATERS
        virtual void update_parameters(lambda::type<size> lambda) = 0; //!< Utility to update all the prameters of the model

                void zero_updater(lambda::type<size> lambda);

                // GCV-COMPUTATION
        virtual Real compute_f(lambda::type<size> lambda) = 0;       //!< Main function, represents the gcv computation

                // OUTPUT MANAGERS
                output_Data<size>  get_output(std::pair<lambda::type<size>, UInt> optimal_pair, const timespec & time_count, const std::vector<Real> & GCV_v, const std::vector<lambda::type<size>> & lambda_v, int termination_);
                void set_output_partial_best(void);
                output_Data<size> get_output_full(void);
                void set_output_partial(void);
                void combine_output_prediction(const VectorXr & f_hat, output_Data<size> & outp, UInt cols);
                //! Virtual Destuctor
        virtual ~GCV_Family(){};
};

//----------------------------------------------------------------------------//
// ** GCV_EXACT **

//! Derived class used for multidimensional lambda exact gcv-based methods
/*!
 This class implements exact methods to identify the correct value of the gcv
 and possibly of its derivatives, by partially solving manually the apply system
 \tparam InputCarrier Carrier-type parameter that contains insight about the problem to be solved
 \tparam size specialization parameter used to characterize the size of the lambda to be used
 \todo 2-D lambda optimization still to be implemented
*/
template<typename InputCarrier, UInt size>
class GCV_Exact: public GCV_Family<InputCarrier, size>
{
/*
        [[ VERSION WITH TIMES STILL TO BE IMPLEMENTED ]]
*/
};

//! Derived class used for unidimensional lambda exact gcv-based methods
/*!
 This class implements exact methods to identify the correct value of the gcv
 and possibly of its derivatives, by partially solving manually the apply system
 \tparam InputCarrier Carrier-type parameter that contains insight about the problem to be solved
*/
template<typename InputCarrier>
class GCV_Exact<InputCarrier, 1>: public GCV_Family<InputCarrier, 1>
{
        private:
                //! An external updater whose purpose is keeping the internal values coherent with the computations to be made from time to time
                GOF_updater<GCV_Exact<InputCarrier, 1>, Real> gu;

                // INTERNAL DATA STRUCTURES
                MatrixXr  R_; 		//!< stores the value of R1^t*R0^{-1}*R1 [size nnodes x nnodes]
                MatrixXr  T_; 		//!< stores the value of Psi^t*Q*Psi+lambda*R [size nnodes x nnodes]
                MatrixXr  V_; 		//!< stores the value of T^{-1}*Psi^t*Q [size nnodes x s]
                MatrixXr  S_;           //!< stores the value of Psi*V [as in Stu-Hunter Sangalli] [size s x s]
                Real      trS_ = 0.0;   //!< stores the value of the trace of S
                MatrixXr  dS_;          //!< stores the derivative of S w.r.t. lambda [size s x s]
                Real      trdS_ = 0.0;  //!< stores the value of the trace of dS
                MatrixXr  ddS_;         //!< stores the second derivative of S w.r.t. lambda [size s x s]
                Real      trddS_ = 0.0; //!< stores the value of the trace of ddS
                Real      lambdaT = -1; //!< stores the lambdaT for parabolic case

                //! Additional utility matrices [just the ones for the specific carrier that is proper of the problem]
                AuxiliaryData<InputCarrier> adt;

                // COMPUTERS and DOF methods
                void compute_z_hat (lambda::type<1> lambda) override;
                void update_dof(lambda::type<1> lambda)     override;
                void update_dor(lambda::type<1> lambda)     override;

                // SETTERS
                void set_R_(void);
                void set_R_(Real lambdaT);
                void set_T_(lambda::type<1> lambda);
                void set_V_(void);
                void set_S_and_trS_(void);
                void set_dS_and_trdS_(void);
                void set_ddS_and_trddS_(void);
                void set_iter_trS_(Real lambdaS);

                // UTILITIES
                void LeftMultiplybyPsiAndTrace(Real & trace, MatrixXr & ret, const MatrixXr & mat);

                // GLOBAL UPDATERS
                void update_matrices(lambda::type<1> lambda);

        public:
                // CONSTRUCTORS
                //! Constructor of the class given the InputCarrier
                /*!
                 \param the_carrier the structure from which to take all the data for the derived classes
                 \sa set_R_()
                */
                GCV_Exact<InputCarrier, 1>(InputCarrier & the_carrier_):
                        GCV_Family<InputCarrier, 1>(the_carrier_)
                        {
                                this->set_R_(); // this matrix is unchanged during the whole procedure, thus it's set once and for all
                        }
                //! Constructor of the class given the InputCarrier and lambdaT (for Temporal and parabolic case)
                /*!
                \param the_carrier the structure from which to take all the data for the derived classes
                \param lambdaT parameter of the PDE
                 \sa set_R_()
                 \sa set_lambdaT()
                */
                GCV_Exact<InputCarrier, 1>(InputCarrier & the_carrier_, Real lambdaT_):
                        GCV_Family<InputCarrier, 1>(the_carrier_)
                        {
                        	/*
					This is the constructor for the parabolic spatio-temporal case.
					In the iterative case, the R setter is the same as for the space-only
					problem. In the monolithic solution, we have an ad-hoc setter
				*/
				if(the_carrier_.get_model()->isIter())
					this->set_R_();
				else
					this->set_R_(lambdaT_);
                                this->lambdaT = lambdaT_;
                        }

                // PUBLIC UPDATERS
                //inline void set_lambdaT(Real lambdaT_){this->set_R_(lambdaT_);} // parabolic case

                void update_parameters(lambda::type<1> lambda) override;

                void first_updater(lambda::type<1> lambda);
                void second_updater(lambda::type<1> lambda);

                // GCV-COMPUTATION
                Real compute_f( lambda::type<1> lambda) override;
                Real compute_fp(lambda::type<1> lambda);
                Real compute_fs(lambda::type<1> lambda);
                //! Virtual Destuctor
        virtual ~GCV_Exact(){};
};


//! Derived class used for bidimensional lambda exact gcv-based methods
/*!
 This class implements exact methods to identify the correct value of the gcv
 and possibly of its derivatives, by partially solving manually the apply system
 \tparam InputCarrier Carrier-type parameter that contains insight about the problem to be solved
*/
template<typename InputCarrier>
class GCV_Exact<InputCarrier, 2>: public GCV_Family<InputCarrier, 2>
{
        private:
                //! An external updater whose purpose is keeping the internal values coherent with the computations to be made from time to time
                GOF_updater<GCV_Exact<InputCarrier, 2>, lambda::type<2>> gu;

                // INTERNAL DATA STRUCTURES
                MatrixXr  R_;           //!< stores the value of R1^t*R0^{-1}*R1 [size nnodes x nnodes]
                MatrixXr  T_;           //!< stores the value of Psi^t*Q*Psi+lambda*R [size nnodes x nnodes]
                MatrixXr  V_;           //!< stores the value of T^{-1}*Psi^t*Q [size nnodes x s]
                MatrixXr  S_;           //!< stores the value of Psi*V [as in Stu-Hunter Sangalli] [size s x s]
                Real      trS_ = 0.0;   //!< stores the value of the trace of S
                MatrixXr  dS_;          //!< stores the derivative of S w.r.t. lambda [size s x s]
                Real      trdS_ = 0.0 ; //!< stores the value of the trace of dS
                MatrixXr  ddS_;         //!< stores the second derivative of S w.r.t. lambdaS [size s x s]
                Real      trddS_ = 0.0; //!< stores the value of the trace of ddS
                MatrixXr  time_dS_;           //!< stores the derivative of S w.r.t. lambdaT [size s x s]
                Real      time_trdS_ = 0.0 ;  //!< stores the value of the trace of time_dS
                MatrixXr  time_ddS_;          //!< stores the second derivative of S w.r.t. lambdaT [size s x s]
                Real      time_trddS_ = 0.0;  //!< stores the value of the trace of time_ddS
                MatrixXr  time_ddS_mxd_ ;          //!< stores the second derivative of S w.r.t. lambdaS and w.r.t lambdaT [size s x s]
                Real      time_trddS_mxd_ = 0.0;   //!< stores the value of the trace of time_ddS_mxd

                //! Additional utility matrices [just the ones for the specific carrier that is proper of the problem]
                AuxiliaryData<InputCarrier> adt;
                AuxiliaryData<InputCarrier> time_adt;

                // COMPUTERS and DOF methods
                void compute_z_hat (lambda::type<2> lambda) override;
                void update_dof(lambda::type<2> lambda)     override;
                void update_dor(lambda::type<2> lambda)     override;

                // SETTERS
                void set_R_(void);                      // separable case
                void set_T_(lambda::type<2> lambda);
                void set_V_(void);
                void set_S_and_trS_(void);
                void set_dS_and_trdS_(void);
                void set_ddS_and_trddS_(void);
                void set_ddS_and_trddS_mxd_(void);
                
                // UTILITIES
                void LeftMultiplybyPsiAndTrace(Real & trace, MatrixXr & ret, const MatrixXr & mat);

                // GLOBAL UPDATERS
                void update_matrices(lambda::type<2> lambda);

        public:
                // CONSTRUCTORS
                //! Constructor of the class given the InputCarrier
                /*!
                 \param the_carrier the structure from which to take all the data for the derived classes
                 \sa set_R_()
                */
                GCV_Exact<InputCarrier, 2>(InputCarrier & the_carrier_):
                        GCV_Family<InputCarrier, 2>(the_carrier_)
                        {
                                this->set_R_();         // this matrix is unchanged during the whole procedure, thus it's set once and for all
                                time_adt.flag_time = true;      // GCV_Exact<InputCarrier, 2> is used only in separable case
                        }

                // PUBLIC UPDATERS
                void update_parameters(lambda::type<2> lambda) override;

                void first_updater(lambda::type<2> lambda);
                void second_updater(lambda::type<2> lambda);

                // GCV-COMPUTATION
                Real compute_f(lambda::type<2> lambda) override;
                lambda::type<2> compute_fp(lambda::type<2> lambda);
                MatrixXr compute_fs(lambda::type<2> lambda);
                //! Virtual Destuctor
        virtual ~GCV_Exact(){};
};

//----------------------------------------------------------------------------//
// ** GCV_STOCHASTIC **

//! Derived class used for multidimensional lambda stochastic gcv-based methods
/*!
 This class implements an efficent stochastic approximation method to identify
 a soound approxiamtion of thee gcv, the method involved is known as Stochastic
 Woodbury decomposition and directly takes advantage of the function apply proper
 of any carrier to perform its computations.
 \tparam InputCarrier Carrier-type parameter that contains insight about the problem to be solved
 \tparam size specialization parameter used to characterize the size of the lambda to be used
*/
template<typename InputCarrier, UInt size>
class GCV_Stochastic: public GCV_Family<InputCarrier, size>
{
        private:
                //! An external updater whose purpose is keeping the internal values coherent with the computations to be made from time to time
                GOF_updater<GCV_Stochastic<InputCarrier, size>, lambda::type<size>> gu;

                // INTERNAL DATA STRUCTURES
                MatrixXr US_;           //!< binary{+1/-1} random matrix used for stochastic gcv computations [size s x #realizations]
                MatrixXr USTpsi;       //!< US^T*Psi
                MatrixXr b;             //! Right hand side o solution
                bool     us = false;    //!< keeps track of US_ matrix being already computed or not

                // COMPUTERS and DOF methods
                void compute_z_hat (lambda::type<size> lambda) override;
                void update_dof(lambda::type<size> lambda)     override;
                void update_dor(lambda::type<size> lambda)     override;

                // SETTERS
                void set_US_(void);
                
                Real lambdaT;

        public:
                // CONSTRUCTORS
                //! Constructor of the class given the InputCarrier, also computes the US_ matrix
                /*!
                 \param the_carrier the structure from which to take all the data for the derived classes
                 \sa set_US_()
                */
                GCV_Stochastic(InputCarrier & the_carrier_, bool flag_used):
                        GCV_Family<InputCarrier, size>(the_carrier_)
                        {
                                MatrixXr m = this->the_carrier.get_opt_data()->get_DOF_matrix();
                                if(m.cols()>0 && m.rows()>0 && flag_used)
                                {
                                        this->set_US_(); // this matrix is unchanged during the whole procedure, thus it's set once and for all
                                }
                        }
                        
                 //template<typename std::enable_if<size==2, bool>::type = true>
                 GCV_Stochastic(InputCarrier & the_carrier_, Real lambdaT_): //parabolic case
                        GCV_Stochastic(the_carrier_, true)
                        {
                        	this->lambdaT = lambdaT_;
                        }

                // PUBLIC UPDATERS
                void update_parameters(lambda::type<size> lambda) override;

                void first_updater(lambda::type<size> lambda)  {; /*Dummy*/} //!< Dummy function needed for consistency of the external updater
                void second_updater(lambda::type<size> lambda) {; /*Dummy*/} //!< Dummy function needed for consistency of the external updater

                // GCV-COMPUTATION
                Real compute_f( lambda::type<size> lambda) override;
                lambda::type<size> compute_fp(lambda::type<size> lambda) {return lambda::init<size>(-1); /*Dummy*/} //!< Dummy function needed for consistency of the external updater
                typename std::conditional<size==1, Real, MatrixXr>::type compute_fs(lambda::type<size> lambda)
                        {return lambda::init<size>(-1); /*Dummy*/} //!< Dummy function needed for consistency of the external updater
                //! Virtual Destuctor
        virtual ~GCV_Stochastic(){};
};

//----------------------------------------------------------------------------//
// *** K-FOLD CV ***

template<typename InputCarrier, UInt size>
class K_fold_CV: public Lambda_optimizer<InputCarrier, size>
{
/*
         [[ TO BE IMPLEMENTED ]]
*/
};

#include "Lambda_Optimizer_imp.h"

#endif

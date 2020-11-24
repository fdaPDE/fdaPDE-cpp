#ifndef __LAMBDA_OPTIMIZER_H__
#define __LAMBDA_OPTIMIZER_H__

// HEADERS
#include "../../FdaPDE.h"
#include "Auxiliary_Optimizer.h"
#include "Carrier.h"
#include "Gof_Updater.h"
#include "../../FE_Assemblers_Solvers/Include/Solver.h"
#include <algorithm>

// CLASSES
// **** GENERAL METHODS ***

//! Father class used for multidiemnsional lambda optimization purposes
/*!
 This virtual class stores the model to be used by all its children,
 i. e. the classes actually instantiating a real optimization method and
 performing evaluations.
 \tparam InputCarrier Carrier-type parameter that contains insight about the problem to be solved
 \tparam size specialization parameter used to characterize the size of the lambda to be used
 \todo 2-D lambda optimization still to be implemented
*/
template <typename InputCarrier, UInt size>
class Lambda_optimizer
{
/*
        [[ VERSION WITH TIMES STILL TO BE IMPLEMENTED ]]
*/
};


//! Father class used for unidimenional lambda optimization purposes
/*!
 This virtual class stores the model to be used by all its children,
 i. e. the classes actually instantiating a real optimization method and
 performing evaluations. Specialized version for unidimensional problems.
 \tparam InputCarrier Carrier-type parameter that contains insight about the problem to be solved
*/
template <typename InputCarrier>
class Lambda_optimizer<InputCarrier, 1>
{
        protected:
                //! Model containing all the information necessary for the computation of the optimal value
                InputCarrier & the_carrier;

                // CONSTRUCTORS
                //! Constructor of the class given the InputCarrier
                /*!
                 \param the_carrier the structure from which to take all the data for the derived classes
                */
                Lambda_optimizer<InputCarrier, 1>(InputCarrier & the_carrier_):
                        the_carrier(the_carrier_) {}

                // UPDATERS
                //! A pure virtual member used in children classes for updates of internal data
                /*!
                 \param lambda the value of lambda with which to perform the update
                */
        virtual void update_parameters(Real lambda) = 0;

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
 \todo 2-D lambda optimization still to be implemented
*/
template <typename InputCarrier, UInt size>
class GCV_Family: public Lambda_optimizer<InputCarrier, size>
{
/*
        [[ VERSION WITH TIMES STILL TO BE IMPLEMENTED ]]
*/
};

//! Father class used for unidimensional lambda gcv-based methods
/*!
 This virtual class inherits from the generic unidimensional optimizer Lambda_optimizer
 and contains the methods used to compute the main statistics
 related to the optimization problem, like predicted output, errors and variances
 also stores a container to retrieve the data and pass them to the user.
 Specialized version for unidimensional problems.
 \tparam InputCarrier Carrier-type parameter that contains insight about the problem to be solved
*/
template <typename InputCarrier>
class GCV_Family<InputCarrier, 1>: Lambda_optimizer<InputCarrier, 1>
{
        protected:
                //! Model containing all the information necessary for the computation of the optimal value
                using  Lambda_optimizer<InputCarrier, 1>::the_carrier;

                // Output data
                VectorXr        z_hat;                  //!< Model predicted values in the locations [size s]
                VectorXr        eps_hat;                //!< Model predicted error in the locations (residuals) [size s]
                Real            SS_res = 0.0;           //!< Model predicted sum of squares of the residuals
                Real            rmse = 0.0;             //!< Model root mean squared error
                Real            sigma_hat_sq = 0.0;     //!< Model estimated variance of error
                UInt            s;                      //!< Model number of observations (i.e. #locations)
                output_Data     output;                 //!< Output, needed to be user-available, necessarily public

                // Degrees of freedom
                Real            dof = 0.0;              //!< tr(S) + q, degrees of freedom of the model
                Real            dor = 0.0;              //!< s - dof, degrees of freedom of the residuals

                UInt            use_index = -1;         //!< Index of the DOF_matrix to be used, if non empty

                // SETTERS of the putput data
        virtual void compute_z_hat(Real lambda) = 0;    //!< Utility to compute the size of predicted value in the locations
                void compute_z_hat_from_f_hat(const VectorXr & f_hat);
                void compute_eps_hat(void);
                void compute_SS_res(void);
                void compute_rmse(void);
                void compute_sigma_hat_sq(void);
                void compute_s(void);

                // UPDATERS
                void update_errors(Real lambda);

                // DOF methods
        virtual void update_dof(Real lambda) = 0;       //!< Utility to compute the degrees of freedom of the model
        virtual void update_dor(Real lambda) = 0;       //!< Utility to compute the degrees of freedom of the residuals

                // CONSTRUCTORS
                //! Constructor of the class given the InputCarrier
                /*!
                 \param the_carrier the structure from which to take all the data for the derived classes
                 \sa compute_s()
                */
                GCV_Family<InputCarrier, 1>(InputCarrier & the_carrier_):
                        Lambda_optimizer<InputCarrier, 1>(the_carrier_)
                        {
                                this->compute_s();      // stores immediately the number of locations
                        }

        public:
                // UTILITY FOR DOF MATRIX
        inline  void set_index(UInt index){this->use_index = index;}

                // PUBLIC UPDATERS
        virtual void update_parameters(Real lambda) = 0; //!< Utility to update all the prameters of the model

                void zero_updater(Real lambda);

                // GCV-COMPUTATION
        virtual Real compute_f( Real lambda) = 0;       //!< Main function, represents the gcv computation

                // OUTPUT MANAGERS
                output_Data  get_output(std::pair<Real, UInt> optimal_pair, const timespec & time_count, const std::vector<Real> & GCV_v, const std::vector<Real> & lambda_v, int termination_);
                void set_output_partial_best(void);
                output_Data get_output_full(void);
                void set_output_partial(void);
                void combine_output_prediction(const VectorXr & f_hat, output_Data & outp, UInt cols);
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

                //! Additional utility matrices [just the ones for the specific carrier that is proper of the problem]
                AuxiliaryData<InputCarrier> adt;

                // COMPUTERS and DOF methods
                void compute_z_hat (Real lambda) override;
                void update_dof(Real lambda)     override;
                void update_dor(Real lambda)     override;

                // SETTERS
                void set_R_(void);
                void set_T_(Real lambda);
                void set_V_(void);
                void set_S_and_trS_(void);
                void set_dS_and_trdS_(void);
                void set_ddS_and_trddS_(void);

                // UTILITIES
                void LeftMultiplybyPsiAndTrace(Real & trace, MatrixXr & ret, const MatrixXr & mat);

                // GLOBAL UPDATERS
                void update_matrices(Real lambda);

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

                // PUBLIC UPDATERS
                void update_parameters(Real lambda) override;

                void first_updater(Real lambda);
                void second_updater(Real lambda);

                // GCV-COMPUTATION
                Real compute_f( Real lambda) override;
                Real compute_fp(Real lambda);
                Real compute_fs(Real lambda);
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
 \todo 2-D lambda optimization still to be implemented
*/
template<typename InputCarrier, UInt size>
class GCV_Stochastic: public GCV_Family<InputCarrier, size>
{
/*
        [[ VERSION WITH TIMES STILL TO BE IMPLEMENTED ]]
*/
};

//! Derived class used for unidimensional lambda stochastic gcv-based methods
/*!
 This class implements an efficent stochastic approximation method to identify
 a sound approximation of the gcv. The algorithm involved is known as Stochastic
 Woodbury decomposition and directly takes advantage of the function apply proper
 of any Carrier to perform its computations. This template is a specialization for
 the unidimensional case.
 \tparam InputCarrier Carrier-type parameter that contains insight about the problem to be solved
*/
template<typename InputCarrier>
class GCV_Stochastic<InputCarrier, 1>: public GCV_Family<InputCarrier, 1>
{
        private:
                //! An external updater whose purpose is keeping the internal values coherent with the computations to be made from time to time
                GOF_updater<GCV_Stochastic<InputCarrier, 1>, Real> gu;

                // INTERNAL DATA STRUCTURES
                MatrixXr US_;           //!< binary{+1/-1} random matrix used for stochastic gcv computations [size s x #realizations]
                MatrixXr USTpsi;       //!< US^T*Psi
                MatrixXr b;             //! Right hand side o solution
                bool     us = false;    //!< keeps track of US_ matrix being already computed or not

                // COMPUTERS and DOF methods
                void compute_z_hat (Real lambda) override;
                void update_dof(Real lambda)     override;
                void update_dor(Real lambda)     override;

                // SETTERS
                void set_US_(void);

        public:
                // CONSTRUCTORS
                //! Constructor of the class given the InputCarrier, also computes the US_ matrix
                /*!
                 \param the_carrier the structure from which to take all the data for the derived classes
                 \sa set_US_()
                */
                GCV_Stochastic<InputCarrier, 1>(InputCarrier & the_carrier_, bool flag_used):
                        GCV_Family<InputCarrier, 1>(the_carrier_)
                        {
                                MatrixXr m = this->the_carrier.get_opt_data()->get_DOF_matrix();
                                if(m.cols()>0 && m.rows()>0 && flag_used)
                                {
                                        this->set_US_(); // this matrix is unchanged during the whole procedure, thus it's set once and for all
                                }
                        }

                // PUBLIC UPDATERS
                void update_parameters(Real lambda) override;

                void first_updater(Real lambda)  {; /*Dummy*/} //!< Dummy function needed for consistency of the external updater
                void second_updater(Real lambda) {; /*Dummy*/} //!< Dummy function needed for consistency of the external updater

                // GCV-COMPUTATION
                Real compute_f( Real lambda) override;
                Real compute_fp(Real lambda) {return 0; /*Dummy*/} //!< Dummy function needed for consistency of the external updater
                Real compute_fs(Real lambda) {return 0; /*Dummy*/} //!< Dummy function needed for consistency of the external updater
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

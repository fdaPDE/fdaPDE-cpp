#ifndef __AUXILIARY_OPTIMIZER_H__
#define __AUXILIARY_OPTIMIZER_H__

// HEADERS
#include <functional>
#include <string>
#include "../../FdaPDE.h"
#include "../../FE_Assemblers_Solvers/Include/Solver.h"
#include "../../Global_Utilities/Include/Solver_Definitions.h"
#include "Carrier.h"
#include "Solution_Builders.h"

// CLASSES
//! Template class for data storing and efficient management.
/*!
 General purpose class storing data useful for fastening computation in
 Lambda_optimizer derived classes and AuxiliaryOptimizer. Its content
 are matrices, vectors and doubles useful for GCV calculations.
 \tparam InputCarrier the type of Carrier used in the optimization.
 \tparam Enable dummy typename for SFINAE instantiation of a more refined version for problems with forcing terms.
 \sa Carrier, AuxiliaryOptimizer, Lambda_optimizer
*/
template<typename InputCarrier, typename Enable = void>
struct AuxiliaryData
{
        MatrixXr K_;                            //!< Stores T^{-1}*R                            [nnodes x nnodes]
        MatrixXr F_;                            //!< Stores K*v                                 [nnodes x nnodes]
        VectorXr t_;                            //!< Stores dS*z;
        Real     a_;                            //!< Stores <eps_hat, dS*z>
        Real     b_;                            //!< Stores <t, Q*t>
        Real     c_;                            //!< Stores <eps_hat, ddS*z>
};

//! Template class for data storing and efficient management in forcing term based problems
/*!
 General purpose class storing data useful for fastening computation in
 Lambda_optimizer derived classes and AuxiliaryOptimizer, spercialized under forcing
 term based problems. Its content are matrices, vectors and doubles
 useful for GCV calculations.
 \tparam InputCarrier the type of Carrier used in the optimization.
 \tparam typename for SFINAE instantiation of this more refined version for problems with forcing terms, distinguishing from base version.
 \sa AuxiliaryData, Carrier, AuxiliaryOptimizer, Lambda_optimizer
*/
template<typename InputCarrier>
struct AuxiliaryData<InputCarrier, typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>, t_type>::value>::type>
{
        MatrixXr K_;                            //!< Stores T^{-1}*R                            [nnodes x nnodes]
        MatrixXr F_;                            //!< Stores K*v                                 [nnodes x nnodes]
        VectorXr t_;                            //!< Stores dS*z;
        Real     a_;                            //!< Stores the value of <eps_hat, dS*z>
        Real     b_;                            //!< Stores <t, Q*t>
        Real     c_;                            //!< Stores <eps_hat, ddS*z>
        VectorXr f_;                            //!< Stores R1^T*R0^{-1}*u
        VectorXr g_;                            //!< Stores T^{-1}*f
        VectorXr h_;                            //!< Stores (lambda*K-I)*g
        VectorXr p_;                            //!< Stores Psi*h-t
        VectorXr r_;                            //!< Stores Q*s

        void left_multiply_by_psi(const InputCarrier & carrier, VectorXr & ret, const VectorXr & vec);
};


//! General purpose class to support efficient case-driven computation of Lambda_Optimizer
/*!
 This struct is a collection of static methods called "universal" in their name
 whose main purpose is to use SFINAE on the template InputCarrier type to provide
 correct implementation for each possible method in Lambda_Optimizer derived classes.
 Since functions are static no object of this class needs to be ever created
 \sa AuxiliaryData,  Lambda_optimizer
*/
struct AuxiliaryOptimizer
{
        static void bc_utility(MatrixXr & mat, const std::vector<UInt> * bc_idxp);
        static void bc_utility(SpMat & mat, const std::vector<UInt> * bc_idxp);
        /* -------------------------------------------------------------------*/

        //! SFINAE based method to compute matrix R in case of Forced problem
        /*!
         \param R a reference to the matrix to be computed
         \param carrier the Carrier-type object containing the data
         \param adt the AuxiliaryData type to store useful byproducts of R matrix computation
         \return an integer signaling the correct ending of the process
         \note AuxiliaryOptimizer f_ is computed at this level
        */
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>, t_type>::value, UInt>::type
                universal_R_setter(MatrixXr & R, const InputCarrier & carrier, AuxiliaryData<InputCarrier> & adt);

        //! SFINAE based method to compute matrix R in case of non-Forced problem
        /*!
         \param R a reference to the matrix to be computed
         \param carrier the Carrier-type object containing the data
         \param adt the AuxiliaryData type to store useful byproducts of R matrix computation
         \return an integer signaling the correct ending of the process
        */
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>, f_type>::value, UInt>::type
                universal_R_setter(MatrixXr & R, const InputCarrier & carrier, AuxiliaryData<InputCarrier> & adt);
        /* -------------------------------------------------------------------*/

        //! SFINAE based method to compute matrix T in case of Areal problem
        /*!
         \param T a reference to the matrix to be computed
         \param carrier the Carrier-type object containing the data
         \return an integer signaling the correct ending of the process
        */
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>, t_type>::value, UInt>::type
                universal_T_setter(MatrixXr & T, InputCarrier & carrier);

        //! SFINAE based method to compute matrix T in case of pointwise problem
        /*!
         \param T a reference to the matrix to be computed
         \param carrier the Carrier-type object containing the data
         \return an integer signaling the correct ending of the process
        */
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>, f_type>::value, UInt>::type
                universal_T_setter(MatrixXr & T, InputCarrier & carrier);
        /* -------------------------------------------------------------------*/

        //! SFINAE based method to compute matrix V in case of Forced problem
        /*!
         \param V a reference to the matrix to be computed
         \param T a const reference to matrix T
         \param R a const reference to matrix R
         \param carrier the Carrier-type object containing the data
         \param adt the AuxiliaryData type to store useful byproducts of V matrix computation
         \return an integer signaling the correct ending of the process
         \note AuxiliaryOptimizer K_ and  g_ are computed at this level
        */
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>, t_type>::value, UInt>::type
                universal_V_setter(MatrixXr & V, const MatrixXr & T, const MatrixXr & R, const InputCarrier & carrier, AuxiliaryData<InputCarrier> & adt);

        //! SFINAE based method to compute matrix V in case of non-Forced problem
        /*!
         \param V a reference to the matrix to be computed
         \param T a const reference to matrix T
         \param R a const reference to matrix R
         \param carrier the Carrier-type object containing the data
         \param adt the AuxiliaryData type to store useful byproducts of V matrix computation
         \return an integer signaling the correct ending of the process
         \note AuxiliaryOptimizer K_ is computed at this level
        */
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>, f_type>::value, UInt>::type
                universal_V_setter(MatrixXr & V, const MatrixXr & T, const MatrixXr & R, const InputCarrier & carrier, AuxiliaryData<InputCarrier> & adt);
        /* -------------------------------------------------------------------*/

        //! SFINAE based method to compute matrix E in case of Areal problem
        /*!
         \param E a reference to the matrix to be computed
         \param carrier the Carrier-type object containing the data
         \return an integer signaling the correct ending of the process
        */
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>, t_type>::value, UInt>::type
                universal_E_setter(MatrixXr & E, const InputCarrier & carrier);

        //! SFINAE based method to compute matrix E in case of pointwise problem
        /*!
         \param E a reference to the matrix to be computed
         \param carrier the Carrier-type object containing the data
         \return an integer signaling the correct ending of the process
        */
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>, f_type>::value, UInt>::type
                universal_E_setter(MatrixXr & E, const InputCarrier & carrier);

        static void set_E_ln_W_ptw(MatrixXr & E, const std::vector<UInt> * kp, const MatrixXr * Qp, UInt nr, UInt s);
        static void set_E_lnn_W_ptw(MatrixXr & E, const SpMat * psi_tp, const MatrixXr * Qp);
        static void set_E_W_a(MatrixXr & E, const SpMat * psi_tp, const MatrixXr * Qp, const VectorXr * Ap);
        static void set_E_nW_a(MatrixXr & E, const SpMat * psi_tp, const VectorXr * Ap);
        /* -------------------------------------------------------------------*/

        //! SFINAE based method to compute predictions in locations in case of Forced problem
        /*!
         \param z_hat a reference to the VectorXr to be computed
         \param carrier the Carrier-type object containing the data
         \param S a const reference to matrix S
         \param adt the AuxiliaryData type to store useful byproducts of z_hat computation
         \param lambda the Real datum used as smoothing parameter
         \return an integer signaling the correct ending of the process
         \note AuxiliaryOptimizer r_ is computed at this level
        */
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>, t_type>::value, UInt>::type
                universal_z_hat_setter(VectorXr & z_hat, InputCarrier & carrier, const MatrixXr & S, AuxiliaryData<InputCarrier> & adt, const Real lambda);

        //! SFINAE based method to compute predictions in locations in case of non-Forced problem
        /*!
         \param z_hat a reference to the VectorXr to be computed
         \param carrier the Carrier-type object containing the data
         \param S a const reference to matrix S
         \param adt the AuxiliaryData type to store useful byproducts of z_hat computation
         \param lambda the Real datum used as smoothing parameter
         \return an integer signaling the correct ending of the process
        */
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>, f_type>::value, UInt>::type
                universal_z_hat_setter(VectorXr & z_hat, InputCarrier & carrier, const MatrixXr & S, AuxiliaryData<InputCarrier> & adt, const Real lambda);

        //! Utility to compute the common part of universal_z_hat_setter among Forced and non-Forced problems
        /*!
         \param z_hat a reference to the VectorXr to be computed
         \param carrier the Carrier-type object containing the data
         \param S a const reference to matrix S
        */
        template<typename InputCarrier>
        static void common_z_hat_part(VectorXr & z_hat, InputCarrier & carrier, const MatrixXr & S);
        /* -------------------------------------------------------------------*/

        //! SFINAE based method to compute right hand term for stochastic dof evaluation, areal type
        /*!
         \param b a reference to the MatrixXr to be computed
         \param carrier the Carrier-type object containing the data
         \param US a stochastic matrix used for purpose of computing stochastic dofs
         \param nnodes number of nodes of the mesh
         \return an integer signaling the correct ending of the process
        */
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>, t_type>::value, UInt>::type
                universal_b_setter(MatrixXr & b, InputCarrier & carrier, const MatrixXr & US, const UInt nnodes);

        //! SFINAE based method to compute right hand term for stochastic dof evaluation, pointwise type
        /*!
         \param b a reference to the MatrixXr to be computed
         \param carrier the Carrier-type object containing the data
         \param US a stochastic matrix used for purpose of computing stochastic dofs
         \param nnodes number of nodes of the mesh
         \return an integer signaling the correct ending of the process
        */
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>, f_type>::value, UInt>::type
                universal_b_setter(MatrixXr & b, InputCarrier & carrier, const MatrixXr & US, const UInt nnodes);
        /* -------------------------------------------------------------------*/

        //! SFINAE based method: general updater of first derivative for forcing term data
        /*!
         \param adt the AuxiliaryData type to store useful terms
         \param carrier the Carrier-type object containing the data
         \param dS const reference of the derivative matrix of S
         \param eps const reference of the error
         \param lambda smoothing parameter for which the update has to be performed
         \return an integer signaling the correct ending of the process
         \note this function updates t_,h_,p_ and a_ of AuxiliaryData
        */
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>, t_type>::value, UInt>::type
                universal_first_updater(AuxiliaryData<InputCarrier> & adt, const InputCarrier & carrier, const MatrixXr & dS, const VectorXr & eps, const Real lambda);

        //! SFINAE based method: general updater of first derivative for non-Forced data
        /*!
         \param adt the AuxiliaryData type to store useful terms
         \param carrier the Carrier-type object containing the data
         \param dS const reference of the derivative matrix of S
         \param eps const reference of the error
         \param lambda smoothing parameter for which the update has to be performed
         \return an integer signaling the correct ending of the process
         \note this function updates t_, and a_ of AuxiliaryData
        */
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>, f_type>::value, UInt>::type
                universal_first_updater(AuxiliaryData<InputCarrier> & adt, const InputCarrier & carrier, const MatrixXr & dS, const VectorXr & eps, const Real lambda);
        /* -------------------------------------------------------------------*/

        //! SFINAE based method: general updater of second derivative for Forced data
        /*!
         \param adt the AuxiliaryData type to store useful terms
         \param carrier the Carrier-type object containing the data
         \param ddS const reference of the second derivative matrix of S
         \param eps const refernce of the error
         \param lambda smoothing parameter for which the update has to be performed
         \return an integer signaling the correct ending of the process
         \note this function updates b_, and c_ of AuxiliaryData
        */
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>, t_type>::value, UInt>::type
                universal_second_updater(AuxiliaryData<InputCarrier> & adt, InputCarrier & carrier, const MatrixXr & ddS, const VectorXr & eps, const Real lambda);

        //! SFINAE based method: general updater of second derivative for non-Forced data
        /*!
         \param adt the AuxiliaryData type to store useful terms
         \param carrier the Carrier-type object containing the data
         \param ddS const reference of the second derivative matrix of S
         \param eps const refernce of the error
         \param lambda smoothing parameter for which the update has to be performed
         \return an integer signaling the correct ending of the process
         \note this function updates b_, and c_ of AuxiliaryData
        */
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>, f_type>::value, UInt>::type
                universal_second_updater(AuxiliaryData<InputCarrier> & adt, InputCarrier & carrier, const MatrixXr & ddS, const VectorXr & eps, const Real lambda);
        /* -------------------------------------------------------------------*/

        //! SFINAE based method for purpose of gcv coputation (NOW FAKE SFINAE)
        /*!
         \param s number of observations
         \param sigma_hat_sq esimated variance of the error
         \param dor (s-dof)
         \return the value of the GCV
         \todo dependence on template might be removed if you can extend symmetry of terms also to temporal data
         \note the SFINAE use in this case is fake since we have created a perfect symmetry between the forcing term and non-forcing term cases, still left for eventual Temporal
        */
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<t_type,t_type>::value, Real>::type
                universal_GCV(const Real s, const Real sigma_hat_sq, const Real dor);
        /* -------------------------------------------------------------------*/

        //! SFINAE based method for purpose of gcv first derivative coputation (NOW FAKE SFINAE)
        /*!
         \param adt the AuxiliaryData type to get useful terms
         \param s number of observations
         \param sigma_hat_sq esimated variance of the error
         \param dor (s-dof)
         \param trdS trace of dS matrix
         \return the value of the GCV first derivative
         \todo dependence on template might be removed if you can extend symmetry of terms also to temporal data
         \note the SFINAE use in this case is fake since we have created a perfect symmetry between the forcing term and non-forcing term cases, still left for eventual Temporal
        */
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<t_type, t_type>::value, Real>::type
                universal_GCV_d(const AuxiliaryData<InputCarrier> & adt, const Real s, const Real sigma_hat_sq, const Real dor, const Real trdS);

        /* -------------------------------------------------------------------*/

        //! SFINAE based method for purpose of gcv second derivative coputation (NOW FAKE SFINAE)
        /*!
         \param adt the AuxiliaryData type to get useful terms
         \param s number of observations
         \param sigma_hat_sq esimated variance of the error
         \param dor (s-dof)
         \param trdS trace of dS matrix
         \param trddS trace of ddS matrix
         \return the value of the GCV second derivative
         \todo dependence on template might be removed if you can extend symmetry of terms also to temporal data
         \note the SFINAE use in this case is fake since we have created a perfect symmetry between the forcing term and non-forcing term cases, still left for eventual Temporal
        */
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<t_type, t_type>::value, Real>::type
                universal_GCV_dd(const AuxiliaryData<InputCarrier> & adt, const Real s, const Real sigma_hat_sq, const Real dor, const Real trdS, const Real trddS);
};

#include "Auxiliary_Optimizer_imp.h"

#endif

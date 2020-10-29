#ifndef __AUXILIARY_OPTIMIZER_IMP_H__
#define __AUXILIARY_OPTIMIZER_IMP_H__

//! Utility used for efficient multiplication of a vector for matrix Psi
/*!
 \param carrier the Carrier-type object containing the data
 \param ret the vector where to store the result, passed by reference
 \param vec the vector to be left multiplied by Psi
*/
template<typename InputCarrier>
void AuxiliaryData<InputCarrier, typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value>::type>::left_multiply_by_psi(const InputCarrier & carrier, VectorXr & ret, const VectorXr & vec)
{
        if (carrier.loc_are_nodes())
        {
                const UInt s = carrier.get_n_obs();
                ret = VectorXr::Zero(s);

                const std::vector<UInt> * kp = carrier.get_obs_indicesp();
                for (UInt i = 0; i < s; i++)
                        ret.coeffRef(i) += vec.coeff((*kp)[i]);
        }
        else
        {
                // Psi is full
                ret = (*carrier.get_psip())*vec;
        }
}

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_R_setter(MatrixXr & R, const InputCarrier & carrier, AuxiliaryData<InputCarrier> & adt)
        {
                SpMat  R1p_= *carrier.get_R1p();         // Get the value of matrix R1
                const std::vector<UInt> * bc_indices = carrier.get_bc_indicesp();
                AuxiliaryOptimizer::bc_utility(R1p_, bc_indices);

                Eigen::SparseLU<SpMat> factorized_R0p(*(carrier.get_R0p()));
                R = (R1p_).transpose()*factorized_R0p.solve(R1p_);     // R == _R1^t*R0^{-1}*R1
                adt.f_ = ((R1p_).transpose())*factorized_R0p.solve((*carrier.get_up()));

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,f_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_R_setter(MatrixXr & R, const InputCarrier & carrier, AuxiliaryData<InputCarrier> & adt)
        {
                SpMat  R1p_= *carrier.get_R1p();         // Get the value of matrix R1
                const std::vector<UInt> * bc_indices = carrier.get_bc_indicesp();
                AuxiliaryOptimizer::bc_utility(R1p_, bc_indices);

                Eigen::SparseLU<SpMat>factorized_R0p(*(carrier.get_R0p()));
                R = (R1p_).transpose()*factorized_R0p.solve(R1p_);     // R == _R1^t*R0^{-1}*R1

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_T_setter(MatrixXr & T, InputCarrier & carrier)
        {
                const VectorXr * Ap = carrier.get_Ap();  // For areal data
                const SpMat * psip = carrier.get_psip();
                const SpMat * psi_tp = carrier.get_psi_tp();
                const std::vector<UInt> * bc_idxp = carrier.get_bc_indicesp();

                MatrixXr aux = (*psi_tp)*(*Ap).asDiagonal()*carrier.lmbQ(*psip);
                AuxiliaryOptimizer::bc_utility(aux, bc_idxp);
                T += aux; // Add correction

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,f_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_T_setter(MatrixXr & T, InputCarrier & carrier)
        {
                const SpMat * psip = carrier.get_psip();
                const SpMat * psi_tp = carrier.get_psi_tp();
                const std::vector<UInt> * bc_idxp = carrier.get_bc_indicesp();

                MatrixXr aux = (*psi_tp)*carrier.lmbQ(*psip);
                AuxiliaryOptimizer::bc_utility(aux, bc_idxp);
                T += aux; // Add correction

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_V_setter(MatrixXr & V, const MatrixXr & T, const MatrixXr & R, const InputCarrier & carrier, AuxiliaryData<InputCarrier> & adt)
        {
                // Based on the type of problem at hand select at compile time the most useful Factorizer
                typedef typename std::conditional<std::is_base_of<Areal, InputCarrier>::value, Eigen::PartialPivLU<MatrixXr>, Eigen::LDLT<MatrixXr>>::type Factorizer;
                Factorizer factorized_T(T);      // define a factorization of the defined type

                if(!carrier.is_areal() && !carrier.has_W())
                {
                        // Q == I
                        const SpMat * psi_tp = carrier.get_psi_tp();
                        V = factorized_T.solve(MatrixXr(*psi_tp));      // find the value of V = T^{-1}*Psi^t
                }
                else
                {
                        MatrixXr E_;    // Declare an empty auxiliary matrix
                        const UInt ret =  AuxiliaryOptimizer::universal_E_setter<InputCarrier>(E_, carrier);
                        V = factorized_T.solve(E_);     // find the value of V = T^{-1}*E
                }
                adt.K_ = factorized_T.solve(R);         // K = T^{-1}*R
                adt.g_ = factorized_T.solve(adt.f_);

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,f_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_V_setter(MatrixXr & V, const MatrixXr & T, const MatrixXr & R, const InputCarrier & carrier, AuxiliaryData<InputCarrier> & adt)
        {
                // Based on the type of problem at hand select at compile time the most useful Factorizer
                typedef typename std::conditional<std::is_base_of<Areal, InputCarrier>::value, Eigen::PartialPivLU<MatrixXr>, Eigen::LDLT<MatrixXr>>::type Factorizer;
                Factorizer factorized_T(T);	        // define a factorization of the defined type

                if(!carrier.is_areal() && !carrier.has_W())
                {
                        // Q == I
                        const SpMat * psi_tp = carrier.get_psi_tp();
                        V = factorized_T.solve(MatrixXr(*psi_tp));      // find the value of V = T^{-1}*Psi^t
                }
                else
                {
                        MatrixXr E_;                // Declare an empty auxiliary matrix
                        const UInt ret =  AuxiliaryOptimizer::universal_E_setter<InputCarrier>(E_, carrier);
                        V = factorized_T.solve(E_);          // find the value of V = T^{-1}*E
                }
                adt.K_ = factorized_T.solve(R);              // K = T^{-1}*R

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_E_setter(MatrixXr & E, const InputCarrier & carrier)
        {
                const VectorXr * Ap = carrier.get_Ap();
                if (carrier.has_W())
                {
                        // Psi is full && Q != I
                        const SpMat * psi_tp = carrier.get_psi_tp();
                        const MatrixXr * Qp = carrier.get_Qp();
                        AuxiliaryOptimizer::set_E_W_a(E, psi_tp, Qp, Ap);

                }
                else
                {
                        // Psi is full && Q == I
                        const SpMat * psi_tp = carrier.get_psi_tp();
                        AuxiliaryOptimizer::set_E_nW_a(E, psi_tp, Ap);
                }
                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,f_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_E_setter(MatrixXr & E, const InputCarrier & carrier)
        {
                const MatrixXr * Qp = carrier.get_Qp();        // Q != I
                if (carrier.loc_are_nodes())
                {
                        // Psi is permutation
                        const UInt nr = carrier.get_psip()->cols();
                        const UInt s = carrier.get_n_obs();
                        const std::vector<UInt> * kp  = carrier.get_obs_indicesp();
                        AuxiliaryOptimizer::set_E_ln_W_ptw(E, kp, Qp, nr, s);
                }
                else
                {
                        // Psi is full
                        const SpMat * psi_tp = carrier.get_psi_tp();
                        AuxiliaryOptimizer::set_E_lnn_W_ptw(E, psi_tp, Qp);
                }
                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_z_hat_setter(VectorXr & z_hat, InputCarrier & carrier, const MatrixXr & S, AuxiliaryData<InputCarrier> & adt, const Real lambda)
        {
                common_z_hat_part(z_hat, carrier, S);

                adt.left_multiply_by_psi(carrier, adt.r_, adt.g_);

                if (carrier.has_W())
                {
                        adt.r_ = lambda*carrier.lmbQ(adt.r_);

                }
                else
                {
                        adt.r_ = lambda * adt.r_;
                }

                z_hat += adt.r_;

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,f_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_z_hat_setter(VectorXr & z_hat, InputCarrier & carrier, const MatrixXr & S, AuxiliaryData<InputCarrier> & adt, const Real lambda)
        {
                common_z_hat_part(z_hat, carrier, S);

                return 0;
        }

template<typename InputCarrier>
void AuxiliaryOptimizer::common_z_hat_part(VectorXr & z_hat, InputCarrier & carrier, const MatrixXr & S)
{
        const VectorXr * zp = carrier.get_zp();
        if(carrier.has_W())
        {
                const MatrixXr * Hp = carrier.get_Hp();
                z_hat = ((*Hp)+carrier.lmbQ(S))*(*zp);
        }
        else
        {
                z_hat = S*(*zp);
        }
}

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_b_setter(MatrixXr & b, InputCarrier & carrier, const MatrixXr & US, const UInt nnodes)
        {
                if (carrier.has_W())
                {
                        b.topRows(nnodes) = (*carrier.get_psi_tp())*(carrier.get_Ap()->asDiagonal())*carrier.lmbQ(US);
                }
                else
                {
                        b.topRows(nnodes) = (*carrier.get_psi_tp())*(carrier.get_Ap()->asDiagonal())*US;
                }
                return 0;
        }


template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,f_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_b_setter(MatrixXr & b, InputCarrier & carrier, const MatrixXr & US, const UInt nnodes)
        {
                if (carrier.has_W())
                {
                        b.topRows(nnodes) = (*carrier.get_psi_tp())*carrier.lmbQ(US);
                }
                else
                {
                        b.topRows(nnodes) = (*carrier.get_psi_tp())*US;
                }
                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_first_updater(AuxiliaryData<InputCarrier> & adt, const InputCarrier & carrier, const MatrixXr & dS, const VectorXr & eps, const Real lambda)
        {
                const VectorXr * zp = carrier.get_zp();
                adt.t_ = dS*(*zp);
                MatrixXr temp = lambda*adt.K_;
                for (UInt i=0; i<temp.cols(); i++)
                {
                        temp.coeffRef(i,i) -= 1;
                }
                adt.h_ = temp*adt.g_;
                adt.left_multiply_by_psi(carrier, adt.p_, adt.h_);
                adt.p_ -= adt.t_;
                adt.a_ = eps.transpose()*adt.p_;  // note different from previous case!!

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,f_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_first_updater(AuxiliaryData<InputCarrier> & adt, const InputCarrier & carrier, const MatrixXr & dS, const VectorXr & eps, const Real lambda)
        {
                const VectorXr * zp = carrier.get_zp();
                adt.t_ = dS*(*zp);
                adt.a_ = -eps.transpose()*adt.t_;

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_second_updater(AuxiliaryData<InputCarrier> & adt, InputCarrier & carrier, const MatrixXr & ddS, const VectorXr & eps, const Real lambda)
        {
                const VectorXr * zp = carrier.get_zp();
                if (carrier.has_W())
                        adt.b_ = adt.p_.transpose()*VectorXr(carrier.lmbQ(adt.p_));
                else
                        adt.b_ = adt.p_.squaredNorm();

                VectorXr aux;
                adt.left_multiply_by_psi(carrier, aux, -2*adt.K_*adt.h_);

                adt.c_ = eps.transpose()*(-ddS*(*zp) + aux);

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,f_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_second_updater(AuxiliaryData<InputCarrier> & adt, InputCarrier & carrier, const MatrixXr & ddS, const VectorXr & eps, const Real lambda)
        {
                const VectorXr * zp = carrier.get_zp();
                if (carrier.has_W())
                        adt.b_ = adt.t_.transpose()*VectorXr(carrier.lmbQ(adt.t_));
                else
                        adt.b_ = adt.t_.squaredNorm();
                adt.c_ = -eps.transpose()*ddS*(*zp);

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<t_type,t_type>::value, Real>::type
        AuxiliaryOptimizer::universal_GCV(const Real s, const Real sigma_hat_sq, const Real dor)
        {
                return s*sigma_hat_sq/Real(dor);
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<t_type,t_type>::value, Real>::type
        AuxiliaryOptimizer::universal_GCV_d(const AuxiliaryData<InputCarrier> & adt, const Real s, const Real sigma_hat_sq, const Real dor, const Real trdS)
        {
                return 2*s*(sigma_hat_sq*trdS + adt.a_)/Real(dor*dor);
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<t_type,t_type>::value, Real>::type
        AuxiliaryOptimizer::universal_GCV_dd(const AuxiliaryData<InputCarrier> & adt, const Real s, const Real sigma_hat_sq, const Real dor, const Real trdS, const Real trddS)
        {
                return 2*s*(trdS*(3*sigma_hat_sq*trdS+4*adt.a_)/dor + sigma_hat_sq*trddS + adt.b_ + adt.c_)/Real(dor*dor);
        }

#endif

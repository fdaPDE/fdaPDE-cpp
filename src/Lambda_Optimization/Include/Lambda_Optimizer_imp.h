#ifndef __LAMBDA_OPTIMIZER_IMP_H__
#define __LAMBDA_OPTIMIZER_IMP_H__

// HEADERS
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include "../../Global_Utilities/Include/Timing.h"

// *** GCV-BASED ***

// -- Output managers --
//! Full output for optimization problem
/*! Set and return all output data, dof included, plus lambda final and
 number of iterations.
 \param optimal_pair output of optimization method consisting of gcv value and number of iterations
 \param time_count counts how much time was consumed by the optimization method
 \return full output_Data struct
*/
template<typename InputCarrier>
output_Data GCV_Family<InputCarrier, 1>::get_output(std::pair<Real,UInt> optimal_pair, const timespec & time_count, const std::vector<Real> & GCV_v, const std::vector<Real> & lambda_v, int termination_)
{
        this->output.content            = "full_optimization";
        this->output.lambda_sol         = optimal_pair.first;
        this->output.n_it               = optimal_pair.second;
        this->output.z_hat              = MatrixXr(this->z_hat);
        (this->output.rmse).push_back(this->rmse);
        this->output.sigma_hat_sq       = this->sigma_hat_sq;
        (this->output.dof).push_back(this->dof);
        this->output.time_partial       = time_count.tv_sec + 1e-9*time_count.tv_nsec;
        this->output.GCV_evals          = GCV_v;
        this->output.GCV_opt            = GCV_v[GCV_v.size()-1];
        this->output.lambda_vec         = lambda_v;
        this->output.lambda_pos         = GCV_v.size();
        this->output.termination        = termination_;
        this->output.betas              = this->the_carrier.get_model()->getBeta();
        return this->output;
}

//! Full output for dof computation
/*! Set and return all output data, dof included.
*/
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::set_output_partial(void)
{
        (this->output.rmse).push_back(this->rmse);
        (this->output.dof).push_back(this->dof);

}

//! Getter of the fill output
/*!
 \return the full output_Data struct
*/
template<typename InputCarrier>
output_Data GCV_Family<InputCarrier, 1>::get_output_full(void)
{
        return this->output;
}


//! Full output for dof computation
/*! Set and return all output data, dof included.
 \return output_Data struct containing predictions and dof
*/
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::set_output_partial_best(void)
{
        this->output.content            = "full_dof_grid";
        this->output.z_hat              = MatrixXr(this->z_hat);
        this->output.sigma_hat_sq       = this->sigma_hat_sq;

}

//! Partial output, for prediction evaluation without dof
/*! Set and return all output data butdof related
 \param f_hat the system soultion from which to evaluate the predictions
 \return output_Data struct containing predictions
*/
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::combine_output_prediction(const VectorXr & f_hat, output_Data & outp, UInt cols)
{
        this->compute_z_hat_from_f_hat(f_hat);
        this->compute_eps_hat();
        this->compute_SS_res();
        this->compute_rmse();

        if(outp.content != "prediction")
                outp.content = "prediction";

        outp.z_hat.col(cols) = this->z_hat;

        outp.rmse.push_back(this->rmse);
}

// -- Setters --
//! Utility to compute the predicted value in the locations given system solution f_hat
/*!
 \param f_hat top block of solution vector produced by apply methods
*/
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::compute_z_hat_from_f_hat(const VectorXr & f_hat)
{
        // z_hat  = H*z+Q*Psi*g_hat

        if (this->the_carrier.has_W())
        {
                this->z_hat = (*this->the_carrier.get_Hp())*(*this->the_carrier.get_zp()) + this->the_carrier.lmbQ((*this->the_carrier.get_psip())*f_hat);
        }
        else
        {
                this->z_hat = (*this->the_carrier.get_psip())*f_hat;
        }

        // Debugging purpose print
        /* Rprintf("z_hat \n");
           for(UInt i = 0; i < this->s-1; i++)
                  Rprintf("%f, ", this->z_hat[i]);
           Rprintf("%f", this->z_hat[s-1]);
           Rprintf("\n");
        */
}

//! Utility to compute the predicted residuals in the locations
/*!
 \pre z_hat must have been computed
*/
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::compute_eps_hat(void)
{
        // eps_hat = z-z_hat
        this->eps_hat = (*this->the_carrier.get_zp())-this->z_hat;

        // Debugging purpose print
        /* Rprintf("Eps_hat \n");
           for(UInt i = 0; i < this->s-1; i++)
                  Rprintf("%f, ", this->eps_hat[i]);
           Rprintf("%f", this->eps_hat[s-1]);
           Rprintf("\n");
        */
}

//! Utility to compute the sum of the squares of the residuals
/*!
 \pre compute_eps_hat() must have been called
 \sa compute_eps_hat()
*/
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::compute_SS_res(void)
{
        // SS_res = ||eps_hat||^2
        this->SS_res = this->eps_hat.squaredNorm();

        // Debugging purpose print
        // Rprintf("SS_res  = %f\n", this->SS_res);
}

//! Utility to compute the root mean square error
/*!
 \pre compute_SS_res() must have been called
 \sa compute_SS_res()
*/
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::compute_rmse(void)
{
        // rmse =std::sqrt(SS_res/#locations)
        this->rmse = std::sqrt(this->SS_res/Real(this->s));

        // Debugging purpose print
        // Rprintf("RMSE  = %f\n", this->rmse);
}

//! Utility to compute the estimated variance of the error
/*!
 \pre compute_SS_res() and compute_dor() must have been called
 \sa compute_SS_res(), compute_dor()
*/
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::compute_sigma_hat_sq(void)
{
        // sigma_hat^2 = SS_res/dor
        this->sigma_hat_sq = this->SS_res/Real(this->dor);

        // Debugging purpose print
        // Rprintf("sigma_hat_sq = %f\n", this->sigma_hat_sq);
}

//! Utility to compute the size of the model
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::compute_s(void)
{
        // s = #locations
        this->s = this->the_carrier.get_n_obs();

        // Debugging purpose print
        // Rprintf("s [# locations]  = %d\n", this->s);
}


// -- Updaters --
//! Utility to update the output-error parameters, fundamental for a correct computation of the gcv
/*!
 \param lambda the actual value of lambda to be used for the update
 \sa compute_eps_hat(), compute_SS_res(), compute_rmse(), update_dof(Real lambda),update_dor(Real lambda) and compute_sigma_hat_sq().
*/
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::update_errors(Real lambda)
{
        // this order must be kept
        this->compute_eps_hat();
        this->compute_SS_res();
        this->compute_rmse();
        this->update_dof(lambda);
        this->update_dor(lambda);
        this->compute_sigma_hat_sq();
}

//! Update all parameters needed to compute the gcv function, depending on lambda
/*!
 \param lambda the actual value of lambda to be used for the update
 \sa update_parameters(Real lambda)
*/
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::zero_updater(Real lambda)
{
        // Virtual update, depends on the gcv computational method [exact or stochastic]
        this->update_parameters(lambda);
}

//----------------------------------------------------------------------------//
// ** GCV_EXACT **

// -- Setters --
//! Method to set the value of member R_
/*!
 \remark R = R1^t * R0^{-1} * R1 therefore is NOT dependent on \lambda
*/
template<typename InputCarrier>
void GCV_Exact<InputCarrier, 1>::set_R_(void)
{
        /* Debugging purpose timer [part I]
         timer Time_partial_n;
	 Time_partial_n.start();
         Rprintf("WARNING: start taking time to build R inverse matrix \n");
        */

        const UInt ret = AuxiliaryOptimizer::universal_R_setter<InputCarrier>(this->R_, this->the_carrier, this->adt);

        /* Debugging purpose timer [part II]
         Rprintf("WARNING: partial time after the building R inverse matrix\n");
	 timespec T_n = Time_partial_n.stop();
        */
}

//! Method to set the value of member T_
/*!
 \remark T = D + \lambda * R where D is the top-left block of the matrix DMat
 \pre set_R_ must be called before set_T_, the matrix D_ [DataMatrix] must be constructed in the model, s must be defined
 \sa set_R_(void), getDataMatrix(SpMat & DMat)
 \param lambda the value for which to perform the optimization
*/
template<typename InputCarrier>
void GCV_Exact<InputCarrier, 1>::set_T_(Real lambda)
{
        this->T_ = lambda*this->R_;
        const UInt ret = AuxiliaryOptimizer::universal_T_setter<InputCarrier>(this->T_, this->the_carrier);
}

//! Method to set the value of member V_
/*!
 \remark V = T^{-1}*Psi^t*Q
 \pre set_T_ must be called before set_V_
 \sa set_T_(Real lambda)
*/
template<typename InputCarrier>
void GCV_Exact<InputCarrier, 1>::set_V_(void)
{
        const UInt ret = AuxiliaryOptimizer::universal_V_setter<InputCarrier>(this->V_, this->T_, this->R_, this->the_carrier, this->adt);
}

//! Method to set the value of member S_ and its trace trS_
/*!
 \remark S = Psi*V
 \pre set_V_ must be called before set_S_
 \sa set_V_(void)
*/
template<typename InputCarrier>
void GCV_Exact<InputCarrier, 1>::set_S_and_trS_(void)
{
        this->trS_ = 0.0;
        this->LeftMultiplybyPsiAndTrace(this->trS_, this->S_, this->V_);
}

//! Method to set the value of member dS_ and its trace trdS_, also computes utility matrix K_
/*!
 \remark dS_ = -Psi*T^{-1}*R*V
 \pre definition of base matrices R_, T_, V_
 \sa set_V_(void)
*/
template<typename InputCarrier>
void GCV_Exact<InputCarrier, 1>::set_dS_and_trdS_(void)
{
        // dS_ = -Psi*(Psi^t*Q*Psi+lambda*R1^t*R0^-1*R1)^{-1}*R1^t*R0^{-1}*R1*(Psi^t*Q*Psi+lambda*R1^t*R0^{-1}*R1)^{-1}*Psi^t*Q
        //     = -Psi*T^{-1}*R*V
        //     =  Psi*(-K*V)
        //    :=  -Psi*F
        this->adt.F_= this->adt.K_*this->V_;  // F = K*V
        this->trdS_ = 0.0;

        this->LeftMultiplybyPsiAndTrace(this->trdS_, this->dS_, -this->adt.F_);
}

//! Method to set the value of member ddS_ and its trace trddS_
/*!
 \remark ddS_ = -2*Psi*K^2*V
 \pre definition of V_ and K_
 \sa set_V_(void), set_dS_and_trdS_(void)
*/
template<typename InputCarrier>
void GCV_Exact<InputCarrier, 1>::set_ddS_and_trddS_(void)
{
        // ddS_ = 2*Psi*K^2*V = 2*K*F
        MatrixXr G_ = 2*this->adt.K_*this->adt.F_; // G = 2*K^2*V
        this->trddS_ = 0.0;

        this->LeftMultiplybyPsiAndTrace(this->trddS_, this->ddS_, G_);
}

// -- Utilities --
//! Utility to left multiply a matrix by Psi_ and compute the trace of the new matrix
/*!
 \param trace real where to store the value of the computed matrix
 \param ret sparse matrix where to store the computed matrix
 \param mat matrix to left multiply by Psi_
 \sa set_S_and_trS_(void), set_dS_and_trdS_(void), set_ddS_and_trddS_(void)
*/
template<typename InputCarrier>
void GCV_Exact<InputCarrier, 1>::LeftMultiplybyPsiAndTrace(Real & trace, MatrixXr & ret, const MatrixXr & mat)
{
        if (this->the_carrier.loc_are_nodes())
        {
                // Psi is permutation

                // THEORETICAL REMARK:
                // Since Psi is a rectangular permutation matrix, if function
                // k: loctions -> nodes s.t. Psi = Indicator(i,k[i]) then
                // Psi*F   == Indicator(i,j)*f_{k[i]j}

                // IMPLEMENTATION OF THE REMARK:
                // the number of non-null entries of E is at most s^2,
                // we reserve a vector containing such entries and
                // we set the final matrix from these triplets

                ret = MatrixXr::Zero(this->s,this->s); // Initialize return matrix as #locations*#locations

                const std::vector<UInt> * kp = this->the_carrier.get_obs_indicesp(); // Get z [observations]
                for (UInt i = 0; i < this->s; i++)
                        for (UInt j = 0; j < this->s; j++)
                        {
                                if (i == j) // diagonal block, also update trace
                                {
                                        Real v = mat.coeff((*kp)[i], j);
                                        trace += v;
                                        ret.coeffRef(i,i) += mat.coeff((*kp)[i], i);
                                }
                                else    // just update return matrix
                                {
                                        ret.coeffRef(i,j) += mat.coeff((*kp)[i], j);
                                }
                        }
        }
        else
        {
                // Psi is full, compute matrix and trace directly
                ret = (*this->the_carrier.get_psip())*mat;
                for (int i = 0; i < this->s; ++i)
                        trace += ret.coeff(i, i);
        }
}

// -- Computers and dof --
//! Utility to compute the predicted values in the locations
/*!
 \param lambda value of the optimization parameter
*/
template<typename InputCarrier>
void GCV_Exact<InputCarrier, 1>::compute_z_hat(Real lambda)
{
        UInt ret;
        if (this->the_carrier.get_bc_indicesp()->size()==0)
                ret = AuxiliaryOptimizer::universal_z_hat_setter<InputCarrier>(this->z_hat, this->the_carrier, this->S_, this->adt, lambda);
        else {

                const UInt nnodes    = this->the_carrier.get_n_nodes();
                const VectorXr f_hat = VectorXr(this->the_carrier.apply(lambda)).head(nnodes);

                // Compute the predicted values in the locations from the f_hat
                this->compute_z_hat_from_f_hat(f_hat);

               }
        // Debugging purpose print
        /* Rprintf("z_hat \n");
           for(UInt i = 0; i < this->s-1; i++)
                  Rprintf("%f, ", this->z_hat[i]);
           Rprintf("%f", this->z_hat[s-1]);
           Rprintf("\n");
        */
}

//! Utility to compute the degrees of freedom of the model
/*!
 \param lambda value of the optimization parameter
*/
template<typename InputCarrier>
void GCV_Exact<InputCarrier, 1>::update_dof(Real lambda)
{
        // dof = tr(S) + #covariates
	this->dof = this->trS_;

        if(this->the_carrier.has_W()) // add number of covariates, if present
                this->dof += (*this->the_carrier.get_Wp()).cols();

        // Debugging purpose
        // Rprintf("DOF: %f\n", this->dof);
}

//! Utility to compute the degrees of freedom of the residuals
/*!
 \param lambda value of the optimization parameter
 \pre update_dof() must have been called
 \sa update_dof()
*/
template<typename InputCarrier>
void GCV_Exact<InputCarrier, 1>::update_dor(Real lambda)
{
        // dor = #locations - dof
        this->dor = this->s-this->dof*this->the_carrier.get_opt_data()->get_tuning();

        if (this->dor < 0)   // Just in case of bad computation
        {
                Rprintf("WARNING: Some values of the trace of the matrix S('lambda') are inconstistent.\n");
                Rprintf("This might be due to ill-conditioning of the linear system.\n");
                Rprintf("Try increasing value of 'lambda'. Value of 'lambda' that produces an error is: %e \n", lambda);
        }

        // Debugging purpose
        // Rprintf("DOR: %f\n", this->dor);
}

// -- Global Updaters --
//! Utility to update the gcv-exact parameters, fundamental for a correct computation of the gcv
/*!
 \param lambda the actual value of lambda to be used for the update
 \sa set_T(Real lambda), set_V(), set_S_and_trS_() and compute_z_hat(Real lambda)
*/
template<typename InputCarrier>
void GCV_Exact<InputCarrier, 1>::update_matrices(Real lambda)
{
        // this order must be kept
        this->set_T_(lambda);
        this->set_V_();
        this->set_S_and_trS_();
        this->compute_z_hat(lambda);
}

// -- Public updaters --
//! Setting all the parameters which are recursively lambda dependent
/*!
 \remark The order in which functions are invoked is essential for the consistency of the procedure
 \sa update_matrices(Real lambda), update_errors(Real lambda)
*/
template<typename InputCarrier>
void GCV_Exact<InputCarrier, 1>::update_parameters(Real lambda)
{
        // this order must be kept
        this->update_matrices(lambda);
        this->update_errors(lambda);
}

//! Update all parameters needed to compute the gcv fist derivative, depending on lambda
/*!
 \param lambda the actual value of lambda to be used for the update
 \sa zero_updater(Real lambda), second_updater(Real lambda)
*/
template<typename InputCarrier>
void GCV_Exact<InputCarrier, 1>::first_updater(Real lambda)
{
        this->set_dS_and_trdS_();       // set first derivative of S and its trace
        UInt ret = AuxiliaryOptimizer::universal_first_updater<InputCarrier>(this->adt, this->the_carrier, this->dS_, this->eps_hat, lambda);
}

//! Update all parameters needed to compute the gcv second derivative, depending on lambda
/*!
 \param lambda the actual value of lambda to be used for the update
 \sa zero_updater(Real lambda), first_updater(Real lambda)
*/
template<typename InputCarrier>
void GCV_Exact<InputCarrier, 1>::second_updater(Real lambda)
{
        this->set_ddS_and_trddS_();     // set second derivative of S and its trace
        UInt ret = AuxiliaryOptimizer::universal_second_updater<InputCarrier>(this->adt, this->the_carrier, this->ddS_, this->eps_hat, lambda);
}

// -- GCV and derivatives --
// GCV function and derivatives
//! Main function computes the gcv in an exact fashion, depending on lambda
/*!
 Compact computation:
 GCV = s*(z-zhat)^t*(z-zhat)/(s-(q+trS))^2
     = SS_res*s/(dor^2)
     = sigma_hat_^2*s/dor
 \param lambda the actual value of lambda to be used for the computation
 \return the value of the gcv
*/
template<typename InputCarrier>
Real GCV_Exact<InputCarrier, 1>::compute_f(Real lambda)
{
        // call external updater to update [if needed] the parameters for gcv calculus
        this->gu.call_to(0, lambda, this);

        // compute the value of the gcv
        Real GCV_val =
                AuxiliaryOptimizer::universal_GCV<InputCarrier>(this->s, this->sigma_hat_sq, this->dor);

        // Debugging purpose print
        //Rprintf("LAMBDA = %f\n",lambda);
	//Rprintf("GCV = %f\n",GCV_val);

	return GCV_val;
}

//! Computes the gcv firt derivative in an exact fashion, depending on lambda
/*!
 Compact computation:
 dGCV(lambda)/dlambda = s * (d(1/dor(lambda)^2)/dlambda * SSres + dSSres(lambda)/dlambda * 1/dor^2)
                                     [1]                                        [2]
 where [1] = 2/dor^3 * d(tr(S))/dlambda = 2/dor^3 * tr(Phi*T^{-1}*R*V)
 and   [2] = 2*eps_hat^t*d(eps_hat)/dlambda = -2*eps^hat*dS
 summing: 2*s/(dor^2) * (sigma_hat_^2*tr(dS/dlambda) - eps_hat*dS/dlambda*z)
 \param lambda the actual value of lambda to be used for the computation
 \return the value of the gcv first derivative
 \sa compute_f(Real lambda)
*/
template<typename InputCarrier>
Real GCV_Exact<InputCarrier, 1>::compute_fp(Real lambda)
{
        // call external updater to update [if needed] the parameters for gcv first derivative
        this->gu.call_to(1, lambda, this);

        // compute the value of gcv first derivative
	Real GCV_der_val =
                AuxiliaryOptimizer::universal_GCV_d<InputCarrier>(this->adt, this->s, this->sigma_hat_sq, this->dor, this->trdS_);

        // Debugging purpose print
	//Rprintf("GCV_derivative = %f\n", GCV_der_val);

	return GCV_der_val;
}

//! Computes the gcv second derivative in an exact fashion, depending on lambda
/*!
 \param lambda the actual value of lambda to be used for the computation
 \return the value of the gcv second derivative
 \sa compute_f(Real lambda), compute_fp(Real lambda)
*/
template<typename InputCarrier>
Real GCV_Exact<InputCarrier, 1>::compute_fs(Real lambda)
{
        // call external updater to update [if needed] the parameters for gcv second derivative
        this->gu.call_to(2, lambda, this);

        // compute the value of gcv second derivative
	Real GCV_sec_der_val =
                AuxiliaryOptimizer::universal_GCV_dd<InputCarrier>(this->adt, this->s, this->sigma_hat_sq, this->dor, this->trdS_, this->trddS_);

        // Debugging purpose print
	//Rprintf("GCV_second_derivative = %f\n", GCV_sec_der_val);

	return GCV_sec_der_val;
}

//----------------------------------------------------------------------------//
// ** GCV_STOCHASTIC **

// -- Setters --
//! Setter of the stochastic binary matrix US_ needed for dof methods
template<typename InputCarrier>
void GCV_Stochastic<InputCarrier, 1>::set_US_(void)
{
        /* Debugging purpose timer [part I]
         timer Time_partial;
         Time_partial.start();
         Rprintf("WARNING: start taking time set_US\n");
        */

        // [[TODO Creation of the random generators ]]
        // set the seed
        UInt seed = this->the_carrier.get_opt_data()->get_seed();
        if(seed == 0)
        {
                seed = std::chrono::system_clock::now().time_since_epoch().count();
        }

        std::default_random_engine generator(seed);
	std::bernoulli_distribution distribution(0.5); // define random Be(p), p = 0.5

        // get number of relaizations to be performed and resize the matrix accordingly
        const UInt nr = this->the_carrier.get_opt_data()->get_nrealizations();
        this->US_ = MatrixXr::Zero(this->s, nr);

        // fill the binary matrix entries accordig to the stochastic algorithm
        for (UInt i=0; i<this->s; ++i)
                for (UInt j=0; j<nr; ++j)
                {
                        if (distribution(generator))
                        {
                                this->US_.coeffRef(i, j) = 1.0;
                        }
                        else
                        {
                                this->US_.coeffRef(i, j) = -1.0;
                        }
                }

        this->USTpsi = this->US_.transpose()*(*this->the_carrier.get_psip());

        // Define the first right hand side : | I  0 |^T * psi^T * Q * u
        UInt nnodes = this->the_carrier.get_n_nodes();
        this->b = MatrixXr::Zero(2*nnodes, this->US_.cols());
        UInt ret = AuxiliaryOptimizer::universal_b_setter(this->b, this->the_carrier, this->US_, nnodes);
        // Validate the completion of the task
        this->us = true;
        //Debugging purpose
        //Rprintf("Random matrix for stochastic method succesfully built.\n");

        /* Debugging purpose timer [part II]
         Rprintf("WARNING: time after the set_US method\n");
         timespec T = Time_partial.stop();
        */
}

// -- Computers and dof --
//! Utility to compute the degrees of freedom of the model, uses Stochastic Woodbury algorithm
/*!
 \param lambda value of the optimization parameter
*/
template<typename InputCarrier>
void GCV_Stochastic<InputCarrier, 1>::update_dof(Real lambda)
{
        MatrixXr m = this->the_carrier.get_opt_data()->get_DOF_matrix();
        if(m.cols()!=1 || m.rows()<=this->use_index)
        {
                /* Debugging purpose timer [part I]
                 timer Time_partial;
                 Time_partial.start();
                 Rprintf("WARNING: start taking time update_dof\n");
                */

        	UInt nnodes = this->the_carrier.get_n_nodes();
                UInt nr     = this->the_carrier.get_opt_data()->get_nrealizations();

                if(this->us == false) // check is US matrix has been defined
                {
                        this->set_US_();
                }


        	// Solve the system
            	MatrixXr x = this->the_carrier.apply_to_b(b, lambda);

        	VectorXr edf_vect(nr);
        	Real q = 0;

        	// Degrees of freedom = q + E[ u^T * psi * | I  0 |* x ]
        	if (this->the_carrier.has_W())
                {
        		q = this->the_carrier.get_Wp()->cols();
        	}

        	// For any realization we calculate the degrees of freedom
        	for (UInt i = 0; i < nr; ++i)
                {
        		edf_vect(i) = this->USTpsi.row(i).dot(x.col(i).head(nnodes)) + q;
        	}

        	// Estimates: sample mean, sample variance
        	this->dof = edf_vect.sum()/nr;

                // Deugging purpose print
                // Rprintf("DOF:%f\n", this->dof);

                /* Debugging purpose timer [part II]
                 Rprintf("WARNING: time after the update_dof method\n");
                 timespec T = Time_partial.stop();
                */
        }
        else
        {
                Rprintf("No DOF computation required\n");
                this->dof = m(this->use_index,0);
                //std::cout<< this->dof << std::endl;
        }
}

//! Utility to compute the degrees of freedom of the residuals
/*!
 \param lambda value of the optimization parameter
 \pre update_dof() must have been called
 \sa update_dof()
*/
template<typename InputCarrier>
void GCV_Stochastic<InputCarrier, 1>::update_dor(Real lambda)
{
        // dor = #locations - dof
        this->dor = this->s-this->dof*this->the_carrier.get_opt_data()->get_tuning();

        if (this->dor < 0)   // Just in case of bad computation
        {
                Rprintf("WARNING: Some values of the trace of the matrix S('lambda') are inconstistent.\n");
                Rprintf("This might be due to ill-conditioning of the linear system.\n");
                Rprintf("Try increasing value of 'lambda'. Value of 'lambda' that produces an error is: %e \n", lambda);
        }

        // Debugging purpose
        // Rprintf("DOR: %f\n", this->dor);
}

//! Utility to compute the predicted values in the locations
/*!
 \param lambda value of the optimization parameter
*/
template<typename InputCarrier>
void GCV_Stochastic<InputCarrier, 1>::compute_z_hat(Real lambda)
{
        /* Debugging purpose timer [part I]
         timer Time_partial;
         Time_partial.start();
         Rprintf("WARNING: start taking time compute_z_hat\n");
        */

        // Solve the system to find the predicted values of the spline coefficients
        const UInt nnodes    = this->the_carrier.get_n_nodes();
        const VectorXr f_hat = VectorXr(this->the_carrier.apply(lambda)).head(nnodes);

        // Compute the predicted values in the locations from the f_hat
        this->compute_z_hat_from_f_hat(f_hat);

        /* Debugging purpose timer [part II]
         Rprintf("WARNING: time after the compute_z_hat method\n");
         timespec T = Time_partial.stop();
        */
}

// -- Updaters --
//! Setting all the parameters which are recursively lambda dependent
/*!
 \remark The order in which functions are invoked is essential for the consistency of the procedure
 \sa compute_z_hat(Real lambda), update_errors(Real lambda)
*/
template<typename InputCarrier>
void GCV_Stochastic<InputCarrier, 1>::update_parameters(Real lambda)
{
        this->compute_z_hat(lambda);
        this->update_errors(lambda);
}

// -- GCV function --
// GCV function and derivatives
//! Main function computes the gcv in an exact fashion, depending on lambda
/*!
 Compact computation:
 GCV = s*(z-zhat)^t*(z-zhat)/(s-(q+trS))^2
     = SS_res*s/(dor^2)
     = sigma_hat_^2*s/dor
 \param lambda the actual value of lambda to be used for the computation
 \return the value of the gcv
*/
template<typename InputCarrier>
Real GCV_Stochastic<InputCarrier, 1>::compute_f(Real lambda)
{
        // call external updater to update [if needed] the parameters for gcv calculus
        this->gu.call_to(0, lambda, this);

        // compute the value of the gcv
        Real GCV_val =
                AuxiliaryOptimizer::universal_GCV<InputCarrier>(this->s, this->sigma_hat_sq, this->dor);

        // Debugging purpose print
        // Rprintf("LAMBDA = %f\n",lambda);
	// Rprintf("GCV = %f\n",GCV_val);

	return GCV_val;
}

//----------------------------------------------------------------------------//

#endif

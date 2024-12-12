// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __ESF_H__
#define __ESF_H__

// questi sono da controllare 
#include <fdaPDE/linear_algebra.h>
#include <fdaPDE/utils.h>

#include "../model_macros.h"
#include "../model_traits.h"
#include "../model_base.h"
#include "srpde.h"
#include "strpde.h"
#include "exact_edf.h"
#include "stochastic_edf.h" 
#include "inference_base.h"
#include "inference.h"

#include <algorithm>


using fdapde::core::is_empty;

namespace fdapde {
namespace models {

template <typename Model, typename Strategy> class ESF: public InferenceBase<Model>{

    private:
        struct ExactInverse{
            DMatrix<double> compute(Model m){
            return inverse(m.E());       
            }
        };

        struct NonExactInverse{
            SpMatrix<double> compute(Model m){
               return Base::invE_approx(m);  
            }
        };

        int n_flip = 1000; //default value
        int set_seed = 0;
        DMatrix<double> Lambda_ {};
        
        bool is_speckman_aux_computed = false;
        DVector<double> Speckman_aux_ranges;                       

        DMatrix<double> Qp_ {};
        int p_l_ = 0;   // number of locations for inference on f
        SpMatrix<double> Psi_p_ {};   // Psi only in the locations for inference
        DMatrix<double> Qp_dec_ {}; // Decomposition of Qp matrix
        DVector<int> mesh_nodes_ {};  // if inference is performed on a subset of mesh nodes

     
    public:
        using Base = InferenceBase<Model>;
        using Base::m_;
        using Base::beta_;  
        using Base::f0_;
        using Base::C_;
        using Base::locations_f_;
        using Base::alpha_f_;
        using Base::beta0_;
        using Solver = typename std::conditional<std::is_same<Strategy, exact>::value, ExactInverse, NonExactInverse>::type;
        Solver s_;

        using Base::V_;

        // constructors
        ESF() = default;                   // deafult constructor
        ESF(const Model& m): Base(m) {};     // constructor      

        DVector<double> p_value(CIType type) override{
            // extract matrix C (in the eigen-sign-flip case we cannot have linear combinations, but we can have at most one 1 for each column of C) 
            fdapde_assert(!is_empty(C_));        

            if(is_empty(beta0_)){ 
                Base::setBeta0(DVector<double>::Zero(m_.beta().size()));
            }

            int p = C_.rows(); 

            DVector<double> result;     // declare the vector that will store the p-values
            
            // compute Lambda
            if(is_empty(Lambda_)){
                V();
            }

            Eigen::SelfAdjointEigenSolver<DMatrix<double>> solver(Lambda_); // compute eigenvectors and eigenvalues of Lambda

            DMatrix<double> eigenvalues = solver.eigenvalues();
            DMatrix<double> eigenvectors = solver.eigenvectors();
            // Store beta_hat
            DVector<double> beta_hat = m_.beta();
            DVector<double> beta_hat_mod = beta_hat;
            
            if(type == simultaneous){  
                // SIMULTANEOUS   
                for(int i = 0; i < p; ++i){
                    for(int j = 0; j < C_.cols(); j++){
                        if(C_(i,j) > 0){
                            beta_hat_mod[j] = beta0_[j];
                        }
                    }
                }
                
                // partial residuals
                DMatrix<double> X = m_.X();
                DMatrix<double> res_H0 = m_.y() - X * beta_hat_mod; 
                // W^t * V * D
                DMatrix<double> Xt = (C_ * X.transpose()) * eigenvectors * eigenvalues.asDiagonal();   
                DVector<double> Tilder = eigenvectors.transpose() * res_H0;   

                // Initialize observed statistic and sign-flipped statistic
                DVector<double> stat = Xt * Tilder;
                DVector<double> stat_flip = stat;

                //Random sign-flips
                std::default_random_engine eng;
                std::uniform_int_distribution<int> distr(0, 1); 

                //if we have a seed 
                if(set_seed != 0) {
                    eng.seed(set_seed);
                } else {
                    std::random_device rd; 
                    eng.seed(rd()); // random seed 
                }

                int up = 0;
                int down = 0;

                DVector<double> Tilder_perm = Tilder;
                
                for(int i = 0; i < n_flip; i++){
                    for(int j = 0; j < Xt.cols(); ++j){
                        int flip = 2 * distr(eng) - 1;
                        Tilder_perm(j) = Tilder(j) * flip;
                    }
                    stat_flip = Xt * Tilder_perm; // Flipped statistic

                    if(is_Unilaterally_Greater(stat_flip, stat)){ 
                        up = up + 1;
                    }
                    else{ 
                    if(is_Unilaterally_Smaller(stat_flip, stat)){ 
                        down = down + 1;
                        }                    
                    }
                }
                
                double pval_up = static_cast<double>(up) / n_flip;
                double pval_down = static_cast<double>(down) / n_flip;

                result.resize(p); 
                result(0) = 2 * std::min(pval_up, pval_down); 
                for(int k = 1; k < p; k++){
                result(k) = 0.0;
                }
            }
            else{
                // ONE AT THE TIME   
                DMatrix<double> res_H0(Lambda_.cols(), p);
                DMatrix<double> X = m_.X();
                for(int i = 0; i < p; ++i){

                beta_hat_mod = beta_hat;

                for(int j = 0; j < C_.cols(); ++j){
                    if(C_(i,j)>0){
                        beta_hat_mod[j] = beta0_[j];
                    }
                }
                // compute the partial residuals
                res_H0.col(i) = m_.y() - X * beta_hat_mod;
                }
                // compute the vectors needed for the statistic
                DMatrix<double> Xt = (C_ * X.transpose()) * eigenvectors * eigenvalues.asDiagonal();   	// W^t * V * D
                DMatrix<double> Tilder = eigenvectors.transpose() * res_H0;   			        		// V^t * partial_res_H0

                // Observed statistic
                DMatrix<double> stat = Xt * Tilder;
                DMatrix<double> stat_flip = stat;

                // Random sign-flips
                std::default_random_engine eng;
                std::uniform_int_distribution<int> distr(0, 1); 

                //if we have a seed 
                if(set_seed != 0) {
                    eng.seed(set_seed);
                } else {
                    std::random_device rd; 
                    eng.seed(rd()); // random seed 
                }
                DVector<double> up = DVector<double>::Zero(p);
                DVector<double> down = DVector<double>::Zero(p);
                
                DMatrix<double> Tilder_perm = Tilder;
                
                for(int i = 0; i < n_flip; ++i){
                    for(int j = 0; j < Xt.cols(); ++j){
                        int flip = 2 * distr(eng) - 1;
                        Tilder_perm.row(j) = Tilder.row(j) * flip;
                    }
                    stat_flip = Xt * Tilder_perm; // Flipped statistic
                
                    for(int k = 0; k < p; ++k){
                        if(stat_flip(k, k) > stat(k, k)){
                            ++up(k);
                        }else{
                            ++down(k);
                        }
                    } 
                }
                
                DVector<double> pval_up = up.array() / static_cast<double>(n_flip);
                DVector<double> pval_down = down.array() / static_cast<double>(n_flip);

                result.resize(p);
                result = 2 * min(pval_up, pval_down); // Obtain the bilateral p_value starting from the unilateral
            } 
            return result;
        }   


        // CI for beta
        DMatrix<double> computeCI(CIType type) override{
            if(is_empty(Lambda_)){
                V();  // compute Lambda
            }
            DVector<double> beta_hat = m_.beta(); // store beta_hat
            DVector<double> beta_hat_mod = beta_hat;
            
            fdapde_assert(!is_empty(C_));     
            int p = C_.rows(); 

            DVector<int> beta_in_test; 
            beta_in_test.resize(p);
            for(int i=0; i<p; i++){
                for(int j=0; j<C_.cols(); j++){
                    if(C_(i,j)>0){
                        beta_in_test[i]=j;
                    }
                }
            }

            Eigen::SelfAdjointEigenSolver<DMatrix<double>> solver(Lambda_); // compute eigenvectors and eigenvalues of Lambda
            DMatrix<double> eigenvalues = solver.eigenvalues();
            DMatrix<double> eigenvectors = solver.eigenvectors();

            // declare the matrix that will store the intervals
            DMatrix<double> result;
            result.resize(p, 2);

            // compute the initial ranges from speckman's CI (initial guess for CI) 
            if(!is_speckman_aux_computed){
                Compute_speckman_aux();
            }

            // this vector will store the tolerance for each interval upper/lower limit
            DVector<double> ESF_bisection_tolerances = 0.1*Speckman_aux_ranges; // 0.1 of the speckman CI as maximum tolerance
            

            // define storage structures for bisection algorithms
            DVector<double> UU; // Upper limit for Upper bound
            UU.resize(p);
            DVector<double> UL; // Lower limit for Upper bound
            UL.resize(p);
            DVector<double> LU; // Upper limit for Lower bound
            LU.resize(p);
            DVector<double> LL; // Lower limit for Lower bound
            LL.resize(p);



            // ESF initialization with Speckman intervals as initial guess 
            for(int i = 0; i < p; ++i){
                double half_range = Speckman_aux_ranges(i);

                // compute the limits of the interval
                result(i,0) = beta_hat(beta_in_test[i]) - half_range;
                LU(i)=result(i,0)+0.5*half_range;
                LL(i)=result(i,0)-0.5*half_range;
                result(i,1) = beta_hat(beta_in_test[i]) + half_range;
                UU(i)=result(i,1) +0.5*half_range;
                UL(i)=result(i,1) -0.5*half_range; 	
            }


            // define booleans used to understand which CI need to be computed forward on
            std::vector<bool> converged_up(p,false);
            std::vector<bool> converged_low(p,false);
            bool all_betas_converged=false;

            // matrix that stores p_values of the bounds at actual step
            DMatrix<double> local_p_values;
            local_p_values.resize(4,p);
            
            // compute the vectors needed for the statistic
            DMatrix<double> TildeX = (C_ * m_.X().transpose()) * eigenvectors * eigenvalues.asDiagonal();   	// W^t * V * D
            DMatrix<double> Tilder_star = eigenvectors.transpose();   			        		// V^t
            // Select eigenvalues that will not be flipped basing on the estimated bias carried
            DVector<double> Tilder_hat = eigenvectors.transpose()* (m_.y() - (m_.X())* beta_hat); // This vector represents Tilder using only beta_hat, needed for bias estimation
            DVector<double> Partial_res_H0_CI;
            Partial_res_H0_CI.resize(Lambda_.cols());
            

            // fill the p_values matrix
            for (int i=0; i<p; i++){
                DMatrix<double> TildeX_loc = TildeX.row(beta_in_test[i]);
                beta_hat_mod = beta_hat;

                // compute the partial residuals and p value
                beta_hat_mod(beta_in_test[i])=UU(i); // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
                Partial_res_H0_CI = m_.y() - (m_.X()) * (beta_hat_mod); // (y-W*beta_hat(non in test)-W*UU[i](in test))
                local_p_values(0,i)=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);

                // compute the partial residuals and p value
                beta_hat_mod(beta_in_test[i])=UL(i); // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
                Partial_res_H0_CI =  m_.y() - (m_.X()) * (beta_hat_mod); // (y-W*beta_hat(non in test)-W*UL[i](in test))
                local_p_values(1,i)=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);

                // compute the partial residuals and p value
                beta_hat_mod(beta_in_test[i])=LU(i); // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
                Partial_res_H0_CI =  m_.y() - (m_.X()) * (beta_hat_mod); // (y-W*beta_hat(non in test)-W*LU[i](in test))
                local_p_values(2,i)=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);

                // compute the partial residuals and p value
                beta_hat_mod(beta_in_test[i])=LL(i); // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
                Partial_res_H0_CI = m_.y() - (m_.X()) * (beta_hat_mod); // (y-W*beta_hat(non in test)-W*LL[i](in test))
                local_p_values(3,i)=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);

            }

            // extract the CI significance (1-confidence)
            double alpha=0.05;

            if(type == one_at_the_time){
                alpha=0.5*alpha;
            }else{
                alpha=0.5/p*alpha;
            }
        
            int Max_Iter=50;
            int Count_Iter=0;
            while((!all_betas_converged) & (Count_Iter<Max_Iter)){
            
                // Compute all p_values (only those needed)
                for (int i=0; i<p; i++){

                DMatrix<double> TildeX_loc= TildeX.row(beta_in_test[i]);
                beta_hat_mod = beta_hat;

                if(!converged_up[i]){
                    if(local_p_values(0,i)>alpha){ // Upper-Upper bound excessively tight

                    UU(i)=UU(i)+0.5*(UU(i)-UL(i));
                
                    // compute the partial residuals
                    beta_hat_mod(beta_in_test[i])=UU(i); // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
                    Partial_res_H0_CI =  m_.y() - (m_.X()) * (beta_hat_mod); // (z-W*beta_hat(non in test)-W*UU[i](in test))
                    local_p_values(0,i)=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);
            
                    }else{
            
                        if(local_p_values(1,i)<alpha){ // Upper-Lower bound excessively tight
                            UL(i)=beta_hat(beta_in_test[i])+0.5*(UL(i)-beta_hat(beta_in_test[i]));
                    
                            // compute the partial residuals
                            beta_hat_mod(beta_in_test[i])=UL(i); // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
                            Partial_res_H0_CI =  m_.y() - (m_.X()) * (beta_hat_mod); // (z-W*beta_hat(non in test)-W*UL[i](in test))
                            local_p_values(1,i)=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);

                        }else{//both the Upper bounds are well defined

                            if(UU(i)-UL(i)<ESF_bisection_tolerances(i)){

                            converged_up[i]=true;

                            }else{

                                double proposal=0.5*(UU(i)+UL(i));
            
                                // compute the partial residuals
                                beta_hat_mod(beta_in_test[i])=proposal; // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
                                Partial_res_H0_CI =  m_.y() - (m_.X()) * (beta_hat_mod); // (z-W*beta_hat(non in test)-W*proposal)
                                double prop_p_val=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);

                                if(prop_p_val<=alpha){UU(i)=proposal; local_p_values(0,i)=prop_p_val;}else{UL(i)=proposal;local_p_values(1,i)=prop_p_val;}
                            }
                        }
                    }
                }  

                if(!converged_low[i]){
                    if(local_p_values(2,i)<alpha){ // Lower-Upper bound excessively tight

                        LU(i)=beta_hat(beta_in_test[i])-0.5*(beta_hat(beta_in_test[i])-LU(i));
    
                        // compute the partial residuals
                        beta_hat_mod(beta_in_test[i])=LU(i); // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
                        Partial_res_H0_CI =  m_.y() - (m_.X()) * (beta_hat_mod); // (z-W*beta_hat(non in test)-W*LU[i](in test))
                        local_p_values(2,i)=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);
                    
                    }else{
    
                        if(local_p_values(3,i)>alpha){ // Lower-Lower bound excessively tight
                            LL(i)=LL(i)-0.5*(LU(i)-LL(i));
                    
                            // compute the partial residuals
                            beta_hat_mod(beta_in_test[i])=LL(i); // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
                            Partial_res_H0_CI =  m_.y() - (m_.X()) * (beta_hat_mod);// (z-W*beta_hat(non in test)-W*LL[i](in test))
                            local_p_values(3,i)=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);

                        }else{//both the Upper bounds are well defined

                            if(LU(i)-LL(i)<ESF_bisection_tolerances(i)){

                                converged_low[i]=true;

                            }else{

                                double proposal=0.5*(LU(i)+LL(i));
    
                                // compute the partial residuals
                                beta_hat_mod(beta_in_test[i])=proposal; // beta_hat_mod(i) = beta_hat(i) if i not in test; beta_HP otherwise
                                Partial_res_H0_CI =  m_.y() - (m_.X()) * (beta_hat_mod); // (z-W*beta_hat(non in test)-W*proposal)
                                double prop_p_val=compute_CI_aux_beta_pvalue(Partial_res_H0_CI, TildeX_loc, Tilder_star);

                                if(prop_p_val<=alpha){LL(i)=proposal; local_p_values(3,i)=prop_p_val;}else{LU(i)=proposal;local_p_values(2,i)=prop_p_val;}
                            }
                        }
                    }
                }  // end if
            }  // end for
            all_betas_converged =true;
            for(int j=0; j<p; j++){

                if(!converged_up[j] || !converged_low[j]){
                    all_betas_converged=false;
                }
            }

            Count_Iter++;

            }
            result.resize(p,2);
            // for each row of C matrix
            for(int i=0; i<p; ++i){
                    
                if(Count_Iter < Max_Iter){ // No discrepancy between beta_hat(i) and ESF, bisection converged
                    // Limits of the interval
                    result(i,0) = 0.5*(LU(i)+LL(i));
                    result(i,1) = 0.5*(UU(i)+UL(i)); 
                    }else{ // Not converged in time
                    // Limits of the interval
                    result(i,0) = 10e20;
                    result(i,1) = 10e20; 
                }
            }
            
            return result;
            
        };

        void Compute_speckman_aux(void){
            
            fdapde_assert(!is_empty(C_));  

            double alpha_=0.05;
            //check if Lambda has been computed
            if(!is_empty(Lambda_)){
                V();
            }

            // extract W and W^T
            DMatrix<double> X = m_.X();
            DMatrix<double> X_t = m_.X().transpose();  
            // Decomposition of [W^t * Lambda^2 * W] 
            Eigen::PartialPivLU<DMatrix<double>> XLX_dec; 
            XLX_dec.compute((X_t)*(Lambda_*Lambda_)*(X));
    
            // get the residuals 
            DVector<double> eps_hat = (m_.y() - m_.fitted());
            // build squared residuals
            DVector<double> Res2=eps_hat.array()*eps_hat.array();
    
            // resize the variance-covariance matrix
            int q = C_.cols();
            DMatrix<double> V;
            V.resize(q,q);
            DMatrix<double> diag = Res2.asDiagonal();
    
            V = (XLX_dec).solve((X_t)*(Lambda_*Lambda_)*Res2.asDiagonal()*(Lambda_*Lambda_)*(X)*(XLX_dec).solve(DMatrix<double>::Identity(q,q))); // V = [(W*Lambda2*W)^-1 * Res2 * (W*Lambda2*W)^-1]

            // Extract the quantile needed for the computation of upper and lower bounds
            double quant = normal_standard_quantile(1 - alpha_/2);            

            // extract matrix C 
            int p = C_.rows(); 
            Speckman_aux_ranges.resize(p);
    
            // for each row of C matrix
            for(int i=0; i<p; ++i){
                DVector<double> col = C_.row(i);
                // compute the standard deviation of the linear combination and half range of the interval
                double sd_comb = std::sqrt(col.adjoint()*V*col);
                double half_range=sd_comb*quant;
                // save the half range
                Speckman_aux_ranges(i)=half_range;
            }

            this->is_speckman_aux_computed = true;

            return;
        }


        double compute_CI_aux_beta_pvalue(const DVector<double> & partial_res_H0_CI, const DMatrix<double> & TildeX,  const  DMatrix<double> & Tilder_star) const {
            // declare the vector that will store the p-values
            double result;
        
            // compute the vectors needed for the statistic 
            DVector<double> Tilder = Tilder_star * partial_res_H0_CI;

            // Initialize observed statistic and sign_flipped statistic
            DMatrix<double> stat_temp = TildeX*Tilder;
            double stat=stat_temp(0);
            double stat_flip=stat;

            // Random sign-flips
            std::default_random_engine eng;
            std::uniform_int_distribution<int> distr(0, 1); 

            //if we have a seed 
            if(set_seed != 0) {
                eng.seed(set_seed);
            } else {
                std::random_device rd; 
                eng.seed(rd()); // random seed 
            }
            double count_Up = 0;   // Counter for the number of flipped statistics that are larger the observed statistic
            double count_Down = 0; // Counter for the number of flipped statistics that are smaller the observed statistic
                
            DVector<double> Tilder_perm=Tilder;
        
            // get the number of flips
            int nflip=n_flip;

            for(int i=0;i<nflip;i++){
                for(int j=0;j<TildeX.cols();j++){
                    int flip;
                    flip=2*distr(eng)-1;
                    Tilder_perm(j)=Tilder(j)*flip;
                }
                    DMatrix<double> stat_flip_temp = TildeX*Tilder_perm; 
                    stat_flip= stat_flip_temp(0);// Flipped statistic
                    if(stat_flip > stat){ ++count_Up;}else{ 
                        if(stat_flip < stat){ ++count_Down;}  
                    }
            }
                
            double pval_Up = count_Up/n_flip;     
            double pval_Down = count_Down/n_flip; 

            result = std::min(pval_Up, pval_Down); // select the correct unilateral p_value 

            return result;
            
        };


        void Psi_p(){
            // case in which the locations are extracted from the observed ones
            if(is_empty(locations_f_) && is_empty(mesh_nodes_)){
                Psi_p_ = m_.Psi();
                if(p_l_ == 0)
                    p_l_ = Psi_p_.rows();
            }
            else if(is_empty(mesh_nodes_)){
                int m = locations_f_.size();
                SpMatrix<double> Psi = m_.Psi().transpose();
                Psi_p_.resize(m, Psi.rows());
                for(int j = 0; j < m; ++j) {
                    int row = locations_f_[j];
                    for(SpMatrix<double>::InnerIterator it(Psi, row); it; ++it) {
                        Psi_p_.insert(j, it.row()) = it.value();
                    }                
                }
                Psi_p_.makeCompressed();
            }
            else {
                // opposite dimensions of the Psi in the other cases
                int m = mesh_nodes_.size();
                SpMatrix<double> Psi = m_.PsiESF().transpose();
                SpMatrix<double> temp_psip;
                temp_psip.resize(m, Psi.rows());
                for(int j = 0; j < m; ++j) {
                    int row = mesh_nodes_[j];
                    for(SpMatrix<double>::InnerIterator it(Psi, row); it; ++it) {
                        temp_psip.insert(j, it.row()) = it.value();
                    }       
                }
                Psi_p_ = temp_psip;
            }
            Psi_p_.makeCompressed();
            if(p_l_ == 0)
                p_l_ = Psi_p_.rows();
        }

        DVector<double> yp(){
            if(is_empty(locations_f_)  && is_empty(mesh_nodes_))
                return m_.y();

            else if (!is_empty(locations_f_)){        
                int m = locations_f_.size();
                DVector<double> y = m_.y();
                DVector<double> yp;
                yp.resize(m);
                for(int j = 0; j < m; ++j) {
                    int row = locations_f_[j];
                    yp.row(j) = y.row(row);
                }
                return yp;
            }
            else{
                int m = mesh_nodes_.size();
                DVector<double> y = m_.y();
                DVector<double> yp;
                yp.resize(m);
                for(int j = 0; j < m; ++j) {
                    int row = mesh_nodes_[j];
                    yp.row(j) = y.row(row);
                }
                return yp;
            }
        }

        DMatrix<double> Xp(){
            // case in which the locations are extracted from the observed ones
            if(is_empty(locations_f_) && is_empty(mesh_nodes_)){
                return m_.X();
            }
            if(!m_.has_covariates())
                return m_.X();
            if(!is_empty(locations_f_)){
                int m = locations_f_.size();
                DMatrix<double> X = m_.X();
                DMatrix<double> Xp;
                Xp.resize(m, X.cols());
                for(int j = 0; j < m; ++j) {
                    int row = locations_f_[j];
                    Xp.row(j) = X.row(row);
                }
                return Xp;
            }   
            if(!is_empty(mesh_nodes_)){
                int m = mesh_nodes_.size();
                DMatrix<double> X = m_.X();
                DMatrix<double> Xp;
                Xp.resize(m, X.cols());
                for(int j = 0; j < m; ++j) {
                    int row = mesh_nodes_[j];
                    Xp.row(j) = X.row(row);
                }
                return Xp;
            }   
            return m_.X();  
        }

        DMatrix<double> Wp(){
            //  set it as a identity matrix
            if(is_empty(Psi_p_)){
                Psi_p();
            }   
            DMatrix<double> id = DMatrix<double>::Identity(p_l_, p_l_);
            return id;
        }

        // computes matrix Q = W(I - X*(X^\top*W*X)^{-1}*X^\top*W)
        void Qp() {
            if (!m_.has_covariates()){
                Qp_ = Wp() * DMatrix<double>::Identity(p_l_, p_l_);
            }   
            else{
                DMatrix<double> Xp_ = Xp();      // X^\top*W,   q x p_l_ 
                DMatrix<double> XpTXp = Xp_.transpose() * Xp_;
                Qp_ = DMatrix<double>::Identity(p_l_, p_l_) - Xp_ * XpTXp.ldlt().solve(Xp_.transpose());
            }

        }

        void Qp_dec(){
            // compute the Q only on the locations in which you want inference
            if(is_empty(Qp_))
                Qp();
            // Q * (y - epislon) = Q * r
            // Eigen decomposition
            Eigen::SelfAdjointEigenSolver<DMatrix<double>> Q_eigen(Qp_);  
            // matrix V is the matrix p_l_ x (p_l_-q) of the nonzero eigenvectors of Q
            DMatrix<double> Q_eigenvectors = Q_eigen.eigenvectors();
            if(m_.has_covariates())
                Qp_dec_ = Q_eigenvectors.rightCols(p_l_ - m_.X().cols());
            else
                Qp_dec_ = Q_eigenvectors;
        }



        double f_p_value(){
            if(is_empty(Psi_p_))
                Psi_p();
            if(is_empty(f0_)){
                Base::setf0(DVector<double>::Zero(Psi_p_.rows()));
            }
            if(is_empty(Qp_dec_))
                Qp_dec();

            // test statistic when Pi = Id
            DVector<double> VQr = Qp_dec_.transpose() * Qp_ * (yp() - f0_);

            DVector<double> Ti = Psi_p_.transpose() * Qp_dec_ * VQr;

            // save the rank of Ti
            double Ti_rank = Ti.array().square().sum();

            // random sign-flips
            // Bernoulli dist (-1, 1) with p = 0.5
            std::default_random_engine eng;
            std::uniform_int_distribution<int> distr(0, 1); 

            //if we have a set seed 
            if(set_seed != 0) {
                eng.seed(set_seed);
            } else {
                std::random_device rd; 
                eng.seed(rd()); // random seed 
            }
            int count = 0;
            DVector<double> tp_vqr = VQr; 
            DVector<double> Tp = Ti;

            for(int i = 0; i < n_flip; ++i){
                for(int j = 0; j < VQr.size(); ++j){
                    int flip = 2 * distr(eng) - 1;
                    tp_vqr(j) = VQr(j) * flip;
                }
                Tp = Psi_p_.transpose() * Qp_dec_ * tp_vqr;
                // flipped statistics
                double Tp_rank = Tp.array().square().sum();
                if(Tp_rank >= Ti_rank){
                    ++count;
                } 
            }
            double p_value = count/static_cast<double>(n_flip);
            return p_value;
        }

        
        double sign_flip_p_value(){
            if(is_empty(Psi_p_))
                Psi_p();
            if(is_empty(f0_)){
                Base::setf0(DVector<double>::Zero(Psi_p_.rows()));
            }
            if(is_empty(Qp_)){
                Qp();
            }
            // test statistic when Pi = Id
            DVector<double> VQr = Qp_ * (yp() - f0_);

            DVector<double> Ti = Psi_p_.transpose() * VQr;

            // save the rank of Ti
            double Ti_rank = Ti.array().square().sum();

            // random sign-flips
            // Bernoulli dist (-1, 1) with p = 0.5
            std::default_random_engine eng;
            std::uniform_int_distribution<int> distr(0, 1); 

            //if we have a set seed 
            if(set_seed != 0) {
                eng.seed(set_seed);
            } else {
                std::random_device rd; 
                eng.seed(rd()); // random seed 
            }
            int count = 0;
            DVector<double> tp_vqr = VQr; 
            DVector<double> Tp = Ti;
        
            for(int i = 0; i < n_flip; ++i){
                for(int j = 0; j < VQr.size(); ++j){
                    int flip = 2 * distr(eng) - 1;
                    tp_vqr(j) = VQr(j) * flip;
                }
                Tp = Psi_p_.transpose() * tp_vqr;
                // flipped statistics
                double Tp_rank = Tp.array().square().sum();
            
                if(Tp_rank >= Ti_rank){
                    ++count;
                } 
            }
            double p_value = static_cast<double>(count)/static_cast<double>(n_flip);

            return p_value;
        }


        double f_CI_p_val(const DVector<double> res, const int curr_index){
            if(is_empty(Qp_dec_))
                Qp_dec();

            DVector<double> Vres = Qp_dec_.transpose() * res;

            DVector<double> Ti = Psi_p_.transpose() * Qp_dec_ * Vres;

            // random sign-flips
            // Bernoulli dist (-1, 1) with p = 0.5
            std::default_random_engine eng;
            std::uniform_int_distribution<int> distr(0, 1); 

            //if we have a set seed 
            if(set_seed != 0) {
                eng.seed(set_seed);
            } else {
                std::random_device rd; 
                eng.seed(rd()); // random seed 
            }
            int up = 0;
            int down = 0;
            DVector<double> Tp = Ti;
            DVector<double> perm_res = Vres;

            for(int i = 0; i < n_flip; ++i){
                    for(int j = 0; j < Vres.size(); ++j){
                        int flip = 2 * distr(eng) - 1;
                        perm_res(j) = Vres(j) * flip;
                    }
                    Tp = Psi_p_.transpose() * Qp_dec_ * perm_res;

                    if(Tp(curr_index) > Ti(curr_index)){
                        ++up;
                    }
                    else if (Tp(curr_index) < Ti(curr_index)){
                        ++down;
                    }
            }

            return std::min(up, down);

        }



        DVector<double> Wald_range(){
            // create Wald object
            Wald<Model, Strategy> Wald_(m_);

            DMatrix<double> id = DMatrix<double>::Identity(m_.Psi().cols(), m_.Psi().cols());
            
            DMatrix<double> Wald_CI = Wald_.f_CI();
            int size = Wald_CI.rows();
            DVector<double> ranges {};
            ranges.resize(size);

            for(int i = 0; i < size; ++i){
                ranges(i) = 0.5 * (Wald_CI(i, 1) - Wald_CI(i, 0));
            }

            return ranges;
        }
        
        DMatrix<double> f_CI(){
            // only if the locations are a subset of the mesh
            // compute the Q only on the locations in which you want inference
            if(is_empty(Qp_))
            Qp();
            
            DVector<double> fp = Psi_p_ * m_.f();
            int nloc = fp.size();

            DMatrix<double> CI;
            CI.resize(nloc, 2);

            DVector<double> Wald_ranges_ = Wald_range();

            // bisection tolerance
            DVector<double> bis_tol = 0.2 * Wald_ranges_;

            DVector<double> UU; // Upper limit for Upper bound
            UU.resize(nloc);
            DVector<double> UL; // Lower limit for Upper bound
            UL.resize(nloc);
            DVector<double> LU; // Upper limit for Lower bound
            LU.resize(nloc);
            DVector<double> LL; // Lower limit for Lower bound
            LL.resize(nloc);

            // initialize the intervals with a Wald estimate
            for(int i = 0; i < nloc; ++i){
                double range = Wald_ranges_(i);
                CI(i, 0) = fp(i) - range; 
                CI(i, 1) = fp(i) + range;
                UU(i) = CI(i, 1) + 0.5 * range;
                UL(i) = CI(i, 1) - 0.5 * range;
                LU(i) = CI(i, 0) + 0.5 * range;
                LL(i) = CI(i, 0) - 0.5 * range;
            }

            // vector of booleans to check which inetrvals have converged
            std::vector<bool> converged_up(nloc, false);
            std::vector<bool> converged_low(nloc, false);
            bool converged = false;
            
            // we need p-values ofr UU, UL, LU, LL
            DMatrix<double> p_values(4, nloc);
            // need this to compute the p-values
            DVector<double> res(nloc);

            // for getting the current index
            if(is_empty(locations_f_)){
                locations_f_ = DVector<int>::LinSpaced(nloc, 0, nloc-1);
            }

            for(int i = 0; i < nloc; ++i){
                DVector<double> fp_UU = fp;
                DVector<double> fp_UL = fp;
                DVector<double> fp_LU = fp;
                DVector<double> fp_LL = fp;

                res = Qp_ * (yp() - fp_UU); 
                p_values(0,i) = f_CI_p_val(res, locations_f_[i]);
                res = Qp_ * (yp() - fp_UL); 
                p_values(1,i) = f_CI_p_val(res, locations_f_[i]);
                res = Qp_ * (yp() - fp_LU); 
                p_values(2,i) = f_CI_p_val(res, locations_f_[i]);
                res = Qp_ * (yp() - fp_LL); 
                p_values(3,i) = f_CI_p_val(res, locations_f_[i]);

            }
            // one at the time confidence intervals
            int max_it = 5;
            int count = 0;
            double ppval = 0;

            while(!converged && count < max_it){

                // compute p-values needed
                for(int i = 0; i < nloc; ++i){
                    DVector<double> fp_UU = fp;
                    DVector<double> fp_UL = fp;
                    DVector<double> fp_LU = fp;
                    DVector<double> fp_LL = fp;

                    if(!converged_up[i]){
                    // check upper upper
                        if(p_values(0, i) > 0.5 * alpha_f_){
                            UU(i) = UU(i) + 0.5 * (UU(i) - UL(i));
                            res = Qp_ * (yp() - fp_UU); 
                            p_values(0,i) = f_CI_p_val(res, locations_f_[i]);
                        }
                        else if(p_values(1, i) < 0.5 * alpha_f_){
                            UL(i) = fp(i) + 0.5 * (UL(i) - fp(i));
                            res = Qp_ * (yp() - fp_UL); 
                            p_values(1,i) = f_CI_p_val(res, locations_f_[i]);
                        }
                        else{
                            if((UU(i) - UL(i)) < bis_tol(i)){
                                converged_up[i] = true;
                            }
                            else{
                                double proposal = 0.5 * (UU(i) + UL(i));
                                DVector<double> fpp = fp;
                                res = Qp_ * (yp() - fpp);
                                ppval = f_CI_p_val(res, locations_f_[i]);
                                if(ppval <= 0.5 * alpha_f_){
                                    UU(i) = proposal;
                                    p_values(0, i) = ppval;
                                }
                                else{
                                    UL(i) = proposal;
                                    p_values(1, i) = ppval;
                                }
                            }
                        }
                    } // end if

                    if(!converged_low[i]){
                    // check lower upper
                        if(p_values(2, i) < 0.5 * alpha_f_){
                            LU(i) = fp(i) - 0.5 * (fp(i) - LU(i));
                            res = Qp_ * (yp() - fp_LU); 
                            p_values(2,i) = f_CI_p_val(res, locations_f_[i]);
                        }
                        else if(p_values(3, i) > 0.5 * alpha_f_){
                            LL(i) = LL(i) + 0.5 * (LU(i) - LL(i));
                            res = Qp_ * (yp() - fp_LL); 
                            p_values(3,i) = f_CI_p_val(res, locations_f_[i]);
                        }
                        else{
                            if((LU(i) - LL(i)) < bis_tol(i)){
                                converged_low[i] = true;
                            }
                            else{
                                double proposal = 0.5 * (LU(i) + LL(i));
                                DVector<double> fpp = fp;
                                res = Qp_ * (yp() - fpp);
                                ppval = f_CI_p_val(res, locations_f_[i]);
                                if(ppval <= 0.5 * alpha_f_){
                                    LL(i) = proposal;
                                    p_values(3, i) = ppval;
                                }
                                else{
                                    LU(i) = proposal;
                                    p_values(2, i) = ppval;
                                }
                            }
                        }

                    } // end if
                } // for end

                converged = true;

                for(int j = 0; j < nloc; ++j){
                    if(!converged_up[j] || !converged_low[j]){
                        converged = false;
                    }
                }

                count++;
                
            }
            if(!converged){
                std::cout << "Not all ESF CI converged" << std::endl;
            }

            for(int k = 0; k < nloc; ++k){
                if(!converged_up[k] || !converged_low[k]){
                    // set a unfeasible number if the interval has not converged
                    CI(k, 0) = 10e20;
                    CI(k, 1) = 10e20;
                }    
                else{
                    CI(k, 0) = 0.5 * (LU(k) + LL(k));
                    CI(k, 1) = 0.5 * (UU(k) + UL(k));
                }        
            }

            return CI;       
        }

        void V() override{
            Lambda_ = DMatrix<double>::Identity(m_.n_obs(), m_.n_obs()) - m_.Psi() * s_.compute(m_) * m_.PsiTD();
        }   

        inline bool is_Unilaterally_Greater (DVector<double> v, DVector<double> u){
            int q = v.size();
            for (int i = 0; i < q; ++i){
                if(v(i) <= u(i)){
                    return false;
                }
            }
            return true;
        }

        inline bool is_Unilaterally_Smaller (DVector<double> v, DVector<double> u){
            int q = v.size();
            for (int i = 0; i < q; ++i){
                if(v(i) >= u(i)){
                    return false;
                }
            }
            return true;
        }

        inline DVector<double> min(const DVector<double> & v, const DVector<double> & u){
            DVector<double> result;
            result.resize(v.size());
            for(int i = 0; i < v.size(); ++i){
                result(i) = std::min(v(i), u(i));
            }
            return result;
        } 
        
        void setNflip(int m){
            n_flip = m;
        }

        void setseed(int k){
            set_seed=k;
        }
            
        void setMesh_loc(DVector<int> m_nodes){
            mesh_nodes_ = m_nodes;
        }

};

} // namespace models
} // namespace fdapde

#endif // __ESF_H__

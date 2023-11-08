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

#ifndef __MSQRPDE_H__
#define __MSQRPDE_H__

#include <fdaPDE/pde.h>
#include <fdaPDE/utils.h>
using fdapde::core::PDEBase;

#include "../model_base.h"
#include "../model_macros.h"
#include "../model_traits.h"
#include "../sampling_design.h"
#include "distributions.h"
#include "fpirls.h"
#include "regression_base.h"

namespace fdapde{
namespace models{
	    
  
template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver>
    class MSQRPDE : public RegressionBase<MSQRPDE<PDE, RegularizationType, SamplingDesign, Solver>> {  
        // compile time checks
        static_assert(std::is_base_of<PDEBase, PDE>::value);
    private:
        typedef RegressionBase<MSQRPDE<PDE, RegularizationType, SamplingDesign, Solver>> Base;

        unsigned int h_;                      // number of quantile orders 
        const std::vector<double> alphas_;    // quantile order 

        // algorithm's parameters 
        double gamma0_ = 1.0;                  // crossing penalty 
        double eps_ = 1e-6;                    // crossing tolerance 
        double C_ = 1.5;                       // crossing penalty factor
        double tolerance_ = 1e-5;              // convergence tolerance 
        double tol_weights_ = 1e-6;            // weights tolerance
        std::size_t max_iter_ = 50;            // inner max number of iterations 
        std::size_t max_iter_global_ = 100;    // outer max number of iterations 
        std::size_t k_ = 0;                    // inner iteration index
        std::size_t iter_ = 0;                 // outer iteration index

        // linear system  
        SparseBlockMatrix<double,2,2> A_{};         // system matrix of non-parametric problem (2hN x 2hN matrix)
        fdapde::SparseLU<SpMatrix<double>> invA_;   // factorization of matrix A
        DVector<double> b_{};                       // system rhs 

        // room for solution 
        DVector<double> f_curr_{};     // current estimate of the spatial field f (1 x h*N vector)
        DVector<double> fn_curr_{};    // current estimate of the spatial field f_n (1 x h*n vector)
        DVector<double> g_curr_{};     // current PDE misfit
        DVector<double> beta_curr_{};  // current estimate of the coefficient vector (1 x h*q vector)

        DVector<double> f_init_{}; 
        DVector<double> fn_init_{}; 
        DVector<double> g_init_{}; 
        DVector<double> beta_init_{}; 

        // room for algorithm quantities 
        SpMatrix<double> Ih_;                 // identity h x h 
        SpMatrix<double> In_;                 // identity n x n
        SpMatrix<double> Iq_;                 // identity q x q 
        SpMatrix<double> Ihn_;                // identity h*n x h*n 
        SpMatrix<double> D_{};                
        SpMatrix<double> D_script_{}; 
        DVector<double> l_hn_{}; 
        SpMatrix<double> Psi_multiple_{}; 
        SpMatrix<double> R0_multiple_{};
        SpMatrix<double> R1_multiple_{}; 
        DVector<double> z_{};
        DVector<double> w_{};
        DiagMatrix<double> W_bar_{};
        SpMatrix<double> W_multiple_{}; 
        DMatrix<double> X_multiple_{};
        DMatrix<double> XtWX_multiple_{};
        Eigen::PartialPivLU<DMatrix<double>> invXtWX_multiple_{}; 
        DMatrix<double> Q_multiple_{}; 
        DMatrix<double> H_multiple_{}; 
        DMatrix<double> U_multiple_{};
        DMatrix<double> V_multiple_{};

    public:
        IMPORT_REGRESSION_SYMBOLS;
        using Base::invXtWX_;
        using Base::lambda_D;   // smoothing parameter in space
        using Base::n_basis;    // number of spatial basis
        using Base::P;          // discretized penalty matrix: P = \lambda_D*(R1^T*R0^{-1}*R1)
        using Base::W_;         // weight matrix
        using Base::XtWX_; 

        DVector<double> lambdas_D;       // smoothing parameters in space  

        // constructor
        MSQRPDE() = default;
        MSQRPDE(const PDE& pde, std::vector<double>& alphas = {0.1, 0.5, 0.9}) : Base(pde), alphas_(alphas) {
            // // Check if the provided quantile orders are an increasing sequence 
            // auto i = std::adjacent_find(alphas.begin(), alphas.end(), std::greater_equal<int>());
            // assert(i == alphas.end(), "Quantile orders must be an increasing sequence"); 
            h_ = alphas_.size();
        }; 

        // getters
        const DVector<double>& f() const { return f_curr_; };            // estimate of spatial field
        const DVector<double>& g() const { return g_curr_; };            // PDE misfit
        const DVector<double>& beta() const { return beta_curr_; };      // estimate of regression coefficients
        // const SparseBlockMatrix<double,2,2>& A_mult() const { return A_; };              // debug
        const SpMatrix<double>& Psi_mult() const { return Psi_multiple_; };              // debug
        // const SpMatrix<double>& W_mult() const { return W_multiple_debug; };             // debug
        // const SpMatrix<double>& D_script() const { return D_script_; };                  // debug  
        // const DiagMatrix<double>& Delta_mult() const { return Delta_debug; };            // debug
        // const DiagMatrix<double>& Wbar_mult() const { return W_bar_debug; };             // debug 
        // const DMatrix<double>& Q_mult() const { return Q_multiple_; };                   // debug
        // const DMatrix<double>& H_mult_debug() const { return H_multiple_debug; };        // debug 
        // const DMatrix<double>& X_mult() const { return X_multiple_; };                   // debug 
        // const SpMatrix<double>&  Dscriptj() const { return Dscriptj_debug; };            // debug 
        // const DMatrix<double>& XtWX_multiple() const { return XtWX_multiple_debug; };    // debug

        // ModelBase implementation
        void init_model();
        void update_to_weights() { return; };   // update model object in case of changes in the weights matrix
        virtual void solve(); // finds a solution to the smoothing problem

        // Utilities 
        const bool crossing_constraints() const;

        void setLambdas_D(DMatrix<double> l){
            lambdas_D.resize(l.rows()); 
            for(auto i = 0; i < l.rows(); ++i)
                lambdas_D(i) = l(i,0);  
        }
            
        void assemble_matrices(){

            // room for solution 
            f_curr_.resize(h_*n_basis());
            fn_curr_.resize(h_*n_obs());
            g_curr_.resize(h_*n_basis());
            f_init_.resize(h_*n_basis());
            fn_init_.resize(h_*n_obs());
            g_init_.resize(h_*n_basis());

            if(has_covariates()){
                beta_curr_.resize(h_*q());
                beta_init_.resize(h_*q());
            }

            // set identity matrices 
            Ih_.resize(h_, h_); 
            Ih_.setIdentity();
            In_.resize(n_obs(), n_obs()); 
            In_.setIdentity();
            Iq_.resize(q(), q()); 
            Iq_.setIdentity();
            Ihn_.resize(h_*n_obs(), h_*n_obs()); 
            Ihn_.setIdentity();

            // assemble FEM, mass and stiffness matrices
            Psi_multiple_ = Kronecker(Ih_, Psi()); 
            R0_multiple_ = Kronecker(SpMatrix<double>(DiagMatrix<double>(lambdas_D)), R0());
            R1_multiple_ = Kronecker(SpMatrix<double>(DiagMatrix<double>(lambdas_D)), R1());

            // assemble
            if(has_covariates()){ 
                X_multiple_ = Kronecker(DMatrix<double>(Ih_), X()); 
                XtWX_multiple_.resize(h_*q(), h_*q()); 
                U_multiple_.resize(2*h_*n_basis(), h_*q());
                V_multiple_.resize(h_*q(), 2*h_*n_basis());
            }

            // assemble l_hn_
            l_hn_.resize(h_*n_obs()); 
            l_hn_ = DVector<double>::Zero(h_*n_obs());
            l_hn_.block(0,0, n_obs(), 1) = -DVector<double>::Ones(n_obs());
            l_hn_.block((h_-1)*n_obs(),0, n_obs(), 1) = DVector<double>::Ones(n_obs());


            SpMatrix<double> E_{};
            E_.resize(h_-1, h_); 
            std::vector<fdapde::Triplet<double>> tripletList;
            tripletList.reserve(2*(h_-1));

            for(std::size_t i = 0; i < h_-1; ++i){
                tripletList.emplace_back(i, i+1, 1.0);
                tripletList.emplace_back(i, i, -1.0);
            }

            E_.setFromTriplets(tripletList.begin(), tripletList.end());
            E_.makeCompressed();

            if(has_covariates()){
                D_ = Kronecker(E_, Iq_); 
            }       
            D_script_ = Kronecker(E_, In_);  

        }

        DVector<double> rho_alpha(const double&, const DVector<double>&) const; 
        DVector<double> fitted(unsigned int j) const; 
        DVector<double> fitted() const; 
        void abs_res_adj(DVector<double>& w); 
        double model_loss() const; 
        double crossing_penalty() const;
        double crossing_penalty_f() const; 
        double crossing_penalty_param() const; 

        const DMatrix<double>& H_multiple(); 
        const DMatrix<double>& Q_multiple(); 

        virtual ~MSQRPDE() = default;
    };

    template <typename PDE_, typename RegularizationType_, typename SamplingDesign_, typename Solver_>
        struct model_traits<MSQRPDE<PDE_, RegularizationType_, SamplingDesign_, Solver_>> {
        typedef PDE_ PDE;
        typedef SpaceOnly regularization;
        typedef SamplingDesign_ sampling;
        typedef MonolithicSolver solver;
        enum { N = PDE::N, M = PDE::M, n_lambda = 1 };
    };


    // msqrpde trait
    template <typename Model>
    struct is_msqrpde { static constexpr bool value = is_instance_of<Model, MSQRPDE>::value; };


    // perform proper initialization and update of model. Computes quantites which can be reused
    // across many calls to solve() and are **not affected by a change in the data**.
    // It is implicitly called by ModelBase::init() as part of the initialization process.
    // NB: a change in the smoothing parameter must trigger a re-initialization of the model
    template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver>
    void MSQRPDE<PDE, RegularizationType, SamplingDesign, Solver>::init_model() {

        // Assemble matrices
        assemble_matrices();  

        // Definition of h SQRPDE models for initialization 
        for(std::size_t j = 0; j < h_; ++j){
            SQRPDE<PDE, fdapde::models::SpaceOnly, SamplingDesign, fdapde::models::MonolithicSolver> 
                model_j(pde(), alphas_[j]);

            // solver initialization
            model_j.data() = data();
            model_j.set_lambda_D(lambdas_D[j]);     
            model_j.set_spatial_locations(this->locs());
            model_j.init_pde();
            model_j.init_regularization();
            model_j.init_sampling();    
            model_j.init_nan();
            model_j.init_model();
            model_j.solve();

            f_curr_.block(j*n_basis(), 0, n_basis(), 1) = model_j.f();
            fn_curr_.block(j*n_obs(), 0, n_obs(), 1) = Psi()*model_j.f();
            g_curr_.block(j*n_basis(), 0, n_basis(), 1) = model_j.g();
            if(has_covariates()){
                beta_curr_.block(j*q(), 0, q(), 1) = model_j.beta();
            }

        }

        // init = curr 
        f_init_ = f_curr_; 
        fn_init_ = fn_curr_; 
        g_init_ = g_curr_; 
        if(has_covariates()){
            beta_init_ = beta_curr_;
            std::cout << "Beta_init_ : " << beta_init_ << std::endl; 

        }

        std::cout << "Range f_init_ : " << f_init_.minCoeff() << " , " << f_init_.maxCoeff() << std::endl ; 
        std::cout << "Range g_init_ : " << g_init_.minCoeff() << " , " << g_init_.maxCoeff() << std::endl ; 
        std::cout << "Model loss of the initializazion: " << model_loss() << " + penalty = " << g_curr_.dot(R0_multiple_*g_curr_) << std::endl;  
        std::cout << "----- crossing global at init = " << crossing_penalty() << std::endl;
        return;
    }

    // finds a solution to the SR-PDE smoothing problem
    template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver>
        void MSQRPDE<PDE, RegularizationType, SamplingDesign, Solver>::solve() {

        w_.resize(h_*n_obs()); 
        W_bar_.resize(h_*n_obs());    

        DiagMatrix<double> Delta_; 
        Delta_.resize((h_-1)*n_obs()); 

        z_.resize(h_*n_obs());

        DVector<double> t{};    
        t.resize(h_*n_obs()); 

        double crossing_penalty_init = crossing_penalty(); 

        unsigned int count_gamma = 0;
        while(crossing_constraints() && iter_ < max_iter_global_){ 

            std::cout << "----------------Gamma = " << gamma0_ << std::endl; 
            // algorithm stops when an enought small difference between two consecutive values of the J is recordered
            double J_old = tolerance_+1; double J_new = 0;
            k_ = 0;

            while(k_ < max_iter_ && std::abs(J_new - J_old) > tolerance_){    // && crossing_constraints()

                std::cout << "--------------------------  k_ = " << k_ << std::endl; 

                // assemble W, Delta, z 
                DVector<double> delta_((h_-1)*n_obs()); 

                for(int j = 0; j < h_; ++j){

                    DVector<double> abs_res_j;
                    DVector<double> delta_j; 
                    DVector<double> z_j;
                    
                    abs_res_j = (y() - fitted(j)).cwiseAbs(); 

                    if(j < h_-1) {
                        delta_j = (2*(eps_*DVector<double>::Ones(n_obs()) - D_script_.block(j*n_obs(), 0, n_obs(), h_*n_obs())*fitted())).cwiseAbs().cwiseInverse(); 
                    }
                    
                    z_j = y() - (1 - 2*alphas_[j])*abs_res_j; 

                    abs_res_adj(abs_res_j);

                    w_.block(j*n_obs(), 0, n_obs(), 1) = 2*n_obs()*abs_res_j;

                    if(j < h_-1) 
                        delta_.block(j*n_obs(), 0, n_obs(), 1) = delta_j; 

                    z_.block(j*n_obs(), 0, n_obs(), 1) = z_j;         
                }

                Delta_.diagonal() = delta_;
                W_bar_.diagonal() = w_.cwiseInverse(); 
                W_multiple_ = SpMatrix<double>(W_bar_) + gamma0_*D_script_.transpose()*Delta_*D_script_; 

                // assemble t 
                t = D_script_.transpose()*Delta_*eps_*DVector<double>::Ones((h_-1)*n_obs()) + 0.5*l_hn_; 

                // assemble system matrix for nonparameteric part
                A_ = SparseBlockMatrix<double,2,2>
                (-Psi_multiple_.transpose()*W_multiple_*Psi_multiple_,    R1_multiple_.transpose(),
                R1_multiple_,                                           R0_multiple_            );
                
                // cache non-parametric matrix factorization for reuse
                invA_.compute(A_);

                // prepare rhs of linear system
                b_.resize(A_.rows());
                b_.block(h_*n_basis(),0, h_*n_basis(),1) = DVector<double>::Zero(h_*n_basis());  // b_g = 0 

                DVector<double> sol; // room for problem' solution      
                
                if(!has_covariates()){ // nonparametric case     

                    // update rhs of SR-PDE linear system
                    b_.block(0,0, h_*n_basis(),1) = -Psi_multiple_.transpose()*(W_bar_*z_ + gamma0_*t);

                    // solve linear system A_*x = b_
                    sol = invA_.solve(b_);

                    f_curr_ = sol.head(h_*n_basis());
                    fn_curr_ = Psi_multiple_*f_curr_; 

                    } else{ // parametric case
    
                        XtWX_multiple_ = X_multiple_.transpose()*W_multiple_*X_multiple_;
                        invXtWX_multiple_ = XtWX_multiple_.partialPivLu(); 

                        // update rhs of SR-PDE linear system
                        b_.block(0,0, h_*n_basis(),1) = -Psi_multiple_.transpose()*(Ihn_ - H_multiple().transpose())*(W_bar_*z_ + gamma0_*t);  

                        // definition of matrices U and V  for application of woodbury formula
                        U_multiple_ = DMatrix<double>::Zero(2*h_*n_basis(), h_*q());
                        U_multiple_.block(0,0, h_*n_basis(), h_*q()) = Psi_multiple_.transpose()*W_multiple_*X_multiple_;
                        V_multiple_ = DMatrix<double>::Zero(h_*q(), 2*h_*n_basis());
                        V_multiple_.block(0,0, h_*q(), h_*n_basis()) = X_multiple_.transpose()*W_multiple_*Psi_multiple_;
                        // solve system (A_ + U_*(X^T*W_*X)*V_)x = b using woodbury formula from NLA module
                        sol = SMW<>().solve(invA_, U_multiple_, XtWX_multiple_, V_multiple_, b_); 
                        // store result of smoothing 
                        f_curr_    = sol.head(h_*n_basis());
                        fn_curr_ = Psi_multiple_*f_curr_; 
                        beta_curr_ = invXtWX_multiple_.solve(X_multiple_.transpose()*(W_bar_*z_ - W_multiple_*fn_curr_ + gamma0_*t));
                        
                        std::cout << "Beta_curr: " << beta_curr_ << std::endl; 

                    }
                    // store PDE misfit
                    g_curr_ = sol.tail(h_*n_basis());
                    
                    std::cout << "Range g_curr: " << g_curr_.minCoeff() << " , " << g_curr_.maxCoeff() << std::endl;
                    std::cout << "Range f_curr: " << f_curr_.minCoeff() << " , " << f_curr_.maxCoeff() << std::endl; 
                    
                    // update J 
                    J_old = J_new; 
                    J_new = model_loss() + g_curr_.dot(R0_multiple_*g_curr_) + gamma0_*crossing_penalty();   // R0 multiple already contains lambdas!
                    

                    std::cout << "----- crossing global at iter k = " << crossing_penalty() << std::endl;
                    std::cout << "----- J_old = " << J_old << std::endl;
                    std::cout << "----- J_new = " << J_new << " = " << model_loss() << " + " << g_curr_.dot(R0_multiple_*g_curr_) << std::endl;  
                    std::cout << "----- J_old - J_new = " << J_old - J_new << std::endl; 

                    k_++;  
            }

            double crossing_penalty_new = crossing_penalty(); 

            std::cout << "#################### cross new:  " << crossing_penalty_new << " , cross init: " << crossing_penalty_init << std::endl; 
            std::cout << "#################### cross new-init: " << std::setprecision(10) << (crossing_penalty_new - crossing_penalty_init) << std::endl; 
            
            gamma0_ *= C_;  
            iter_++;     
            
        }

    return;
    }   

    // Utilities 

    template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver>
        DVector<double> MSQRPDE<PDE, RegularizationType, SamplingDesign, Solver>::fitted() const{

        DVector<double> fit = fn_curr_; 
        if(has_covariates())
            fit += X_multiple_*beta_curr_;
        return fit; 
    }

    template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver>
        DVector<double> MSQRPDE<PDE, RegularizationType, SamplingDesign, Solver>::fitted(unsigned int j) const{
        return fitted().block(j*n_obs(), 0, n_obs(), 1); 
    }

    template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver>
        const bool MSQRPDE<PDE, RegularizationType, SamplingDesign, Solver>::crossing_constraints() const {
        // Return true if the current estimate of quantiles is crossing, false otherwise 

        bool crossed = false; 
        
        unsigned int j = 0;
        while(!crossed && j < h_-1){ 
            DVector<double> crossing_value = fitted(j+1)-fitted(j); 
            unsigned int i = 0; 
            while(!crossed && i < n_obs()){          
                if(crossing_value(i) < eps_)
                    crossed = true;  
                i++;
            }
            j++; 
        }

        if(crossed)
            std::cout << "crossed = true" << std::endl;
        else{
            std::cout << "crossed = false" << std::endl;
        }

        return crossed;
    }

    template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver>
        double MSQRPDE<PDE, RegularizationType, SamplingDesign, Solver>::model_loss() const{

        double loss = 0.; 
        for(auto j = 0; j < h_; ++j)
            loss += (rho_alpha(alphas_[j], y() - fitted(j))).sum(); 
        return loss/n_obs();

    }

    template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver>
        double MSQRPDE<PDE, RegularizationType, SamplingDesign, Solver>::crossing_penalty() const{
  
        // compute value of the unpenalized unconstrained functional J: 
        double pen = 0.; 
        for(int j = 0; j < h_-1; ++j){
            pen += (eps_*DVector<double>::Ones(n_obs()) - (fitted(j+1) - fitted(j))).cwiseMax(0.).sum(); 
        }   
        return pen; 

    }

    template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver>
        double MSQRPDE<PDE, RegularizationType, SamplingDesign, Solver>::crossing_penalty_f() const{
  
        // compute value of the unpenalized unconstrained functional J: 
        double pen = 0.; 
        for(int j = 0; j < h_-1; ++j){
            pen += (eps_*DVector<double>::Ones(n_obs()) - (fn_curr_.block((j+1)*n_obs(), 0, n_obs(), 1) - fn_curr_.block(j*n_obs(), 0, n_obs(), 1))).cwiseMax(0.).sum(); 
        }
        
        return pen; 
    }

    template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver>
        double MSQRPDE<PDE, RegularizationType, SamplingDesign, Solver>::crossing_penalty_param() const{
  
        // compute value of the unpenalized unconstrained functional J: 
        double pen = 0.; 
        for(int j = 0; j < h_-1; ++j){
            pen += (eps_*DVector<double>::Ones(n_obs()) - (X()*(beta_curr_.block((j+1)*q(), 0, q(), 1) - beta_curr_.block(j*q(), 0, q(), 1)))).cwiseMax(0.).sum(); 

        }
    
        return pen; 
    }

    template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver>
        const DMatrix<double>& MSQRPDE<PDE, RegularizationType, SamplingDesign, Solver>::H_multiple() {
        // compute H = X*(X^T*W*X)^{-1}*X^T*W
        H_multiple_ = X_multiple_*(invXtWX_multiple_.solve(X_multiple_.transpose()*W_multiple_));

        return H_multiple_;
    }

    template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver>
        const DMatrix<double>& MSQRPDE<PDE, RegularizationType, SamplingDesign, Solver>::Q_multiple() {
        // compute Q = W(I - H) = W ( I - X*(X^T*W*X)^{-1}*X^T*W ) 
        Q_multiple_ = W_multiple_*(DMatrix<double>::Identity(n_obs()*h_, n_obs()*h_) - H_multiple());

        return Q_multiple_;
    }

    // returns the pinball loss at a specific x 
    template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver>
    DVector<double> MSQRPDE<PDE, RegularizationType, SamplingDesign, Solver>::rho_alpha(const double& alpha, const DVector<double>& x) const{ 
        return 0.5*x.cwiseAbs() + (alpha - 0.5)*x; 
    }

    template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver>
        void MSQRPDE<PDE, RegularizationType, SamplingDesign, Solver>::abs_res_adj(DVector<double>& res) {
            unsigned int count_debug = 1; 
            for(int i = 0; i < res.size(); ++i) {
                if(res(i) < tol_weights_) {
                    count_debug++; 
                    res(i) += tol_weights_;  
                }            
            }
    }



} // namespace models
} // namespace fdapde
    
#endif // __MSQRPDE_H__
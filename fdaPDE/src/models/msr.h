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


#ifndef __MIXED_EFFECTS_SPATIAL_REGRESSION_H__
#define __MIXED_EFFECTS_SPATIAL_REGRESSION_H__

#include "header_check.h"

namespace fdapde {

template <typename VariationalSolver> class MSRPDE {

    private:
        using solver_t = std::decay_t<VariationalSolver>;
        using vector_t = Eigen::Matrix<double, Dynamic, 1>;
        using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
        using sparse_matrix_t = Eigen::SparseMatrix<double>;
        static constexpr int n_lambda = solver_t::n_lambda;

   public:
    MSRPDE() noexcept = default;
    template <typename GeoFrame, typename Penalty>
    MSRPDE(const std::string& formula, const GeoFrame& gf, Penalty&& penalty) noexcept :
        solver_() {
        discretize(penalty.get());
        analyze_data(formula, gf);
    }

    // modifiers
    template <typename... Args> void discretize(Args&&... args) { solver_.discretize(std::forward<Args>(args)...); }
    template <typename GeoFrame, typename WeightMatrix>
    void analyze_data(const std::string& formula, const GeoFrame& gf, const WeightMatrix& W) {
        fdapde_assert(gf.n_layers() == 1);
        Formula formula_(formula);
        n_obs_ = gf[0].rows();
        n_covs_ = 0;
        for (const std::string& token : formula_.covs()) {
            if (gf.contains(token)) { n_covs_++; }
        } 

        // check there is only one level, stop execution if not
        std::set<std::string> levels;
        for(const auto& token : formula_.efxs()) levels.insert(token.efx());
        fdapde_assert(std::cmp_equal(levels.size() FDAPDE_COMMA 1));
        // extract effects
        std::vector<std::string> efxs;
        for (const auto& token : formula_.efxs()) {
            if (token.cov() == std::string("1") || gf.contains(token.cov())) {
                efxs.push_back(gf.contains(token.cov()) ? token.cov() : std::string("1"));
            }
        }
        n_random_covs_ = efxs.size();

        // assemble Z matrix
        Z_.resize(n_obs_, n_random_covs_);    // note: at this step, n_obs_ = n_locs_
        std::cout << "Number of random effects = " << n_random_covs_ << std::endl;
        for (int i = 0; i < n_random_covs_; ++i) {
            if (efxs[i] == std::string("1")) {   // random intercept
                Z_.col(i).setOnes();
            } else {
                gf[0].data().template col<double>(efxs[i]).assign_to(Z_.col(i));
            }
        }

        // extract grouping structure
        group_ids_.resize(n_obs_);
        gf[0].template col<double>("g").assign_to(group_ids_);

        std::unordered_set<unsigned int> unique_ids(group_ids_.begin(), group_ids_.end());
        n_groups_ = unique_ids.size();
        std::cout << "Number of groups = " << n_groups_ << std::endl; 

        // Extract the size of each group
        group_sizes_.resize(n_groups_);
        loc_to_glob_map_.resize(n_groups_);

        for(int i=0; i < group_ids_.size(); ++i){  
            int i_loc = group_ids_(i)-1;           // ATT -1 because group ids start from 1          
            group_sizes_[i_loc] += 1;	           // update the counter
            loc_to_glob_map_[i_loc].push_back(i);  // map the local index to the global one
        }

        // // debug 
        // std::cout << "printing loc_to_glob_map_" << std::endl; 
        // for(int k=0; k<loc_to_glob_map_.size(); ++k){
        //     std::cout << std::endl; 
        //     for(int j=0; j<loc_to_glob_map_[k].size(); ++j){
        //         std::cout << loc_to_glob_map_[k][j] << "; ";  
        //     }    
        // }

        // save NAN pattern before correction with zeros 
        na_pattern_ = na_matrix(gf[0].data().template col<double>(formula_.lhs()));
        
        // // debug 
        // std::cout << "printing na_pattern_" << std::endl; 
        // for(int k=0; k<na_pattern_.size(); ++k){
        //     std::cout << na_pattern_[k] << "; "; 
        // }
        // std::cout << std::endl;

        solver_.analyze_data(formula, gf, W);   // M qui avviene la normalizzazione di W_ del solver e vengono corrette per i NA la Psi e la y del solver 
        y_ = solver_.response();                // M corretta per NA
        
        n_obs_ = solver_.n_obs(); // M added
        std::cout << "n_obs_ = " << n_obs_ << std::endl;
        std::cout << "number of data = " << gf[0].rows() << std::endl;
    }
    template <typename GeoFrame> void analyze_data(const std::string& formula, const GeoFrame& gf) {
        return analyze_data(formula, gf, vector_t::Ones(gf[0].rows()).asDiagonal());  // ATT ficticious weights initialization
    }
    void set_fpirls_tolerance(double tol) { tol_ = tol; }  // set FPIRLS convergence tolerance

    // fitting
    // Functional penalized iterative reweighted least squares
    template <typename... Args>
        requires(sizeof...(Args) > 0)
    auto fit(Args&&... args) {
        vector_t lambda(n_lambda);
        internals::for_each_index_and_args<sizeof...(Args)>(
          [&]<int Ns_, typename Ts_>(const Ts_& ts) {
              if (Ns_ < n_lambda) {
                  fdapde_static_assert(
                    std::is_convertible_v<Ts_ FDAPDE_COMMA double>, INVALID_SMOOTHING_PARAMETER_TYPE);
                  lambda[Ns_] = ts;
              }
          },
          args...);
        matrix_t y = y_;

        // initialization (nota M: eseguito due volte la prima volta che viene chiamato il fit, ma è necessario farlo se viene chiamato più volte il fit)      
        initial_weights_();   // smart Delta initialization + correction of Z for missing values
        solver_.update_response_and_weights(y, sparse_mat_weights_);   // restore solver state (qui i pesi vengono anche normalizzati)

        // nota: no scale lambda here (only for quantile regression)
        // mu_ = solver_.Psi() * solver_.f();   // M: no need of mu_ since does not enter in the abs residuals computation.... 

        double Jold = std::numeric_limits<double>::max(), Jnew = 0;
        n_iter_ = 0;
        std::cout << "Start FPIRLS with max_iter_=" << max_iter_ << " and tolerance=" << tol_ << std::endl;
        while (n_iter_ < max_iter_ && std::abs(Jnew - Jold) > tol_) {
            
            // std::cout << "FPIRLS iteration " << n_iter_+1 << std::endl;

            // compute pseudo observations
            py_ = y; 	

            // compute ZtildeTZtilde_ with current Delta_ and invert it 
            update_ZtildeTZtilde_(); 
              
            // compute weights
            update_pW_(); 

            // set weights and pseudo-observations to zero where there are missing values
            for(std::size_t i=0; i < y.size(); ++i){   // qui voglio loopare su tutto il vettore, non solo su quelli osservati 
            
                // Set to zeros the weights and pseudo-observations where there are missing values 
                if(na_pattern_[i]){
                    
                    // py_(i)=0.;  ATT commentato rispetto al codice vecchio, altrimenti perdo dopo nella chiamata solver_.update_response_and_weights non ha il giusto na_pattern, e n_obs e' errato... ---> comunque il risultato esce identico, anche nel caso con i missing data, quindi era una correzione inutile
                    
                    int block_idx = group_ids_(i)-1;  // ATT -1 because group ids start from 1
                    std::vector<unsigned int> glob_idxs_of_block = loc_to_glob_map_[block_idx]; 
                    unsigned int block_row_idx; 
                    for(int idx=0; idx < glob_idxs_of_block.size(); ++idx){
                        if(glob_idxs_of_block[idx] == i){
                            block_row_idx = idx; 
                        }
                    }
                    
                    for(int col_idx=0; col_idx<group_sizes_[block_idx]; ++col_idx){
                        pW_(block_idx)(block_row_idx, col_idx) = 0.;   // NOTA: se ho già settato le righe di Z a zero, questo set è inutile perchè questi sono già zeri
                    }

                }
            
            }
            update_sparse_mat_weights_(); // update sparse_mat_weights_ with current NA pattern of pW_

            // \argmin_{\beta, f} [ 1/n * \norm(W^{1/2} * (y - X * \beta - f_n))^2 + P_{\lambda}(f) ]
	        solver_.update_response_and_weights(py_, sparse_mat_weights_);           
            solver_.fit(std::forward<Args>(args)...);
            mu_ = fitted();   // fn + X%*%beta (no random part here!) 

            compute_bhat_();          // mu_ is needed here
            compute_sigma_sq_hat_();  // note: default value is false => metodo Melchionda (calcolo senza edf nelle fpirls iterations)  --> so this value is not stochastic since there are no stochastic edf
            build_LTL_();
            compute_C_();

            // update Delta_
            for(auto k=0; k<n_random_covs_; ++k){
                Delta_(k) = C_(k,k) * std::sqrt(n_groups_);
            }           

            // prepare for next iteration
            double data_loss = data_loss_();
            Jold = Jnew;
            Jnew = data_loss + solver_.ftPf(lambda);
            n_iter_++;

            std::cout << "data_loss=" << data_loss << std::endl; 
            std::cout << "penalty=" << solver_.ftPf(lambda) << std::endl; 

            std::cout << "|DeltaJ| at iter" << n_iter_ << " =" << std::abs(Jnew - Jold) << std::endl;
        }

        std::cout << "FPIRLS terminated after " << n_iter_ << " iterations with |DeltaJ|=" << std::abs(Jnew - Jold) << std::endl;
        std::cout << "Computing variance estimates at convergence..." << std::endl; 

        // compute sigma_sq_hat (with edf) at convergence 
        compute_sigma_sq_hat_(true);   
        // compute Sigma_b_ matrix at convergence 
        Sigma_b_ = Delta_;
        for(auto k=0; k < n_random_covs_; ++k){
            Sigma_b_(k) *= Delta_(k);
            Sigma_b_(k) = sigma_sq_hat_/Sigma_b_(k);   // ATT: assumes independence between random components
        }
        
        std::cout << "Final Delta_ = " << Delta_ << std::endl;
        std::cout << "Final sigma_sq_hat_ = " << sigma_sq_hat_ << std::endl;
        std::cout << "Final Sigma_b = " << Sigma_b_ << std::endl;

	return std::make_pair(solver_.f(), solver_.beta());
    }
    template <typename... Args> auto fit(Args&&... args) { return fit(std::forward<Args>(args)...); }

    // observers
    const vector_t& f() const { return solver_.f(); }
    vector_t fn() const { return solver_.fn(); }   // no const because the solver's getter is not const 
    const vector_t& beta() const { return solver_.beta(); }
    const std::vector<vector_t>& b_hat() const { return b_hat_; }
    const vector_t& misfit() const { return solver_.misfit(); }   // M: senza i random effects
    int n_covs() const { return n_covs_; }
    int n_random_covs() const { return n_random_covs_; }
    int n_obs() const { return n_obs_; }
    double edf(int r = 100, int seed = random_seed) { return solver_.edf(r, seed); }
    const vector_t& response() const { return solver_.response(); }
    double sigma_sq_hat() const {return sigma_sq_hat_;}
    const vector_t& Sigma_b() const {return Sigma_b_;}
    vector_t fitted() const {     // M: senza i random effects. Nota: qui NON devono esserci gli zeri in corrispondenza dei NA! 
        vector_t fitted_ = solver_.Psi() * f();
        if (n_covs_ != 0) { fitted_ += solver_.design_matrix() * beta(); }
        return fitted_;
    }
    vector_t random_effects() const {   // M: solo la parte dei random effects del fit 
        vector_t random_effects = vector_t::Zero(y_.size()); 
        for(int i=0; i<n_groups_; ++i){
            random_effects(loc_to_glob_map_[i]) += Z_by_group_(i) * b_hat_[i];
        }
        return random_effects;
    }
    const BinaryMatrix<-1, 1>& na_pattern() const {return na_pattern_;}
    unsigned int n_iter() const {return n_iter_;}

    // modifiers
    void set_fpirls_max_iter(int max_iter) { 
        std::cout << "setting max_iter fpirls to " << max_iter << std::endl; 
        max_iter_ = max_iter; 
    }    

    void set_likelihood_dataloss_type(bool type) {
        std::cout << "setting likelihood_dataloss_type to " << type << std::endl;
        likelihood_dataloss_type_ = type;
    }

    void set_compute_sigma_with_edf(bool type) {
        std::cout << "setting compute_sigma_with_edf to " << type << std::endl;
        compute_sigma_with_edf_ = type;
    }


    // Generalized Cross Validation index
    struct gcv_t : public ScalarFieldBase<n_lambda, gcv_t> {
        using Base = ScalarFieldBase<1, gcv_t>;
        static constexpr int StaticInputSize = n_lambda;
        static constexpr int NestAsRef = 0;
        static constexpr int XprBits = 0;
        using Scalar = double;
        using InputType = Vector<Scalar, StaticInputSize>;
        using edf_cache_t = std::unordered_map<
          std::array<double, StaticInputSize>, double, internals::std_array_hash<double, StaticInputSize>>;

        gcv_t() noexcept = default;
        gcv_t(MSRPDE* model, const edf_cache_t& edf_cache) :
            model_(model),
            n_(model->n_obs()),
            q_(model->n_covs()),
            p_(model->n_random_covs()),
            edf_cache_(edf_cache),
            r_(100),
            seed_(random_seed) { }
        gcv_t(MSRPDE* model, const edf_cache_t& edf_cache, int r, int seed) :
            model_(model), n_(model->n_obs()), q_(model->n_covs()), p_(model->n_random_covs()), edf_cache_(edf_cache), r_(r), seed_(seed) { }
        gcv_t(MSRPDE* model) : gcv_t(model, edf_cache_t()) { }
        gcv_t(MSRPDE* model, int r, int seed) : gcv_t(model, edf_cache_t(), r, seed) { }

        template <typename InputType_>
            requires(internals::is_subscriptable<InputType_, int>)
        constexpr double operator()(const InputType_& lambda) {
            return internals::apply_index_pack<n_lambda>([&]<int... Ns_>() { return operator()(lambda[Ns_]...); });
        }
        template <typename... LambdaT>
            requires(std::is_convertible_v<LambdaT, double> && ...)
        constexpr double operator()(LambdaT... lambda) {
            model_->fit(static_cast<double>(lambda)...);
            std::array<double, StaticInputSize> lambda_vec {lambda...};
            if (edf_cache_.find(lambda_vec) == edf_cache_.end()) {   // cache Tr[S]
                edf_cache_[lambda_vec] = model_->edf(r_, seed_);
            }
            double dor = n_ - (q_ + edf_cache_.at(lambda_vec));   // residual degrees of freedom
            
	        double norm = 0.;
            vector_t op1 = model_->response();            
            vector_t op2 = model_->fitted() + model_->random_effects();  // M: va aggiunta la random part perché non è contenuta in fitted() 
            
            // compute norm only on observed data
            for (int i = 0; i < op1.size(); ++i) {
                if (!model_->na_pattern()[i]) norm += (op2.coeff(i, 0) - op1.coeff(i, 0))*(op2.coeff(i, 0) - op1.coeff(i, 0));
            }
            // return (n_ / std::pow(dor, 2)) * (model_->fitted() - model_->response()).squaredNorm();  --> M: non tiene conto dei NA!
            return (n_ / std::pow(dor, 2)) * norm;    
        }
        // observers
        const edf_cache_t& edf_cache() const { return edf_cache_; }
        edf_cache_t& edf_cache() { return edf_cache_; }
       private:
        MSRPDE* model_;
        int n_ = 0, q_ = 0, p_ = 0;
        edf_cache_t edf_cache_;
        // stochastic edf approximation parameter
        int r_, seed_;
    };
    friend gcv_t;
    gcv_t gcv() { return gcv_t(this); }
    gcv_t gcv(const typename gcv_t::edf_cache_t& edf_cache) { return gcv_t(this, edf_cache); }
    gcv_t gcv(int r, int seed) { return gcv_t(this, r, seed); }
    gcv_t gcv(const typename gcv_t::edf_cache_t& edf_cache, int r, int seed) { return gcv_t(this, edf_cache, r, seed); }

    private:
        vector_t y_;
        BinaryMatrix<-1, 1> na_pattern_; 

        sparse_matrix_t sparse_mat_weights_ {}; 
        Eigen::Matrix<matrix_t, Dynamic, 1> pW_ {};   // diagonal blocks of W^k (one block for each group)
   
        vector_t py_;       // pseudo observations 

        vector_t mu_;      // \mu^k = [ \mu^k_1, ..., \mu^k_n ] : fitted vector at step k
        Eigen::SparseLU<sparse_matrix_t> invA_;

        unsigned int n_groups_; 
        std::vector<unsigned int> group_sizes_; 
        std::vector<std::vector<unsigned int>> loc_to_glob_map_; 
        Eigen::Matrix<unsigned int, Dynamic, 1> group_ids_; 

        matrix_t Z_;                                // random effect design matrix
        Eigen::Matrix<matrix_t, Dynamic, 1> Z_by_group_;  // vectors of blocks storing the random effect matrices Z_i, i=1,...,n_groups_ (by groups)
        Eigen::Matrix<matrix_t, Dynamic, 1> ZTZ_;  // vectors of blocks storing the random effect matrices Z_i^T * Z_i
        Eigen::Matrix<Eigen::LLT<matrix_t>, Dynamic, 1> ZtildeTZtilde_; 
        std::vector<vector_t> b_hat_;  // vector containing the estimates of the random effects for each group
        Eigen::LLT<matrix_t> LTL_; 
        matrix_t C_;       // Cholesky factor of LTL_
        vector_t Delta_;   // ATT: we are assuming independent random effects  
        vector_t Sigma_b_; // ATT: we are assuming independent random effects  
        double sigma_sq_hat_;     // estimate of the variance of the errors 

        int max_iter_ = 200;    // fpirls maximum iteration number
        double tol_ = 1e-6;     // fpirls convergence tolerance

        solver_t solver_;
        int n_obs_ = 0, n_covs_ = 0, n_random_covs_ = 0; 
        int n_iter_ = 0;

        bool likelihood_dataloss_type_ = false;   // default: loss "stile FPIRLS" (sempre positiva)
        bool compute_sigma_with_edf_ = true;      // default: calcolo sigma_sq_hat_ CON edf anche nelle fpirls iterations
        // nota: la scelta dei default è basato su quanto osservato in test 6 

    private:

        // Smart Delta initialization
        void Delta_init_(){

            // // for debug: print all the vector na_pattern_
            // for(int k=0; k<na_pattern_.size(); ++k){
            //     std::cout << na_pattern_[k] << "; "; 
            // }
            
            // Compute Z_by_group_ and ZTZ_ for each group
            for(int i=0; i < n_groups_; ++i){

                Z_by_group_(i) = matrix_indexing_(Z(), loc_to_glob_map_[i]);  
                
                // metto a zero le righe di Z che hanno missing data -> questo serve per calcolo di Delta_, Ztilde e quindi b_i, sigma_sq_hat_ etc..
                for(int glob_idx : loc_to_glob_map_[i]){
                    
                    if(na_pattern_[glob_idx]){ 
                        
                        std::vector<unsigned int> glob_idxs_of_block = loc_to_glob_map_[i]; 
                        unsigned int block_row_idx; 
                        for(int idx=0; idx < glob_idxs_of_block.size(); ++idx){
                            if(glob_idxs_of_block[idx] == glob_idx){
                                block_row_idx = idx; 
                            }
                        }

                        // correct Z_by_group_
                        (Z_by_group_(i).row(block_row_idx)).setZero();

                        // correct Z_
                        (Z_.row(glob_idx)).setZero();
                    }
                }  
                
                ZTZ_(i) = Z_by_group_(i).transpose() * Z_by_group_(i);
            }

            // initialize Delta_
            for(int k=0; k < n_random_covs_; ++k){
                Delta_(k) = 0.; 
                for(int i=0; i < n_groups_; i++){
                    for(int j=0; j < group_sizes_[i]; j++){
                        Delta_(k) += Z_by_group_(i)(j,k) * Z_by_group_(i)(j,k); 
                    }
                }
                Delta_(k) = std::sqrt( Delta_(k)/n_groups_ ) * 3 / 8;  // ATT 3/8 fa zero (o metti il punto o metti 3/8 dopo)
            }
        
        }

        // Update sparse_mat_weights_ with current pW_
        void update_sparse_mat_weights_(){

            sparse_mat_weights_.setZero();  // clean the object 

            // Store room for triplets
            std::vector<Eigen::Triplet<double>> triplet_list;  
            unsigned int size_triplet = 0; 
            for(const auto& pw_block : pW_) {
                size_triplet += pw_block.rows()*pw_block.cols();
            }
            triplet_list.reserve(size_triplet);

            // Fill the sparse matrix using the local to global map
            for(auto k=0; k<n_groups_; ++k){
                const auto pw_block = pW_(k); 
                for(auto i=0; i<group_sizes_[k]; ++i){
                    for(auto j=0; j<group_sizes_[k]; ++j){
                        triplet_list.emplace_back(loc_to_glob_map_[k][i], loc_to_glob_map_[k][j], pw_block(i, j)); 
                    }
                }
            }

            // finalize construction
            sparse_mat_weights_.setFromTriplets(triplet_list.begin(), triplet_list.end());
            sparse_mat_weights_.makeCompressed();

        }

        // Pre-allocate memory for all quatities
        void resize_matrices_(){
            std::cout << "n_groups_ = " << n_groups_ << std::endl;
            Z_by_group_.resize(n_groups_);
            for(int i=0; i<n_groups_; ++i){
                Z_by_group_(i).resize(group_sizes_[i], n_random_covs_); 
            }
            ZTZ_.resize(n_groups_);
            for(int i=0; i<n_groups_; ++i){
                ZTZ_(i).resize(n_random_covs_, n_random_covs_); 
            }
            ZtildeTZtilde_.resize(n_groups_);
            Delta_.resize(n_random_covs_);
            Sigma_b_.resize(n_random_covs_);

            pW_.resize(n_groups_); 
            for(int i=0; i<n_groups_; ++i){
                pW_(i).resize(group_sizes_[i], group_sizes_[i]);
            }

            int total_size = 0;  // total size of the sparse matrix
            for(const auto& mat : pW_) {
                total_size += mat.rows();   
            }
            sparse_mat_weights_.resize(total_size, total_size); // resize sparse_mat_weights_
        }

        // Weights initialization
        void initial_weights_() {
            resize_matrices_(); // allocate memory
            Delta_init_();      // smart initialization of Delta_
            update_ZtildeTZtilde_(); 
            update_pW_(); 
        }

        void update_ZtildeTZtilde_(){
            for(auto i=0; i < n_groups_; ++i){

                // For each group, consider Z^T Z
                matrix_t ZtildeTZtilde_temp = ZTZ_(i);
                
                // Then add the diagonal elements stored in D_
                for(auto k=0; k < n_random_covs_; ++k){
                    ZtildeTZtilde_temp(k,k) += Delta_(k) * Delta_(k);
                }	
                
                ZtildeTZtilde_(i).compute(ZtildeTZtilde_temp);
            }
        }

        void update_pW_(){

            for(int i=0; i < n_groups_; ++i){

                // Compute the current block with the Woodbury identity
                pW_(i) = - Z_by_group_(i) * ZtildeTZtilde_(i).solve( Z_by_group_(i).transpose() );  // I_ni - Z_i * (ZtildeTZtilde_i_)^(-1) * Z_i^T

                // Add the identity matrix (1 on the diagonal)
                for(int k=0; k < group_sizes_[i]; ++k){
                    pW_(i)(k,k) = 1 + pW_(i)(k,k);	
                }
            } 
            
            update_sparse_mat_weights_();  // update sparse_mat_weights_ with current pW_

        }

        // returns the data loss (J_parametric)
        double data_loss_() const { 
            
            if(likelihood_dataloss_type_){

                // J_parametric = -0.5*(n_groups*p-n)log(sigma^2) - 0.5*(|| Delta*b_i ||/ sigma)^2 + n_groups_*log(det(Delta))
                std::cout << "!!! --versione MELCHIONDA data loss-- !!!" << std::endl;
                
                double data_loss_value = 0.;

                // cast to int to avoid overflow
                int signed_int = n_groups_*n_random_covs_ - n_obs();  
                
                data_loss_value -= signed_int * std::log(sigma_sq_hat_);
                
                for(auto i=0; i < n_groups_; ++i){
                    // log-likelihood of random effects	(completed outside the for cycle)
                    vector_t Deltab_i = Delta_.asDiagonal() * b_hat_[i];
                    data_loss_value -= ( Deltab_i ).dot( Deltab_i ) / sigma_sq_hat_;
                }
                
                // Compute the determinant of Delta (NOTE: Delta is a diagonal matrix stored in a vector!)
                double detDelta = 1.;
                for(auto k=0; k < n_random_covs_; ++k){
                    detDelta *= Delta_(k);
                }
                data_loss_value += 2 * n_groups_ * std::log(detDelta);

                return data_loss_value/2;  

            } else{

                std::cout << "!!! --versione FPIRLS data loss-- !!!" << std::endl;
                double data_loss_value = 0.;
                // Compute the square root of the weights matrix with Cholosky
                Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> chol(sparse_mat_weights_);
                Eigen::SparseMatrix<double> sqrtW = chol.matrixL();                
                vector_t data_loss_vector = sqrtW * (py_ - (mu_ + random_effects()) ); 
                
                for(int i = 0; i < data_loss_vector.size(); ++i) {
                    if(!na_pattern_[i]) data_loss_value += (data_loss_vector.coeff(i, 0))*(data_loss_vector.coeff(i, 0));
                }

                return data_loss_value / n_obs_;


            }
  




        }

        // helper functions
        vector_t vector_indexing_(const vector_t& big_vector, const std::vector<unsigned int> ids){

            // ids: vector of global indexes 

            vector_t small_vector(ids.size());
            for(int k=0; k < ids.size(); ++k){
                small_vector(k) = big_vector(ids[k]);
            }
            
            return small_vector;
        }

        matrix_t matrix_indexing_(const matrix_t& big_matrix, const std::vector<unsigned int> row_ids){

            // row_ids: vector of global row indexes
            
            matrix_t small_matrix(row_ids.size(), big_matrix.cols());

            // Reconstruct the block using loc_to_glob_map
            for (unsigned int i = 0; i < row_ids.size(); ++i) {
                small_matrix.row(i) = big_matrix.row(row_ids[i]);
            }
            return small_matrix;

        }

        matrix_t X() const { return solver_.design_matrix(); }
        matrix_t Z() const { return Z_; } 

        // compute b_hat_
        void compute_bhat_() { 
            b_hat_.resize(n_groups_); 
            for(auto i=0; i < n_groups_; ++i){
                vector_t res_i = vector_indexing_(y_, loc_to_glob_map_[i]);
                res_i -= vector_indexing_(mu_, loc_to_glob_map_[i]);   // here mu_ does NOT contain the random effect
                
                // M: for missing (ATT aggiunto rispetto a Melchionda)
                for(int j=0; j<res_i.size(); ++j){
                    if(na_pattern_[loc_to_glob_map_[i][j]]){
                        res_i(j) = 0.; 
                    }
                }

                
                b_hat_[i] = ZtildeTZtilde_(i).solve( Z_by_group_(i).transpose() * res_i );
            }
        }

        // compute sigma_sq_hat_
        void compute_sigma_sq_hat_(bool edf_flag=false, int seed=1234) {

            sigma_sq_hat_ = 0.;	

            // set to fit vector zeros where there are missing values, so that the residual is computed correctly (necessary since mu_ does not have zeros in NA indexes)
            vector_t fit_adj = mu_;    
            for(int j=0; j<fit_adj.size(); ++j){
                if(na_pattern_[j]){
                    fit_adj(j) = 0.; 
                }
            }

            for(auto i=0; i < n_groups_; i++){ 

                vector_t mu_i = vector_indexing_(fit_adj, loc_to_glob_map_[i]);
                vector_t res_i = vector_indexing_(y_, loc_to_glob_map_[i]);  // note: y_ has zeros in correspondence of missing values
                
                // res_i -= ( mu_i + Z_by_group_(i)*b_hat_[i] );   //  here the residual contains the random effect too 
                // // note: sum only over observed data (funziona se nell'init abbiamo settato in Z_by_group_ zero rows in correspondence of missing values)

                // M: set zeros in random part for missing values case (when Z_by_group_ is not zeroed)
                vector_t Zb_i = Z_by_group_(i)*b_hat_[i];
                for(int j=0; j<res_i.size(); ++j){
                    if(na_pattern_[loc_to_glob_map_[i][j]]){
                        Zb_i(j) = 0.; 
                    } 
                }
                res_i -= ( mu_i + Zb_i );   //  here the residual contains the random effect too 
                // note: sum only over observed data

                sigma_sq_hat_ += res_i.dot(res_i);
            }


            if(compute_sigma_with_edf_){
                // Versione Pigani 
                std::cout << "--sigma2 computation with edf at each iteration--" << std::endl;
                double edf_value = edf(100, seed);   // here we set a seed for the edf stochastic computation for reproducibility
                if(n_covs_ != 0){
                    edf_value += n_covs_;   // n_random_covs_? Pigani non lo mette  
                }
                sigma_sq_hat_ /= (n_obs()-edf_value); 
            } else{
                // Versione Melchionda 
                std::cout << "--sigma2 computation with edf only at convergence--" << std::endl;
                if(edf_flag){   

                    double edf_ = edf(100, seed);   // here we set a seed for the edf stochastic computation for reproducibility
                    if(n_covs_ != 0){
                        edf_ += n_covs_;   // +m*n_random_covs_?
                    }
                    sigma_sq_hat_ /= (n_obs()-edf_); 

                    std::cout << "edf = " << std::setprecision(16) << edf_ << std::endl;

                } else{
                    sigma_sq_hat_ /= n_obs();  
                }
    

            }
  
        }

        // compute LTL_ 
        void build_LTL_(){

            matrix_t LTL_temp = matrix_t::Zero(n_random_covs_, n_random_covs_);

            for(auto i=0; i < n_groups_; ++i){
                LTL_temp += b_hat_[i]*(b_hat_[i]).transpose() / sigma_sq_hat_;
                LTL_temp += ZtildeTZtilde_(i).solve(matrix_t::Identity(n_random_covs_, n_random_covs_));
            }
            
            LTL_.compute(LTL_temp); // ATT: rispetto alla formula della teoria, non manca un /n_groups_? no, perchè inserito dopo il calcolo di C 
        }

        // compute C_ 
        void compute_C_(){
            matrix_t C_temp = LTL_.matrixL();
            C_ = C_temp.triangularView<Eigen::Lower>().solve( matrix_t::Identity(n_random_covs_, n_random_covs_) );
        }



}; 

// deduction guide
template <typename GeoFrame, typename Penalty>
MSRPDE(const std::string& formula, const GeoFrame& gf, Penalty&& solver)
  -> MSRPDE<typename Penalty::solver_t>;

}   // namespace fdapde

#endif   // __MIXED_EFFECTS_SPATIAL_REGRESSION_H__
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

#ifndef __GENERALIZED_SPATIAL_REGRESSION_H__
#define __GENERALIZED_SPATIAL_REGRESSION_H__

#include "header_check.h"

namespace fdapde {

template <typename VariationalSolver, typename Distribution> class GSRPDE {
   private:
    using solver_t = std::decay_t<VariationalSolver>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    static constexpr int n_lambda = solver_t::n_lambda;
   public:
    GSRPDE() noexcept = default;
    template <typename GeoFrame, typename Penalty>
    GSRPDE(const std::string& formula, const GeoFrame& gf, const Distribution& distr, Penalty&& penalty) noexcept :
        distr_(distr), solver_() {
        fdapde_assert(gf.n_layers() == 1);
        Formula formula_(formula);
        n_obs_ = gf[0].rows();
        n_covs_ = 0;
        for (const std::string& token : formula_.rhs()) {
            if (gf.contains(token)) { n_covs_++; }
        }
        // discretize
        if constexpr (requires(Penalty p) { p.get(); }) {
            solver_ = solver_t(formula, gf, penalty.get());
        } else {
            solver_ = solver_t(formula, gf, penalty(gf.template triangulation<0>()).get());
        }
	y_ = solver_.response();
    }
    
    // Functional penalized iterative reweighted least squares
    template <typename... LambdaT>
        requires(std::is_convertible_v<LambdaT, double> && ...)
    void fit(LambdaT... lambda) {
        // initialize mean vector
        vector_t y = y_;
        solver_.update_response_and_weights(y, vector_t::Ones(n_obs_).asDiagonal());   // restore solver state
        if constexpr (requires(Distribution d, vector_t v) { d.transform(v); }) {
            mu_ = distr_.transform(y);
        } else {
            mu_ = y;
        }
        double Jold = std::numeric_limits<double>::max(), Jnew = 0;
        n_iter_ = 0;
        while (n_iter_ < max_iter_ && std::abs(Jnew - Jold) > tol_) {
            vector_t G = distr_.der_link(mu_);   // G^(k) = diag(g'(\mu^(k)_1), ..., g'(\mu^(k)_n))
            pW_ = ((G.array().pow(2) * distr_.variance(mu_).array()).inverse()).matrix();
            py_ = G.asDiagonal() * (y - mu_) + distr_.link(mu_);
            // \argmin_{\beta, f} [ \norm(W^{1/2} * (y - X * \beta - f_n))^2 + P_{\lambda}(f) ]
	    solver_.update_response_and_weights(py_, pW_.asDiagonal());
            solver_.fit(lambda...);
            mu_ = distr_.inv_link(fitted());    
            // prepare for next iteration
            double data_loss =
              (distr_.variance(mu_).array().sqrt().inverse().matrix().asDiagonal() * (y - mu_)).squaredNorm() / n_obs_;
            Jold = Jnew;
            Jnew = data_loss + solver_.ftPf(lambda...);
	    n_iter_++;
        }
	return;
    }
    // observers
    const vector_t& f() const { return solver_.f(); }
    const vector_t& beta() const { return solver_.beta(); }
    int n_covs() const { return n_covs_; }
    int n_obs() const { return n_obs_; }
    double edf(int r = 100, int seed = random_seed) { return solver_.edf(r, seed); }
    const vector_t& response() const { return solver_.response(); }
    vector_t fitted() const {
        matrix_t fitted_ = solver_.Psi() * f();
        if (n_covs_ != 0) { fitted_ += solver_.design_matrix() * beta(); }
        return fitted_;
    }

    // Generalized Cross Validation index
    struct gcv_t : public ScalarFieldBase<n_lambda, gcv_t> {
        using Base = ScalarFieldBase<1, gcv_t>;
        static constexpr int StaticInputSize = n_lambda;
        static constexpr int NestAsRef = 0;
        static constexpr int XprBits = 0;
        using Scalar = double;
        using InputType = Vector<Scalar, StaticInputSize>;

        gcv_t() noexcept = default;
        gcv_t(GSRPDE* model) : model_(model), n_(model->n_obs()), q_(model->n_covs()), r_(100), seed_(random_seed) { }
        gcv_t(GSRPDE* model, int r, int seed) :
            model_(model), n_(model->n_obs()), q_(model->n_covs()), r_(r), seed_(seed) { }

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
            if (edf_map_.find(lambda_vec) == edf_map_.end()) {   // cache Tr[S]
                edf_map_[lambda_vec] = model_->edf(r_, seed_);
            }
            double dor = n_ - (q_ + edf_map_.at(lambda_vec));   // residual degrees of freedom
	    // compute total deviance
            vector_t mu = model_->distr_.inv_link(model_->fitted());
            double total_deviance = 0;
            for (int i = 0; i < n_; ++i) { total_deviance += model_->distr_.deviance(mu[i], model_->y_(i, 0)); }
            return (n_ / std::pow(dor, 2)) * total_deviance;
        }
       private:
        GSRPDE* model_;
        int n_ = 0, q_ = 0;
        std::unordered_map<
          std::array<double, StaticInputSize>, double, internals::std_array_hash<double, StaticInputSize>>
          edf_map_;
        // stochastic edf approximation parameter
        int r_, seed_;
    };
    gcv_t gcv() { return gcv_t(this); }
    gcv_t gcv(int r, int seed) { return gcv_t(this, r, seed); }

    // inference
  
   private:
    vector_t y_;
    vector_t mu_;          // \mu^k = [ \mu^k_1, ..., \mu^k_n ] : mean vector at step k
    vector_t py_;          // \tilde y^k = G^k(y-u^k) + \theta^k
    vector_t pW_;          // diagonal of W^k = ((G^k)^{-2})*((V^k)^{-1})
    int max_iter_ = 200;   // fpirls maximum iteration number
    double tol_ = 1e-6;    // fprils convergence tolerance

    Distribution distr_;
    solver_t solver_;
    int n_obs_ = 0, n_covs_ = 0;
    int n_iter_ = 0;
};

// deduction guide
template <typename GeoFrame, typename Distribution, typename Penalty>
GSRPDE(const std::string& formula, const GeoFrame& gf, const Distribution& distr, Penalty&& solver)
  -> GSRPDE<typename Penalty::solver_t, Distribution>;

}   // namespace fdapde

#endif   // __GENERALIZED_SPATIAL_REGRESSION_H__

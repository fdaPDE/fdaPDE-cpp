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

#ifndef __SPATIAL_REGRESSION_H__
#define __SPATIAL_REGRESSION_H__

#include "header_check.h"

namespace fdapde {

template <typename VariationalSolver> class SRPDE {
   private:
    using solver_t = std::decay_t<VariationalSolver>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    static constexpr int n_lambda = solver_t::n_lambda;
   public:
    SRPDE() noexcept = default;
    template <typename GeoFrame, typename Penalty>
    SRPDE(const std::string& formula, const GeoFrame& gf, Penalty&& penalty) noexcept : solver_() {
        fdapde_assert(gf.n_layers() == 1);	
        Formula formula_(formula);
	n_obs_  = gf[0].rows();
	n_covs_ = 0;
        for (const std::string& token : formula_.rhs()) {
            if (gf.contains(token)) { n_covs_++; }
        }
        if constexpr (requires(Penalty p) { p.get(); }) {
            solver_ = solver_t(formula, gf, penalty.get());
        } else {
            solver_ = solver_t(formula, gf, penalty(gf.template triangulation<0>()).get());
        }
    }
    template <typename... LambdaT>
        requires(std::is_convertible_v<LambdaT, double> && ...) ||
                (sizeof...(LambdaT) == 1 && (internals::is_vector_like_v<LambdaT> && ...))
    void fit(LambdaT... lambda) {
        solver_.fit(lambda...);
    }
    // observers
    const vector_t& f() const { return solver_.f(); }
    const vector_t& beta() const { return solver_.beta(); }
    int n_covs() const { return n_covs_; }
    int n_obs() const { return n_obs_; }
    double edf(int r = 100, int seed = random_seed) { return solver_.edf(r, seed); }
    const vector_t& response() const { return solver_.response(); }
    vector_t fitted() const {
        vector_t fitted_ = solver_.fn();
        if constexpr (requires(solver_t s) { s.design_matrix(); }) {
            if (n_covs_ != 0) { fitted_ += solver_.design_matrix() * beta(); }
        }
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
        using edf_cache_t = std::unordered_map<
          std::array<double, StaticInputSize>, double, internals::std_array_hash<double, StaticInputSize>>;

        gcv_t() noexcept = default;
        gcv_t(SRPDE* model, const edf_cache_t& edf_cache) :
            model_(model), n_(model->n_obs()), q_(model->n_covs()), edf_cache_(edf_cache), r_(100), seed_(random_seed) { }
        gcv_t(SRPDE* model, const edf_cache_t& edf_cache, int r, int seed) :
            model_(model), n_(model->n_obs()), q_(model->n_covs()), edf_cache_(edf_cache), r_(r), seed_(seed) { }
        gcv_t(SRPDE* model) : gcv_t(model, edf_cache_t()) { }
        gcv_t(SRPDE* model, int r, int seed) : gcv_t(model, edf_cache_t(), r, seed) { }

        template <typename InputType_>
            requires(internals::is_subscriptable<InputType_, int>)
        constexpr double operator()(const InputType_& lambda) {
            return internals::apply_index_pack<n_lambda>([&]<int... Ns_>() { return operator()(lambda[Ns_]...); });
        }
        template <typename... LambdaT>
            requires(std::is_convertible_v<LambdaT, double> && ...) && (sizeof...(LambdaT) == StaticInputSize)
        constexpr double operator()(LambdaT... lambda) {
            model_->fit(static_cast<double>(lambda)...);
            std::array<double, StaticInputSize> lambda_vec {lambda...};
            if (edf_cache_.find(lambda_vec) == edf_cache_.end()) {   // cache Tr[S]
                edf_cache_[lambda_vec] = model_->edf(r_, seed_);
            }
            double dor = n_ - (q_ + edf_cache_.at(lambda_vec));   // residual degrees of freedom
            return (n_ / std::pow(dor, 2)) * (model_->fitted() - model_->response()).squaredNorm();
        }
        // observers
        const edf_cache_t& edf_cache() const { return edf_cache_; }
        edf_cache_t& edf_cache() { return edf_cache_; }
       private:
        SRPDE* model_;
        int n_ = 0, q_ = 0;
        edf_cache_t edf_cache_;
        // stochastic edf approximation parameter
        int r_, seed_;
    };
    gcv_t gcv() { return gcv_t(this); }
    gcv_t gcv(const typename gcv_t::edf_cache_t& edf_cache) { return gcv_t(this, edf_cache); }
    gcv_t gcv(int r, int seed) { return gcv_t(this, r, seed); }
    gcv_t gcv(const typename gcv_t::edf_cache_t& edf_cache, int r, int seed) { return gcv_t(this, edf_cache, r, seed); }

    // inference
  
   private:
    solver_t solver_;
    int n_obs_ = 0, n_covs_ = 0;
};

// deduction guide
template <typename GeoFrame, typename Penalty>
SRPDE(const std::string& formula, const GeoFrame& gf, Penalty&& solver) -> SRPDE<typename Penalty::solver_t>;

}   // namespace fdapde

#endif //  __SPATIAL_REGRESSION_H__

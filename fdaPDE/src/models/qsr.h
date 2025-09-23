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

#ifndef __QUANTILE_SPATIAL_REGRESSION_H__
#define __QUANTILE_SPATIAL_REGRESSION_H__

#include "header_check.h"

namespace fdapde {

template <typename VariationalSolver> class QSRPDE {
   private:
    using solver_t = std::decay_t<VariationalSolver>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    static constexpr int n_lambda = solver_t::n_lambda;
   public:
    QSRPDE() noexcept = default;
    template <typename GeoFrame, typename Penalty>
    QSRPDE(const std::string& formula, const GeoFrame& gf, double alpha, Penalty&& penalty) noexcept :
        solver_(), alpha_(alpha) {
        discretize(penalty.get());
        analyze_data(formula, gf);
    }
    template <typename GeoFrame, typename Penalty>   // default to median fitting
    QSRPDE(const std::string& formula, const GeoFrame& gf, Penalty&& penalty) noexcept :
        QSRPDE(formula, gf, 0.5, penalty) { }

    // modifiers
    void set_level(double alpha) { alpha_ = alpha; }
    template <typename... Args> void discretize(Args&&... args) { solver_.discretize(std::forward<Args>(args)...); }
    template <typename GeoFrame, typename WeightMatrix>
    void analyze_data(const std::string& formula, const GeoFrame& gf, const WeightMatrix& W) {
        fdapde_assert(gf.n_layers() == 1);
        Formula formula_(formula);
        n_obs_ = gf[0].rows();
        n_covs_ = 0;
        for (const std::string& token : formula_.rhs()) {
            if (gf.contains(token)) { n_covs_++; }
        }
        solver_.analyze_data(formula, gf, W);
        y_ = solver_.response();
    }
    template <typename GeoFrame> void analyze_data(const std::string& formula, const GeoFrame& gf) {
        return analyze_data(formula, gf, vector_t::Ones(gf[0].rows()).asDiagonal());
    }
    // fitting
    // Functional penalized iterative reweighted least squares
    template <typename... Args>
        requires(sizeof...(Args) > 0)
    auto fit(double alpha, Args&&... args) {
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
        // initialization
        solver_.update_response_and_weights(y, vector_t::Ones(n_obs_).asDiagonal());   // restore solver state
        internals::apply_index_pack_and_args<sizeof...(Args)>(
          [&]<typename... Ts_>(Ts_... ts) {
              solver_.nonparametric_fit([&]() {
                  if constexpr (Ts_::index < n_lambda) {
                      return 2. * ts.value;   // scale smoothing parameter
                  } else {
                      return ts.value;
                  }
              }()...);
          },
          args...);
        mu_ = solver_.Psi() * solver_.f();	
        double Jold = std::numeric_limits<double>::max(), Jnew = 0;
        n_iter_ = 0;
        while (n_iter_ < max_iter_ && std::abs(Jnew - Jold) > tol_) {
            vector_t abs_res = (y - mu_).array().abs();
            // W_i = 0.5 * (abs_res[i] + tol_weights_) if abs_res[i] < tol_weights, W_i = 0.5 * abs_res[i] otherwise
            pW_ = (abs_res.array() < tol_weights_)
                    .select((2. * (abs_res.array() + tol_weights_)).inverse(), (2. * abs_res.array()).inverse());
            py_ = y - (1 - 2. * alpha) * abs_res;	  
            // \argmin_{\beta, f} [ 1/n * \norm(W^{1/2} * (y - X * \beta - f_n))^2 + P_{\lambda}(f) ]
	    solver_.update_response_and_weights(py_, pW_.asDiagonal());
            solver_.fit(std::forward<Args>(args)...);
            mu_ = fitted();
            // prepare for next iteration
            double data_loss = (pW_.cwiseSqrt().matrix().asDiagonal() * (py_ - mu_)).squaredNorm() / n_obs_;
            Jold = Jnew;
            Jnew = data_loss + solver_.ftPf(lambda);
            n_iter_++;
        }
	return std::make_pair(solver_.f(), solver_.beta());
    }
    template <typename... Args> auto fit(Args&&... args) { return fit(alpha_, std::forward<Args>(args)...); }
    // observers
    const vector_t& f() const { return solver_.f(); }
    const vector_t& beta() const { return solver_.beta(); }
    const vector_t& misfit() const { return solver_.misfit(); }
    int n_covs() const { return n_covs_; }
    int n_obs() const { return n_obs_; }
    double edf(int r = 100, int seed = random_seed) { return solver_.edf(r, seed); }
    const vector_t& response() const { return solver_.response(); }
    vector_t fitted() const {
        vector_t fitted_ = solver_.Psi() * f();
        if (n_covs_ != 0) { fitted_ += solver_.design_matrix() * beta(); }
        return fitted_;
    }
    // modifiers
    void set_pinball_smoothing_factor(double eps) { eps_ = eps; }

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
        gcv_t(QSRPDE* model, const edf_cache_t& edf_cache) :
            model_(model),
            n_(model->n_obs()),
            q_(model->n_covs()),
            edf_cache_(edf_cache),
            r_(100),
            seed_(random_seed) { }
        gcv_t(QSRPDE* model, const edf_cache_t& edf_cache, int r, int seed) :
            model_(model), n_(model->n_obs()), q_(model->n_covs()), edf_cache_(edf_cache), r_(r), seed_(seed) { }
        gcv_t(QSRPDE* model) : gcv_t(model, edf_cache_t()) { }
        gcv_t(QSRPDE* model, int r, int seed) : gcv_t(model, edf_cache_t(), r, seed) { }

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
            double pinball = 0;
            for (int i = 0; i < n_; ++i) {
                pinball += model_->pinball_loss(model_->y_[i] - model_->mu_[i], std::pow(10, model_->eps_));
            }
	    return (std::pow(pinball, 2) / std::pow(dor, 2));
        }
        // observers
        const edf_cache_t& edf_cache() const { return edf_cache_; }
        edf_cache_t& edf_cache() { return edf_cache_; }
       private:
        QSRPDE* model_;
        int n_ = 0, q_ = 0;
        edf_cache_t edf_cache_;
        // stochastic edf approximation parameter
        int r_, seed_;
    };
    friend gcv_t;
    gcv_t gcv() { return gcv_t(this); }
    gcv_t gcv(const typename gcv_t::edf_cache_t& edf_cache) { return gcv_t(this, edf_cache); }
    gcv_t gcv(int r, int seed) { return gcv_t(this, r, seed); }
    gcv_t gcv(const typename gcv_t::edf_cache_t& edf_cache, int r, int seed) { return gcv_t(this, edf_cache, r, seed); }

    // inference
  
   private:
    vector_t y_;
    double alpha_ = 0.5;   // quantile order (default to median)
    vector_t py_;          // y - (1 - 2 * alpha) * |y - X * beta - f|
    vector_t pW_;          // diagonal of W^k = 1 / (2 * n * |y - X * beta - f|)
    vector_t mu_;          // \mu^k = [ \mu^k_1, ..., \mu^k_n ] : quantile vector at step k
    double eps_ = -3;      // pinball loss smoothing factor (defaulted to non-smoothed pinball)
    int max_iter_ = 200;   // fpirls maximum iteration number
    double tol_ = 1e-6;    // fprils convergence tolerance
    double tol_weights_ = 1e-6;

    solver_t solver_;
    int n_obs_ = 0, n_covs_ = 0;
    int n_iter_ = 0;

    double pinball_loss(double x, double eps) const {   // quantile check function
        return (alpha_ - 1) * x + eps * fdapde::log1pexp(x / eps);
    };
    // non-smoothed pinball
    double pinball_loss(double x) const { return 0.5 * std::abs(x) + (alpha_ - 0.5) * x; }
};

// deduction guide
template <typename GeoFrame, typename Penalty>
QSRPDE(const std::string& formula, const GeoFrame& gf, double alpha, Penalty&& solver)
  -> QSRPDE<typename Penalty::solver_t>;

}   // namespace fdapde

#endif   // __QUANTILE_SPATIAL_REGRESSION_H__

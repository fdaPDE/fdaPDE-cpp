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
    QSRPDE(const std::string& formula, const GeoFrame& gf, Penalty&& penalty) noexcept :
        solver_() {
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
    }

    // Functional penalized iterative reweighted least squares
    template <typename... LambdaT>
        requires(std::is_convertible_v<LambdaT, double> && ...)
    void fit(LambdaT... lambda) {
        matrix_t y = solver_.response();
        // initialization
        internals::apply_index_pack<n_lambda>([&]<int... Ns_>() { solver_((2 * lambda[Ns_])...); });
        mu_ = fitted();
        double Jold = std::numeric_limits<double>::max(), Jnew = 0;
        n_iter_ = 0;
        while (n_iter_ < max_iter_ && std::abs(Jnew - Jold) > tol_) {
            vector_t abs_res = (y - mu_).array().abs();
            // W_i = 0.5 * (abs_res[i] + tol_weights_) if abs_res[i] < tol_weights, W_i = 0.5 * abs_res[i] otherwise
            pW_ = (abs_res.array() < tol_weights_)
                    .select((2. * (abs_res.array() + tol_weights_)).inverse(), (2. * abs_res.array()).inverse());
            py_ = y - (1. - 2. * alpha_) * abs_res;	  
            // \argmin_{\beta, f} [ 1/n * \norm(W^{1/2} * (y - X * \beta - f_n))^2 + P_{\lambda}(f) ]
	    solver_.update_response_and_weights(py_, pW_.asDiagonal());
            solver_.fit(lambda...);
            mu_ = fitted();
            // prepare for next iteration
            double data_loss = (pW_.cwiseSqrt().matrix().asDiagonal() * (py_ - mu_)).squaredNorm();
            Jold = Jnew;
            Jnew = data_loss + solver_.ftPf()[0];
	    n_iter_++;
        }
	return;
    }
    // observers
    const matrix_t& f() const { return solver_.f(); }
    const matrix_t& beta() const { return solver_.beta(); }
    int n_covs() const { return n_covs_; }
    int n_obs() const { return n_obs_; }
    double edf() { return solver_.edf(); }
    const matrix_t& response() const { return solver_.response(); }
    matrix_t fitted() const {
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
        gcv_t(QSRPDE* model) : model_(model), n_(model->n_obs()), q_(model->n_covs()) { }

        template <typename InputType_>
            requires(internals::is_subscriptable<InputType_, int>)
        constexpr double operator()(const InputType_& lambda) {
            return internals::apply_index_pack<n_lambda>([&]<int... Ns_>() { return operator()(lambda[Ns_]...); });
        }
        template <typename... LambdaT>
            requires(std::is_convertible_v<LambdaT, double> && ...)
        constexpr double operator()(LambdaT... lambda) {
            model_->fit(lambda...);
            int dor = n_ - (q_ + model_->edf());   // residual degrees of freedom
            double pinball = 0;
            for (int i = 0; i < n_; ++i) {
                pinball += model_->pinball_loss(model_->response()[i] - mu_[i], std::pow(10, eps_));
            }
            return std::pow(pinball, 2) / std::pow(dor, 2);
        }
       private:
        QSRPDE* model_;
        int n_ = 0, q_ = 0;
    };
    gcv_t gcv() { return gcv_t(this); }

    // inference
  
   private:
    double alpha_ = 0.5;   // quantile order (default to median)
    vector_t py_;          // y - (1 - 2 * alpha) * |y - X * beta - f|
    vector_t pW_;          // diagonal of W^k = 1 / (2 * n * |y - X * beta - f|)
    vector_t mu_;          // \mu^k = [ \mu^k_1, ..., \mu^k_n ] : quantile vector at step k
    double eps_ = -1e-1;   // pinball loss smoothing factor
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
QSRPDE(const std::string& formula, const GeoFrame& gf, Penalty&& solver) -> QSRPDE<typename Penalty::solver_t>;

}   // namespace fdapde

#endif   // __QUANTILE_SPATIAL_REGRESSION_H__

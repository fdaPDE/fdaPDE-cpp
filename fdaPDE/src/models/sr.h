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
        if constexpr (requires(Penalty p) { p.get(); }) {
            solver_ = solver_t(formula, gf, penalty.get());
        } else {
            solver_ = solver_t(formula, gf, penalty(gf.template triangulation<0>()).get());
        }

        Formula formula_(formula);
	n_obs_  = gf[0].rows();
	n_covs_ = 0;
        for (const std::string& token : formula_.rhs()) {
            if (gf.contains(token)) { n_covs_++; }
        }
    }
    template <typename... LambdaT>
        requires(std::is_convertible_v<LambdaT, double> && ...)
    void fit(LambdaT... lambda) {
        solver_(lambda...);
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
        gcv_t(SRPDE* model) : model_(model), n_(model->n_obs()), q_(model->n_covs()) { }

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
            return (n_ / std::pow(dor, 2)) * (model_->fitted() - model_->response()).squaredNorm();
        }
       private:
        SRPDE* model_;
        int n_ = 0, q_ = 0;
    };
    gcv_t gcv() { return gcv_t(this); }

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

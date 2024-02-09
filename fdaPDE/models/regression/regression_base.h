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

#ifndef __REGRESSION_BASE_H__
#define __REGRESSION_BASE_H__

#include <fdaPDE/utils.h>
#include <fdaPDE/linear_algebra.h>
#include "../model_macros.h"
#include "../model_traits.h"
#include "../space_only_base.h"
#include "../space_time_base.h"
#include "../space_time_separable_base.h"
#include "../space_time_parabolic_base.h"
#include "../sampling_design.h"
#include "gcv.h"
#include "stochastic_edf.h"
using fdapde::core::BinaryVector;

namespace fdapde {
namespace models {

// base class for any *regression* model
template <typename Model, typename RegularizationType>
class RegressionBase :
    public select_regularization_base<Model, RegularizationType>::type,
    public SamplingBase<Model> {
   protected:
    DiagMatrix<double> W_ {};   // diagonal matrix of weights (implements possible heteroscedasticity)
    DMatrix<double> XtWX_ {};   // q x q dense matrix X^\top*W*X
    DMatrix<double> T_ {};      // T = \Psi^\top*Q*\Psi + P (required by GCV)
    Eigen::PartialPivLU<DMatrix<double>> invXtWX_ {};   // factorization of the dense q x q matrix XtWX_.
    // missing data and masking logic
    BinaryVector<fdapde::Dynamic> nan_mask_;     // indicator function over missing observations
    BinaryVector<fdapde::Dynamic> y_mask_;       // discards i-th observation from the fitting if y_mask_[i] == true
    std::size_t n_nan_ = 0;                      // number of missing entries in observation vector
    SpMatrix<double> B_;                         // matrix \Psi corrected for NaN and masked observations

    // matrices required for Woodbury decomposition
    DMatrix<double> U_;   // [\Psi^\top*D*W*X, 0]
    DMatrix<double> V_;   // [X^\top*W*\Psi,   0]

    // room for problem solution
    DVector<double> f_ {};      // estimate of the spatial field (1 x N vector)
    DVector<double> g_ {};      // PDE misfit
    DVector<double> beta_ {};   // estimate of the coefficient vector (1 x q vector)
   public:
    using Base = typename select_regularization_base<Model, RegularizationType>::type;
    using Base::df_;                    // BlockFrame for problem's data storage
    using Base::idx;                    // indices of observations
    using Base::n_basis;                // number of basis function over domain D
    using Base::P;                      // discretized penalty matrix
    using Base::R0;                     // mass matrix
    using SamplingBase<Model>::D;       // for areal sampling, matrix of subdomains measures, identity matrix otherwise
    using SamplingBase<Model>::Psi;     // matrix of spatial basis evaluation at locations p_1 ... p_n
    using SamplingBase<Model>::PsiTD;   // block \Psi^\top*D (not nan-corrected)
    using Base::model;

    RegressionBase() = default;
    // space-only and space-time parabolic constructor (they require only one PDE)
    fdapde_enable_constructor_if(has_single_penalty, Model) RegressionBase(const pde_ptr& pde, Sampling s) :
        Base(pde), SamplingBase<Model>(s) {};
    // space-time separable constructor
    fdapde_enable_constructor_if(has_double_penalty, Model)
      RegressionBase(const pde_ptr& space_penalty, const pde_ptr& time_penalty, Sampling s) :
        Base(space_penalty, time_penalty), SamplingBase<Model>(s) {};

    // getters
    const DMatrix<double>& y() const { return df_.template get<double>(OBSERVATIONS_BLK); }   // observation vector y
    std::size_t q() const {
        return df_.has_block(DESIGN_MATRIX_BLK) ? df_.template get<double>(DESIGN_MATRIX_BLK).cols() : 0;
    }
    const DMatrix<double>& X() const { return df_.template get<double>(DESIGN_MATRIX_BLK); }   // covariates
    const DiagMatrix<double>& W() const { return W_; }                                         // observations' weights
    const DMatrix<double>& XtWX() const { return XtWX_; }
    const Eigen::PartialPivLU<DMatrix<double>>& invXtWX() const { return invXtWX_; }
    const DVector<double>& f() const { return f_; };         // estimate of spatial field
    const DVector<double>& g() const { return g_; };         // PDE misfit
    const DVector<double>& beta() const { return beta_; };   // estimate of regression coefficients
    const BinaryVector<fdapde::Dynamic>& nan_mask() const { return nan_mask_; }
    BinaryVector<fdapde::Dynamic> masked_obs() const { return y_mask_ | nan_mask_; }
    std::size_t n_obs() const { return y().rows() - masked_obs().count(); }   // number of (active) observations
    // getters to Woodbury decomposition matrices
    const DMatrix<double>& U() const { return U_; }
    const DMatrix<double>& V() const { return V_; }
    // access to NaN corrected \Psi and \Psi^\top*D matrices
    const SpMatrix<double>& Psi() const { return !is_empty(B_) ? B_ : Psi(not_nan()); }
    auto PsiTD() const { return !is_empty(B_) ? B_.transpose() * D() : Psi(not_nan()).transpose() * D(); }
    bool has_covariates() const { return q() != 0; }                 // true if the model has a parametric part
    bool has_weights() const { return df_.has_block(WEIGHTS_BLK); }  // true if heteroskedastic observation are provided
    bool has_nan() const { return n_nan_ != 0; }                     // true if there are missing data
    // setters
    void set_mask(const BinaryVector<fdapde::Dynamic>& mask) {
        fdapde_assert(mask.size() == Base::n_locs());
        if (mask.any()) {
            model().runtime().set(runtime_status::require_psi_correction);
            y_mask_ = mask;   // mask[i] == true, removes the contribution of the i-th observation from fit
        }
    }
    // efficient left multiplication by matrix Q = W(I - X*(X^\top*W*X)^{-1}*X^\top*W)
    DMatrix<double> lmbQ(const DMatrix<double>& x) const {
        if (!has_covariates()) return W_ * x;
        DMatrix<double> v = X().transpose() * W_ * x;   // X^\top*W*x
        DMatrix<double> z = invXtWX_.solve(v);          // (X^\top*W*X)^{-1}*X^\top*W*x
        // compute W*x - W*X*z = W*x - (W*X*(X^\top*W*X)^{-1}*X^\top*W)*x = W(I - H)*x = Q*x
        return W_ * x - W_ * X() * z;
    }
    // computes fitted values \hat y = \Psi*f_ + X*beta_
    DMatrix<double> fitted() const {
        fdapde_assert(!is_empty(f_));
        DMatrix<double> hat_y = Psi(not_nan()) * f_;
        if (has_covariates()) hat_y += X() * beta_;
        return hat_y;
    }
    // efficient evaluation of the term \lambda*f^\top*P*f, for a supplied PDE misfit g and solution estimate f
    double ftPf(const SVector<Base::n_lambda>& lambda, const DVector<double>& f, const DVector<double>& g) const {
        fdapde_assert(f.rows() == R0().rows() && g.rows() == R0().rows());
        if constexpr (!is_space_time_separable<Model>::value) {
            return lambda[0] * g.dot(R0() * g);   // \int_D (Lf-u)^2 = g^\top*R_0*g = f^\top*(R_1^\top*(R_0)^{-1}*R_1)*f
        } else {
            return lambda[0] * g.dot(R0() * g) + lambda[1] * f.dot(Base::PT() * f);
        }
    }
    double ftPf(const SVector<Base::n_lambda>& lambda) const {
        if (is_empty(g_)) return f().dot(Base::P(lambda) * f());   // fallback to explicit f^\top*P*f
        return ftPf(lambda, f(), g());
    }
    double ftPf() const { return ftPf(Base::lambda()); }
    // GCV support
    template <typename EDFStrategy_, typename... Args> GCV gcv(Args&&... args) {
        return GCV(Base::model(), EDFStrategy_(std::forward<Args>(args)...));
    }
    const DMatrix<double>& T() {   // T = \Psi^\top*Q*\Psi + P
        T_ = PsiTD() * lmbQ(Psi()) + P();
        return T_;
    }
    // data dependent regression models' initialization logic
    void analyze_data() {
        // initialize empty masks
        if (!y_mask_.size()) y_mask_.resize(Base::n_locs());
        if (!nan_mask_.size()) nan_mask_.resize(Base::n_locs());
        // compute q x q dense matrix X^\top*W*X and its factorization
        if (has_weights() && df_.is_dirty(WEIGHTS_BLK)) {
            W_ = df_.template get<double>(WEIGHTS_BLK).col(0).asDiagonal();
            model().runtime().set(runtime_status::require_W_update);
        } else if (is_empty(W_)) {
            // default to homoskedastic observations
            W_ = DVector<double>::Ones(Base::n_locs()).asDiagonal();
        }
        // compute q x q dense matrix X^\top*W*X and its factorization
        if (has_covariates() && (df_.is_dirty(DESIGN_MATRIX_BLK) || df_.is_dirty(WEIGHTS_BLK))) {
            XtWX_ = X().transpose() * W_ * X();
            invXtWX_ = XtWX_.partialPivLu();
        }
        // derive missingness pattern from observations vector (if changed)
        if (df_.is_dirty(OBSERVATIONS_BLK)) {
            n_nan_ = 0;
            for (std::size_t i = 0; i < df_.template get<double>(OBSERVATIONS_BLK).size(); ++i) {
                if (std::isnan(y()(i, 0))) {   // requires -ffast-math compiler flag to be disabled
                    nan_mask_.set(i);
                    n_nan_++;
                    df_.template get<double>(OBSERVATIONS_BLK)(i, 0) = 0.0;   // zero out NaN
                }
            }
            if (has_nan()) model().runtime().set(runtime_status::require_psi_correction);
        }
        return;
    }
    // correct \Psi setting to zero rows corresponding to masked observations
    void correct_psi() {
        if (masked_obs().any()) B_ = (~masked_obs().blk_repeat(1, n_basis())).select(Psi(not_nan()));
    }
};

}   // namespace models
}   // namespace fdapde

#endif   // __REGRESSION_BASE_H__

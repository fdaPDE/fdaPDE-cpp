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
    DMatrix<double> XtWX_ {};   // q x q dense matrix X^T*W*X
    DMatrix<double> T_ {};      // T = \Psi^T*Q*\Psi + P (required by GCV)
    Eigen::PartialPivLU<DMatrix<double>> invXtWX_ {};   // factorization of the dense q x q matrix XtWX_.
    // missing data and masking logic
    BinaryVector<fdapde::Dynamic> nan_mask_;   // indicator function over missing observations
    BinaryVector<fdapde::Dynamic> y_mask_;     // discards i-th observation from the fitting if y_mask_[i] == true
    std::size_t n_nan_ = 0;                    // number of missing entries in observation vector
    SpMatrix<double> B_;                       // matrix \Psi corrected for NaN and masked observations

    // matrices required for Woodbury decomposition
    DMatrix<double> U_;   // [\Psi^T*D*W*X, 0]
    DMatrix<double> V_;   // [X^T*W*\Psi,   0]

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
    using SamplingBase<Model>::D;       // for areal sampling, matrix of subdomains measures, identity matrix otherwise
    using SamplingBase<Model>::Psi;     // matrix of spatial basis evaluation at locations p_1 ... p_n
    using SamplingBase<Model>::PsiTD;   // block \Psi^T*D (not nan-corrected)
    using Base::model;

    RegressionBase() = default;
    // space-only constructor
    fdapde_enable_constructor_if(is_space_only, Model) RegressionBase(const pde_ptr& pde, Sampling s) :
        Base(pde), SamplingBase<Model>(s) {};
    // space-time constructors
    fdapde_enable_constructor_if(is_space_time_separable, Model)
      RegressionBase(const pde_ptr& space_penalty, const pde_ptr& time_penalty, Sampling s) :
        Base(space_penalty, time_penalty), SamplingBase<Model>(s) {};
    fdapde_enable_constructor_if(is_space_time_parabolic, Model)
      RegressionBase(const pde_ptr& pde, Sampling s, const DVector<double>& time) :
        Base(pde, time), SamplingBase<Model>(s) {};

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
    std::size_t n_obs() const { return y().rows() - n_nan_; }   // number of observations (nan corrected)
    // getters to Woodbury decomposition matrices
    const DMatrix<double>& U() const { return U_; }
    const DMatrix<double>& V() const { return V_; }
    // access to NaN corrected \Psi and \Psi^T*D matrices
    const SpMatrix<double>& Psi() const { return !is_empty(B_) ? B_ : Psi(not_nan()); }
    auto PsiTD() const { return !is_empty(B_) ? B_.transpose() * D() : Psi(not_nan()).transpose() * D(); }

    // utilities
    bool has_covariates() const { return q() != 0; }                 // true if the model has a parametric part
    bool has_weights() const { return df_.has_block(WEIGHTS_BLK); }  // true if heteroscedastic observation are provided
    bool has_nan() const { return n_nan_ != 0; }                     // true if there are missing data

    // an efficient way to perform a left multiplication by Q implementing the following
    //  given the design matrix X, the weight matrix W and x
    //    compute v = X^T*W*x
    //    solve Yz = v
    //    return Wx - WXz = W(I-H)x = Qx
    // it is required to having assigned a design matrix X to the model before calling this method
    DMatrix<double> lmbQ(const DMatrix<double>& x) const {
        if (!has_covariates()) return W_ * x;
        DMatrix<double> v = X().transpose() * W_ * x;   // X^T*W*x
        DMatrix<double> z = invXtWX_.solve(v);          // (X^T*W*X)^{-1}*X^T*W*x
        // compute W*x - W*X*z = W*x - (W*X*(X^T*W*X)^{-1}*X^T*W)*x = W(I - H)*x = Q*x
        return W_ * x - W_ * X() * z;
    }
    DMatrix<double> fitted() const {   // computes fitted values \hat y = \Psi*f_ + X*beta_
        DMatrix<double> hat_y = Psi(not_nan()) * f_;
        if (has_covariates()) hat_y += X() * beta_;
        return hat_y;
    }
    // GCV support
    template <typename EDFStrategy_, typename... Args> GCV gcv(Args&&... args) {
        return GCV(Base::model(), EDFStrategy_(std::forward<Args>(args)...));
    }
    const DMatrix<double>& T() {   // T = \Psi^T*Q*\Psi + P
        T_ = PsiTD() * lmbQ(Psi()) + P();
        return T_;
    }
    // data dependent regression models' initialization logic
    void analyze_data() {
        // compute q x q dense matrix X^T*W*X and its factorization
        if (has_weights() && df_.is_dirty(WEIGHTS_BLK)) {
            W_ = df_.template get<double>(WEIGHTS_BLK).col(0).asDiagonal();
            model().runtime().set(runtime_status::require_W_update);
        } else if (is_empty(W_)) {
            // default to homoskedastic observations
            W_ = DVector<double>::Ones(Base::n_locs()).asDiagonal();
        }
        // compute q x q dense matrix X^T*W*X and its factorization if covariates are supplied
        if (has_covariates() && (df_.is_dirty(DESIGN_MATRIX_BLK) || df_.is_dirty(WEIGHTS_BLK))) {
            XtWX_ = X().transpose() * W_ * X();
            invXtWX_ = XtWX_.partialPivLu();
        }
        // derive missingness pattern from observations vector (if changed)
        if (df_.is_dirty(OBSERVATIONS_BLK)) {
            nan_mask_.resize(y().rows());
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

    // assembles matrix B (sets to zero \Psi rows corresponding to missing or masked observations)
    void correct_psi() {
        BinaryVector<fdapde::Dynamic> mask(y().rows());
        if (has_nan()) mask = mask | nan_mask_;
        if (y_mask_.size() != 0) mask = mask | y_mask_;
        // assemble matrix B_ according to given mask
        if (mask.any()) { B_ = mask_matrix(Psi(not_nan()), mask); }
        return;
    }
    // if mask[i] is true, removes the contribution of the i-th observation from the model's fitting
    void mask_obs(const BinaryVector<fdapde::Dynamic>& y_mask) {
        fdapde_assert(y_mask.size() == y().rows());
        model().runtime().set(runtime_status::require_psi_correction);
        y_mask_ = y_mask;
	return;
    }
};

}   // namespace models
}   // namespace fdapde

#endif   // __REGRESSION_BASE_H__

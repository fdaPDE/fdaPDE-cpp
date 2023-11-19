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
#include "../model_macros.h"
#include "../model_traits.h"
#include "../space_only_base.h"
#include "../space_time_base.h"
#include "../space_time_separable_base.h"
#include "../space_time_parabolic_base.h"
#include "../sampling_design.h"
#include "gcv.h"
#include "stochastic_edf.h"

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
    std::unordered_set<std::size_t> nan_idxs_;   // indexes of missing observations
    SpMatrix<double> B_;                         // matrix \Psi corrected for NaN observations

    // matrices required for Woodbury decomposition
    DMatrix<double> U_;   // [\Psi^T*D*W*X, 0]
    DMatrix<double> V_;   // [X^T*W*\Psi,   0]

    // room for problem solution
    DVector<double> f_ {};      // estimate of the spatial field (1 x N vector)
    DVector<double> g_ {};      // PDE misfit
    DVector<double> beta_ {};   // estimate of the coefficient vector (1 x q vector)
   public:
    using Base = typename select_regularization_base<Model, RegularizationType>::type;
    using Base::df_;           // BlockFrame for problem's data storage
    using Base::idx;           // indices of observations
    using Base::n_basis;       // number of basis function over domain D
    using Base::pde_;          // differential operator L
    using Base::P;             // discretized penalty matrix
    using SamplingBase<Model>::D;     // for areal sampling, matrix of subdomains measures, identity matrix otherwise
    using SamplingBase<Model>::Psi;   // matrix of spatial basis evaluation at locations p_1 ... p_n
    using Base::model;
  
    RegressionBase() = default;
    // space-only constructor
    template <
      typename U = Model,   // fake type to enable substitution in SFINAE
      typename std::enable_if<std::is_same<typename U::RegularizationType, SpaceOnly>::value, int>::type = 0>
    RegressionBase(const pde_ptr& pde, Sampling s) : Base(pde), SamplingBase<Model>(s) {};
    // space-time constructor
    template <
      typename U = Model,   // fake type to enable substitution in SFINAE
      typename std::enable_if<!std::is_same<typename U::RegularizationType, SpaceOnly>::value, int>::type = 0>
    RegressionBase(const pde_ptr& pde, Sampling s, const DVector<double>& time) :
        Base(pde, time), SamplingBase<Model>(s) {};
    // copy constructor, copy only pde object (as a consequence also the problem domain)
    RegressionBase(const RegressionBase& rhs) : Base(rhs), SamplingBase<Model>(rhs) { pde_ = rhs.pde_; }

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
    const std::unordered_set<std::size_t>& nan_idxs() const { return nan_idxs_; }   // missing data indexes
    std::size_t n_obs() const { return y().rows(); }   // number of observations
    // getters to Woodbury decomposition matrices
    const DMatrix<double>& U() const { return U_; }
    const DMatrix<double>& V() const { return V_; }
    // access to NaN corrected \Psi and \Psi^T*D matrices
    const SpMatrix<double>& Psi() const { return has_nan() ? B_ : Psi(not_nan()); }
    auto PsiTD() const { return has_nan() ? B_.transpose() * D() : Psi(not_nan()).transpose() * D(); }
  
    // utilities
    bool has_covariates() const { return q() != 0; }                 // true if the model has a parametric part
    bool has_weights() const { return df_.has_block(WEIGHTS_BLK); }  // true if heteroscedastic observation are assumed
    bool has_nan() const { return nan_idxs_.size() != 0; }           // true if there are missing data

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
    template <template <typename> typename edf_evaluation_strategy, typename... Args>
    GCV<Model, edf_evaluation_strategy> gcv(Args&&... args) {
      return GCV<Model, edf_evaluation_strategy>(Base::model(), std::forward<Args>(args)...);
    }
    const DMatrix<double>& T() {   // T = \Psi^T*Q*\Psi + P
        T_ = PsiTD() * lmbQ(Psi()) + P();
        return T_;
    }

    void init_data() {
        if (has_weights() && df_.is_dirty(WEIGHTS_BLK)) {   // update observations' weights if provided
            W_ = df_.template get<double>(WEIGHTS_BLK).col(0).asDiagonal();
	    model().runtime().set(runtime_status::require_W_update);
        } else if (is_empty(W_)) {
            // default to homoskedastic observations
            W_ = DVector<double>::Ones(Base::n_locs()).asDiagonal();
        }
        if (has_covariates() && (df_.is_dirty(DESIGN_MATRIX_BLK) || df_.is_dirty(WEIGHTS_BLK))) {
            // compute q x q dense matrix X^T*W*X and its factorization
            XtWX_ = X().transpose() * W_ * X();
            invXtWX_ = XtWX_.partialPivLu();
        }
	// clear dirty bits
	df_.clear_dirty_bit(WEIGHTS_BLK);
	df_.clear_dirty_bit(DESIGN_MATRIX_BLK);
        return;
    }
    void init_nan() {   // regression models' missing data logic (called by SamplingBase::init_sampling())
        // derive missingness pattern
        nan_idxs_.clear();   // empty nan indexes set
        for (std::size_t i = 0; i < df_.template get<double>(OBSERVATIONS_BLK).size(); ++i) {
            if (std::isnan(y()(i, 0))) {   // requires -ffast-math compiler flag to be disabled
                nan_idxs_.insert(i);
                df_.template get<double>(OBSERVATIONS_BLK)(i, 0) = 0.0;   // zero out NaN
            }
        }
        // matrix B assembly logic (set to zero rows corresponding to missing observations)
        if (has_nan()) {
            // reserve space
            std::size_t n = Psi(not_nan()).rows();
            std::size_t N = Psi(not_nan()).cols();
            B_.resize(n, N);
            // triplet list to fill sparse matrix
            std::vector<fdapde::Triplet<double>> triplet_list;
            triplet_list.reserve(n * N);
            for (int k = 0; k < Psi(not_nan()).outerSize(); ++k)
                for (SpMatrix<double>::InnerIterator it(Psi(not_nan()), k); it; ++it) {
                    if (nan_idxs_.find(it.row()) == nan_idxs_.end()) {
                        // no missing data at this location
                        triplet_list.emplace_back(it.row(), it.col(), it.value());
                    }
                }
            // finalize construction
            B_.setFromTriplets(triplet_list.begin(), triplet_list.end());
            B_.makeCompressed();
        }
        return;
    }
};

}   // namespace models
}   // namespace fdapde

#endif   // __REGRESSION_BASE_H__

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

namespace fdapde {
namespace models {

// base class for any *regression* fdaPDE model
template <typename Model>
class RegressionBase :
    public select_regularization_type<Model>::type,
    public SamplingDesign<Model, typename model_traits<Model>::sampling> {
   protected:
    DiagMatrix<double> W_ {};   // diagonal matrix of weights (implements possible heteroscedasticity)
    DMatrix<double> XtWX_ {};   // q x q dense matrix X^T*W*X
    Eigen::PartialPivLU<DMatrix<double>> invXtWX_ {};   // factorization of the dense q x q matrix XtWX_.

    // matrices required for Woodbury decomposition
    DMatrix<double> U_;   // [\Psi^T*D*W*X, 0]
    DMatrix<double> V_;   // [X^T*W*\Psi,   0]

    // quantites related to missing data setting
    std::unordered_set<std::size_t> nan_idxs_;   // indexes of missing observations
    SpMatrix<double> B_;   // matrix \Psi where rows corresponding to NaN observations are zeroed

    // room for problem solution
    DVector<double> f_ {};      // estimate of the spatial field (1 x N vector)
    DVector<double> g_ {};      // PDE misfit
    DVector<double> beta_ {};   // estimate of the coefficient vector (1 x q vector)
   public:
    typedef typename model_traits<Model>::PDE PDE;   // PDE used for regularization in space
    typedef typename select_regularization_type<Model>::type Base;
    typedef SamplingDesign<Model, typename model_traits<Model>::sampling> SamplingBase;
    using Base::df_;           // BlockFrame for problem's data storage
    using Base::idx;           // indices of observations
    using Base::n_basis;       // number of basis function over domain D
    using Base::pde_;          // differential operator L
    using SamplingBase::D;     // matrix of subdomains measures (for areal sampling)
    using SamplingBase::Psi;   // matrix of spatial basis evaluation at locations p_1 ... p_n

    RegressionBase() = default;
    // space-only constructor
    template <
      typename U = Model,   // fake type to enable substitution in SFINAE
      typename std::enable_if<std::is_same<typename model_traits<U>::regularization, SpaceOnly>::value, int>::type = 0>
    RegressionBase(const PDE& pde) :
        select_regularization_type<Model>::type(pde),
        SamplingDesign<Model, typename model_traits<Model>::sampling>() {};
    // space-time constructor
    template <
      typename U = Model,   // fake type to enable substitution in SFINAE
      typename std::enable_if<!std::is_same<typename model_traits<U>::regularization, SpaceOnly>::value, int>::type = 0>
    RegressionBase(const PDE& pde, const DVector<double>& time) :
        select_regularization_type<Model>::type(pde, time),
        SamplingDesign<Model, typename model_traits<Model>::sampling>() {};
    // copy constructor, copy only pde object (as a consequence also the problem domain)
    RegressionBase(const RegressionBase& rhs) { pde_ = rhs.pde_; }

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
    std::size_t n_obs() const { return y().rows() - nan_idxs_.size(); }   // number of observations
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
    DMatrix<double> lmbQ(const DMatrix<double>& x) const;            // efficient multiplication by matrix Q
    DMatrix<double> fitted() const;   // computes fitted values \hat y = \Psi*f_ + X*beta_

    // initialization methods
    void update_data();   // update model's status to data (called by ModelBase::setData())
    void init_nan();      // regression models' missing data logic (called by SamplingBase::init_sampling())
    void init_nan(const std::unordered_set<std::size_t>& nan_idxs_input);
};

// implementative details

// an efficient way to perform a left multiplication by Q implementing the following
//  given the design matrix X, the weight matrix W and x
//    compute v = X^T*W*x
//    solve Yz = v
//    return Wx - WXz = W(I-H)x = Qx
// it is required to having assigned a design matrix X to the model before calling this method
template <typename Model> DMatrix<double> RegressionBase<Model>::lmbQ(const DMatrix<double>& x) const {
    if (!has_covariates()) return W_ * x;
    DMatrix<double> v = X().transpose() * W_ * x;   // X^T*W*x
    DMatrix<double> z = invXtWX_.solve(v);          // (X^T*W*X)^{-1}*X^T*W*x
    // compute W*x - W*X*z = W*x - (W*X*(X^T*W*X)^{-1}*X^T*W)*x = W(I - H)*x = Q*x
    return W_ * x - W_ * X() * z;
}

// initialization stuffs depending on the data
template <typename Model> void RegressionBase<Model>::update_data() {
    // default to homoskedastic observations
    DVector<double> W = DVector<double>::Ones(Base::n_locs());
    if (has_weights()) {   // update observations' weights if provided
        W = df_.template get<double>(WEIGHTS_BLK).col(0);
    }
    W_ = W.asDiagonal();
    // model is semi-parametric
    if (has_covariates()) {
        // compute q x q dense matrix X^T*W*X and its factorization
        XtWX_ = X().transpose() * W_ * X();
        invXtWX_ = XtWX_.partialPivLu();
    }
}

// computes fitted values \hat y = \Psi*f_ + X*beta_
template <typename Model> DMatrix<double> RegressionBase<Model>::fitted() const {
    DMatrix<double> hat_y = Psi(not_nan()) * f_;
    if (has_covariates()) hat_y += X() * beta_;
    return hat_y;
}

// missing data logic
template <typename Model> void RegressionBase<Model>::init_nan() {
    // derive missingness pattern
    nan_idxs_.clear();   // empty nan indexes set
    for (std::size_t i = 0; i < df_.template get<double>(OBSERVATIONS_BLK).size(); ++i) { 
        if(std::isnan(y()(i, 0))) {   // requires -ffast-math compiler flag to be disabled
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

// missing data logic correction used in fpirls
template <typename Model> void RegressionBase<Model>::init_nan(const std::unordered_set<std::size_t>& nan_idxs_input) {
    // derive missingness pattern
    nan_idxs_.clear();   // empty nan indexes set
    nan_idxs_ = nan_idxs_input ; 
        
    // matrix B assembly logic (set to zero rows corresponding to missing observations)
    if (has_nan()) {
        // reserve space
        std::size_t n = Psi(not_nan()).rows();
        std::size_t N = Psi(not_nan()).cols();
        B_.resize(n, N);
        // triplet list to fill sparse matrix
        std::vector<fdapde::Triplet<double>> triplet_list;
        std::cout << "regression base init nan pt1" << std::endl;
        std::cout << "regression base init nan: ndof = " << model_traits<Model>::PDE::SolverType::n_dof_per_element << std::endl;
        std::cout << "regression base init nan: n*ndof = " << n * model_traits<Model>::PDE::SolverType::n_dof_per_element << std::endl;
        triplet_list.reserve(n * model_traits<Model>::PDE::SolverType::n_dof_per_element);    // M  --> prima: n*N
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


// trait to detect if a type is a regression model
template <typename T> struct is_regression_model {
    static constexpr bool value = fdapde::is_base_of_template<RegressionBase, T>::value;
};

}   // namespace models
}   // namespace fdapde

#endif   // __REGRESSION_BASE_H__

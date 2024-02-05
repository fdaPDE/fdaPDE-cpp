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

#ifndef __FPLS_H__
#define __FPLS_H__

#include <fdaPDE/utils.h>
#include <Eigen/SVD>

#include "../../calibration/calibration_base.h"
#include "../../calibration/off.h"
#include "../../calibration/gcv.h"
using fdapde::calibration::Calibrator;
#include "functional_base.h"
#include "regularized_svd.h"
#include "center.h"

namespace fdapde {
namespace models {

// FPLS (Functional Partial Least Square regression) model signature
template <typename RegularizationType_>
class FPLS : public FunctionalBase<FPLS<RegularizationType_>, RegularizationType_> {
   public:
    using RegularizationType = std::decay_t<RegularizationType_>;
    using This = FPLS<RegularizationType>;
    using Base = FunctionalBase<This, RegularizationType>;
    using SmootherType = std::conditional_t<is_space_only<This>::value, SRPDE, STRPDE<RegularizationType, monolithic>>;
    IMPORT_MODEL_SYMBOLS;
    using Base::df_;
    using Base::n_basis;
    using Base::n_obs;
    using Base::n_stat_units;
    using Base::X;   // n_stat_units \times n_locs data matrix

    // constructors
    FPLS() = default;
    fdapde_enable_constructor_if(is_space_only, This)
      FPLS(const pde_ptr& pde, Sampling s, RegularizedSVD<sequential> rsvd) :
        Base(pde, s), rsvd_(rsvd) {};
    fdapde_enable_constructor_if(is_space_time_separable, This)
      FPLS(const pde_ptr& space_penalty, const pde_ptr& time_penalty, Sampling s, RegularizedSVD<sequential> rsvd) :
        Base(space_penalty, time_penalty, s), rsvd_(rsvd) {};

    void init_model() {
        // initialize smoothing solver for regression step
        if constexpr (is_space_only<SmootherType>::value) { smoother_ = SmootherType(Base::pde(), Base::sampling()); }
	else {
            smoother_ = SmootherType(Base::pde(), Base::time_pde(), Base::sampling());
            smoother_.set_temporal_locations(Base::time_locs());
        }
        smoother_.set_spatial_locations(Base::locs());
        if (!calibrator_) {   // smoothing solver's calibration strategy fallback
            if (rsvd_.calibration() == Calibration::off) {
                calibrator_ = calibration::Off {Base::lambda()};
            } else {
                calibrator_ = calibration::GCV {core::Grid<Dynamic> {}, StochasticEDF(100)}(rsvd_.lambda_grid());
	    }
        }
        return;
    }
    void solve() {
        // allocate space
        W_.resize(n_basis(), n_comp_);        // optimal direction in X space
        C_.resize(n_basis(), n_comp_);        // optimal X loadings
        V_.resize(Y().cols(), n_comp_);       // optimal direction in Y space
        D_.resize(Y().cols(), n_comp_);       // optimal Y loadings
        T_.resize(n_stat_units(), n_comp_);   // X latent component

        // copy original data to avoid side effects
        DMatrix<double> X_h = X(), Y_h = Y();

        for (std::size_t h = 0; h < n_comp_; ++h) {
            // correlation maximization
            // solves \argmin_{v,w} \norm_F{Y_h^\top*X_h - v^\top*w}^2 + (v^\top*v)*P_{\lambda}(w)
	  rsvd_.compute(Y_h.transpose() * X_h, *this, 1);
            W_.col(h) = rsvd_.loadings();
            V_.col(h) = rsvd_.scores() / rsvd_.loadings_norm()[0];
            T_.col(h) = X_h * Psi() * W_.col(h);   // X latent component

            // regression: solves \argmin_{c} \norm_F{X_h - t*c^\top}^2 + P_{\lambda}(c)
            C_.col(h) = smooth_mean(X_h, T_.col(h), smoother_, calibrator_);
            D_.col(h) = Y_h.transpose() * T_.col(h) / T_.col(h).squaredNorm();

            // deflation
            X_h -= T_.col(h) * (Psi() * C_.col(h)).transpose();
            Y_h -= T_.col(h) * D_.col(h).transpose();
        }

        B_ = W_ * (C_.transpose() * Psi().transpose() * Psi() * W_).partialPivLu().solve(D_.transpose());
        return;
    }

    // getters
    const DMatrix<double>& Y() const { return df_.template get<double>(OBSERVATIONS_BLK); }
    const DMatrix<double>& X() const { return df_.template get<double>(DESIGN_MATRIX_BLK); }
    const DMatrix<double>& X_space_directions() const { return W_; }
    const DMatrix<double>& Y_space_directions() const { return V_; }
    const DMatrix<double>& X_latent() const { return T_; }
    const DMatrix<double>& X_loadings() const { return C_; }
    const DMatrix<double>& Y_loadings() const { return D_; }
    DMatrix<double> fitted() const { return X_latent() * Y_loadings().transpose(); }
    DMatrix<double> reconstructed() const { return X_latent() * X_loadings().transpose(); }
    const DMatrix<double>& B() const { return B_; }
    // setters
    void set_ncomp(std::size_t n_comp) { n_comp_ = n_comp; }
    void set_rsvd(const RegularizedSVD<sequential>& rsvd) { rsvd_ = rsvd; }
    template <typename CalibratorType_> void set_smoothing_step_calibrator(CalibratorType_&& calibrator) {
        calibrator_ = calibrator;
    }
   private:
    SmootherType smoother_;                 // smoothing algorithm used in regression step
    Calibrator<SmootherType> calibrator_;   // calibration strategy used in regression step
    RegularizedSVD<sequential> rsvd_;       // RSVD solver employed in correlation maximization step
    std::size_t n_comp_ = 3;                // number of latent components

    // problem solution
    DMatrix<double> W_;   // optimal directions in X space
    DMatrix<double> V_;   // optimal directions in Y space
    DMatrix<double> T_;   // latent components
    DMatrix<double> C_;   // optimal X loadings
    DMatrix<double> D_;   // optimal Y loadings
    DMatrix<double> B_;
};

}   // namespace models
}   // namespace fdapde

#endif   // __FPLS_H__

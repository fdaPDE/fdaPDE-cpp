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

#ifndef __GSRPDE_H__
#define __GSRPDE_H__

#include <fdaPDE/pde.h>
#include <fdaPDE/utils.h>

#include "../model_macros.h"
#include "regression_base.h"
#include "fpirls.h"

namespace fdapde {
namespace models {
  
// base class for GSRPDE model
template <typename RegularizationType_>
class GSRPDE : public RegressionBase<GSRPDE<RegularizationType_>, RegularizationType_> {
   public:
    using RegularizationType = std::decay_t<RegularizationType_>;
    using This = GSRPDE<RegularizationType>;
    using Base = RegressionBase<GSRPDE<RegularizationType>, RegularizationType>;
    // import commonly defined symbols from base
    IMPORT_REGRESSION_SYMBOLS;
    using Base::invXtWX_;   // LU factorization of X^T*W*X
    using Base::lambda_D;   // smoothing parameter in space
    using Base::P;          // discretized penalty
    using Base::W_;         // weight matrix
    using Base::XtWX_;      // q x q matrix X^T*W*X
    // constructor
    GSRPDE() = default;
    // space-only and space-time parabolic constructor
    fdapde_enable_constructor_if(has_single_penalty, This)
      GSRPDE(const pde_ptr& pde, Sampling s, const Distribution& distr) :
        Base(pde, s), distr_(distr) {
        fpirls_ = FPIRLS<This>(this, tol_, max_iter_);
    };
    // space-time separable constructor
    fdapde_enable_constructor_if(has_double_penalty, This)
      GSRPDE(const pde_ptr& space_penalty, const pde_ptr& time_penalty, Sampling s, const Distribution& distr) :
        Base(space_penalty, time_penalty, s), distr_(distr) {
        fpirls_ = FPIRLS<This>(this, tol_, max_iter_);
    };

    // setters
    void set_fpirls_tolerance(double tol) { tol_ = tol; }
    void set_fpirls_max_iter(std::size_t max_iter) { max_iter_ = max_iter; }
    void init_model() { fpirls_.init(); };
    void solve() {
        fdapde_assert(y().rows() != 0);
        // execute FPIRLS for minimization of functional \norm{V^{-1/2}(y - \mu)}^2 + \lambda \int_D (Lf - u)^2
        fpirls_.compute();
        // fpirls_ converged: extract matrix W and solution estimates
        W_ = fpirls_.solver().W();
        f_ = fpirls_.solver().f();
        g_ = fpirls_.solver().g();
        // store parametric part
        if (has_covariates()) {
            beta_ = fpirls_.solver().beta();
            XtWX_ = fpirls_.solver().XtWX();
            invXtWX_ = fpirls_.solver().invXtWX();
            U_ = fpirls_.solver().U();
            V_ = fpirls_.solver().V();
        }
        invA_ = fpirls_.solver().invA();
        return;
    }
    // required by FPIRLS (see fpirls.h for details)
    // initalizes mean vector \mu
    void fpirls_init() { mu_ = distr_.preprocess(y()); };
    // computes W^k = ((G^k)^{-2})*((V^k)^{-1}) and y^k = G^k(y-u^k) + \theta^k
    void fpirls_compute_step() {
        DVector<double> theta_ = distr_.link(mu_);   // \theta^k = (g(\mu^k_1), ..., g(\mu^k_n))
        DVector<double> G_ = distr_.der_link(mu_);   // G^k = diag(g'(\mu^k_1), ..., g'(\mu^k_n))
        DVector<double> V_ = distr_.variance(mu_);   // V^k = diag(v(\mu^k_1), ..., v(\mu^k_n))
        pW_ = ((G_.array().pow(2) * V_.array()).inverse()).matrix();
        py_ = G_.asDiagonal() * (y() - mu_) + theta_;
    }
    // updates mean vector \mu after WLS solution
    void fpirls_update_step(const DMatrix<double>& hat_f, const DMatrix<double>& hat_beta) {
        mu_ = distr_.inv_link(hat_f);
    }
    // returns the data loss \norm{V^{-1/2}(y - \mu)}^2
    double data_loss() const {
        DVector<double> V = distr_.variance(mu_).array().sqrt().inverse().matrix();
        return (V.asDiagonal() * (y() - mu_)).squaredNorm();
    }
    const DVector<double>& py() const { return py_; }
    const DVector<double>& pW() const { return pW_; }
    const fdapde::SparseLU<SpMatrix<double>>& invA() const { return invA_; }
    // GCV support
    double norm(const DMatrix<double>& op1, const DMatrix<double>& op2) const {   // total deviance \sum dev(\hat y - y)
        DMatrix<double> mu = distr_.inv_link(op1);
        double result = 0;
        for (std::size_t i = 0; i < n_locs(); ++i) {
            if (!Base::masked_obs()[i]) result += distr_.deviance(mu.coeff(i, 0), op2.coeff(i, 0));
        }
        return result;
    }
   private:
    Distribution distr_ {};
    DVector<double> py_;   // \tilde y^k = G^k(y-u^k) + \theta^k
    DVector<double> pW_;   // diagonal of W^k = ((G^k)^{-2})*((V^k)^{-1})
    DVector<double> mu_;   // \mu^k = [ \mu^k_1, ..., \mu^k_n ] : mean vector at step k
    fdapde::SparseLU<SpMatrix<double>> invA_;

    // FPIRLS parameters (set to default)
    FPIRLS<This> fpirls_;
    std::size_t max_iter_ = 200;
    double tol_ = 1e-4;
};
  
}   // namespace models
}   // namespace fdapde

#endif   // __GSRPDE_H__

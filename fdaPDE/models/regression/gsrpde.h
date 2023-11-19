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
   private:
    Distribution distr_ {};
    DVector<double> py_;   // \tilde y^k = G^k(y-u^k) + \theta^k
    DVector<double> pW_;   // diagonal of W^k = ((G^k)^{-2})*((V^k)^{-1})
    DVector<double> mu_;   // \mu^k = [ \mu^k_1, ..., \mu^k_n ] : mean vector at step k
    DMatrix<double> T_;    // T = \Psi^T*Q*\Psi + P

    // FPIRLS parameters (set to default)
    std::size_t max_iter_ = 200;
    double tol_ = 1e-4;
   public:
    using RegularizationType = RegularizationType_;
    using Base = RegressionBase<GSRPDE<RegularizationType>, RegularizationType>;
    // import commonly defined symbols from base
    IMPORT_REGRESSION_SYMBOLS;
    using Base::lambda_D;   // smoothing parameter in space
    using Base::P;          // discretized penalty
    using Base::W_;         // weight matrix
    // constructor
    GSRPDE() = default;
    // space-only constructor
    template <
      typename U = RegularizationType,
      typename std::enable_if<std::is_same<U, SpaceOnly>::value, int>::type = 0>
    GSRPDE(const pde_ptr& pde, Sampling s, const Distribution& distr) : Base(pde, s), distr_(distr) {};
    // space-time constructor
    template <
      typename U = RegularizationType,
      typename std::enable_if<!std::is_same<U, SpaceOnly>::value, int>::type = 0>
    GSRPDE(const pde_ptr& pde, const DVector<double>& time, Sampling s, const Distribution& distr) :
        Base(pde, s, time), distr_(distr) {};

    // setters
    void set_fpirls_tolerance(double tol) { tol_ = tol; }
    void set_fpirls_max_iter(std::size_t max_iter) { max_iter_ = max_iter; }

    void init_data()  { return; };
    void init_model() { return; };
    void solve() {   // finds a solution to the smoothing problem
        // execute FPIRLS for minimization of functional \norm{V^{-1/2}(y - \mu)}^2 + \lambda \int_D (Lf - u)^2
        FPIRLS<decltype(*this)> fpirls(*this, tol_, max_iter_);   // FPIRLS engine
        fpirls.compute();

        // fpirls converged: extract matrix W and solution estimates
        W_ = fpirls.solver().W();
        f_ = fpirls.solver().f();
        if (has_covariates()) { beta_ = fpirls.solver().beta(); }
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
    double norm(const DMatrix<double>& op1, const DMatrix<double>& op2) const {   // total deviance \sum dev(\hat y - y)
        DMatrix<double> mu = distr_.inv_link(op1);
        double result = 0;
        for (std::size_t i = 0; i < op2.rows(); ++i) { result += distr_.deviance(mu.coeff(i, 0), op2.coeff(i, 0)); }
        return result;
    }

    virtual ~GSRPDE() = default;
};
  
}   // namespace models
}   // namespace fdapde

#endif   // __GSRPDE_H__

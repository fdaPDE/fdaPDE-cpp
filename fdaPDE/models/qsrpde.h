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

#ifndef __QSRPDE_H__
#define __QSRPDE_H__

#include <fdaPDE/pde.h>
#include <fdaPDE/utils.h>

#include "../model_macros.h"
#include "fpirls.h"
#include "regression_base.h"

namespace fdapde {
namespace models {

template <typename RegularizationType_>
class QSRPDE : public RegressionBase<QSRPDE<RegularizationType_>, RegularizationType_> {
   public:
    using RegularizationType = std::decay_t<RegularizationType_>;
    using This = QSRPDE<RegularizationType>;
    using Base = RegressionBase<QSRPDE<RegularizationType>, RegularizationType>;
    // import commonly defined symbols from base
    IMPORT_REGRESSION_SYMBOLS;
    using Base::invXtWX_;   // LU factorization of X^T*W*X
    using Base::lambda_D;   // smoothing parameter in space
    using Base::n_basis;    // number of spatial basis functions
    using Base::P;          // discretized penalty matrix
    using Base::W_;         // weight matrix
    using Base::XtWX_;      // q x q matrix X^T*W*X


    // constructor
    QSRPDE() = default;
    // space-only and space-time parabolic constructor
    fdapde_enable_constructor_if(has_single_penalty, This)
    QSRPDE(const pde_ptr& pde, Sampling s, double alpha = 0.5) : Base(pde, s), alpha_(alpha) {
        fpirls_ = FPIRLS<This>(this, tol_, max_iter_);
    };
    // space-time separable constructor
    fdapde_enable_constructor_if(has_double_penalty, This)
    QSRPDE(const pde_ptr& space_penalty, const pde_ptr& time_penalty, Sampling s, double alpha = 0.5) :
        Base(space_penalty, time_penalty, s), alpha_(alpha) {
        fpirls_ = FPIRLS<This>(this, tol_, max_iter_);
    };
    // setter
    void set_fpirls_tolerance(double tol) { tol_ = tol; }
    void set_fpirls_max_iter(std::size_t max_iter) { max_iter_ = max_iter; }
    void set_alpha(double alpha) { alpha_ = alpha; }
    void set_eps_power(double eps) { eps_ = eps; }
    void set_weights_tolerance(double tol_weights) { tol_weights_ = tol_weights; }

    void init_model() { 
        fpirls_.init();
    }
    void solve() {

        // std::cout << "qsrpde solve: n tot " << y().size() << std::endl;
        // std::cout << "qsrpde solve: n " << n_obs() << std::endl;
        // std::cout << "qsrpde solve: n missing " << Base::masked_obs().count() << std::endl;

        // execute FPIRLS_ for minimization of functional \norm{V^{-1/2}(y - \mu)}^2 + \lambda \int_D (Lf - u)^2
        fpirls_.compute();
        // fpirls_ converged: store solution estimates
        W_ = fpirls_.solver().W();
        f_ = fpirls_.solver().f();
        g_ = fpirls_.solver().g();
        // parametric part, if problem was semi-parametric
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

    // required by FPIRLS_ (see fpirls_.h for details)
    // initalizes mean vector \mu
    void fpirls_init() {
        // TODO: use a standardized solver (we can exploit the fpirls solver)
        // non-parametric and semi-parametric cases coincide here, since beta^(0) = 0
        SparseBlockMatrix<double, 2, 2> A(
          -PsiTD() * Psi() / n_obs(), 2 * lambda_D() * R1().transpose(),   // NB: observe the 2 * here
          lambda_D() * R1(),          lambda_D() * R0()                );
        if constexpr (is_space_time_separable<This>::value) {
            A.block(0, 0) -= 2*Base::lambda_T() * Kronecker(Base::P1(), Base::pde().mass());
        }
        fdapde::SparseLU<SpMatrix<double>> invA;
        invA.compute(A);
        // assemble rhs of srpde problem
        DVector<double> b(A.rows());
        b.block(n_basis(), 0, n_basis(), 1) = lambda_D() * u();
        b.block(0, 0, n_basis(), 1) = -PsiTD() * y() / n_obs();
        mu_ = Psi(not_nan()) * (invA.solve(b)).head(n_basis());
    }
    // computes W^k = diag(1/(2*n*|y - X*beta - f|)) and y^k = y - (1-2*alpha)|y - X*beta - f|
    void fpirls_compute_step() {
        DVector<double> abs_res = (y() - mu_).array().abs();
        // W_i = 1/(2*n*(abs_res[i] + tol_weights_)) if abs_res[i] < tol_weights, W_i = 1/(2*n*abs_res[i]) otherwise
        pW_ =
          (abs_res.array() < tol_weights_)
            .select(
              (2 * n_obs() * (abs_res.array() + tol_weights_)).inverse(), (2 * n_obs() * abs_res.array()).inverse());
        py_ = y() - (1 - 2. * alpha_) * abs_res;


        // ATT: riaggiunto  
        for(std::size_t i=0; i<n_locs(); ++i){
            if(Base::nan_mask()[i]){
                py_(i)=0.; 
                pW_(i)=0.; 
            }
        }

    }
    // updates mean vector \mu after WLS solution
    void fpirls_update_step(const DMatrix<double>& hat_f, const DMatrix<double>& hat_beta) { mu_ = hat_f; }
    // returns the data loss \norm{diag(W)^{-1/2}(y - \mu)}^2
    double data_loss() const { return (pW_.cwiseSqrt().matrix().asDiagonal() * (py_ - mu_)).squaredNorm(); }
    const DVector<double>& py() const { return py_; }
    const DVector<double>& pW() const { return pW_; }
    const fdapde::SparseLU<SpMatrix<double>>& invA() const { return invA_; }
    // GCV support
    double norm(const DMatrix<double>& op1, const DMatrix<double>& op2) const {

        double result = 0;
        for (std::size_t i = 0; i < n_locs(); ++i) {
            if (!Base::masked_obs()[i]) result += pinball_loss(op2.coeff(i, 0) - op1.coeff(i, 0), std::pow(10, eps_));
        }
        return std::pow(result, 2) / n_obs();
    }
   private:
    double alpha_ = 0.5;      // quantile order (default to median)
    DVector<double> py_ {};   // y - (1-2*alpha)|y - X*beta - f|
    DVector<double> pW_ {};   // diagonal of W^k = 1/(2*n*|y - X*beta - f|)
    DVector<double> mu_;      // \mu^k = [ \mu^k_1, ..., \mu^k_n ] : quantile vector at step k
    fdapde::SparseLU<SpMatrix<double>> invA_;

    FPIRLS<This> fpirls_;          // fpirls algorithm
    std::size_t max_iter_ = 200;   // maximum number of iterations in fpirls before forced stop
    double tol_ = 1e-6;            // fprils convergence tolerance
    double tol_weights_ = 1e-6;

    double eps_ = -1.0;   // pinball loss smoothing factor
    double pinball_loss(double x, double eps) const {   // quantile check function
        return (alpha_ - 1) * x + eps * fdapde::log1pexp(x / eps);
    };
    double pinball_loss(double x) const { return 0.5 * std::abs(x) + (alpha_ - 0.5) * x; }   // non-smoothed pinball
};

}   // namespace models
}   // namespace fdapde

#endif   // __QSRPDE_H__

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

#include <memory>
#include <type_traits>

#include <fdaPDE/pde.h>
#include <fdaPDE/utils.h>
using fdapde::core::PDEBase;

#include "../model_base.h"
#include "../model_macros.h"
#include "../model_traits.h"
#include "../sampling_design.h"
#include "regression_base.h"
#include "fpirls.h"

namespace fdapde {
namespace models {

// base class for GSRPDE model
template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver, typename Distribution>
class GSRPDE : public RegressionBase<GSRPDE<PDE, RegularizationType, SamplingDesign, Solver, Distribution>> {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
   private:
    typedef RegressionBase<GSRPDE<PDE, RegularizationType, SamplingDesign, Solver, Distribution>> Base;
    Distribution distr_ {};
    DVector<double> py_;   // \tilde y^k = G^k(y-u^k) + \theta^k
    DVector<double> pW_;   // diagonal of W^k = ((G^k)^{-2})*((V^k)^{-1})
    DVector<double> mu_;   // \mu^k = [ \mu^k_1, ..., \mu^k_n ] : mean vector at step k
    DMatrix<double> T_;    // T = \Psi^T*Q*\Psi + P

    // FPIRLS parameters (set to default)
    std::size_t max_iter_ = 200;
    double tol_ = 1e-4;
   public:
    IMPORT_REGRESSION_SYMBOLS;
    using Base::lambda_D;   // smoothing parameter in space
    using Base::P;          // discretized penalty
    using Base::W_;         // weight matrix
    // constructor
    GSRPDE() = default;
    // space-only constructor
    template <
      typename U = RegularizationType, typename std::enable_if<std::is_same<U, SpaceOnly>::value, int>::type = 0>
    GSRPDE(const PDE& pde) : Base(pde) {};
    // space-time constructor
    template <
      typename U = RegularizationType, typename std::enable_if<!std::is_same<U, SpaceOnly>::value, int>::type = 0>
    GSRPDE(const PDE& pde, const DVector<double>& time) : Base(pde, time) {};

    // setter
    void set_fpirls_tolerance(double tol) { tol_ = tol; }
    void set_fpirls_max_iter(std::size_t max_iter) { max_iter_ = max_iter; }

    void init_model() { return; };          // update model object in case of **structural** changes in its definition
    void update_to_weights() { return; };   // update model object in case of changes in the weights matrix
    virtual void solve();                   // finds a solution to the smoothing problem

    // required by FPIRLS (see fpirls.h for details)
    // initalizes mean vector \mu
    void fpirls_init() {
        mu_ = y();
        distr_.preprocess(mu_);
    };
    // computes W^k = ((G^k)^{-2})*((V^k)^{-1}) and y^k = G^k(y-u^k) + \theta^k
    void fpirls_pre_solve_step() {
        DVector<double> theta_ = distr_.link(mu_);   // \theta^k = (g(\mu^k_1), ..., g(\mu^k_n))
        DVector<double> G_ = distr_.der_link(mu_);   // G^k = diag(g'(\mu^k_1), ..., g'(\mu^k_n))
        DVector<double> V_ = distr_.variance(mu_);   // V^k = diag(v(\mu^k_1), ..., v(\mu^k_n))
        pW_ = ((G_.array().pow(2) * V_.array()).inverse()).matrix();
        py_ = G_.asDiagonal() * (y() - mu_) + theta_;
    }
    // updates mean vector \mu after WLS solution
    void fpirls_post_solve_step(const DMatrix<double>& hat_fitted) {
        mu_ = distr_.inv_link(hat_fitted);
    }
    // returns the data loss \norm{V^{-1/2}(y - \mu)}^2
    double data_loss() const {
        DVector<double> V = distr_.variance(mu_).array().sqrt().inverse().matrix();
        return (V.asDiagonal() * (y() - mu_)).squaredNorm();
    }
    const DVector<double>& py() const { return py_; }
    const DVector<double>& pW() const { return pW_; }

    // GCV support
    const DMatrix<double>& T();                                                  // T = \Psi^T*Q*\Psi + P
    double norm(const DMatrix<double>& op1, const DMatrix<double>& op2) const;   // total deviance \sum dev(op1 - op2)

    virtual ~GSRPDE() = default;
};

// implementative details
  
// finds a solution to the GSRPDE smoothing problem
template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver, typename Distribution>
void GSRPDE<PDE, RegularizationType, SamplingDesign, Solver, Distribution>::solve() {
    // execute FPIRLS for minimization of functional \norm{V^{-1/2}(y - \mu)}^2 + \lambda \int_D (Lf - u)^2
    FPIRLS<decltype(*this)> fpirls(*this, tol_, max_iter_);   // FPIRLS engine
    fpirls.compute();

    // fpirls converged: extract matrix W and solution estimates
    W_ = fpirls.solver().W();
    f_ = fpirls.solver().f();
    if (has_covariates()) { beta_ = fpirls.solver().beta(); }
    return;
}

template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver, typename Distribution>
const DMatrix<double>& GSRPDE<PDE, RegularizationType, SamplingDesign, Solver, Distribution>::T() {
    if (!has_covariates())   // case without covariates, Q is the identity matrix
        T_ = PsiTD() * W() * Psi() + P();
    else   // general case with covariates
        T_ = PsiTD() * lmbQ(Psi()) + P();
    return T_;
}

// returns the deviance of y - \hat y induced by the considered distribution
template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver, typename Distribution>
double GSRPDE<PDE, RegularizationType, SamplingDesign, Solver, Distribution>::norm(
  const DMatrix<double>& obs, const DMatrix<double>& fitted) const {
    DMatrix<double> mu = distr_.inv_link(fitted);
    double result = 0;
    for (std::size_t i = 0; i < obs.rows(); ++i) { result += distr_.deviance(mu.coeff(i, 0), obs.coeff(i, 0)); }
    return result;
}

template <
  typename PDE_, typename RegularizationType_, typename SamplingDesign_, typename Solver_, typename DistributionType_>
struct model_traits<GSRPDE<PDE_, RegularizationType_, SamplingDesign_, Solver_, DistributionType_>> {
    typedef PDE_ PDE;
    typedef RegularizationType_ regularization;
    typedef SamplingDesign_ sampling;
    typedef Solver_ solver;
    typedef DistributionType_ DistributionType;
    enum { N = PDE::N, M = PDE::M, n_lambda = n_smoothing_parameters<RegularizationType_>::value };
};
// specialization for separable regularization
template <typename PDE_, typename SamplingDesign_, typename Solver_, typename DistributionType_>
struct model_traits<GSRPDE<PDE_, fdapde::models::SpaceTimeSeparable, SamplingDesign_, Solver_, DistributionType_>> {
    typedef PDE_ PDE;
    typedef fdapde::models::SpaceTimeSeparable regularization;
    typedef SplineBasis<3> TimeBasis;   // use cubic B-splines
    typedef SamplingDesign_ sampling;
    typedef Solver_ solver;
    typedef DistributionType_ DistributionType;
    enum { N = PDE::N, M = PDE::M, n_lambda = 2 };
};

// gsrpde trait
template <typename Model> struct is_gsrpde {
    static constexpr bool value = is_instance_of<Model, GSRPDE>::value;
};
  
}   // namespace models
}   // namespace fdapde

#endif   // __GSRPDE_H__

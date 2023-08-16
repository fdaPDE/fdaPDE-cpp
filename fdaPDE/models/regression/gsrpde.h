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
    DiagMatrix<double> W_;
    Distribution distribution_ {};
    DVector<double> py_;   // \tilde y^k = G^k(y-u^k) + \theta^k
    DVector<double> pW_;   // diagonal of W^k = ((G^k)^{-2})*((V^k)^{-1})
    DMatrix<double> T_;    // T = \Psi^T*Q*\Psi + P

    // FPIRLS parameters (set to default)
    std::size_t max_iter_ = 15;
    double tol_ = 0.0002020;
   public:
    IMPORT_REGRESSION_SYMBOLS;
    using Base::lambda_D;   // smoothing parameter in space
    using Base::P;          // discretized penalty
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

    // required by FPIRLS: computes W^k = ((G^k)^{-2})*((V^k)^{-1}) and \tilde y^k = G^k(y-u^k) + \theta^k
    std::tuple<DVector<double>&, DVector<double>&> compute(const DVector<double>& mu);
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
    FPIRLS<decltype(*this), Distribution> fpirls(*this, tol_, max_iter_);   // FPIRLS engine
    fpirls.compute();

    // fpirls converged: extract matrix P and solution estimates
    W_ = fpirls.weights().asDiagonal();
    f_ = fpirls.f();
    if (has_covariates()) { beta_ = fpirls.beta(); }
    return;
}

// required by FPIRLS: computes W^k = ((G^k)^{-2})*((V^k)^{-1}) and \tilde y^k = G^k(y-u^k) + \theta^k
template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver, typename Distribution>
std::tuple<DVector<double>&, DVector<double>&>
GSRPDE<PDE, RegularizationType, SamplingDesign, Solver, Distribution>::compute(const DVector<double>& mu) {
    DVector<double> theta_ = distribution_.link(mu);   // \theta^k = [ g(\mu^k_1), ..., g(\mu^k_n) ]
    DVector<double> G_ = distribution_.der_link(mu);   // G^k = diag(g'(\mu^k_1), ..., g'(\mu^k_n))
    DVector<double> V_ = distribution_.variance(mu);   // V^k = diag(v(\mu^k_1), ..., v(\mu^k_n))
    // compute weight matrix and pseudo-observation vector
    pW_ = ((G_.array().pow(2) * V_.array()).inverse()).matrix();
    py_ = G_.asDiagonal() * (y() - mu) + theta_;
    return std::tie(pW_, py_);
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
    double result = 0;
    for (std::size_t i = 0; i < obs.rows(); ++i) {
        result += distribution_.deviance(obs.coeff(i, 0), fitted.coeff(i, 0));
    }
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
    enum { N = PDE::N, M = PDE::M, R = PDE::R, n_lambda = n_smoothing_parameters<RegularizationType_>::value };
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
    enum { N = PDE::N, M = PDE::M, R = PDE::R, n_lambda = 2 };
};

// gsrpde trait
template <typename Model> struct is_gsrpde {
    static constexpr bool value = is_instance_of<Model, GSRPDE>::value;
};
  
}   // namespace models
}   // namespace fdapde

#endif   // __GSRPDE_H__

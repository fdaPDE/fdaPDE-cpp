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

#ifndef __DEPDE_H__
#define __DEPDE_H__

#include <fdaPDE/optimization.h>
#include "../model_base.h"
#include "../model_macros.h"
#include "../sampling_design.h"
#include "density_estimation_base.h"

namespace fdapde {
namespace models {

template <typename RegularizationType_>
class DEPDE : public DensityEstimationBase<DEPDE<RegularizationType_>, RegularizationType_> {
   private:
    using Base = DensityEstimationBase<DEPDE<RegularizationType_>, RegularizationType_>;
    double tol_ = 1e-5;        // tolerance on custom stopping criterion
    DVector<double> g_init_;   // density initialization basis expansion coefficients
    core::Optimizer<DEPDE<RegularizationType_>> opt_;
   public:
    using RegularizationType = std::decay_t<RegularizationType_>;
    using This = DEPDE<RegularizationType>;
    using Base::grad_int_exp;   // \nabla_g (\int_{\Omega} \exp(g))
    using Base::int_exp;        // \int_{\Omega} \exp(g)
    using Base::n_locs;         // overall number of data locations
    using Base::n_obs;          // number of observations (!= n_locs if masked)
    using Base::P;              // discretized penalty matrix R_1^\top * (R_0)^{-1} * R_1
    using Base::PsiQuad;        // reference basis evaluation at quadrature nodes
    using Base::Upsilon;        // \Upsilon_(i,:) = is_space_only<Model> ? \Psi : \Phi(i,:) \kron S(p_i) \Psi

    // space-only constructor
    template <typename PDE_>
    DEPDE(PDE_&& pde) requires(is_space_only<This>::value) : Base(pde) { }
    // space-time separable constructor
    template <typename SpacePDE_, typename TimePDE_>
    DEPDE(SpacePDE_&& space_penalty, TimePDE_&& time_penalty) requires(is_space_time_separable<This>::value)
        : Base(space_penalty, time_penalty) { }

    // K-fold cross validation index
    struct CVScore {
        CVScore(DEPDE& model) : model_(model) { }
        double operator()(
          const DVector<double>& lambda, [[maybe_unused]] const BinaryVector<fdapde::Dynamic>& train_mask,
          const BinaryVector<fdapde::Dynamic>& test_mask) {
            model_.set_lambda(lambda);
            // fit model on train set
            model_.set_mask(test_mask);   // discard test set from training phase
            model_.init();
            model_.solve();
            double test_err = 0;
            for (int i = 0; i < test_mask.size(); ++i) {
	      if (test_mask[i]) { test_err += std::exp(model_.Psi().row(i) * model_.g()); }
            }
            return model_.int_exp(2. * model_.g()) - 2. / test_mask.count() * test_err;
        }
       private:
        DEPDE& model_;
    };

    // evaluates penalized negative log-likelihood at point
    // L(g) = - 1^\top*\Upsilon*g + \sum_{e \in mesh} w^\top*exp[\Psi_q*g_e] + \lambda_S*g^\top*P*g
    double operator()(const DVector<double>& g) {
        return -(Upsilon() * g).sum() + n_obs() * int_exp(g) + g.dot(P() * g);
    }
    // log-likelihood gradient functor
    // \nabla_g(L(g)) = -\Upsilon^\top*1 + n*\sum_{e \in mesh} w*exp[\Psi_q*g_e]*\Psi_q^\top + 2*P*g
    std::function<DVector<double>(const DVector<double>&)> derive() {
        return [this, dllik = DVector<double>(-Upsilon().transpose() * DVector<double>::Ones(n_locs()))](
                 const DVector<double>& g) {
	    return DVector<double>(dllik + n_obs() * grad_int_exp(g) + 2 * P() * g);
	};
    }
    // optimization algorithm custom stopping criterion
    template <typename OptimizerType> bool opt_stopping_criterion(OptimizerType& opt) {
        // opt.x_old: density function before update step, opt.x_new: density function after update step
        bool stop = true;
        double llik_old = -(Upsilon() * opt.x_old).sum() + n_obs() * int_exp(opt.x_old);
        double llik_new = -(Upsilon() * opt.x_new).sum() + n_obs() * int_exp(opt.x_new);
	stop &= std::abs((llik_new - llik_old) / llik_old) < tol_;
	if(!stop) return false;	
        double penD_old = Base::xtPDx(opt.x_old), penD_new = Base::xtPDx(opt.x_new); //opt.x_old.dot(Base::PD() * opt.x_old);
        stop &= std::abs((penD_new - penD_old) / penD_old) < tol_;
        if (!stop) return false;
        double loss_old = llik_old + Base::lambda_D() * penD_old, loss_new = llik_new + Base::lambda_D() * penD_new;
        // space-time case
        if constexpr (is_space_time_separable<This>::value) {
            double penT_old = Base::xtPTx(opt.x_old), penT_new = Base::xtPTx(opt.x_new);
            loss_old += Base::lambda_T() * penT_old;
            loss_new += Base::lambda_T() * penT_new;
            stop &= std::abs((penT_new - penT_old) / penT_old) < tol_;
        }
	stop &= std::abs((loss_new - loss_old) / loss_old) < tol_;
	return stop;
    }
    // call optimization algorithm for log-likelihood minimization
    void solve() { Base::g_ = opt_.optimize(*this, g_init_); }
    void init_model() { return; }
    // setters
    void set_tolerance(double tol) { tol_ = tol; }
    void set_g_init(const DVector<double>& g_init) { g_init_ = g_init; }
    void set_f_init(const DVector<double>& f_init) { g_init_ = f_init.array().log(); }
    template <typename OptimizerType> void set_optimizer(OptimizerType&& opt) { opt_ = opt; }
    // getters
    const DVector<double>& g_init() const { return g_init_; }
};

}   // namespace models
}   // namespace fdapde

#endif   // __DEPDE_H__

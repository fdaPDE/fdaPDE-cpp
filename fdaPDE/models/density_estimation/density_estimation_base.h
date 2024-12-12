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

#ifndef __DENSITY_ESTIMATION_BASE_H__
#define __DENSITY_ESTIMATION_BASE_H__

#include <fdaPDE/utils.h>
#include "../model_macros.h"
#include "../model_traits.h"
#include "../space_only_base.h"
#include "../space_time_separable_base.h"
#include "../sampling_design.h"

namespace fdapde {
namespace models {

#define SPACE_LOCS "SPACE_LOCS"
#define TIME_LOCS "TIME_LOCS"
  
namespace internal {
consteval int de_select_quadrature(int local_dim) { return local_dim == 1 ? 4 : local_dim == 2 ? 6 : 14; }
};   // namespace internal

template <typename Model, typename RegularizationType>
class DensityEstimationBase :
    public select_regularization_base<Model, RegularizationType>::type,
    public SamplingBase<Model> {
   protected:
    std::function<double(const DVector<double>&)> int_exp_;   // \int_{\mathcal{D}} exp(g(x)), x \in \mathcal{D}
    std::function<DVector<double>(const DVector<double>&)> grad_int_exp_;   // gradient of int_exp
    DVector<double> g_;          // expansion coefficients vector of estimated density function
    DMatrix<double> PsiQuad_;    // reference element basis evaluation at quadrature nodes
    DVector<double> w_;          // quadrature weights
    SpMatrix<double> Upsilon_;   // \Upsilon_(i,:) = \Phi(i,:) \kron S(p_i) \Psi
    SpMatrix<double> B_;         // \Upsilon matrix corrected for masked observations
    BinaryVector<Dynamic> mask_;

    template <typename IntegratorType_, typename BasisType_>
    DMatrix<double> eval_basis_at_quadrature_(BasisType_&& basis, IntegratorType_&& integrator) const {
        using IntegratorType = std::decay_t<IntegratorType_>;
        DMatrix<double> result(IntegratorType::num_nodes, basis.size());
        for (int i = 0; i < IntegratorType::num_nodes; ++i) {
            for (int j = 0; j < basis.size(); ++j) { result(i, j) = basis[j](integrator.nodes[i]); }
        }
        return result;
    }
   public:
    static constexpr int DomainDimension = fdapde::Dynamic;
    using Base = typename select_regularization_base<Model, RegularizationType>::type;
    using Base::model;
    using SamplingBase<Model>::Psi;     // matrix of basis evaluations at data locations
    using SamplingBase<Model>::n_spatial_locs;

    DensityEstimationBase() = default;
    // space-only constructor
    template <typename PDE_>
    DensityEstimationBase(PDE_&& pde)
        requires(is_space_only<Model>::value)
        : Base(pde), SamplingBase<Model>(Sampling::pointwise) {
        pde.init();   // early PDE initialization
        constexpr int local_dim = std::decay_t<PDE_>::SpaceDomainType::local_dim;
        constexpr int n_quadrature_nodes = internal::de_select_quadrature(local_dim);
        // compute reference basis evaluation at quadrature nodes
        core::IntegratorTable<local_dim, n_quadrature_nodes> integrator {};
        PsiQuad_ = eval_basis_at_quadrature_(pde.reference_basis(), integrator);
        w_ = Eigen::Map<SVector<n_quadrature_nodes>>(integrator.weights.data());
        // store functor for approximation of \int_{\mathcal{D}} \exp(g(p)). Computes
        // \sum_{e \in mesh} {e.measure() * \sum_j {w_j * exp[\sum_i {g_i * \psi_i(q_j)}]}}
        int_exp_ = [&](const DVector<double>& g) {
            double result = 0;
            for (auto e = pde.domain().cells_begin(); e != pde.domain().cells_end(); ++e) {
                result += w_.dot((PsiQuad_ * g(pde.dofs().row(e->id()))).array().exp().matrix()) * e->measure();
            }
            return result;
        };
	// store functor for computation of \nabla_g(\int_{\mathcal{D}} \exp(g(p)))
        grad_int_exp_ = [&](const DVector<double>& g) {
            DVector<double> grad = DVector<double>::Zero(g.rows());
            for (auto e = pde.domain().cells_begin(); e != pde.domain().cells_end(); ++e) {
                grad(pde.dofs().row(e->id())) +=
                  PsiQuad_.transpose() *
                  (PsiQuad_ * g(pde.dofs().row(e->id()))).array().exp().cwiseProduct(w_.array()).matrix() *
                  e->measure();
            }
            return grad;
        };
    }
    // space-time separable constructor
    template <typename SpacePDE_, typename TimePDE_>
    DensityEstimationBase(SpacePDE_&& s_pen, TimePDE_&& t_pen)
        requires(is_space_time_separable<Model>::value)
        : Base(s_pen, t_pen), SamplingBase<Model>(Sampling::pointwise) {
        // early penalties initialization
        s_pen.init();
        t_pen.init();
        constexpr int local_dim = std::decay_t<SpacePDE_>::SpaceDomainType::local_dim;
        constexpr int n_quadrature_s = internal::de_select_quadrature(local_dim);
        constexpr int n_quadrature_t = 5;
        // compute reference basis evaluation at quadrature nodes
        core::IntegratorTable<local_dim, n_quadrature_s> s_integrator {};
        core::IntegratorTable<1, n_quadrature_t, core::GaussLegendre> t_integrator {};
        DMatrix<double> Psi = eval_basis_at_quadrature_(s_pen.reference_basis(), s_integrator);
        DMatrix<double> Phi = eval_basis_at_quadrature_(t_pen.reference_basis(), t_integrator);
        PsiQuad_ = Kronecker(Phi, Psi);
        w_ = Kronecker(
          Eigen::Map<SVector<n_quadrature_t>>(t_integrator.weights.data()),
          Eigen::Map<SVector<n_quadrature_s>>(s_integrator.weights.data()));
        // store functor for approximation of double integral \int_T \int_{\mathcal{D}} \exp(g(p, t))
        int_exp_ =
          [&, n = s_pen.reference_basis().size(), m = t_pen.reference_basis().size()](const DVector<double>& g) {
              double result = 0;
              DVector<int> active_dofs(n * m);
              for (auto e = s_pen.domain().cells_begin(); e != s_pen.domain().cells_end(); ++e) {       // space loop
                  for (auto i = t_pen.domain().cells_begin(); i != t_pen.domain().cells_end(); ++i) {   // time loop
                      for (int j = 0; j < m; ++j) {   // compute active dofs
                          active_dofs.middleRows(j * n, n) =
                            s_pen.dofs().row(e->id()).transpose().array() + j * Base::n_spatial_basis();
                      }
                      result += w_.dot((PsiQuad_ * g(active_dofs)).array().exp().matrix()) * e->measure() *
                                (0.5 * i->measure());
                  }
              }
              return result;
          };
        // store functor for computation of \nabla_g(\int_T \int_{\mathcal{D}} \exp(g(p, t)))
        grad_int_exp_ =
	  [&, n = s_pen.reference_basis().size(), m = t_pen.reference_basis().size()](const DVector<double>& g) {
            DVector<double> grad = DVector<double>::Zero(g.rows());
            DVector<int> active_dofs(n * m);
            for (auto e = s_pen.domain().cells_begin(); e != s_pen.domain().cells_end(); ++e) {       // space loop
                for (auto i = t_pen.domain().cells_begin(); i != t_pen.domain().cells_end(); ++i) {   // time loop
                    for (int j = 0; j < m; ++j) {   // compute active dofs
                        active_dofs.middleRows(j * n, n) =
                          s_pen.dofs().row(e->id()).transpose().array() + j * Base::n_spatial_basis();
                    }
                    grad(active_dofs) += PsiQuad_.transpose() *
                                         ((PsiQuad_ * g(active_dofs)).array().exp()).cwiseProduct(w_.array()).matrix() *
                                         e->measure() * (0.5 * i->measure());
                }
            }
            return grad;
        };
    }

    // setters
    // redefine set_data() as it forwards data as locations in density estimation
    void set_data(const BlockFrame<double, int>& df, bool reindex = false) {
        fdapde_assert(
          df.cols() > 0 &&
          ((is_space_only<Model>::value && df.has_block(SPACE_LOCS)) ||
            (is_space_time_separable<Model>::value && df.has_block(SPACE_LOCS) && df.has_block(TIME_LOCS))));
        Base::set_data(df, reindex);
        if (df.has_block(SPACE_LOCS)) SamplingBase<Model>::set_spatial_locations(df.get<double>(SPACE_LOCS));
        if constexpr (is_space_time_separable<Model>::value) Base::set_temporal_locations(df.get<double>(TIME_LOCS));
    }
    void set_mask(const BinaryVector<fdapde::Dynamic>& mask) {
        fdapde_assert(mask.size() == n_locs());
        if (mask.any()) {
            model().runtime().set(runtime_status::require_psi_correction);
            mask_ = mask;   // mask[i] == true, removes the contribution of the i-th observation from fit
        }
    }
    // getters
    int n_obs() const { return Base::df_.rows() - mask_.count(); };
    int n_locs() const { return Base::df_.rows(); }
    const SpMatrix<double>& Psi() const { return is_space_only<Model>::value ? Psi(not_nan()) : Upsilon_; }
    const SpMatrix<double>& Upsilon() const { return !is_empty(B_) ? B_ : Psi(); }
    const DMatrix<double>& PsiQuad() const { return PsiQuad_; }
    const DMatrix<double>& w() const { return w_; }
    double int_exp(const DVector<double>& g) const { return int_exp_(g); }
    double int_exp() const { return int_exp_(g_); }
    DVector<double> grad_int_exp(const DVector<double>& g) const { return grad_int_exp_(g); }
    DVector<double> grad_int_exp() const { return grad_int_exp_(g_); }
    const DVector<double>& g() const { return g_; }          // expansion coefficient vector of log density field
    DVector<double> f() const { return g_.array().exp(); }   // expansion coefficient vector of density field
    const BinaryVector<Dynamic>& masked_obs() const { return mask_; }

    // initialization methods
    void analyze_data() {
        if constexpr (is_space_time_separable<Model>::value) {
            if (is_empty(Upsilon_)) {
                Upsilon_.resize(Base::n_temporal_locs(), Base::n_basis());
                std::vector<fdapde::Triplet<double>> triplet_list;
                // reserve with some bound
                for (int i = 0; i < Base::n_temporal_locs(); ++i) {
                    // kronecker product between Phi i-th row and Psi i-th row
                    SpMatrix<double> tmp =
                      Kronecker(SpMatrix<double>(Base::Phi().row(i)), SpMatrix<double>(Psi().row(i)));
                    for (int j = 0; j < tmp.outerSize(); ++j)
                        for (SpMatrix<double>::InnerIterator it(tmp, j); it; ++it) {
                            triplet_list.emplace_back(i, it.col(), it.value());
                        }
                }
                Upsilon_.setFromTriplets(triplet_list.begin(), triplet_list.end());
                Upsilon_.makeCompressed();
            }
        }
        return;
    }
    void tensorize_psi() { return; }   // avoid tensorization for space-time problems
    void correct_psi() {
      if (masked_obs().any()) B_ = (~masked_obs().repeat(1, Base::n_basis())).select(Psi());
    }
};

}   // namespace models
}   // namespace fdapde

#endif   // __DENSITY_ESTIMATION_BASE_H__

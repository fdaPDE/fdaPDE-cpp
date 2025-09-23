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

#ifndef __FE_DE_ELLIPTIC_SOLVER_H__
#define __FE_DE_ELLIPTIC_SOLVER_H__

#include "header_check.h"

namespace fdapde {
namespace internals {

struct fe_de_elliptic {
   private:
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using binary_t = BinaryMatrix<Dynamic, Dynamic>;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using sparse_solver_t = eigen_sparse_solver_movable_wrap<Eigen::SparseLU<sparse_matrix_t>>;
    template <typename DataLocs>
    static constexpr bool is_valid_data_locs_descriptor_v = std::is_same_v<DataLocs, matrix_t>;
    template <typename Penalty> struct is_valid_penalty {
        static constexpr bool value = requires(Penalty penalty) {
            penalty.bilinear_form();
            penalty.linear_form();
        };
    };
    template <typename Penalty> static constexpr bool is_valid_penalty_v = is_valid_penalty<Penalty>::value;
    // high-order quadrature for integration of constraint \int_D (e^g)
    template <int EmbedDim> struct de_quadrature {
        using type = std::conditional_t<EmbedDim == 1, QS1DP7_, std::conditional_t<EmbedDim == 2, QS2DP4_, QS3DP5_>>;
    };
    template <int EmbedDim> using de_quadrature_t = de_quadrature<EmbedDim>::type;
   public:
    static constexpr int n_lambda = 1;
    using solver_category = de_solver;

    // penalized negative log-likelihood objective functor
    struct llik_t {
        llik_t(fe_de_elliptic& m, double lambda) : m_(std::addressof(m)), lambda_(lambda) { }
        llik_t(fe_de_elliptic& m, double lambda, double tol) : m_(std::addressof(m)), lambda_(lambda), tol_(tol) { }
        // penalized negative log-likelihood at point
        double operator()(const vector_t& g) {
            return -(m_->Psi_ * g).sum() + m_->n_obs_ * m_->int_exp_(g) + lambda_ * g.dot(m_->P_ * g);
        }
        // gradient functor
        std::function<vector_t(const vector_t&)> gradient() {
            return [this, dllik = vector_t(-m_->Psi_.transpose() * vector_t::Ones(m_->n_obs_))](const vector_t& g) {
                return vector_t(dllik + m_->n_obs_ * m_->grad_int_exp_(g) + 2 * lambda_ * m_->P_ * g);
            };
        }
        // injected optimization stopping criterion
        template <typename Optimizer> bool stop_if(Optimizer& opt) {
            double llik_old = -(m_->Psi_ * opt.x_old).sum() + m_->n_obs_ * m_->int_exp_(opt.x_old);
            double llik_new = -(m_->Psi_ * opt.x_new).sum() + m_->n_obs_ * m_->int_exp_(opt.x_new);
            if (std::abs((llik_new - llik_old) / llik_old) > tol_) { return false; }
            double penD_old = opt.x_old.dot(m_->P_ * opt.x_old);
            double penD_new = opt.x_new.dot(m_->P_ * opt.x_new);
            if (std::abs((penD_new - penD_old) / penD_old) > tol_) { return false; }
            double loss_old = llik_old + lambda_ * penD_old;
            double loss_new = llik_new + lambda_ * penD_new;
            return std::abs((loss_new - loss_old) / loss_old) < tol_;
        }
       private:
        fe_de_elliptic* m_;
        double lambda_;
        double tol_ = 1e-5;
    };

    fe_de_elliptic() noexcept = default;
    template <typename GeoFrame, typename Penalty>
        requires(is_valid_penalty_v<Penalty>)
    fe_de_elliptic(const GeoFrame& gf, Penalty&& penalty) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        discretize(penalty);
        analyze_data(gf);
    }

    // perform finite element based numerical discretization
    template <typename Penalty> void discretize(Penalty&& penalty) {
        using BilinearForm = typename std::decay_t<Penalty>::BilinearForm;
        using LinearForm = typename std::decay_t<Penalty>::LinearForm;
        fdapde_static_assert(
          internals::is_valid_penalty_pair_v<BilinearForm FDAPDE_COMMA LinearForm>, INVALID_PENALTY_DESCRIPTION);
        using FeSpace = typename BilinearForm::TrialSpace;
	using DofHandler = typename FeSpace::DofHandlerType;
	using Triangulation = typename FeSpace::Triangulation;
        constexpr int embed_dim = Triangulation::embed_dim;

        // discretization
        const BilinearForm& bilinear_form = penalty.bilinear_form();
        const LinearForm& linear_form = penalty.linear_form();
        const FeSpace& fe_space = bilinear_form.trial_space();
        const DofHandler& dof_handler = fe_space.dof_handler();
        n_dofs_ = bilinear_form.n_dofs();   // number of basis functions over physical domain
        internals::fe_mass_assembly_loop<FeSpace> mass_assembler(bilinear_form.trial_space());
        R0_ = mass_assembler.assemble();
        R1_ = bilinear_form.assemble();
        u_ = linear_form.assemble();

        // penalty matrix
        sparse_solver_t invR0;
        invR0.compute(R0_);
        P_ = R1_.transpose() * invR0.solve(R1_);
        // store handles for basis system evaluation at locations
        point_eval_ = [fe_space = bilinear_form.trial_space()](const matrix_t& locs) -> decltype(auto) {
            return internals::point_basis_eval(fe_space, locs);
        };

	// geometry
        const Triangulation& triangulation = fe_space.triangulation();
        // eval reference basis at quadrature nodes, store de_quadrature weights
        de_quadrature_t<embed_dim> quad_rule;
        int n_quad_nodes = quad_rule.order;
        int n_shape_functions = fe_space.n_shape_functions();
        matrix_t PsiQuad(n_quad_nodes, n_shape_functions);
	vector_t w(n_quad_nodes);
        for (int i = 0; i < n_quad_nodes; ++i) {
            for (int j = 0; j < n_shape_functions; ++j) {
                PsiQuad(i, j) = fe_space.eval_shape_value(j, quad_rule.nodes.row(i).transpose());
            }
            w[i] = quad_rule.weights[i];
        }
        // store handle for approximation of \int_D (e^g)
        int_exp_ = [&, dof_handler, PsiQuad, w](const vector_t& g) {
            double val_ = 0;
            for (auto it = triangulation.cells_begin(); it != triangulation.cells_end(); ++it) {
                val_ += w.dot((PsiQuad * g(dof_handler.dofs().row(it->id()))).array().exp().matrix()) * it->measure();
            }
	    return val_;
        };
        // store handle for approximation of \nabla_g(\int_D (e^g))
        grad_int_exp_ = [&, dof_handler, PsiQuad, w](const vector_t& g) {
            vector_t grad = vector_t::Zero(g.rows());
            for (auto it = triangulation.cells_begin(); it != triangulation.cells_end(); ++it) {
                grad(dof_handler.dofs().row(it->id())) +=
                  PsiQuad.transpose() *
                  (PsiQuad * g(dof_handler.dofs().row(it->id()))).array().exp().cwiseProduct(w.array()).matrix() *
                  it->measure();
            }
            return grad;
        };
        return;
    }
    // fit from geoframe
    template <typename GeoFrame> void analyze_data(const GeoFrame& gf) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        fdapde_assert(gf.n_layers() == 1 && gf[0].category()[0] == ltype::point);
        n_obs_ = gf[0].rows();
        // eval physical basis at spatial locations
        const auto& spatial_index = geo_index_cast<0, POINT>(gf[0]);
        if (spatial_index.points_at_dofs()) {
            Psi_.resize(n_obs_, n_dofs_);
            Psi_.setIdentity();
        } else {
            Psi_ = point_eval_(spatial_index.coordinates());
        }
	return;
    }
    // main fit entry point
    template <typename Optimizer, typename... Callbacks>
    const vector_t& fit(double lambda, const vector_t& g_init, Optimizer&& opt, Callbacks&&... callbacks) {
        g_ = opt.optimize(llik_t(*this, lambda, tol_), g_init, std::forward<Callbacks>(callbacks)...);
        return g_;
    }
    template <typename Optimizer, typename LambdaT, typename... Callbacks>
        requires(internals::is_vector_like_v<LambdaT>)
    const vector_t& fit(LambdaT&& lambda, const vector_t& g_init, Optimizer&& opt, Callbacks&&... callbacks) {
        fdapde_assert(lambda.size() == n_lambda);
        return fit(lambda[0], g_init, opt, std::forward<Callbacks>(callbacks)...);
    }
    // modifiers
    void set_llik_tolerance(double tol) { tol_ = tol; }

    // observers
    const sparse_matrix_t& mass() const { return R0_; }
    const sparse_matrix_t& stiff() const { return R1_; }
    const sparse_matrix_t& Psi() const { return Psi_; }
    double int_exp(const vector_t& g) const { return int_exp_(g); }
    double int_exp() const { return int_exp_(g_); }
    vector_t grad_int_exp(const vector_t& g) const { return grad_int_exp_(g); }
    vector_t grad_int_exp() const { return grad_int_exp_(g_); }
    const vector_t& log_density() const { return g_; }
    vector_t density() const { return g_.array().exp(); }
    vector_t gn() const { return Psi_ * g_; }
    vector_t fn() const { return Psi_ * g_.array().exp().matrix(); }
   private:
    int n_dofs_ = 0, n_obs_ = 0;
    sparse_matrix_t R0_;    // n_dofs x n_dofs matrix [R0]_{ij} = \int_D \psi_i * \psi_j
    sparse_matrix_t R1_;    // n_dofs x n_dofs matrix [R1]_{ij} = \int_D a(\psi_i, \psi_j)
    sparse_matrix_t Psi_;   // n_obs x n_dofs matrix [Psi]_{ij} = \psi_j(p_i)
    vector_t u_;            // n_dofs x 1 vector u_i = \int_D u * \psi_i

    matrix_t P_;   // n_dofs x n_dofs penalty matrix P_ = R1^\top * (R0)^{-1} * R1
    std::function<double(const vector_t&)> int_exp_;          // \int exp(g)
    std::function<vector_t(const vector_t&)> grad_int_exp_;   // \nabla_g \int exp(g)
    vector_t g_;
    // basis system evaluation handle
    std::function<sparse_matrix_t(const matrix_t& locs)> point_eval_;
    double tol_ = 1e-5;
};

}   // namespace internals

// elliptic solver API
template <typename BilinearForm_, typename LinearForm_> struct fe_de_elliptic {
    using solver_t = internals::fe_de_elliptic;
   private:
    struct penalty_packet {
        using BilinearForm = std::decay_t<BilinearForm_>;
        using LinearForm = std::decay_t<LinearForm_>;
       private:
        BilinearForm bilinear_form_;
        LinearForm linear_form_;
       public:
        penalty_packet(const BilinearForm_& bilinear_form, const LinearForm_& linear_form) :
            bilinear_form_(bilinear_form), linear_form_(linear_form) { }
        // observers
        const BilinearForm& bilinear_form() const { return bilinear_form_; }
        const LinearForm& linear_form() const { return linear_form_; }
    };
   public:
    fe_de_elliptic(const BilinearForm_& bilinear_form, const LinearForm_& linear_form) :
        penalty_(bilinear_form, linear_form) { }
    const penalty_packet& get() const { return penalty_; }
   private:
    penalty_packet penalty_;
};

}   // namespace fdapde

#endif   // __FE_DE_ELLIPTIC_SOLVER_H__

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

#ifndef __FE_DE_SEPARABLE_SOLVER_H__
#define __FE_DE_SEPARABLE_SOLVER_H__

#include "header_check.h"

namespace fdapde {
namespace internals {

struct fe_de_separable {
   private:
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using binary_t = BinaryMatrix<Dynamic, Dynamic>;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using sparse_solver_t = eigen_sparse_solver_movable_wrap<Eigen::SparseLU<sparse_matrix_t>>;
    template <typename DataLocs>
    static constexpr bool is_valid_data_locs_descriptor_v = std::is_same_v<DataLocs, matrix_t>;
    template <typename Penalty> struct is_valid_penalty {
       private:
        using LhsPenalty = typename Penalty::LhsPenalty;
        using RhsPenalty = typename Penalty::RhsPenalty;
        template <typename Penalty_> struct is_valid_penalty_impl {
            static constexpr bool value = requires(Penalty_ penalty) {
                penalty.bilinear_form();
                penalty.linear_form();
            };
        };
       public:
        static constexpr bool value = requires(Penalty penalty) {
            penalty.lhs_penalty();
            penalty.rhs_penalty();
        } && is_valid_penalty_impl<LhsPenalty>::value && is_valid_penalty_impl<RhsPenalty>::value;
    };
    template <typename Penalty> static constexpr bool is_valid_penalty_v = is_valid_penalty<Penalty>::value;

    template <typename Tuple> struct function_space_tuple {
        using type = decltype([]<size_t... Is_>(std::index_sequence<Is_...>) {
            return std::make_tuple(typename std::tuple_element_t<Is_, Tuple>::TrialSpace {}...);
        }(std::make_index_sequence<std::tuple_size_v<Tuple>>()));
    };  
    template <typename Penalty1, typename Penalty2>
    const auto& fe_penalty_(const Penalty1& penalty1, const Penalty2& penalty2) const {
        return select_one_between(
          penalty1, penalty2, []() { return  is_fe_space_v<typename Penalty1::BilinearForm::TrialSpace>; });
    }
    template <typename Penalty1, typename Penalty2>
    const auto& bs_penalty_(const Penalty1& penalty1, const Penalty2& penalty2) const {
        return select_one_between(
          penalty1, penalty2, []() { return !is_fe_space_v<typename Penalty1::BilinearForm::TrialSpace>; });
    }
  
    // high-order quadrature for integration of constraint \int_D (e^g)
    template <int EmbedDim> struct de_fe_quadrature {
        using type = std::conditional_t<EmbedDim == 1, QS1DP7_, std::conditional_t<EmbedDim == 2, QS2DP4_, QS3DP5_>>;
    };
    template <int EmbedDim> using de_fe_quadrature_t = de_fe_quadrature<EmbedDim>::type;
    using de_bs_quadrature_t = QGL1DP9_;
    template <typename InfoT> struct is_valid_info_t {
        static constexpr bool value = requires(InfoT info) { info.penalty; };
    };
    template <typename GeoFrame, typename Penalty> auto tuplify_penalty_(const GeoFrame& gf, Penalty&& penalty) {
        using Penalty_ = std::decay_t<Penalty>;
        return internals::apply_index_pack<n_lambda>([&]<int... Ns_>() {
            return std::make_tuple([&]() {
                using T = std::tuple_element_t<Ns_, Penalty_>;
                if constexpr (requires(T t) { t.get(); }) {
                    return std::get<Ns_>(penalty).get();
                } else {
                    if constexpr (internals::is_pair_v<T>) {
                        return std::get<Ns_>(penalty);   // user supplied penalty pair
                    } else {
                        return std::get<Ns_>(penalty)(gf.template triangulation<Ns_>()).get();
                    }
                }
            }()...);
        });
    }
    // evaluates reference basis system at quadrature nodes (only active dofs considered)
    template <typename FuncSpace, typename Quadrature>
    matrix_t eval_fe_shape_values_at_quadrature_(const FuncSpace& func_space, const Quadrature& quad) const {
        int n_quad_nodes = quad.order;
        int n_shape_functions = func_space.n_shape_functions();
        matrix_t m(n_quad_nodes, n_shape_functions);
        for (int i = 0; i < n_quad_nodes; ++i) {
            for (int j = 0; j < n_shape_functions; ++j) {
                m(i, j) = func_space.eval_shape_value(j, quad.nodes.row(i).transpose());
            }
        }
        return m;
    }
    template <typename FuncSpace, typename Quadrature, typename CellIterator>
    matrix_t
    eval_bs_shape_values_at_quadrature_(const FuncSpace& func_space, const Quadrature& quad, CellIterator it) const {
        int n_quad_nodes = quad.order;
	const auto& dof_handler = func_space.dof_handler();
	std::vector<int> active_dofs = dof_handler.active_dofs(it->id());
        matrix_t m(n_quad_nodes, active_dofs.size());
        double a = it->nodes()[0], b = it->nodes()[1];   // cell range
        for (int i = 0; i < n_quad_nodes; ++i) {
            for (int j = 0; j < active_dofs.size(); ++j) {
                m(i, j) =
                  func_space.eval_cell_value(active_dofs[j], it->id(), (b - a) / 2 * quad.nodes[i] + (b + a) / 2);
            }
        }
        return m;
    }
   public:
    static constexpr int n_lambda = 2;
    using solver_category = de_solver;

    // penalized negative log-likelihood objective functor
    struct llik_t {
        llik_t(fe_de_separable& m, const std::array<double, n_lambda>& lambda) :
            m_(std::addressof(m)), lambda_(lambda) { }
        llik_t(fe_de_separable& m, const std::array<double, n_lambda>& lambda, double tol) :
            m_(std::addressof(m)), lambda_(lambda), tol_(tol) { }
        // penalized negative log-likelihood at point
        double operator()(const vector_t& g) {
            return -(m_->Psi_ * g).sum() + m_->n_obs_ * m_->int_exp_(g) + lambda_[0] * g.dot(m_->PD_ * g) +
                   lambda_[1] * g.dot(m_->PT_ * g);
        }
        // gradient functor
        std::function<vector_t(const vector_t&)> gradient() {
            return [this, dllik = vector_t(-m_->Psi_.transpose() * vector_t::Ones(m_->n_obs_))](const vector_t& g) {
                return vector_t(
                  dllik + m_->n_obs_ * m_->grad_int_exp_(g) + 2 * (lambda_[0] * m_->PD_ + lambda_[1] * m_->PT_) * g);
            };
        }
        // injected optimization stopping criterion
        template <typename Optimizer> bool stop_if(Optimizer& opt) {
            double llik_old = -(m_->Psi_ * opt.x_old).sum() + m_->n_obs_ * m_->int_exp_(opt.x_old);
            double llik_new = -(m_->Psi_ * opt.x_new).sum() + m_->n_obs_ * m_->int_exp_(opt.x_new);
            if (std::abs((llik_new - llik_old) / llik_old) > tol_) { return false; }
            double penD_old = opt.x_old.dot(m_->PD_ * opt.x_old), penT_old = opt.x_old.dot(m_->PT_ * opt.x_old);
            double penD_new = opt.x_new.dot(m_->PD_ * opt.x_new), penT_new = opt.x_new.dot(m_->PT_ * opt.x_new);
            if (
              std::abs((penD_new - penD_old) / penD_old) > tol_ || std::abs((penT_new - penT_old) / penT_old) > tol_) {
                return false;
            }
            double loss_old = llik_old + lambda_[0] * penD_old + lambda_[1] * penT_old;
            double loss_new = llik_new + lambda_[0] * penD_new + lambda_[1] * penT_new;
            return std::abs((loss_new - loss_old) / loss_old) < tol_;
        }
       private:
        fe_de_separable* m_;
        std::array<double, n_lambda> lambda_;
        double tol_ = 1e-5;
    };

    fe_de_separable() noexcept = default;
    template <typename GeoFrame, typename Penalty>
        requires(is_valid_penalty_v<Penalty>)
    fe_de_separable(const GeoFrame& gf, Penalty&& penalty) {
        fdapde_static_assert(GeoFrame::Order == 2, THIS_CLASS_IS_FOR_ORDER_TWO_GEOFRAMES_ONLY);
        discretize(penalty);
        analyze_data(gf);
    }

    template <typename Penalty> void discretize(Penalty&& penalty) {
        using LhsPenalty = typename std::decay_t<Penalty>::LhsPenalty;
        using LhsBilinearForm = typename LhsPenalty::BilinearForm;
        using LhsLinearForm = typename LhsPenalty::LinearForm;
        fdapde_static_assert(
          internals::is_valid_penalty_pair_v<LhsBilinearForm FDAPDE_COMMA LhsLinearForm>, INVALID_PENALTY_DESCRIPTION);
        const LhsPenalty& penalty1 = penalty.lhs_penalty();
        using RhsPenalty = typename std::decay_t<Penalty>::RhsPenalty;
        using RhsBilinearForm = typename RhsPenalty::BilinearForm;
        using RhsLinearForm = typename RhsPenalty::LinearForm;
        fdapde_static_assert(
          internals::is_valid_penalty_pair_v<RhsBilinearForm FDAPDE_COMMA RhsLinearForm>, INVALID_PENALTY_DESCRIPTION);
        const RhsPenalty& penalty2 = penalty.rhs_penalty();

        using FunctionSpaces = typename function_space_tuple<std::tuple<LhsBilinearForm, RhsBilinearForm>>::type;
        using FS1 = std::tuple_element_t<0, FunctionSpaces>;
        using FS2 = std::tuple_element_t<1, FunctionSpaces>;
        // one penalty must be on a FeSpace
        fdapde_static_assert(is_fe_space_v<FS1> || is_fe_space_v<FS2>, NO_FINITE_ELEMENT_SPACE_DETECTED);
        constexpr int fe_space_index = is_fe_space_v<FS1> ? 0 : 1;
        using FeSpace = std::tuple_element_t<fe_space_index, FunctionSpaces>;
        constexpr int bs_space_index = is_fe_space_v<FS1> ? 1 : 0;
        using BsSpace = std::tuple_element_t<bs_space_index, FunctionSpaces>;
        // enforce a space-time (or SpaceMajor) expansion: index 0 refer to the spatial finite element discretization
        // we enforce a space-time (or SpaceMajor) expansion of the field by reordering the forms so that, index 0
        // always refer to the spatial finite element discretization
        const auto& fe_penalty = fe_penalty_(penalty1, penalty2);
        const auto& bs_penalty = bs_penalty_(penalty1, penalty2);
        // get references to bilinear and linear forms
        auto bilinear_form = std::tie(fe_penalty.bilinear_form(), bs_penalty.bilinear_form());
        auto linear_form = std::tie(fe_penalty.linear_form(), bs_penalty.linear_form());
	// function spaces
        const FeSpace& fe_space = std::get<0>(bilinear_form).trial_space();
        const BsSpace& bs_space = std::get<1>(bilinear_form).trial_space();
        fdapde_assert(bs_space.sobolev_regularity() > 1);
        // discretization
        auto assemble_ = [&, this]<int Index>() {
            auto& func_space = std::get<Index>(bilinear_form).trial_space();
            // assemble mass matrix
            TrialFunction u(func_space);
            TestFunction  v(func_space);
            R0__[Index] = integral(func_space.triangulation())(u * v).assemble();
            R1__[Index] = std::get<Index>(bilinear_form).assemble();
        };
        assemble_.template operator()<0>();
        assemble_.template operator()<1>();
	// number of basis functions on physical domain
        n_dofs__[0] = std::get<0>(bilinear_form).trial_space().n_dofs();
        n_dofs__[1] = std::get<1>(bilinear_form).trial_space().n_dofs();
        n_dofs_ = n_dofs__[0] * n_dofs__[1];
        u_.resize(n_dofs_);
        {
            vector_t u = std::get<0>(linear_form).assemble();
            for (int i = 0; i < n_dofs__[1]; ++i) { u_.segment(i * n_dofs__[0], n_dofs__[0]) = u; }
        }
	// penalty matrix
        PT_ = kronecker(R1__[1], R0__[0]);
	sparse_solver_t invR0;
	invR0.compute(R0__[0]);
        PD_ = kronecker(R0__[1], R1__[0].transpose() * invR0.solve(R1__[0]));
	
        // store handle for basis system evaluation at locations
        internals::for_each_index_in_pack<2>([&]<int Ns>() {
            point_eval_[Ns] = [fe_space =
                                 std::get<Ns>(bilinear_form).trial_space()](const matrix_t& locs) -> decltype(auto) {
                return internals::point_basis_eval(fe_space, locs);
            };
        });

        // geometry
        using FeTriangulation = typename FeSpace::Triangulation;
        constexpr int fe_embed_dim = FeTriangulation::embed_dim;
        const auto& D = fe_space.triangulation();
        const auto& T = bs_space.triangulation();
        // eval reference basis at quadrature nodes, store de_quadrature weights
        de_fe_quadrature_t<fe_embed_dim> fe_quad_rule;
        de_bs_quadrature_t bs_quad_rule;
        std::vector<Eigen::Matrix<double, Dynamic, Dynamic>> PsiQuad;   // \psi_i(q_p) \kron \phi_j(q_t), t = 1, ..., m
        Eigen::Matrix<double, Dynamic, 1> w;                            // quadrature weights
        {
            PsiQuad.resize(T.n_cells());
            matrix_t PsiQuad_ = eval_fe_shape_values_at_quadrature_(fe_space, fe_quad_rule);
	    // integration in time
            for (auto it = T.cells_begin(); it != T.cells_end(); ++it) {
                matrix_t PhiQuad = eval_bs_shape_values_at_quadrature_(bs_space, bs_quad_rule, it);
                PsiQuad[it->id()] = kronecker(PhiQuad, PsiQuad_);   // tensorize
            }
            w = kronecker(bs_quad_rule.weights, fe_quad_rule.weights).as_eigen_matrix();
        }
        // store handle for approximation of \int_T \int_D (e^g)
        int_exp_ = [&, PsiQuad, w, Vh = TpSpace(fe_space, bs_space)](const vector_t& g) {
            double result = 0;
            for (auto jt = T.cells_begin(); jt != T.cells_end(); ++jt) {
                for (auto it = D.cells_begin(); it != D.cells_end(); ++it) {
                    result += w.dot((PsiQuad[jt->id()] * g(Vh.dof_handler().active_dofs(it->id(), jt->id())))
                                      .array().exp().matrix()) *
                              it->measure() * (0.5 * jt->measure());
                }	
            }
            return result;
        };
        // store handle for computation of \nabla_g(\int_T \int_D (e^g))
        grad_int_exp_ = [&, PsiQuad, w, Vh = TpSpace(fe_space, bs_space)](const vector_t& g) {
            vector_t grad = vector_t::Zero(g.rows());
            for (auto jt = T.cells_begin(); jt != T.cells_end(); ++jt) {
                for (auto it = D.cells_begin(); it != D.cells_end(); ++it) {
                    std::vector<int> dofs = Vh.dof_handler().active_dofs(it->id(), jt->id());
                    grad(dofs) += PsiQuad[jt->id()].transpose() *
                                  ((PsiQuad[jt->id()] * g(dofs)).array().exp()).cwiseProduct(w.array()).matrix() *
                                  it->measure() * (0.5 * jt->measure());
                }
            }
            return grad;
        };
        return;
    }
    // fit from geoframe
    template <typename GeoFrame> void analyze_data(const GeoFrame& gf) {
        fdapde_static_assert(GeoFrame::Order == 2, THIS_CLASS_IS_FOR_ORDER_TWO_GEOFRAMES_ONLY);
        fdapde_assert(gf.n_layers() == 1 && gf[0].category()[0] == ltype::point && gf[0].category()[1] == ltype::point);
        n_obs_ = gf[0].rows();

        std::array<sparse_matrix_t, 2> Psi__;
        internals::for_each_index_in_pack<n_lambda>([&]<int Ns>() {
            // eval physical basis at spatial locations
            const auto& spatial_index = geo_index_cast<Ns, POINT>(gf[0]);
            n_locs__[Ns] = spatial_index.rows();
            if (spatial_index.points_at_dofs()) {
                Psi__[Ns].resize(n_locs__[Ns], n_dofs__[Ns]);
                Psi__[Ns].setIdentity();
            } else {
                Psi__[Ns] = point_eval_[Ns](spatial_index.coordinates());
            }
        });
        // tensorize physical basis
        Psi_.resize(n_locs__[1], n_dofs_);
        std::vector<Triplet<double>> triplet_list;
        for (int i = 0; i < n_locs__[1]; ++i) {
            // kronecker product between Psi__[1] i-th row and Psi__[0] i-th row
            sparse_matrix_t tmp = kronecker(sparse_matrix_t(Psi__[1].row(i)), sparse_matrix_t(Psi__[0].row(i)));
            for (int j = 0; j < tmp.outerSize(); ++j) {
                for (sparse_matrix_t::InnerIterator it(tmp, j); it; ++it) {
                    triplet_list.emplace_back(i, it.col(), it.value());
                }
            }
        }
        Psi_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        Psi_.makeCompressed();
        return;
    }
    // main fit entry point
    template <typename Optimizer>
    const vector_t& fit(double lambda_D, double lambda_T, const vector_t& g_init, Optimizer&& opt) {
        g_ = opt.optimize(llik_t(*this, std::array {lambda_D, lambda_T}, tol_), g_init);
        return g_;
    }
    template <typename Optimizer, typename LambdaT>
        requires(internals::is_vector_like_v<LambdaT>)
    const vector_t& fit(LambdaT&& lambda, const vector_t& g_init, Optimizer&& opt) {
        fdapde_assert(lambda.size() == n_lambda);
        return fit(lambda[0], lambda[1]);
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
    // not tensorized quantites
    std::array<int, 2> n_dofs__;           // number of spatial and temporal degrees of freedom {n_dofs_D, n_dofs_T}
    std::array<int, 2> n_locs__;           // number of spatial and temporal data locations {n_locs_D, n_locs_T}
    std::array<sparse_matrix_t, 2> R0__;   // {R0_D, R0_T} = { \int_D \psi_i * \psi_j, \int_T \phi_i * \phi_j }
    std::array<sparse_matrix_t, 2> R1__;   // {R1_D, R1_T} = { \int_D a_D(\psi_i, \psi_j), \int_T a_T(\phi_T, \phi_D) }

    sparse_matrix_t R0_;    // n_dofs x n_dofs matrix R0 = R0_T \kron R0_D
    sparse_matrix_t R1_;    // n_dofs x n_dofs matrix R1 = R0_T \kron R1_D
    sparse_matrix_t Psi_;   // n_obs x n_dofs matrix Psi = Psi_T \kron Psi_D
    vector_t u_;            // (n_dofs_D * n_dofs_T) x 1 vector u = [u_1 \ldots u_n, \ldots, u_1 \ldots u_n]

    matrix_t PD_;                                             // matrix PD = R0_T \kron (R1_D^\top * R0_D^{-1} * R1_D)
    matrix_t PT_;                                             // matrix PT = R1_T \kron R0_D
    std::function<double(const vector_t&)> int_exp_;          // functor computing \int exp(g)
    std::function<vector_t(const vector_t&)> grad_int_exp_;   // functor computing \nabla_g \int exp(g)
    vector_t g_;
    // basis system evaluation handle
    std::array<std::function<sparse_matrix_t(const matrix_t& locs)>, 2> point_eval_;
    double tol_ = 1e-5;                                              // tolerance for custom stopping criterion
};

}   // namespace internals

// separable solver API
template <typename... Penalty>
    requires(sizeof...(Penalty) == 2 && (internals::is_pair_v<Penalty> && ...))
struct fe_de_separable {
    using solver_t = internals::fe_de_separable;
   private:
    struct penalty_packet {
        template <typename Penalty_> struct penalty_bit {
            using BilinearForm = std::tuple_element_t<0, std::decay_t<Penalty_>>;
            using LinearForm = std::tuple_element_t<1, std::decay_t<Penalty_>>;
           private:
            BilinearForm bilinear_form_;
            LinearForm linear_form_;
           public:
            penalty_bit(const Penalty_& penalty) :
                bilinear_form_(std::get<0>(penalty)), linear_form_(std::get<1>(penalty)) { }
            // observers
            const BilinearForm& bilinear_form() const { return bilinear_form_; }
            const LinearForm& linear_form() const { return linear_form_; }
        };
        using LhsPenalty = penalty_bit<std::tuple_element_t<0, std::tuple<Penalty...>>>;
        using RhsPenalty = penalty_bit<std::tuple_element_t<1, std::tuple<Penalty...>>>;

        penalty_packet(const Penalty&... penalty) : penalty_(penalty...) { }
        // observers
        const auto& lhs_penalty() const { return std::get<0>(penalty_); }
        const auto& rhs_penalty() const { return std::get<1>(penalty_); }
       private:
        std::tuple<penalty_bit<Penalty>...> penalty_;
    };
   public:
    fe_de_separable(const Penalty&... penalty) : penalty_(penalty...) { }
    const penalty_packet& get() const { return penalty_; }
   private:
    penalty_packet penalty_;
};

}   // namespace fdapde

#endif   // __FE_DE_SEPARABLE_SOLVER_H__

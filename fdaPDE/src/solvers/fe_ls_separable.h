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

#ifndef __FE_LS_SEPARABLE_SOLVER_H__
#define __FE_LS_SEPARABLE_SOLVER_H__

#include "header_check.h"

namespace fdapde {
namespace internals {

// solves \min_{f, \beta} \| W^{1/2} * (y_i - x_i^\top * \beta - f(p_i, t_j)) \|_2^2 + \int_D \int_T (L_D(f) - u_D)^2 +
// \int_T \int_D (L_T(f) - u_T)^2
class fe_ls_separable_mono {
   private:
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using binary_t = BinaryMatrix<Dynamic, Dynamic>;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using diag_matrix_t   = Eigen::DiagonalMatrix<double, Dynamic, Dynamic>;
    using sparse_solver_t = eigen_sparse_solver_movable_wrap<Eigen::SparseLU<sparse_matrix_t>>;
    using dense_solver_t  = Eigen::PartialPivLU<matrix_t>;
    template <typename DataLocs>
    static constexpr bool is_valid_data_locs_descriptor_v =
      std::is_same_v<DataLocs, matrix_t> || std::is_same_v<DataLocs, binary_t>;
    template <typename InfoT> struct is_valid_info_t {
        static constexpr bool value = requires(InfoT info) { info.penalty; };
    };

    template <typename Tuple> struct function_space_tuple {
        using type = decltype([]<size_t... Is_>(std::index_sequence<Is_...>) {
            return std::make_tuple(typename std::tuple_element_t<Is_, Tuple>::TrialSpace {}...);
        }(std::make_index_sequence<std::tuple_size_v<Tuple>>()));
    };
    template <typename Penalty1, typename Penalty2>
    const auto& fe_penalty_(const Penalty1& penalty1, const Penalty2& penalty2) const {
        return select_one_between(
          penalty1, penalty2, []() { return  is_fe_space_v<typename std::tuple_element_t<0, Penalty1>::TrialSpace>; });
    }
    template <typename Penalty1, typename Penalty2>
    const auto& bs_penalty_(const Penalty1& penalty1, const Penalty2& penalty2) const {
        return select_one_between(
          penalty1, penalty2, []() { return !is_fe_space_v<typename std::tuple_element_t<0, Penalty1>::TrialSpace>; });
    }
   public:
    static constexpr int n_lambda = 2;
    using solver_category = ls_solver;
   private:  
    // evaluation of basis system at spatial locations
    template <typename... DataLocs>
        requires((is_valid_data_locs_descriptor_v<DataLocs> && ...) && (sizeof...(DataLocs) == n_lambda))
    void eval_basis_at_(const DataLocs&... locs) {
        std::array<sparse_matrix_t, n_lambda> Psi__;
        internals::for_each_index_and_args<n_lambda>(
          [&]<int Ns, typename Ts>(const Ts& locs) {
              if constexpr (std::is_same_v<Ts, matrix_t>) {   // pointwise sampling
                  Psi__[Ns] = point_eval_[Ns](locs);
              } else {   // areal sampling
                  const auto& [psi, measure_vect] = areal_eval_[Ns](locs);
                  Psi__[Ns] = psi;
              }
          },
          locs...);
        Psi_ = kronecker(Psi__[1], Psi__[0]);
        D_ = vector_t::Ones(n_locs_).asDiagonal();
        return;
    }
    // optimized basis evaluation at geoframe
    template <typename GeoFrame> void eval_basis_at_(const GeoFrame& gf) {
        std::array<sparse_matrix_t, 2> Psi__;
        internals::for_each_index_in_pack<n_lambda>([&]<int Ns>() {
            switch (gf.category(0)[Ns]) {
            case ltype::point: {
                const auto& spatial_index = geo_index_cast<Ns, POINT>(gf[0]);
                if (spatial_index.points_at_dofs()) {
                    Psi__[Ns].resize(n_locs_, n_dofs_);
                    Psi__[Ns].setIdentity();
                } else {
                    Psi__[Ns] = point_eval_[Ns](spatial_index.coordinates());
                }
                D_ = vector_t::Ones(n_locs_).asDiagonal();
                break;
            }
            case ltype::areal: {
                const auto& spatial_index = geo_index_cast<Ns, POLYGON>(gf[0]);
                const auto& [psi, measure_vec] = areal_eval_[Ns](spatial_index.incidence_matrix());
                Psi__[Ns] = psi;
                vector_t D(n_locs_);
		int m = n_locs_ / spatial_index.rows();
                for (int i = 0; i < m; ++i) { D.segment(i * measure_vec.rows(), measure_vec.rows()) = measure_vec; }
                D_ = D.asDiagonal();
                break;
            }
            }
        });
        Psi_ = kronecker(Psi__[1], Psi__[0]);
        return;
    }
    // unrolls the penalty tuple and injects them into discretize()
    template <typename Penalty> void discretize_loop_(Penalty&& penalty) {
        internals::apply_index_pack<n_lambda>(
          [&]<int... Ns_>() { discretize([&]() { return std::get<Ns_>(penalty); }()...); });
        return;
    }
   public:
    fe_ls_separable_mono() noexcept = default;
    // construct from formula + geoframe
    template <typename GeoFrame, typename WeightMatrix, typename InfoT>
        requires(is_valid_info_t<InfoT>::value)
    fe_ls_separable_mono(const std::string& formula, const GeoFrame& gf, InfoT&& info, const WeightMatrix& W) {
        fdapde_static_assert(GeoFrame::Order == 2, THIS_CLASS_IS_FOR_ORDER_TWO_GEOFRAMES_ONLY);
	fdapde_assert(gf.n_layers() == 1);
        n_obs_  = gf[0].rows();
        n_locs_ = n_obs_;

        discretize_loop_(info.penalty);
        analyze_data(formula, gf, W);
    }
    template <typename GeoFrame, typename InfoT>
        requires(is_valid_info_t<InfoT>::value)
    fe_ls_separable_mono(const std::string& formula, const GeoFrame& gf, InfoT&& info) :
        fe_ls_separable_mono(formula, gf, info, vector_t::Ones(gf[0].rows()).asDiagonal()) { }
    // construct with no data
    template <typename GeoFrame, typename InfoT, typename WeightMatrix>
        requires(is_valid_info_t<InfoT>::value)
    fe_ls_separable_mono(const GeoFrame& gf, InfoT&& info, const WeightMatrix& W) : W_(W) {
        fdapde_static_assert(GeoFrame::Order == 2, THIS_CLASS_IS_FOR_ORDER_TWO_GEOFRAMES_ONLY);
	fdapde_assert(gf.n_layers() == 1);
        n_obs_  = gf[0].rows();
	n_locs_ = n_obs_;

	discretize_loop_(info.penalty);
	eval_basis_at_(gf);
    }
    template <typename GeoFrame, typename InfoT>
        requires(is_valid_info_t<InfoT>::value)
    fe_ls_separable_mono(const GeoFrame& gf, InfoT&& info) :
        fe_ls_separable_mono(gf, info.penalty, vector_t::Ones(gf[0].rows()).asDiagonal()) { }

    // numerical discretization
    template <typename Penalty1, typename Penalty2> void discretize(Penalty1&& penalty1, Penalty2&& penalty2) {
        fdapde_static_assert(
          internals::is_valid_penalty_pair_v<Penalty1> && internals::is_valid_penalty_pair_v<Penalty2>,
          INVALID_PENALTY_DESCRIPTION);
        using BilinearForms =
          std::tuple<std::tuple_element_t<0, std::decay_t<Penalty1>>, std::tuple_element_t<0, std::decay_t<Penalty2>>>;
        using FunctionSpaces = typename function_space_tuple<BilinearForms>::type;
        using FS1 = std::tuple_element_t<0, FunctionSpaces>;
        using FS2 = std::tuple_element_t<1, FunctionSpaces>;
        // one penalty must be on a FeSpace
        fdapde_static_assert(is_fe_space_v<FS1> || is_fe_space_v<FS2>, NO_FINITE_ELEMENT_SPACE_DETECTED);
        constexpr int fe_space_index = is_fe_space_v<FS1> ? 0 : 1;
        constexpr int bs_space_index = is_fe_space_v<FS1> ? 1 : 0;
        using BsSpace = std::tuple_element_t<bs_space_index, FunctionSpaces>;
        // we enforce a space-time (or SpaceMajor) expansion of the field by reordering the forms so that, index 0
        // always refer to the spatial finite element discretization
        const auto& fe_penalty = fe_penalty_(penalty1, penalty2);
        const auto& bs_penalty = bs_penalty_(penalty1, penalty2);
        // get references to bilinear and linear forms
        auto bilinear_form = std::tie(std::get<0>(fe_penalty), std::get<0>(bs_penalty));
        auto linear_form   = std::tie(std::get<1>(fe_penalty), std::get<1>(bs_penalty));
        {
            const BsSpace& bs_space = std::get<bs_space_index>(bilinear_form).trial_space();
            fdapde_assert(bs_space.sobolev_regularity() > 1);
        }
        // discretization
        auto assemble_ = [&, this]<int Index>() {
            auto& space = std::get<Index>(bilinear_form).trial_space();
            // assemble mass matrix
            TrialFunction u(space);
            TestFunction  v(space);
            R0__[Index] = integral(space.triangulation())(u * v).assemble();
            R1__[Index] = std::get<Index>(bilinear_form).assemble();
        };
        assemble_.template operator()<0>();
        assemble_.template operator()<1>();	
        // tensorization
        R0_ = kronecker(R0__[1], R0__[0]);   // R0_T \kron R0_D
        R1_ = kronecker(R0__[1], R1__[0]);   // R0_T \kron R1_D
        K_  = kronecker(R1__[1], R0__[0]);   // R1_T \kron R0_D	
        // number of basis functions on physical domain
        n_dofs__[0] = std::get<0>(bilinear_form).trial_space().n_dofs();
        n_dofs__[1] = std::get<1>(bilinear_form).trial_space().n_dofs();
        n_dofs_ = n_dofs__[0] * n_dofs__[1];
        // forcing discretization
        u_.resize(n_dofs_);
        {
            vector_t u = std::get<0>(linear_form).assemble();
            for (int i = 0; i < n_dofs__[1]; ++i) { u_.segment(i * n_dofs__[0], n_dofs__[0]) = u; }
        }
        // store handlers for basis system evaluation at locations
        internals::for_each_index_in_pack<2>(
          [&]<int Ns>() {
              point_eval_[Ns] = [fe_space =
                                   std::get<Ns>(bilinear_form).trial_space()](const matrix_t& locs) -> decltype(auto) {
                  return internals::point_basis_eval(fe_space, locs);
              };
              areal_eval_[Ns] = [fe_space =
                                   std::get<Ns>(bilinear_form).trial_space()](const binary_t& locs) -> decltype(auto) {
                  return internals::areal_basis_eval(fe_space, locs);
              };
          });
        b_.resize(2 * n_dofs_, 1);
        return;
    }
    // non-parametric fit
    // \sum_i w_i * (y_i - f(p_i))^2 + \int_D (Lf - u)^2
    template <typename DataLocs1, typename DataLocs2, typename WeightMatrix>
        requires(is_valid_data_locs_descriptor_v<DataLocs1> && is_valid_data_locs_descriptor_v<DataLocs2>)
    void analyze_data(const DataLocs1& locs1, const DataLocs2& locs2, const matrix_t& y, const WeightMatrix& W) {
        fdapde_assert(
          locs1.rows() > 0 && locs2.rows() > 0 && y.rows() == locs1.rows() * locs2.rows() && y.cols() == 1 &&
          W.rows() == locs1.rows() * locs2.rows() && W.rows() == W.cols());
        n_obs_  = y.rows();
	n_locs_ = n_obs_;
        n_covs_ = 0;
        eval_basis_at_(locs1, locs2);   // update \Psi matrix
        update_response_and_weights(y, W);
        return;
    }
    // semi-parametric fit
    // \sum_i w_i * (y_i - x_i^\top * \beta - f(p_i))^2 + \int_D (Lf - u)^2
    template <typename DataLocs1, typename DataLocs2, typename WeightMatrix>
        requires(is_valid_data_locs_descriptor_v<DataLocs1> && is_valid_data_locs_descriptor_v<DataLocs2>)
    void analyze_data(
      const DataLocs1& locs1, const DataLocs2& locs2, const matrix_t& y, const matrix_t& X, const WeightMatrix& W) {
        fdapde_assert(
          locs1.rows() > 0 && locs2.rows() > 0 && y.rows() == locs1.rows() * locs2.rows() && y.cols() == 1 &&
          X.rows() == locs1.rows() * locs2.rows() && W.rows() == locs1.rows() * locs2.rows() && W.rows() == W.cols());
        n_obs_  = y.rows();
	n_locs_ = n_obs_;
        bool require_woodbury_realloc = std::cmp_not_equal(n_covs_, X.cols());
        n_covs_ = X.cols();
        eval_basis_at_(locs1, locs2);   // update \Psi matrix
        if (require_woodbury_realloc) { U_ = matrix_t::Zero(2 * n_dofs_, n_covs_); }
        if (require_woodbury_realloc) { V_ = matrix_t::Zero(n_covs_, 2 * n_dofs_); }
        update_response_and_weights(y, X, W);
        return;
    }
    // fit from formula
    template <typename GeoFrame, typename WeightMatrix>
    void analyze_data(const std::string& formula, const GeoFrame& gf, const WeightMatrix& W) {
        fdapde_static_assert(GeoFrame::Order == 2, THIS_CLASS_IS_FOR_ORDER_TWO_GEOFRAMES_ONLY);
        fdapde_assert(gf.n_layers() == 1);
        n_obs_  = gf[0].rows();
        n_locs_ = n_obs_;
        eval_basis_at_(gf);   // update \Psi matrix

        // parse formula, extract response vector and design matrix
        Formula formula_(formula);
        std::vector<std::string> covs;
        for (const std::string& token : formula_.rhs()) {
            if (gf.contains(token)) { covs.push_back(token); }
        }
        bool require_woodbury_realloc = std::cmp_not_equal(n_covs_, covs.size());
        n_covs_ = covs.size();
        const auto& y_data = gf[0].data().template col<double>(formula_.lhs());
        y_.resize(n_locs_, y_data.blk_sz());
        y_data.assign_to(y_);

        if (b_.cols() != y_.cols()) { b_.resize(2 * n_dofs_, y_.cols()); }
        if (n_covs_ != 0) {
            if (require_woodbury_realloc) { U_ = matrix_t::Zero(2 * n_dofs_, n_covs_); }
            if (require_woodbury_realloc) { V_ = matrix_t::Zero(n_covs_, 2 * n_dofs_); }
            X_.resize(n_locs_, n_covs_);   // assemble design matrix
            for (int i = 0; i < n_covs_; ++i) { gf[0].data().template col<double>(covs[i]).assign_to(X_.col(i)); }
        }
        update_response_and_weights(y_, W);   // updates design_matrix releated matrices as well
        return;
    }

    // modifiers
    void update_response(const vector_t& y) {
        fdapde_assert(Psi_.rows() > 0 && y.rows() == Psi_.rows() && y.cols() == 1);
        y_ = y;
        // correct \Psi for missing observations
        auto nan_pattern = na_matrix(y);
        int old_n_obs = n_obs_;
        if (nan_pattern.any()) {
            n_obs_ = n_locs_ - nan_pattern.count();
            B_ = (~nan_pattern).repeat(1, n_dofs_).select(Psi_, 0);
            y_ = (~nan_pattern).select(y_, 0);
        }
        if (old_n_obs != n_obs_) { W_ *= (double)old_n_obs / n_obs_; }
        b_.block(0, 0, n_dofs_, 1) = -PsiNA().transpose() * D_ * W_ * y;
        return;
    }
    template <typename WeightMatrix> void update_weights(const WeightMatrix& W) {
        fdapde_assert(Psi_.rows() > 0 && W.rows() == n_locs_ && W.rows() == W.cols());
        W_ = W;
	W_ /= n_obs_;
        if (n_covs_ == 0) {
            b_.block(0, 0, n_dofs_, 1) = -PsiNA().transpose() * D_ * W_ * y_;
        } else {
            XtWX_ = X_.transpose() * W_ * X_;
            invXtWX_ = XtWX_.partialPivLu();
            invXtWXXtW_ = invXtWX_.solve(X_.transpose() * W_);   // (X^\top * W * X)^{-1} * (X^\top * W)
            // woodbury decomposition matrices
            U_.block(0, 0, n_dofs_, n_covs_) = PsiNA().transpose() * D_ * W_ * X_;
            V_.block(0, 0, n_covs_, n_dofs_) = X_.transpose() * W_ * PsiNA();
            b_.block(0, 0, n_dofs_, 1) = -PsiNA().transpose() * D_ * internals::lmbQ(W_, X_, invXtWX_, y_);
        }
	W_changed_ = true;
	return;
    }
    template <typename WeightMatrix> void update_response_and_weights(const vector_t& y, const WeightMatrix& W) {
        fdapde_assert(
          Psi_.rows() > 0 && y.rows() == n_locs_ && y.cols() == 1 && W.rows() == W.cols() && W.rows() == n_locs_);
        y_ = y;
        // correct \Psi for missing observations
        auto nan_pattern = na_matrix(y);
        if (nan_pattern.any()) {
            n_obs_ = n_locs_ - nan_pattern.count();
            B_ = (~nan_pattern).repeat(1, n_dofs_).select(Psi_, 0);
            y_ = (~nan_pattern).select(y_, 0);
        }
        update_weights(W);
        return;
    }

    // main fit entry point
    std::pair<vector_t, vector_t> fit(double lambda_D, double lambda_T) {
        fdapde_assert(lambda_D > 0 && lambda_T > 0 && n_dofs_ > 0 && n_obs_ > 0);
	std::array<double, n_lambda> lambda {lambda_D, lambda_T};
        if (lambda_saved_.value() != lambda || W_changed_) {
            // assemble and factorize system matrix for nonparameteric part
            SparseBlockMatrix<double, 2, 2> A(
              -PsiNA().transpose() * D_ * W_ * PsiNA() - lambda_T * K_, lambda_D * R1_.transpose(), lambda_D * R1_,
              lambda_D * R0_);
            invA_.compute(A);
            W_changed_ = false;
        }
        if (lambda_saved_.value() != lambda) {
            // linear system rhs
            b_.block(n_dofs_, 0, n_dofs_, 1) = lambda_D * u_;
        }
	lambda_saved_ = lambda;
        vector_t x;
        if (n_covs_ == 0) { 
            x = invA_.solve(b_);
            f_ = x.topRows(n_dofs_);
        } else {   // parametric case
            x = woodbury_system_solve(invA_, U_, XtWX_, V_, b_);
            f_ = x.topRows(n_dofs_);
            beta_ = invXtWXXtW_ * (y_ - Psi_ * f_);
        } 
        g_ = x.bottomRows(n_dofs_);
        return std::make_pair(f_, beta_);
    }
    template <typename LambdaT>
        requires(internals::is_vector_like_v<LambdaT>)
    std::pair<vector_t, vector_t> fit(LambdaT&& lambda) {
        fdapde_assert(lambda.size() == n_lambda);
        return fit(lambda[0], lambda[1]);
    }
    // perform a nonparametric_fit, e.g. discarding possible covariates
    vector_t nonparametric_fit(double lambda_D, double lambda_T) {
        fdapde_assert(lambda_D > 0 && lambda_T > 0 && n_dofs_ > 0 && n_obs_ > 0);
        std::array<double, n_lambda> lambda {lambda_D, lambda_T};
        if (lambda_saved_.value() != lambda) {
            // assemble and factorize system matrix for nonparameteric part
            SparseBlockMatrix<double, 2, 2> A(
              -PsiNA().transpose() * D_ * W_ * PsiNA() - lambda_T * K_, lambda_D * R1_.transpose(), lambda_D * R1_,
              lambda_D * R0_);
            invA_.compute(A);
        }
        vector_t x;
        if (n_covs_ == 0) {   // equivalent to calling fit(lambda)
            if (lambda_saved_.value() != lambda) { b_.block(n_dofs_, 0, n_dofs_, 1) = lambda_D * u_; }
            x = invA_.solve(b_);
        } else {
            vector_t b(2 * n_dofs_);
            // assemble nonparametric linear system rhs
            b.block(0, 0, n_dofs_, 1) = -PsiNA().transpose() * D_ * W_ * y_;
            b.block(n_dofs_, 0, n_dofs_, 1) = lambda_D * u_;   
            x = invA_.solve(b);
        }
        lambda_saved_ = lambda;
        f_ = x.topRows(n_dofs_);
        g_ = x.bottomRows(n_dofs_);
        return f_;
    }

    // hutchinson approximation for Tr[S]
    double edf(int r = 100, int seed = random_seed) {
        fdapde_assert(lambda_saved_.has_value());
        if (!Ys_.has_value() || !Bs_.has_value()) {
            int seed_ = (seed == random_seed) ? std::random_device()() : seed;
            std::mt19937 rng(seed_);
            rademacher_distribution rademacher;
            Us_->resize(n_locs_, r);
            for (int i = 0; i < n_locs_; ++i) {
                for (int j = 0; j < r; ++j) { Us_->operator()(i, j) = rademacher(rng); }
            }
            Ys_ = Us_->transpose() * Psi_;
            Bs_ = matrix_t::Zero(2 * n_dofs_, r);   // implicitly enforce homogeneous forcing
        }
        if (n_covs_ == 0) {
            Bs_->topRows(n_dofs_) = -PsiNA().transpose() * D_ * W_ * (*Us_);
        } else {
            Bs_->topRows(n_dofs_) = -PsiNA().transpose() * D_ * internals::lmbQ(W_, X_, invXtWX_, *Us_);
        }
        matrix_t x = n_covs_ == 0 ? invA_.solve(*Bs_) : woodbury_system_solve(invA_, U_, XtWX_, V_, *Bs_);
        double trS = 0;   // monte carlo Tr[S] approximation
        for (int i = 0; i < r; ++i) { trS += Ys_->row(i).dot(x.col(i).head(n_dofs_)); }
        return trS / r;
    }
    template <typename... LambdaT>
        requires(
          (sizeof...(LambdaT) == 1 && (internals::is_vector_like_v<LambdaT> && ...)) ||
          (sizeof...(LambdaT) == n_lambda && (std::is_floating_point_v<LambdaT> && ...)))
    double edf(const LambdaT&... lambda, int r = 100, int seed = random_seed) {
        std::array<double, n_lambda> lambda_;
        if constexpr (sizeof...(LambdaT) == 1) {
            internals::for_each_index_and_args<sizeof...(LambdaT)>([&]<int Ns_, typename Ts_>(const Ts_& ts) {
                fdapde_assert(ts.size() == n_lambda && ts[0] > 0 && ts[1] > 0);
                lambda_[0] = ts[0];
                lambda_[1] = ts[1];
            });
        } else {
            std::array<double, n_lambda> lambda__ {static_cast<double>(lambda)...};
            fdapde_assert(lambda__[0] > 0 && lambda__[1] > 0);
            lambda_[0] = lambda__[0];
	    lambda_[1] = lambda__[1];
        }
        if (lambda_saved_.value() != lambda_) {
            SparseBlockMatrix<double, 2, 2> A(
              -PsiNA().transpose() * D_ * W_ * PsiNA() - lambda_[1] * K_, lambda_[0] * R1_.transpose(),
              lambda_[0] * R1_, lambda_[0] * R0_);
            invA_.compute(A);
            lambda_saved_ = lambda_;
        }
        return edf(r, seed);
    }
    // penalty matrix: \lambda_D * R0_T \kron (R1_D^\top * R0_D^{-1} * R1_D) + \lambda_T * R1_T \kron R0_D
    matrix_t P(double lambda_D, double lambda_T) const {
        if (!PT_.has_value()) { PT_ = kronecker(R1__[1], R0__[0]); }
        if (!PD_.has_value()) {
            sparse_solver_t invR0;
            invR0.compute(R0__[0]);
            PD_ = kronecker(R0__[1], R1__[0].transpose() * invR0.solve(R1__[0]));
        }
        return lambda_D * (*PD_) + lambda_T * (*PT_);
    }
    template <typename LambdaT>
        requires(internals::is_vector_like_v<LambdaT>)
    matrix_t P(const LambdaT& lambda) const {
        fdapde_assert(lambda.size() == n_lambda);
        return P(lambda[0], lambda[1]);
    }
    matrix_t P() const { return P(1.0, 1.0); }
    double ftPf(double lambda_D, double lambda_T) {
        if (std::array<double, n_lambda> {lambda_D, lambda_T} != lambda_saved_ || W_changed_) {
            fit(lambda_D, lambda_T);
        }
        return f_.dot(P(lambda_D, lambda_T) * f_);
    }
    template <typename LambdaT>
        requires(internals::is_vector_like_v<LambdaT>)
    double ftPf(const LambdaT& lambda) {
        fdapde_assert(lambda.size() == n_lambda);
        return internals::apply_index_pack<n_lambda>([&]<int... Ns>() { return ftPf(lambda[Ns]...); });
    }
    vector_t lmbPsi(const vector_t& rhs) const { return Psi_ * rhs; }
    vector_t fn() const { return Psi_ * f_; }

    // observers
    int n_dofs() const { return n_dofs_; }
    const sparse_matrix_t& mass() const { return R0_; }
    const sparse_matrix_t& stiff() const { return R1_; }
    const sparse_matrix_t& Psi() const { return Psi_; }
    const sparse_matrix_t& PsiNA() const { return B_.has_value() ? *B_ : Psi_; }
    const vector_t& force() const { return u_; }
    const vector_t& f() const { return f_; }
    const vector_t& beta() const { return beta_; }
    const vector_t& misfit() const { return g_; }
    const matrix_t& design_matrix() const { return X_; }
    const vector_t& response() const { return y_; }
   protected:
    std::optional<std::array<double, n_lambda>> lambda_saved_ = std::array<double, n_lambda> {-1, -1};
    sparse_solver_t invA_;
    matrix_t b_;
    // matrices for hutchinson stochastic estimation of Tr[S]
    std::optional<matrix_t> Ys_;
    std::optional<matrix_t> Bs_;
    std::optional<matrix_t> Us_;

    int n_dofs_ = 0, n_locs_ = 0, n_obs_ = 0, n_covs_ = 0;
    // not tensorized quantities
    std::array<int, 2> n_dofs__;           // number of spatial and temporal degrees of freedom {n_dofs_D, n_dofs_T}
    std::array<sparse_matrix_t, 2> R0__;   // {R0_D, R0_T} = { \int_D \psi_i * \psi_j, \int_T \phi_i * \phi_j }
    std::array<sparse_matrix_t, 2> R1__;   // {R1_D, R1_T} = { \int_D a_D(\psi_i, \psi_j), \int_T a_T(\phi_T, \phi_D) }

    sparse_matrix_t R0_;    // n_dofs x n_dofs matrix R0 = R0_T \kron R0_D
    sparse_matrix_t R1_;    // n_dofs x n_dofs matrix R1 = R0_T \kron R1_D
    sparse_matrix_t K_;     // n_dofs x n_dofs matrix  K = R1_T \kron R0_D
    sparse_matrix_t Psi_;   // n_obs x n_dofs matrix Psi = Psi_T \kron Psi_D
    vector_t u_;            // (n_dofs_D * n_dofs_T) x 1 vector u = [u_1 \ldots u_n, \ldots, u_1 \ldots u_n]
    diag_matrix_t D_;       // vector of regions' measures (areal sampling)
    mutable sparse_solver_t invR0_;
    std::optional<sparse_matrix_t> B_;            // \Psi matrix corrected for missing observations
    mutable std::optional<sparse_matrix_t> PD_;   // matrix PD = R0_T \kron (R1_D^\top * R0_D^{-1} * R1_D)
    mutable std::optional<sparse_matrix_t> PT_;   // matrix PT = R1_T \kron R0_D
    vector_t f_, beta_, g_;
    // basis system evaluation handles
    std::array<std::function<sparse_matrix_t(const matrix_t& locs)>, n_lambda> point_eval_;
    std::array<std::function<std::pair<sparse_matrix_t, vector_t>(const binary_t& locs)>, n_lambda> areal_eval_;

    matrix_t X_;               // n_obs x n_covs design matrix
    vector_t y_;               // n_obs x 1 observation vector
    sparse_matrix_t W_;        // n_obs x n_obs matrix of observation weights
    matrix_t U_, V_;           // (2 * n_dofs) x n_covs matrices [\Psi^\top * D * W * y, 0] and [X^\top * W * \Psi, 0]
    matrix_t XtWX_;            // n_covs x n_covs matrix X^\top * W * X
    dense_solver_t invXtWX_;   // factorization of n_covs x n_covs matrix X^\top * W * X
    matrix_t invXtWXXtW_;      // n_covs x n_obs matrix (X^\top * X)^{-1} * (X^\top * W)
    bool W_changed_;
};

}   // namespace internals

// separable monolithic solver method
template <typename... Penalty>
    requires(sizeof...(Penalty) == 2 && (internals::is_valid_penalty_pair_v<Penalty> && ...))
struct fe_ls_separable_mono {
    using solver_t = internals::fe_ls_separable_mono;
   private:
    struct info_t {
        std::tuple<Penalty...> penalty;
    };
   public:
    fe_ls_separable_mono(const Penalty&... penalty) : info_(std::make_tuple(penalty...)) { }
    const info_t& get() const { return info_; }
   private:
    info_t info_;
};

namespace internals {

// central difference time integration loop
class fe_ls_separable_cdti {
   private:
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using binary_t = BinaryMatrix<Dynamic, Dynamic>;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using diag_matrix_t   = Eigen::DiagonalMatrix<double, Dynamic, Dynamic>;
    using sparse_solver_t = eigen_sparse_solver_movable_wrap<Eigen::SparseLU<sparse_matrix_t>>;
    using dense_solver_t  = Eigen::PartialPivLU<matrix_t>;
    template <typename DataLocs>
    static constexpr bool is_valid_data_locs_descriptor_v =
      std::is_same_v<DataLocs, matrix_t> || std::is_same_v<DataLocs, binary_t>;
    template <typename InfoT> struct is_valid_info_t {
        static constexpr bool value = requires(InfoT info) {
            info.penalty;
	    info.max_iter;
	    info.tol;
        };
    };
  
    class block_map_t {
        static constexpr int Order = 3;
        using Scalar = double;
        using storage_t = MdMap<Scalar, full_dynamic_extent_t<Order>, internals::layout_left>;

        storage_t data_;
        int rows_ = 0, cols_ = 0;
        int blk_rows_ = 0, blk_cols_ = 0;   // single block size
       public:
        block_map_t() noexcept = default;
        template <typename DataT>
            requires(internals::is_eigen_dense_xpr_v<DataT> &&
                     std::is_same_v<typename std::decay_t<DataT>::Scalar, double>)
        block_map_t(DataT&& data, int rows, int blk_rows, int blk_cols) :
            data_(data.data(), rows, data.cols(), (data.rows() / rows)),
            rows_(rows),
            cols_(data.cols()),
            blk_rows_(blk_rows),
            blk_cols_(blk_cols) {
            fdapde_assert(data.rows() % rows == 0 && rows % blk_rows == 0 && data.cols() % blk_cols == 0);
        }
        template <typename DataT>
            requires(internals::is_eigen_dense_xpr_v<DataT> &&
                     std::is_same_v<typename std::decay_t<DataT>::Scalar, double>)
        block_map_t(DataT&& data, int rows) :   // divide data in (data.rows() / rows) blocks of size rows x data.cols()
            data_(data.data(), rows, data.cols(), (data.rows() / rows)),
            rows_(rows),
	    cols_(data.cols()),
            blk_rows_(rows),
            blk_cols_(data.cols()) {
            fdapde_assert(data.rows() % rows == 0);
        }
        block_map_t(const block_map_t& other) :
            rows_(other.rows_), cols_(other.cols_), blk_rows_(other.blk_rows_), blk_cols_(other.blk_cols_) {
            for (size_t i = 0; i < data_.size(); ++i) { data_.data()[i] = other.data_.data()[i]; }
        }
        block_map_t& operator=(const block_map_t& other) {
            for (size_t i = 0; i < data_.size(); ++i) { data_.data()[i] = other.data_.data()[i]; }
            rows_ = other.rows_;
            cols_ = other.cols_;
            blk_rows_ = other.blk_rows_;
            blk_cols_ = other.blk_cols_;
            return *this;
        }
        // observers
        auto operator()(int i, int k) const {   // get i-th row-block of k-th time instant
            auto slice_ = data_.template slice<2>(k);
            return slice_.as_eigen_map().block(i * blk_rows_, 0, blk_rows_, cols_);
        }
        auto topRows(int j, int k) const {   // get the first j row-blocks of the k-th time instant
            return data_.template slice<2>(k).as_eigen_map().block(0, 0, j * blk_rows_, cols_);
        }
        auto operator()(int k) const { return data_.template slice<2>(k).as_eigen_map(); }
        int size() const { return data_.size(); }
        // modifiers
        auto operator()(int i, int k) {
            auto slice_ = data_.template slice<2>(k);
            return slice_.as_eigen_map().block(i * blk_rows_, 0, blk_rows_, cols_);
        }
        auto topRows(int j, int k) {
            return data_.template slice<2>(k).as_eigen_map().block(0, 0, j * blk_rows_, cols_);
        }
        auto operator()(int k) { return data_.template slice<2>(k).as_eigen_map(); }
    };

    template <typename DataLocs>
        requires(is_valid_data_locs_descriptor_v<DataLocs>)
    void eval_spatial_basis_at_(const DataLocs& locs) {
        if constexpr (std::is_same_v<DataLocs, matrix_t>) {   // pointwise sampling
            Psi_ = point_eval_(locs);
            D_ = vector_t::Ones(n_locs_).asDiagonal();
        } else {   // areal sampling
            const auto& [psi, measure_vect] = areal_eval_(locs);
            Psi_ = psi;
            D_ = measure_vect.asDiagonal();
        }
        fdapde_assert(n_locs_ == Psi_.rows());
        return;
    }
    template <typename GeoFrame> void eval_spatial_basis_at_(const GeoFrame& gf) {
        switch (gf.category(0)[0]) {
        case ltype::point: {
            const auto& spatial_index = geo_index_cast<0, POINT>(gf[0]);
            n_ = spatial_index.rows();
            if (spatial_index.points_at_dofs()) {
                Psi_.resize(n_, n_);
                Psi_.setIdentity();
            } else {
                Psi_ = point_eval_(spatial_index.coordinates());
            }
            D_ = vector_t::Ones(n_).asDiagonal();
            break;
        }
        case ltype::areal: {
            const auto& spatial_index = geo_index_cast<0, POLYGON>(gf[0]);
            n_ = spatial_index.rows();
            const auto& [psi, measure_vect] = areal_eval_(spatial_index.incidence_matrix());
            Psi_ = psi;
            D_ = measure_vect.asDiagonal();
            break;
        }
        }
        return;
    }
    // J(f,g) = \sum_{k=1}^m (y^k - \Psi*f^k)^T*(y^k - \Psi*f^k) + \lambda_S*(g^k)^T*(g^k) + \lambda_T*(l^k)^T*(l^k)
    double J_(const block_map_t& y, const block_map_t& x, double lambda_D, double lambda_T) const {
        double sse = 0;
        for (int t = 0; t < m_; ++t) {
            sse += (y(t) - Psi_ * x(0, t)).squaredNorm() + lambda_D * x(1, t).squaredNorm() +
                   lambda_T * x(2, t).squaredNorm();
        }
        return sse;
    }
   public:
    static constexpr int n_lambda = 2;
    using solver_category = ls_solver;

    fe_ls_separable_cdti() noexcept = default;
    template <typename GeoFrame, typename InfoT, typename WeightMatrix>
        requires(is_valid_info_t<InfoT>::value)
    fe_ls_separable_cdti(const std::string& formula, const GeoFrame& gf, InfoT&& info, const WeightMatrix& W) :
        max_iter_(info.max_iter), tol_(info.tol) {
        fdapde_static_assert(GeoFrame::Order == 2, THIS_CLASS_IS_FOR_ORDER_TWO_GEOFRAMES_ONLY);
	fdapde_assert(gf.n_layers() == 1);
        n_obs_  = gf[0].rows();
        n_locs_ = n_obs_;

        discretize(info.penalty);
        analyze_data(formula, gf, W);
    }
    template <typename GeoFrame, typename InfoT>
        requires(is_valid_info_t<InfoT>::value)
    fe_ls_separable_cdti(const std::string& formula, const GeoFrame& gf, InfoT&& info) :
        fe_ls_separable_cdti(formula, gf, info, vector_t::Ones(gf[0].rows()).asDiagonal()) { }
    // construct with no data
    template <typename GeoFrame, typename InfoT, typename WeightMatrix>
        requires(is_valid_info_t<InfoT>::value)
    fe_ls_separable_cdti(const GeoFrame& gf, InfoT&& info, const WeightMatrix& W) :
        W_(W), max_iter_(info.max_iter), tol_(info.tol) {
        fdapde_static_assert(GeoFrame::Order == 2, THIS_CLASS_IS_FOR_ORDER_TWO_GEOFRAMES_ONLY);
	fdapde_assert(gf.n_layers() == 1);
	n_obs_  = gf[0].rows();
	n_locs_ = n_obs_;

	// extract temporal mesh
        const auto& time_index = geo_index_cast<1, POINT>(gf[0]);
        const auto& time_coords = time_index.coordinates();
        m_ = time_coords.rows();
        fdapde_assert(m_ > 0 && time_coords.cols() == 1);
        DeltaT_ = time_coords(1, 0) - time_coords(0, 0);
        for (int i = 1; i < m_ - 1; ++i) {
            double lag_i = time_coords(i + 1, 0) - time_coords(i, 0);
            fdapde_assert(DeltaT_ > 0 && lag_i > 0 && almost_equal(DeltaT_ FDAPDE_COMMA lag_i));
        }

        discretize(info.penalty);
        u_.resize(n_dofs_ * m_);
        for (int i = 0; i < m_; ++i) { u_.segment(i * n_dofs_, n_dofs_) = u_space_; }
        eval_spatial_basis_at_(gf);
    }
    template <typename GeoFrame, typename InfoT>
        requires(is_valid_info_t<InfoT>::value)
    fe_ls_separable_cdti(const GeoFrame& gf, InfoT&& info) :
        fe_ls_separable_cdti(gf, info, vector_t::Ones(gf[0].rows()).asDiagonal()) { }

    template <typename Penalty> void discretize(Penalty&& penalty) {
        // fdapde_static_assert(internals::is_valid_penalty_pair_v<Penalty>, INVALID_PENALTY_DESCRIPTION);
        using BilinearForm = std::tuple_element_t<0, std::decay_t<Penalty>>;
        using LinearForm = std::tuple_element_t<1, std::decay_t<Penalty>>;
        using FeSpace = typename BilinearForm::TrialSpace;
        fdapde_static_assert(
          std::is_same_v<typename FeSpace::discretization_category FDAPDE_COMMA finite_element_tag>,
          NO_FINITE_ELEMENT_SPACE_DETECTED);
        // discretization
        const BilinearForm& bilinear_form = std::get<0>(penalty);
        const LinearForm& linear_form = std::get<1>(penalty);
        n_dofs_ = bilinear_form.n_dofs();
        internals::fe_mass_assembly_loop<FeSpace> mass_assembler(bilinear_form.trial_space());
        R0_ = mass_assembler.assemble();
        R1_ = bilinear_form.assemble();
	u_space_ = linear_form.assemble();
	// store handles for basis system evaluation at locations
        point_eval_ = [fe_space = bilinear_form.trial_space()](const matrix_t& locs) -> decltype(auto) {
            return internals::point_basis_eval(fe_space, locs);
        };
        areal_eval_ = [fe_space = bilinear_form.trial_space()](const binary_t& locs) -> decltype(auto) {
            return internals::areal_basis_eval(fe_space, locs);
        };
	return;
    }
    // non-parametric fit
    // \sum_i w_i * (y_i - f(p_i))^2 + \int_D (Lf - u)^2
    template <typename DataLocs1, typename DataLocs2, typename WeightMatrix>
        requires(is_valid_data_locs_descriptor_v<DataLocs1> && is_valid_data_locs_descriptor_v<DataLocs2>)
    void analyze_data(const DataLocs1& locs1, const DataLocs2& locs2, const matrix_t& y, const WeightMatrix& W) {
        fdapde_static_assert(
          std::is_same_v<DataLocs2 FDAPDE_COMMA matrix_t>,
          ITERATIVE_SEPARABLE_PENALIZATION_REQUIRES_POINT_TEMPORAL_EVALUATIONS_ONLY);
        fdapde_assert(
          locs1.rows() > 0 && locs2.rows() > 0 && y.rows() == locs1.rows() * locs2.rows() && y.cols() == 1 &&
          W.rows() == locs1.rows() * locs2.rows() && W.rows() == W.cols());
        n_obs_ = y.rows();
        n_locs_ = n_obs_;
        n_covs_ = 0;
        eval_spatial_basis_at_(locs1);   // update \Psi matrix

        m_ = locs2.rows();
        fdapde_assert(m_ > 0 && locs2.cols() == 1);
        DeltaT_ = locs2(1, 0) - locs2(0, 0);
        for (int i = 1; i < m_ - 1; ++i) {
            double lag_i = locs2(i + 1, 0) - locs2(i, 0);
            fdapde_assert(DeltaT_ > 0 && lag_i > 0 && almost_equal(DeltaT_ FDAPDE_COMMA lag_i));
        }
	// update forcing
        if (u_.size() == 0) {
            u_.resize(n_dofs_ * m_);
            for (int i = 0; i < m_; ++i) { u_.segment(i * n_dofs_, n_dofs_) = u_space_; }
        }
        update_response_and_weights(y, W);
        return;
    }
    // fit from formula
    template <typename GeoFrame, typename WeightMatrix>
    void analyze_data(const std::string& formula, const GeoFrame& gf, const WeightMatrix&) {
        fdapde_static_assert(GeoFrame::Order == 2, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        fdapde_assert(gf.n_layers() == 1 && gf[0].category()[1] == ltype::point);
        n_obs_ = gf[0].rows();
        n_locs_ = n_obs_;

        // extract time step and number of time instants
        const auto& time_index = geo_index_cast<1, POINT>(gf[0]);
        const auto& time_coords = time_index.coordinates();
        m_ = time_coords.rows();
        fdapde_assert(m_ > 0 && time_coords.cols() == 1);
        DeltaT_ = time_coords(1, 0) - time_coords(0, 0);
        for (int i = 1; i < m_ - 1; ++i) {
            double lag_i = time_coords(i + 1, 0) - time_coords(i, 0);
            fdapde_assert(DeltaT_ > 0 && lag_i > 0 && almost_equal(DeltaT_ FDAPDE_COMMA lag_i));
        }
        if (u_.size() == 0) {
            u_.resize(n_dofs_ * m_);
            for (int i = 0; i < m_; ++i) { u_.segment(i * n_dofs_, n_dofs_) = u_space_; }
        }
        eval_spatial_basis_at_(gf);   // update \Psi matrix
        // parse formula, extract response vector
        Formula formula_(formula);
        const auto& y_data = gf[0].data().template col<double>(formula_.lhs());
        y_.resize(n_locs_, y_data.blk_sz());
        y_data.assign_to(y_);

        // update_weights(vector_t::Ones(n_).asDiagonal());
    }

    // modifiers
    // void update_response(const vector_t& y) {
    //     fdapde_assert(Psi_.rows() > 0 && y.rows() == n_locs_ && y.cols() == 1);
    //     y_ = y;
    // 	b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * W_ * y;
    //     return;
    // }
    // template <typename WeightMatrix> void update_weights(const WeightMatrix& W) {
    //     fdapde_assert(Psi_.rows() > 0 && W.rows() == n_locs_ && W.rows() == W.cols());
    //     W_ = W;
    //     // W_ /= n_obs_; ---------------------- divide by n_obs_ or just by n_ ????
    //     b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * W_ * y_;
    //     W_changed_ = true;
    //     return;
    // }
    // template <typename WeightMatrix> void update_response_and_weights(const vector_t& y, const WeightMatrix& W) {
    //     fdapde_assert(
    //       Psi_.rows() > 0 && y.rows() == n_locs_ && y.cols() == 1 && W.rows() == W.cols() && W.rows() == n_);
    //     y_ = y;
    // 	update_weights(W);
    // 	return;
    // }
   private:
    // iterative scheme implementation
    template <typename ResponseT> auto fit_(ResponseT&& response, double lambda_D, double lambda_T) {
        std::array<double, n_lambda> lambda {lambda_D, lambda_T};
        // define auxiliary structures
        block_map_t y(response, n_);
        block_map_t u(u_, n_dofs_);
        matrix_t x_old_buff(3 * n_dofs_ * m_, response.cols()), x_new_buff(3 * n_dofs_ * m_, response.cols());
        block_map_t x_old(x_old_buff, 3 * n_dofs_, n_dofs_, response.cols());
        block_map_t x_new(x_new_buff, 3 * n_dofs_, n_dofs_, response.cols());
        double alpha = lambda_T / std::pow(DeltaT_, 2);
        {   // compute starting point (f^(k, 0), g^(k, 0), l^(k, 0)) k = 1 ... m
            if (lambda_saved_.value() != lambda) {
                SparseBlockMatrix<double, 2, 2> A_(
                  Psi_.transpose() * D_ * Psi_, lambda_D * R1_.transpose(), lambda_D * R1_, -lambda_D * R0_);
                invAs_.compute(A_);
            }
            vector_t b_(2 * n_dofs_);
            for (int t = 0; t < m_; ++t) {
                b_ << Psi_.transpose() * D_ * y(t), lambda_D * u(t);
                x_old.topRows(2, t) = invAs_.solve(b_);
            }
            x_old(2, 0).setZero();
            x_old(2, m_ - 1).setZero();
            for (int t = 1; t < m_ - 1; ++t) {
                x_old(2, t) = (x_old(0, t + 1) - 2 * x_old(0, t) + x_old(0, t - 1)) / std::pow(DeltaT_, 2);
            }
        }
        // iterative_tag scheme initialization
        double Jold = std::numeric_limits<double>::max();
        double Jnew = J_(y, x_old, lambda_D, lambda_T);
        if (lambda_saved_.value() != lambda) {
            sparse_matrix_t Zero(n_dofs_, n_dofs_);
            SparseBlockMatrix<double, 3, 3> A_(
              Psi_.transpose() * D_ * Psi_, lambda_D * R1_.transpose(), -2 * alpha * R0_, lambda_D * R1_,
              -lambda_D * R0_, Zero, -2 * alpha * R0_, Zero, -lambda_T * R0_);
            invA_.compute(A_);
        }
        vector_t b_(3 * n_dofs_);
        // iterative loop
        x_new(2, 0).setZero();
        x_new(2, m_ - 1).setZero();
        int i = 1;
        while (i < max_iter_ && std::abs((Jnew - Jold) / Jnew) > tol_) {
            // at step 0: f^(-1, i-1) = l^(-1, i-1) = 0
            b_ << Psi_.transpose() * D_ * y(0) - alpha * R0_ * x_old(2, 1), lambda_D * u(0),
              -alpha * R0_ * x_old(0, 1);
            x_new.topRows(2, 0) = invA_.solve(b_).topRows(2 * n_dofs_);   // l^(0) = 0
            // general step
            for (int t = 1; t < m_ - 1; ++t) {
                b_ << Psi_.transpose() * D_ * y(t) - alpha * R0_ * (x_old(2, t + 1) + x_old(2, t - 1)),
                  lambda_D * u(t), -alpha * R0_ * (x_old(0, t + 1) + x_old(0, t - 1));
                x_new(t) = invA_.solve(b_);
            }
            // at step m_ - 1: f^(m+1, i-1) = l^(m+1, i-1) = 0
            b_ << Psi_.transpose() * D_ * y(m_ - 1) - alpha * R0_ * x_old(2, m_ - 1), lambda_D * u(m_ - 1),
              -alpha * R0_ * x_old(0, m_ - 1);
            x_new.topRows(2, m_ - 1) = invA_.solve(b_).topRows(2 * n_dofs_);   // l^(m_ - 1) = 0
            // prepare for next iteration
            Jold = Jnew;
            Jnew = J_(y, x_new, lambda_D, lambda_T);

            x_old = x_new;
            i++;
	}
        // return result
        vector_t f(n_dofs_ * m_, response.cols());
        vector_t g(n_dofs_ * m_, response.cols());
        for (int i = 0; i < m_; ++i) {
            f.middleRows(i * n_dofs_, n_dofs_) = x_old(0, i);
            g.middleRows(i * n_dofs_, n_dofs_) = x_old(1, i);
        }
        lambda_saved_ = std::array<double, n_lambda> {lambda_D, lambda_T};
        return std::make_pair(f, g);
    }
   public:
    // main fit entry point
    const vector_t& fit(double lambda_D, double lambda_T) {
        const auto& [f, g] = fit_(y_, lambda_D, lambda_T);
        f_ = f;
        g_ = g;
        return f_;
    }
    template <typename LambdaT>
        requires(internals::is_vector_like_v<LambdaT>)
    const vector_t& fit(LambdaT&& lambda) {
        fdapde_assert(lambda.size() == n_lambda);
        return fit(lambda[0], lambda[1]);
    }

    // hutchinson approximation for Tr[S]
    double edf(int r = 100, int seed = random_seed) {
        fdapde_assert(lambda_saved_.has_value());
        if (!Us_.has_value()) {
            int seed_ = (seed == random_seed) ? std::random_device()() : seed;
            std::mt19937 rng(seed_);
            rademacher_distribution rademacher;
            Us_->resize(n_locs_, r);
            for (int i = 0; i < n_locs_; ++i) {
                for (int j = 0; j < r; ++j) { Us_->operator()(i, j) = rademacher(rng); }
            }
        }
        // Tr[S] \approx \sum_{i=1}^r (u_i^\top * S * u_i)
        double trS = 0;
        for (int i = 0; i < r; ++i) {
            const auto& [f, g] = fit_(Us_->col(i), lambda_saved_->operator[](0), lambda_saved_->operator[](1));
            vector_t fn(n_ * m_);
            for (int i = 0; i < m_; ++i) { fn.middleRows(i * n_, n_) = Psi_ * f.middleRows(i * n_dofs_, n_dofs_); }
            trS += Us_->col(i).dot(fn);
        }
        return trS / r;
    }
    template <typename... LambdaT>
        requires(
          (sizeof...(LambdaT) == 1 && (internals::is_vector_like_v<LambdaT> && ...)) ||
          (sizeof...(LambdaT) == n_lambda && (std::is_floating_point_v<LambdaT> && ...)))
    double edf(const LambdaT&... lambda, int r = 100, int seed = random_seed) {
        std::array<double, n_lambda> lambda_;
        if constexpr (sizeof...(LambdaT) == 1) {
            internals::for_each_index_and_args<sizeof...(LambdaT)>([&]<int Ns_, typename Ts_>(const Ts_& ts) {
                fdapde_assert(ts.size() == n_lambda && ts[0] > 0 && ts[1] > 0);
                lambda_[0] = ts[0];
                lambda_[1] = ts[1];
            });
        } else {
            std::array<double, n_lambda> lambda__ {static_cast<double>(lambda)...};
            fdapde_assert(lambda__[0] > 0 && lambda__[1] > 0);
	    lambda_ = lambda__;
        }
        return edf(r, seed);
    }
    vector_t fn() const {
        vector_t fn_(n_ * m_);
        for (int i = 0; i < m_; ++i) { fn_.middleRows(i * n_, n_) = Psi_ * f_.middleRows(i * n_dofs_, n_dofs_); }
        return fn_;
    }

    // observers
    int n_dofs() const { return n_dofs_; }
    const sparse_matrix_t& mass() const { return R0_; }
    const sparse_matrix_t& stiff() const { return R1_; }
    const sparse_matrix_t& Psi() const { return Psi_; }
    const vector_t& force() const { return u_; }
    const vector_t& f() const { return f_; }
    const vector_t& beta() const { return beta_; }
    const vector_t& misfit() const { return g_; }
    const vector_t& response() const { return y_; }
   protected:
    std::optional<std::array<double, n_lambda>> lambda_saved_ = std::array<double, n_lambda> {-1, -1};
    sparse_solver_t invA_, invAs_;
    // matrices for hutchinson stochastic estimation of Tr[S]
    std::optional<matrix_t> Us_;

    int n_dofs_ = 0, n_obs_ = 0, n_covs_ = 0;
    int n_locs_ = 0, n_ = 0, m_ = 0;   // n_: number of spatial locations, m_: number of time instants

    sparse_matrix_t R0_;    // n_dofs x n_dofs matrix [R0]_{ij} = \int_D \psi_i * \psi_j
    sparse_matrix_t R1_;    // n_dofs x n_dofs matrix [R1]_{ij} = \int_D a(\psi_i, \psi_j)
    sparse_matrix_t Psi_;   // n_obs x n_dofs matrix [Psi]_{ij} = \psi_j(p_i)
    vector_t u_;            // (n_dofs * m) x 1 vector u = [u_1 + + R0_*s / DeltaT, u_2, \ldots, u_n]
    vector_t u_space_;
    diag_matrix_t D_;       // vector of regions' measures (areal sampling)
    mutable sparse_solver_t invR0_;
    vector_t f_, beta_, g_;
    // basis system evaluation handles
    std::function<sparse_matrix_t(const matrix_t& locs)> point_eval_;
    std::function<std::pair<sparse_matrix_t, vector_t>(const binary_t& locs)> areal_eval_;

    vector_t y_;          // n_obs x 1 observation vector
    sparse_matrix_t W_;   // n_obs x n_obs matrix of observation weights
    bool W_changed_;

    int max_iter_;   // maximum number of iterations
    double tol_;     // convergence tolerance
    double DeltaT_;
};

}   // namespace internals

// separable central finite differences time stepping method
template <typename Penalty>
    requires(internals::is_valid_penalty_pair_v<Penalty>)
struct fe_ls_separable_cdti {
    using solver_t = internals::fe_ls_separable_cdti;
   private:
    struct info_t {
        Penalty penalty;
        int max_iter = 50;
        double tol = 1e-4;
    };
   public:
    fe_ls_separable_cdti(const Penalty& penalty, int max_iter = 50, double tol = 1e-4) :
        info_(penalty, max_iter, tol) { }
    const info_t& get() const { return info_; }
   private:
    info_t info_;
};

}   // namespace fdapde

#endif // __FE_LS_SEPARABLE_SOLVER_H__

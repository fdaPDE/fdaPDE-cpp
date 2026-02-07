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

#ifndef __FE_LS_PARABOLIC_SOLVER_H__
#define __FE_LS_PARABOLIC_SOLVER_H__

#include "header_check.h"

namespace fdapde {
namespace internals {

// solves \min_{f, \beta} \| W^{1/2} * (y_i - x_i^\top * \beta - f(p_i, t_j)) \|_2^2 +
// \int_D \int_T (\frac{\partial f}{\partial t} + L(f) - u)^2
class fe_ls_parabolic_mono {
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
    template <typename Penalty> struct is_valid_penalty {
        static constexpr bool value = requires(Penalty penalty) {
            penalty.bilinear_form();
            penalty.linear_form();
            penalty.ic();
        };
    };
    template <typename Penalty> static constexpr bool is_valid_penalty_v = is_valid_penalty<Penalty>::value;

    // basis evaluation at locations
    template <typename DataLocs>
        requires(is_valid_data_locs_descriptor_v<DataLocs>)
    void eval_basis_at_(const DataLocs& locs) {
        if constexpr (std::is_same_v<DataLocs, matrix_t>) {   // pointwise sampling
            Psi_ = point_eval_(locs);
            D_ = vector_t::Ones(n_locs_).asDiagonal();
	    fdapde_assert(n_locs_ == Psi_.rows());
        } else {   // areal sampling
            const auto& [psi, measure_vec] = areal_eval_(locs);
            Psi_ = psi;
            D_ = measure_vec.asDiagonal();
	    fdapde_assert(n_locs_ == Psi_.rows());
        }
        return;
    }
    // optimized basis evaluation at geoframe
    template <typename GeoFrame> void eval_basis_at_(const GeoFrame& gf) {
        switch (gf.category(0)[0]) {
        case ltype::point: {
            const auto& spatial_index = geo_index_cast<0, POINT>(gf[0]);
            n_ = spatial_index.rows();
            if (spatial_index.points_at_dofs()) {
                Psi__.resize(n_, n_);
                Psi__.setIdentity();
            } else {
                Psi__ = point_eval_(spatial_index.coordinates());
            }
            D_ = vector_t::Ones(n_locs_).asDiagonal();
            break;
        }
        case ltype::areal: {
            const auto& spatial_index = geo_index_cast<0, POLYGON>(gf[0]);
            const auto& [psi, measure_vec] = areal_eval_(spatial_index.incidence_matrix());
            Psi__ = psi;
	    n_ = spatial_index.rows();
            vector_t D(n_locs_);
            for (int i = 0; i < m_; ++i) { D.segment(i * measure_vec.rows(), measure_vec.rows()) = measure_vec; }
            D_ = D.asDiagonal();
            break;
        }
        }
	return;
    }
    void tensorize_(int m) {
        sparse_matrix_t Im(m, m);   // m x m identity matrix
        Im.setIdentity();
        R0_ = kronecker(Im, R0__);
        R1_ = kronecker(Im, R1__);
        {
            // assemble matrix associated with derivation in time L_
            // [L]_{ii} = 1/DeltaT for i \in {1 ... m} and [L_]_{i,i-1} = -1/DeltaT for i \in {1 ... m-1}
            std::vector<Eigen::Triplet<double>> triplet_list;
            triplet_list.reserve(2 * m);
            // start assembly loop
            triplet_list.emplace_back(0, 0, 1.0 / DeltaT_);
            for (int i = 1; i < m; ++i) {
                triplet_list.emplace_back(i, i, 1.0 / DeltaT_);
                triplet_list.emplace_back(i, i - 1, -1.0 / DeltaT_);
            }
            L__.resize(m, m);
            L__.setFromTriplets(triplet_list.begin(), triplet_list.end());
            L__.makeCompressed();
            L_ = kronecker(L__, R0__);
        }
        u_.resize(n_dofs_ * m);
        for (int i = 0; i < m; ++i) { u_.segment(i * n_dofs_, n_dofs_) = u__; }
        // correct first n discretized force rows as (u_1 + (R0 * s) / DeltaT);
        u_.segment(0, n_dofs_) += (1.0 / DeltaT_) * (R0__ * s_);
        Psi_ = kronecker(Im, Psi__);
        return;
    }
   public:
    static constexpr int n_lambda = 2;
    using solver_category = ls_solver;

    fe_ls_parabolic_mono() noexcept = default;
    // construct from formula + geoframe
    template <typename GeoFrame, typename WeightMatrix, typename Penalty>
        requires(is_valid_penalty_v<Penalty>)
    fe_ls_parabolic_mono(const std::string& formula, const GeoFrame& gf, Penalty&& penalty, const WeightMatrix& W) {
        fdapde_static_assert(GeoFrame::Order == 2, THIS_CLASS_IS_FOR_ORDER_TWO_GEOFRAMES_ONLY);
	fdapde_assert(gf.n_layers() == 1);
        n_obs_ = gf[0].rows();
        n_locs_ = n_obs_;

        discretize(penalty);
        analyze_data(formula, gf, W);
    }
    template <typename GeoFrame, typename Penalty>
        requires(is_valid_penalty_v<Penalty>)
    fe_ls_parabolic_mono(const std::string& formula, const GeoFrame& gf, Penalty&& penalty) :
        fe_ls_parabolic_mono(formula, gf, penalty, vector_t::Ones(gf[0].rows()).asDiagonal()) { }
    // construct with no data
    template <typename GeoFrame, typename WeightMatrix, typename Penalty>
        requires(is_valid_penalty_v<Penalty>)
    fe_ls_parabolic_mono(const GeoFrame& gf, Penalty&& penalty, const WeightMatrix& W) : W_(W) {
        fdapde_static_assert(GeoFrame::Order == 2, THIS_CLASS_IS_FOR_ORDER_TWO_GEOFRAMES_ONLY);
	fdapde_assert(gf.n_layers() == 1);
        n_obs_ = gf[0].rows();
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
        discretize(penalty);
        // basis system evaluation
	eval_basis_at_(gf);
	tensorize_(m_);
    }

    template <typename GeoFrame, typename Penalty>
        requires(is_valid_penalty_v<Penalty>)
    fe_ls_parabolic_mono(const GeoFrame& gf, Penalty&& penalty) :
        fe_ls_parabolic_mono(gf, penalty, vector_t::Ones(gf[0].rows()).asDiagonal()) { }

    template <typename Penalty> void discretize(Penalty&& penalty) {
        using BilinearForm = typename std::decay_t<Penalty>::BilinearForm;
        using LinearForm = typename std::decay_t<Penalty>::LinearForm;
        fdapde_static_assert(
          internals::is_valid_penalty_pair_v<BilinearForm FDAPDE_COMMA LinearForm>, INVALID_PENALTY_DESCRIPTION);
        using FeSpace = typename BilinearForm::TrialSpace;
	// discretization
	const BilinearForm& bilinear_form = penalty.bilinear_form();
        const LinearForm& linear_form = penalty.linear_form();
        n_dofs_ = bilinear_form.n_dofs();   // number of basis functions over physical domain
	internals::fe_mass_assembly_loop<FeSpace> mass_assembler(bilinear_form.trial_space());
        R0__ = mass_assembler.assemble();
        R1__ = bilinear_form.assemble();
	u__ = linear_form.assemble();	
	// store handles for basis system evaluation at locations
        point_eval_ = [fe_space = bilinear_form.trial_space()](const matrix_t& locs) -> decltype(auto) {
            return internals::point_basis_eval(fe_space, locs);
        };
        areal_eval_ = [fe_space = bilinear_form.trial_space()](const binary_t& locs) -> decltype(auto) {
            return internals::areal_basis_eval(fe_space, locs);
        };
	// store initial condition
	s_ = penalty.ic();
	return;
    }
    // non-parametric fit
    // \sum_i w_i * (y_i - f(p_i))^2 + \int_D (df/dt + Lf - u)^2
    template <typename DataLocs1, typename DataLocs2, typename WeightMatrix>
        requires(is_valid_data_locs_descriptor_v<DataLocs1> && is_valid_data_locs_descriptor_v<DataLocs2>)
    void analyze_data(const DataLocs1& locs1, const DataLocs2& locs2, const matrix_t& y, const WeightMatrix& W) {
        fdapde_assert(
          locs1.rows() > 0 && locs2.rows() > 0 && y.rows() == locs1.rows() * locs2.rows() && y.cols() == 1 &&
          W.rows() == locs1.rows() * locs2.rows() && W.rows() == W.cols());
        n_obs_ = y.rows();
        n_locs_ = n_obs_;
        n_covs_ = 0;
        eval_basis_at_(locs1);   // update \Psi matrix

        int m = locs2.rows();
        if (m != m_) {   // re-tensorize all if number of time steps chanded
	  m_ = m;
	  tensorize_(m_);
        } else {
            sparse_matrix_t Im(m, m);   // m x m identity matrix
            Im.setIdentity();
            Psi_ = kronecker(Im, Psi__);
        }
        update_response_and_weights(y, W);
        return;
    }
    // semi-parametric fit
    // \sum_i w_i * (y_i - x_i^\top * \beta - f(p_i))^2 + \int_D (df/dt + Lf - u)^2
    template <typename DataLocs1, typename DataLocs2, typename WeightMatrix>
        requires(is_valid_data_locs_descriptor_v<DataLocs1> && is_valid_data_locs_descriptor_v<DataLocs2>)
    void analyze_data(
      const DataLocs1& locs1, const DataLocs2& locs2, const matrix_t& y, const matrix_t& X, const WeightMatrix& W) {
        fdapde_assert(
          locs1.rows() > 0 && locs2.rows() > 0 && y.rows() == locs1.rows() * locs2.rows() && y.cols() == 1 &&
          X.rows() == locs1.rows() * locs2.rows() && W.rows() == locs1.rows() * locs2.rows() && W.rows() == W.cols());

        n_obs_ = y.rows();
        n_locs_ = n_obs_;
        bool require_woodbury_realloc = std::cmp_not_equal(n_covs_, X.cols());
        n_covs_ = X.cols();
	eval_basis_at_(locs1);   // update \Psi matrix

	int m = locs2.rows();
        if (m != m_) {   // re-tensorize all if number of time steps chanded
	  m_ = m;
	  tensorize_(m_);
        } else {
            sparse_matrix_t Im(m, m);   // m x m identity matrix
            Im.setIdentity();
            Psi_ = kronecker(Im, Psi__);
        }
        if (require_woodbury_realloc) { U_ = matrix_t::Zero(2 * m_ * n_dofs_, n_covs_); }
        if (require_woodbury_realloc) { V_ = matrix_t::Zero(n_covs_, 2 * m_ * n_dofs_); }
        update_response_and_weights(y, X, W);
        return;
    }
    // fit from formula
    template <typename GeoFrame, typename WeightMatrix>
    void analyze_data(const std::string& formula, const GeoFrame& gf, const WeightMatrix& W) {
        fdapde_static_assert(GeoFrame::Order == 2, THIS_CLASS_IS_FOR_ORDER_TWO_GEOFRAMES_ONLY);
        fdapde_assert(gf.n_layers() == 1);
	n_obs_ = gf[0].rows();
	n_locs_ = n_obs_;
	// extract temporal mesh
        const auto& time_index = geo_index_cast<1, POINT>(gf[0]);
        const auto& time_coords = time_index.coordinates();
	bool require_full_tensorization = std::cmp_not_equal(m_, time_coords.rows());
        m_ = time_coords.rows();
        fdapde_assert(m_ > 0 && time_coords.cols() == 1);
        DeltaT_ = time_coords(1, 0) - time_coords(0, 0);
        for (int i = 1; i < m_ - 1; ++i) {
            double lag_i = time_coords(i + 1, 0) - time_coords(i, 0);
            fdapde_assert(DeltaT_ > 0 && lag_i > 0 && almost_equal(DeltaT_ FDAPDE_COMMA lag_i));
        }
        // basis system evaluation
	eval_basis_at_(gf);
        if (require_full_tensorization) {
            tensorize_(m_);
        } else {
            sparse_matrix_t Im(m_, m_);   // m x m identity matrix
            Im.setIdentity();
            Psi_ = kronecker(Im, Psi__);
        }
        // parse formula, extract response vector and design matrix
        Formula formula_(formula);
        std::vector<std::string> covs;
        for (const std::string& token : formula_.covs()) {
            if (gf.contains(token)) { covs.push_back(token); }
        }
	bool require_woodbury_realloc = std::cmp_not_equal(n_covs_, covs.size());
        n_covs_ = covs.size();
	const auto& y_data = gf[0].data().template col<double>(formula_.lhs());
        y_.resize(n_locs_, y_data.blk_sz());
        y_data.assign_to(y_);

	if (b_.cols() != y_.cols()) { b_.resize(2 * m_ * n_dofs_, y_.cols()); }
        if (n_covs_ != 0) {
            if (require_woodbury_realloc) { U_ = matrix_t::Zero(2 * m_ * n_dofs_, n_covs_); }
            if (require_woodbury_realloc) { V_ = matrix_t::Zero(n_covs_, 2 * m_ * n_dofs_); }
            X_.resize(n_obs_, n_covs_);   // assemble design matrix
            for (int i = 0; i < n_covs_; ++i) { gf[0].data().template col<double>(covs[i]).assign_to(X_.col(i)); }
        }
        update_response_and_weights(y_, W); 
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
            B_ = (~nan_pattern).repeat(1, m_ * n_dofs_).select(Psi_, 0);
            y_ = (~nan_pattern).select(y_, 0);
        }
        if (old_n_obs != n_obs_) { W_ *= (double)old_n_obs / n_obs_; } // re-normalize
        b_.block(0, 0, m_ * n_dofs_, 1) = -PsiNA().transpose() * D_ * W_ * y;
        return;
    }
    template <typename WeightMatrix> void update_weights(const WeightMatrix& W) {
        fdapde_assert(Psi_.rows() > 0 && W.rows() == n_locs_ && W.rows() == W.cols());
        W_ = W;
        W_ /= n_obs_;
        if (n_covs_ == 0) {
            b_.block(0, 0, m_ * n_dofs_, 1) = -PsiNA().transpose() * D_ * W_ * y_;
        } else {
            XtWX_ = X_.transpose() * W_ * X_;
            invXtWX_ = XtWX_.partialPivLu();
            invXtWXXtW_ = invXtWX_.solve(X_.transpose() * W_);   // (X^\top * W * X)^{-1} * (X^\top * W)
            // woodbury decomposition matrices
            U_.block(0, 0, m_ * n_dofs_, n_covs_) = PsiNA().transpose() * D_ * W_ * X_;
            V_.block(0, 0, n_covs_, m_ * n_dofs_) = X_.transpose() * W_ * PsiNA();
            b_.block(0, 0, m_ * n_dofs_, 1) = -PsiNA().transpose() * D_ * internals::lmbQ(W_, X_, invXtWX_, y_);
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
            B_ = (~nan_pattern).repeat(1, m_ * n_dofs_).select(Psi_, 0);
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
            // assemble system matrix for the nonparameteric part
            SparseBlockMatrix<double, 2, 2> A(
              -PsiNA().transpose() * D_ * W_ * PsiNA(), lambda_D * (R1_ + lambda_T * L_).transpose(),
              lambda_D * (R1_ + lambda_T * L_), lambda_D * R0_);
            invA_.compute(A);
            W_changed_ = false;
        }
        if (lambda_saved_.value() != lambda) {
            // linear system rhs
            b_.block(n_dofs_ * m_, 0, n_dofs_ * m_, 1) = lambda_D * u_;
        }
        lambda_saved_ = lambda;
        vector_t x;
        if (n_covs_ == 0) {   // nonparametric case
            x = invA_.solve(b_);
            f_ = x.head(n_dofs_ * m_);
        } else {   // parametric case
            x = woodbury_system_solve(invA_, U_, XtWX_, V_, b_);
            f_ = x.head(n_dofs_ * m_);
            beta_ = invXtWXXtW_ * (y_ - Psi_ * f_);
        }
        g_ = x.tail(n_dofs_ * m_);
        return std::make_pair(f_, beta_);
    }
    template <typename LambdaT>
        requires(internals::is_vector_like_v<LambdaT>)
    std::pair<vector_t, vector_t> fit(LambdaT&& lambda) {
        fdapde_assert(lambda.size() == n_lambda);
        return fit(lambda[0], lambda[1]);
    }
    // perform a nonparametric fit, e.g. discarding possible covariates
    vector_t nonparametric_fit(double lambda_D, double lambda_T) {
        fdapde_assert(lambda_D > 0 && lambda_T > 0 && n_dofs_ > 0 && n_obs_ > 0);
        std::array<double, n_lambda> lambda {lambda_D, lambda_T};
        if (lambda_saved_.value() != lambda) {
            // assemble system matrix for the nonparameteric part
            SparseBlockMatrix<double, 2, 2> A(
              -PsiNA().transpose() * D_ * W_ * PsiNA(), lambda_D * (R1_ + lambda_T * L_).transpose(),
              lambda_D * (R1_ + lambda_T * L_), lambda_D * R0_);
            invA_.compute(A);
        }
	vector_t x;
        if (n_covs_ == 0) {   // equivalent to calling fit(lambda)
            if (lambda_saved_.value() != lambda) { b_.block(m_ * n_dofs_, 0, m_ * n_dofs_, 1) = lambda_D * u_; }
            x = invA_.solve(b_);
        } else {
            vector_t b(2 * m_ * n_dofs_);
            // assemble nonparametric linear system rhs
            b.block(0, 0, m_ * n_dofs_, 1) = -PsiNA().transpose() * D_ * W_ * y_;
            b.block(m_ * n_dofs_, 0, m_ * n_dofs_, 1) = lambda_D * u_;   
            x = invA_.solve(b);
        }
        lambda_saved_ = lambda;
        f_ = x.topRows(m_ * n_dofs_);
        g_ = x.bottomRows(m_ * n_dofs_);
        return f_;
    }

    // hutchinson approximation for Tr[S]
    double edf(int r = 100, int seed = random_seed) {
        fdapde_assert(lambda_saved_.has_value());
        if (!Ys_.has_value() || !Bs_.has_value() || r != Us_->rows()) {   // force reconstruction if r differs from old
            int seed_ = (seed == random_seed) ? std::random_device()() : seed;
            std::mt19937 rng(seed_);
            rademacher_distribution rademacher;
            Us_->resize(n_locs_, r);
            for (int i = 0; i < n_locs_; ++i) {
                for (int j = 0; j < r; ++j) { Us_->operator()(i, j) = rademacher(rng); }
            }
            Ys_ = Us_->transpose() * Psi_;
            Bs_ = matrix_t::Zero(2 * m_ * n_dofs_, r);   // implicitly enforce homogeneous forcing
        }
        if (n_covs_ == 0) {
            Bs_->topRows(m_ * n_dofs_) = -PsiNA().transpose() * D_ * W_ * (*Us_);
        } else {
            Bs_->topRows(m_ * n_dofs_) = -PsiNA().transpose() * D_ * internals::lmbQ(W_, X_, invXtWX_, *Us_);
        }
        matrix_t x = n_covs_ == 0 ? invA_.solve(*Bs_) : woodbury_system_solve(invA_, U_, XtWX_, V_, *Bs_);
        double trS = 0;   // monte carlo Tr[S] approximation
        for (int i = 0; i < r; ++i) { trS += Ys_->row(i).dot(x.col(i).head(m_ * n_dofs_)); }
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
              -PsiNA().transpose() * D_ * W_ * PsiNA(), lambda_[0] * (R1_ + lambda_[1] * L_).transpose(),
              lambda_[0] * (R1_ + lambda_[1] * L_), lambda_[0] * R0_);
            invA_.compute(A);
            lambda_saved_ = lambda_;
        }
        return edf(r, seed);
    }

    // penalty matrix
    matrix_t P(double lambda_D, double lambda_T) const {
        if (!invR0_.has_value()) { invR0_->compute(R0__); }
        if (!PT_.has_value()) { PT_ = kronecker(L__, R0__); }
        return lambda_D * (R1__ + lambda_T * (*PT_)).transpose() * invR0_->solve(R1__ + lambda_T * (*PT_));
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
        if (!PT_.has_value()) { PT_ = kronecker(L__, R0__); }
        return lambda_D * g_.dot(R0_ * g_) + lambda_T * f_.dot((*PT_) * f_);
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
    double DeltaT_ = 0;   // time step
    // not tensorized quantities
    sparse_matrix_t R0__;    // n_dofs x n_dofs matrix [R0]_{ij} = \int_D \psi_i * \psi_j
    sparse_matrix_t R1__;    // n_dofs x n_dofs matrix [R1]_{ij} = \int_D a(\psi_i, \psi_j)
    vector_t u__;            // n_dofs x 1 vector [u]_i = \int_D u * \psi_i
    vector_t s_;             // initial condition vector
    sparse_matrix_t Psi__;   // n_locs x n_dofs matrix [Psi]_{ij} = \psi_j(p_i)
    sparse_matrix_t L__;

    sparse_matrix_t R0_;    // (n_dofs * m) x (n_dofs * m) matrix R0 = Im \kron R0__
    sparse_matrix_t R1_;    // (n_dofs * m) x (n_dofs * m) matrix R1 = Im \kron R1__
    sparse_matrix_t L_;     // (n_dofs * m) x (n_dofs * m) matrix L = L__ \kron R0__
    vector_t u_;            // (n_dofs * m) x 1 vector u = [u_1 + + R0__*s / DeltaT, u_2, \ldots, u_n]
    sparse_matrix_t Psi_;   // (n_obs * m) x (n_dofs * m) matrix Psi = Im \kron Psi__
    diag_matrix_t D_;       // vector of regions' measures (areal sampling)
    mutable std::optional<sparse_solver_t> invR0_;
    std::optional<sparse_matrix_t> B_;            // \Psi matrix corrected for missing observations
    mutable std::optional<sparse_matrix_t> PT_;   // (n_dofs * m) x (n_dofs * m) matrix PT = L__ \kron R0__
    vector_t f_, beta_, g_;
    // basis system evaluation handles
    std::function<sparse_matrix_t(const matrix_t& locs)> point_eval_;
    std::function<std::pair<sparse_matrix_t, vector_t>(const binary_t& locs)> areal_eval_;

    int n_, m_;                // number of spatial and temporal locations
    matrix_t X_;               // n_obs x n_covs design matrix
    vector_t y_;               // n_obs x 1 observation vector
    sparse_matrix_t W_;        // n_obs x n_obs matrix of observation weights
    matrix_t U_, V_;           // (2*n_dofs*m) x n_covs matrices [\Psi^\top * D * W * y, 0] and [X^\top * W * \Psi, 0]
    matrix_t XtWX_;            // n_covs x n_covs matrix X^\top * W * X
    matrix_t invXtWXXtW_;      // n_covs x n_obs matrix (X^\top * X)^{-1} * (X^\top W)
    dense_solver_t invXtWX_;   // factorization of n_covs x n_covs matrix X^\top * W * X
    bool W_changed_;
};

// implicit euler time stepping scheme
struct fe_ls_parabolic_ieul {
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
    template <typename Penalty> struct is_valid_penalty {
        static constexpr bool value = requires(Penalty penalty) {
            penalty.ic();
            penalty.bilinear_form();
            penalty.linear_form();
            penalty.max_iter();
            penalty.tol();
        };
    };
    template <typename Penalty> static constexpr bool is_valid_penalty_v = is_valid_penalty<Penalty>::value;
  
    // auxiliary time-mapping data structure
    class block_map_t {
        static constexpr int Order = 3;
        using Scalar = double;
        using storage_t = MdMap<Scalar, full_dynamic_extent_t<Order>, internals::layout_left>;

        storage_t data_;
        int rows_ = 0, cols_ = 0;
        int blk_rows_ = 0, blk_cols_ = 0;
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
            for (std::size_t i = 0; i < data_.size(); ++i) { data_.data()[i] = other.data_.data()[i]; }
        }
        block_map_t& operator=(const block_map_t& other) {
            for (std::size_t i = 0; i < data_.size(); ++i) { data_.data()[i] = other.data_.data()[i]; }
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
        auto operator()(int k) const { return data_.template slice<2>(k).as_eigen_map(); }
        int size() const { return data_.size(); }
        // modifiers
        auto operator()(int i, int k) {
            auto slice_ = data_.template slice<2>(k);
            return slice_.as_eigen_map().block(i * blk_rows_, 0, blk_rows_, cols_);
        }
        auto operator()(int k) { return data_.template slice<2>(k).as_eigen_map(); }
    };

    template <typename DataLocs>
        requires(is_valid_data_locs_descriptor_v<DataLocs>)
    void eval_spatial_basis_at_(const DataLocs& locs) {
        if constexpr (std::is_same_v<DataLocs, matrix_t>) {   // pointwise sampling
            Psi_ = point_eval_(locs);
            D_ = vector_t::Ones(n_).asDiagonal();
        } else {   // areal sampling
            const auto& [psi, measure_vect] = areal_eval_(locs);
            Psi_ = psi;
            D_ = measure_vect.asDiagonal();
        }
        fdapde_assert(n_ == Psi_.rows());
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
            const auto& [psi, measure_vec] = areal_eval_(spatial_index.incidence_matrix());
            Psi_ = psi;
            D_ = measure_vec.asDiagonal();
            break;
        }
        }
	return;
    }
    // J(f,g) = \sum_{k=1}^m (y^(k) - \Psi * f^(k))^\top * (y^(k) - \Psi * f^(k)) + \lambda_D * (g^(k))^\top * (g^(k))
    double J_(const block_map_t& y, const block_map_t& x, double lambda) const {
        double sse = 0;
        for (int t = 0; t < m_; ++t) {
            sse += ((y(t) - Psi_ * x(0, t)).squaredNorm() / n_obs_ + lambda * x(1, t).squaredNorm());
        }
        return sse;
    }
   public:
    static constexpr int n_lambda = 2;
    using solver_category = ls_solver;

    fe_ls_parabolic_ieul() noexcept = default;
    template <typename GeoFrame, typename Penalty, typename WeightMatrix>
        requires(is_valid_penalty_v<Penalty>)
    fe_ls_parabolic_ieul(const std::string& formula, const GeoFrame& gf, Penalty&& penalty, const WeightMatrix& W) {
        fdapde_static_assert(GeoFrame::Order == 2, THIS_CLASS_IS_FOR_ORDER_TWO_GEOFRAMES_ONLY);
	fdapde_assert(gf.n_layers() == 1);
        n_obs_ = gf[0].rows();
	n_locs_ = n_obs_;

	discretize(penalty);
	analyze_data(formula, gf, W);
    }
    template <typename GeoFrame, typename Penalty>
        requires(is_valid_penalty_v<Penalty>)
    fe_ls_parabolic_ieul(const std::string& formula, const GeoFrame& gf, Penalty&& penalty) :
        fe_ls_parabolic_ieul(formula, gf, penalty, vector_t::Ones(gf[0].rows()).asDiagonal()) { }
    // construct with no data
    template <typename GeoFrame, typename Penalty, typename WeightMatrix>
        requires(is_valid_penalty_v<Penalty>)
    fe_ls_parabolic_ieul(const GeoFrame& gf, Penalty&& penalty, const WeightMatrix& W) {
        fdapde_static_assert(GeoFrame::Order == 2, THIS_CLASS_IS_FOR_ORDER_TWO_GEOFRAMES_ONLY);
	fdapde_assert(gf.n_layers() == 1);
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
	discretize(penalty);
	eval_spatial_basis_at_(gf);
    }
    template <typename GeoFrame, typename Penalty>
        requires(is_valid_penalty_v<Penalty>)
    fe_ls_parabolic_ieul(const GeoFrame& gf, Penalty&& penalty) :
        fe_ls_parabolic_ieul(gf, penalty, vector_t::Ones(gf[0].rows()).asDiagonal()) { }

    // finite element discretization of spatial dimension
    template <typename Penalty> void discretize(Penalty&& penalty) {
        using BilinearForm = typename std::decay_t<Penalty>::BilinearForm;
        using LinearForm = typename std::decay_t<Penalty>::LinearForm;
        fdapde_static_assert(
          internals::is_valid_penalty_pair_v<BilinearForm FDAPDE_COMMA LinearForm>, INVALID_PENALTY_DESCRIPTION);
        using FeSpace = typename BilinearForm::TrialSpace;
	// discretization
        const BilinearForm& bilinear_form = penalty.bilinear_form();
        const LinearForm& linear_form = penalty.linear_form();
        n_dofs_ = bilinear_form.n_dofs();   // number of basis functions over physical domain
        internals::fe_mass_assembly_loop<FeSpace> mass_assembler(bilinear_form.trial_space());
        R0_ = mass_assembler.assemble();
        R1_ = bilinear_form.assemble();
        u_  = linear_form.assemble();
	// store handles for basis system evaluation at locations
        point_eval_ = [fe_space = bilinear_form.trial_space()](const matrix_t& locs) -> decltype(auto) {
            return internals::point_basis_eval(fe_space, locs);
        };
        areal_eval_ = [fe_space = bilinear_form.trial_space()](const binary_t& locs) -> decltype(auto) {
            return internals::areal_basis_eval(fe_space, locs);
        };
	// store initial conditions, numerical scheme parameters
	s_ = penalty.ic();
	max_iter_ = penalty.max_iter();
	tol_ = penalty.tol();
        return;
    }
    // non-parametric fit
    // \sum_i w_i * (y_i - f(p_i))^2 + \int_D \int_T (df/dt + Lf - u)^2
    template <typename DataLocs1, typename DataLocs2, typename WeightMatrix>
        requires(is_valid_data_locs_descriptor_v<DataLocs1> && is_valid_data_locs_descriptor_v<DataLocs2>)
    void analyze_data(const DataLocs1& locs1, const DataLocs2 locs2, const matrix_t& y, const WeightMatrix& W) {
        fdapde_static_assert(
          std::is_same_v<DataLocs2 FDAPDE_COMMA matrix_t>,
          ITERATIVE_SEPARABLE_PENALIZATION_REQUIRES_POINT_TEMPORAL_EVALUATIONS_ONLY);
        fdapde_assert(
          locs1.rows() > 0 && locs2.rows() > 0 && y.rows() == locs1.rows() * locs2.rows() && y.cols() == 1 &&
          W.rows() == locs1.rows() * locs2.rows() && W.rows() == W.cols());
        n_obs_ = y.rows();
        n_locs_ = n_obs_;
        n_covs_ = 0;
        eval_spatial_basis_at_(locs1);

        m_ = locs2.rows();
        fdapde_assert(m_ > 0 && locs2.cols() == 1);
        DeltaT_ = locs2(1, 0) - locs2(0, 0);
        for (int i = 1; i < m_ - 1; ++i) {
            double lag_i = locs2(i + 1, 0) - locs2(i, 0);
            fdapde_assert(DeltaT_ > 0 && lag_i > 0 && almost_equal(DeltaT_ FDAPDE_COMMA lag_i));
        }
        // update forcing
	u0_ = u_ + (1.0 / DeltaT_) * (R0_ * s_);
	update_response_and_weights(y, W);
        return;
    }
    // fit from formula
    template <typename GeoFrame, typename WeightMatrix>
    void analyze_data(const std::string& formula, const GeoFrame& gf, const WeightMatrix& W) {
        fdapde_static_assert(GeoFrame::Order == 2, THIS_CLASS_IS_FOR_ORDER_TWO_GEOFRAMES_ONLY);
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
        // update forcing
	u0_ = u_ + (1.0 / DeltaT_) * (R0_ * s_);
	eval_spatial_basis_at_(gf);
    	// parse formula, extract response vector
        Formula formula_(formula);
        const auto& y_data = gf[0].data().template col<double>(formula_.lhs());
        y_.resize(n_locs_, y_data.blk_sz());
        y_data.assign_to(y_);
	
	update_response_and_weights(y_, W);
	return;
    }

    // modifiers
    void update_response(const vector_t& y) {
        fdapde_assert(Psi_.rows() > 0 && y.rows() == n_locs_ && y.cols() == 1);
        y_ = y;
	// correct \Psi for missing observations
        auto nan_pattern = na_matrix(y);
        if (nan_pattern.any()) {
            n_obs_ = n_locs_ - nan_pattern.count();
            B_.resize(m_);
            for (int i = 0; i < m_; ++i) {
                B_[i] = (~nan_pattern.middleRows(i * n_, n_)).repeat(1, n_dofs_).select(Psi_, 0);
            }
            y_ = (~nan_pattern).select(y_, 0);
        }	
    }
    template <typename WeightMatrix> void update_weights(const WeightMatrix& W) {
        fdapde_assert(Psi_.rows() > 0 && W.rows() == n_locs_ && W.rows() == W.cols());
        W_ = W;
        // check if W_ is time-wise block-constant
        W_const_ = true;
        for (int i = 1; i < m_; ++i) {
            if (sparse_matrix_t(W_.block(i * n_, i * n_, n_, n_) - W_.block(0, 0, n_, n_)).sum() != 0) {
                W_const_ = false;
                break;
            }
        }
        W_ /= n_obs_;
        W_changed_ = true;
        return;
    }
    template <typename WeightMatrix> void update_response_and_weights(const vector_t& y, const WeightMatrix& W) {
        update_response(y);
        update_weights (W);
        return;
    }
   private:
    // iterative scheme implementation
    template <typename ResponseT> auto fit_(ResponseT&& response, double lambda_D, double lambda_T) {
        std::array<double, n_lambda> lambda {lambda_D, lambda_T};
        // define auxiliary structures
        block_map_t y(response, n_);
        matrix_t x_old_buff(2 * n_dofs_ * m_, response.cols()), x_new_buff(2 * n_dofs_ * m_, response.cols());
        block_map_t x_old(x_old_buff, 2 * n_dofs_, n_dofs_, response.cols());
        block_map_t x_new(x_new_buff, 2 * n_dofs_, n_dofs_, response.cols());
        double alpha = lambda_D * lambda_T / DeltaT_;
        int n_fact = (B_.size() != 0 || !W_const_) ? m_ : 1;   // != 1 if missing or heteroschedastic observations
        auto PsiNA = [&](int t) -> const sparse_matrix_t& { return B_.size() != 0 ? B_[t] : Psi_; };
        auto invAs = [&](int t) -> sparse_solver_t& { return n_fact != 1 ? invAs_[t] : invAs_[0]; };
        auto invA  = [&](int t) -> sparse_solver_t& { return n_fact != 1 ? invA_ [t] : invA_ [0]; };
        auto W = [&](int i) { return W_.block(i * n_, i * n_, n_, n_); };
        {   // compute starting point (f^(k,0), g^(k,0)) k = 1 ... m
            if (lambda_saved_.value() != lambda || W_changed_) {
                invAs_.resize(n_fact);
                for (int t = 0; t < n_fact; ++t) {
                    SparseBlockMatrix<double, 2, 2> A_(
                      PsiNA(t).transpose() * D_ * W(t) * PsiNA(t), lambda_D * R1_.transpose(), lambda_D * R1_,
                      -lambda_D * R0_);
                    invAs(t).compute(A_);
                }
                sparse_matrix_t G0 = alpha * R0_.transpose() + lambda_D * R1_.transpose();
                invG0_.compute(G0);
            }
            vector_t b_(2 * n_dofs_);
            for (int t = 0; t < m_; ++t) {
                b_ << PsiNA(t).transpose() * D_ * W(t) * y(t), lambda_D * lambda_T * (t == 0 ? u0_ : u_);
                x_old(0, t) = invAs(t).solve(b_).head(n_dofs_);
            }
            b_ = PsiNA(m_ - 1).transpose() * D_ * W(m_ - 1) * (y(m_ - 1) - PsiNA(m_ - 1) * x_old(0, m_ - 1));
            x_old(1, m_ - 1) = invG0_.solve(b_);
            // general step
            for (int t = m_ - 2; t >= 0; --t) {
                b_ << PsiNA(t).transpose() * D_ * W(t) * (y(t) - PsiNA(t) * x_old(0, t)) +
                        alpha * R0_ * x_old(1, t + 1);
                x_old(1, t) = invG0_.solve(b_);
            }
        }
        // iterative scheme initialization
        double Jold = std::numeric_limits<double>::max();
        double Jnew = J_(y, x_old, lambda_D);
        int i = 1;
        if (lambda_saved_.value() != lambda || W_changed_) {
            invA_.resize(n_fact);
            for (int t = 0; t < n_fact; ++t) {
                SparseBlockMatrix<double, 2, 2> A_(
                  PsiNA(t).transpose() * D_ * W(t) * PsiNA(t), lambda_D * R1_.transpose() + alpha * R0_.transpose(),
                  lambda_D * R1_ + alpha * R0_, -lambda_D * R0_);
                invA(t).compute(A_);
            }
        }
        vector_t b_(2 * n_dofs_);
        // iterative loop
        while (i < max_iter_ && std::abs((Jnew - Jold) / Jnew) > tol_) {
            // at step 0, f^(k-1,i-1) is zero
            b_ << PsiNA(0).transpose() * D_ * W(0) * y(0) + alpha * R0_ * x_old(1, 1), lambda_D * u0_;
            x_new(0) = invA(0).solve(b_);
            // general step
            for (int t = 1; t < m_ - 1; ++t) {
                b_ << PsiNA(t).transpose() * D_ * W(t) * y(t) + alpha * R0_ * x_old(1, t + 1),
                  alpha * R0_ * x_old(0, t - 1) + lambda_D * u_;
                x_new(t) = invA(t).solve(b_);
            }
            // at step m_ - 1, g^(k+1,i-1) is zero
            b_ << PsiNA(m_ - 1).transpose() * D_ * W(m_ - 1) * y(m_ - 1),
              alpha * R0_ * x_old(0, m_ - 2) + lambda_D * u_;
            x_new(m_ - 1) = invA(m_ - 1).solve(b_);
            // prepare for next iteration
            Jold = Jnew;
            x_old = x_new;
            Jnew = J_(y, x_new, lambda_D);
            i++;
        }
        // return result
        f_.resize(n_dofs_ * m_, response.cols());
        g_.resize(n_dofs_ * m_, response.cols());
        for (int i = 0; i < m_; ++i) {
            f_.middleRows(i * n_dofs_, n_dofs_) = x_new(0, i);
            g_.middleRows(i * n_dofs_, n_dofs_) = x_new(1, i);
        }
        lambda_saved_ = std::array<double, n_lambda> {lambda_D, lambda_T};
        return std::make_pair(f_, g_);
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
    const vector_t& initial_condition() const { return s_; }
   protected:
    std::optional<std::array<double, n_lambda>> lambda_saved_ = std::array<double, n_lambda> {-1, -1};
    std::vector<sparse_solver_t> invA_, invAs_;
    sparse_solver_t invG0_;
    // matrices for hutchinson stochastic estimation of Tr[S]
    std::optional<matrix_t> Us_;

    int n_dofs_ = 0, n_obs_ = 0, n_covs_ = 0;
    int n_locs_ = 0, n_ = 0, m_ = 0;   // n_: number of spatial locations, m_: number of time instants

    sparse_matrix_t R0_;               // n_dofs x n_dofs matrix [R0]_{ij} = \int_D \psi_i * \psi_j
    sparse_matrix_t R1_;               // n_dofs x n_dofs matrix [R1]_{ij} = \int_D a(\psi_i, \psi_j)
    sparse_matrix_t Psi_;              // n_obs x n_dofs matrix [Psi]_{ij} = \psi_j(p_i)
    std::vector<sparse_matrix_t> B_;   // m x (n_obs x n_obs) vector of na-corrected \Psi matrices
    vector_t u_, u0_;                  // n_dofs x 1 vector [u]_i = \int_D u * \psi_i, u0_ = u_ + (R0 * s) / DeltaT
    diag_matrix_t D_;                  // vector of regions' measures (areal sampling)
    mutable std::optional<sparse_solver_t> invR0_;
    vector_t f_, beta_, g_;
    vector_t s_;   // initial condition vector
    // basis system evaluation handles
    std::function<sparse_matrix_t(const matrix_t& locs)> point_eval_;
    std::function<std::pair<sparse_matrix_t, vector_t>(const binary_t& locs)> areal_eval_;

    vector_t y_;          // n_obs x 1 observation vector
    sparse_matrix_t W_;   // n_obs x n_obs matrix of observation weights
    bool W_changed_, W_const_;   // W_const_ == true \iff W_ is time-wise constant

    double tol_;      // convergence tolerance
    int max_iter_;    // maximum number of iterations
    double DeltaT_;   // time step
};

}   // namespace internals

// parabolic solver API
// monolithic method
template <typename Penalty>
    requires(internals::is_pair_v<Penalty>)
struct fe_ls_parabolic_mono {
    using solver_t = internals::fe_ls_parabolic_mono;
   private:
    struct penalty_packet {
        using BilinearForm = std::tuple_element_t<0, std::decay_t<Penalty>>;
        using LinearForm = std::tuple_element_t<1, std::decay_t<Penalty>>;
       private:
        BilinearForm bilinear_form_;
        LinearForm linear_form_;
        Eigen::Matrix<double, Dynamic, 1> ic_;
       public:
        penalty_packet(
          const BilinearForm& bilinear_form, const LinearForm& linear_form,
          const Eigen::Matrix<double, Dynamic, 1>& ic) :
            bilinear_form_(bilinear_form), linear_form_(linear_form), ic_(ic) { }
        // observers
        const BilinearForm& bilinear_form() const { return bilinear_form_; }
        const LinearForm& linear_form() const { return linear_form_; }
        const Eigen::Matrix<double, Dynamic, 1>& ic() const { return ic_; }
    };
   public:
    template <typename InitialCondition>
    fe_ls_parabolic_mono(const Penalty& penalty, const InitialCondition& ic) :
        penalty_(std::get<0>(penalty), std::get<1>(penalty), ic) { }
    const penalty_packet& get() const { return penalty_; }
   private:
    penalty_packet penalty_;
};
// implicit euler time integration method
template <typename Penalty>
    requires(internals::is_pair_v<Penalty>)
struct fe_ls_parabolic_ieul {
    using solver_t = internals::fe_ls_parabolic_ieul;
   private:
    struct penalty_packet {
        using BilinearForm = std::tuple_element_t<0, std::decay_t<Penalty>>;
        using LinearForm = std::tuple_element_t<1, std::decay_t<Penalty>>;
       private:
        BilinearForm bilinear_form_;
        LinearForm linear_form_;
        Eigen::Matrix<double, Dynamic, 1> ic_;
        int max_iter_ = 50;
        double tol_ = 1e-4;
       public:
        penalty_packet(
          const BilinearForm& bilinear_form, const LinearForm& linear_form, const Eigen::Matrix<double, Dynamic, 1>& ic,
          int max_iter, double tol) :
            bilinear_form_(bilinear_form), linear_form_(linear_form), ic_(ic), max_iter_(max_iter), tol_(tol) { }
        // observers
        const BilinearForm& bilinear_form() const { return bilinear_form_; }
        const LinearForm& linear_form() const { return linear_form_; }
        const Eigen::Matrix<double, Dynamic, 1>& ic() const { return ic_; }
        int max_iter() const { return max_iter_; }
        double tol() const { return tol_; }
    };
   public:
    template <typename InitialCondition>
    fe_ls_parabolic_ieul(const Penalty& penalty, const InitialCondition& ic, int max_iter = 50, double tol = 1e-4) :
        penalty_(std::get<0>(penalty), std::get<1>(penalty), ic, max_iter, tol) { }
    const penalty_packet& get() const { return penalty_; }
   private:
    penalty_packet penalty_;
};

}   // namespace fdapde

#endif // __FE_LS_PARABOLIC_SOLVER_H__

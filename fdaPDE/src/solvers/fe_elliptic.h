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

#ifndef __FE_ELLIPTIC_SOLVER_H__
#define __FE_ELLIPTIC_SOLVER_H__

#include "header_check.h"

namespace fdapde {
namespace internals {

// solves \min_{f, \beta} \| W^{1/2} * (y_i - x_i^\top * \beta - f(p_i)) \|_2^2 + \int_D (Lf - u)^2, L elliptic operator
struct fe_elliptic_solver {
   private:
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using binary_t = BinaryMatrix<Dynamic, Dynamic>;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using diag_matrix_t   = Eigen::DiagonalMatrix<double, Dynamic, Dynamic>;
    using sparse_solver_t = eigen_sparse_solver_movable_wrap<Eigen::SparseLU<sparse_matrix_t>>;
    using dense_solver_t  = Eigen::PartialPivLU<matrix_t>;

    // evaluation of basis system at spatial locations
    template <typename DataLocs>
        requires(std::is_same_v<DataLocs, matrix_t> || std::is_same_v<DataLocs, binary_t>)
    void eval_basis_at_(const DataLocs& locs) {
        if constexpr (std::is_same_v<DataLocs, matrix_t>) {   // pointwise sampling
            Psi_ = point_eval_(locs);
            D_ = vector_t::Ones(n_locs_).asDiagonal();
	    fdapde_assert(n_locs_ == Psi_.rows());
        } else {   // areal sampling
            const auto& [psi, measure_vect] = areal_eval_(locs);
            Psi_ = psi;
            D_ = measure_vect.asDiagonal();
	    fdapde_assert(n_locs_ == Psi_.rows());
        }
        return;
    }
   public:
    static constexpr int n_lambda = 1;

    fe_elliptic_solver() noexcept = default;
    // construct from formula + geoframe
    template <typename GeoFrame, typename Penalty, typename WeightMatrix>
        requires(internals::is_pair_v<Penalty>)
    fe_elliptic_solver(const std::string& formula, const GeoFrame& gf, Penalty&& penalty, const WeightMatrix& W) :
        W_(W) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        fdapde_assert(gf.n_layers() == 1);
        n_obs_  = gf[0].rows();
	n_locs_ = n_obs_;
	
        discretize(penalty);
	analyze_data(formula, gf, W);
    }
    template <typename GeoFrame, typename Penalty>
        requires(internals::is_pair_v<Penalty>)
    fe_elliptic_solver(const std::string& formula, const GeoFrame& gf, Penalty&& penalty) :
        fe_elliptic_solver(formula, gf, penalty, vector_t::Ones(gf[0].rows()).asDiagonal()) { }
    // construct with no data
    template <typename GeoFrame, typename Penalty, typename WeightMatrix>
        requires(internals::is_pair_v<Penalty>)
    fe_elliptic_solver(const GeoFrame& gf, Penalty&& penalty, const WeightMatrix& W) : W_(W) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        fdapde_assert(gf.n_layers() == 1);
        n_obs_  = gf[0].rows();
	n_locs_ = n_obs_;
	
        discretize(penalty);
        // basis system evaluation
        switch (gf.category(0)[0]) {
        case ltype::point: {
            const auto& spatial_index = geo_index_cast<0, POINT>(gf[0]);
            if (spatial_index.points_at_dofs()) {
                Psi_.resize(n_locs_, n_dofs_);
                Psi_.setIdentity();
            } else {
                Psi_ = point_eval_(spatial_index.coordinates());
            }
            D_ = vector_t::Ones(n_locs_).asDiagonal();
            break;
        }
        case ltype::areal: {
            const auto& spatial_index = geo_index_cast<0, POLYGON>(gf[0]);
            const auto& [psi, measure_vect] = areal_eval_(spatial_index.incidence_matrix());
            Psi_ = psi;
            D_ = measure_vect.asDiagonal();
            break;
        }
        }
    }
    template <typename GeoFrame, typename Penalty>
        requires(internals::is_pair_v<Penalty>)
    fe_elliptic_solver(const GeoFrame& gf, Penalty&& penalty) :
        fe_elliptic_solver(gf, penalty, vector_t::Ones(gf[0].rows()).asDiagonal()) { }

    // perform finite element based numerical discretization
    template <typename Penalty> void discretize(Penalty&& penalty) {
        fdapde_static_assert(internals::is_valid_penalty_pair_v<Penalty>, INVALID_PENALTY_DESCRIPTION);
        using BilinearForm = std::tuple_element_t<0, std::decay_t<Penalty>>;
        using LinearForm = std::tuple_element_t<1, std::decay_t<Penalty>>;
        using FeSpace = typename BilinearForm::TrialSpace;
	// discretization
        const BilinearForm& bilinear_form = std::get<0>(penalty);
        const LinearForm& linear_form = std::get<1>(penalty);
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
	b_.resize(2 * n_dofs_, 1);
        return;
    }
    // non-parametric fit
    // \sum_i w_i * (y_i - f(p_i))^2 + \int_D (Lf - u)^2
    template <typename DataLocs, typename WeightMatrix>
        requires(std::is_same_v<DataLocs, matrix_t> || std::is_same_v<DataLocs, binary_t>)
    void analyze_data(const DataLocs& locs, const matrix_t& y, const WeightMatrix& W) {
        fdapde_assert(
          locs.rows() > 0 && y.rows() == locs.rows() && y.cols() == 1 && W.rows() == locs.rows() &&
          W.rows() == W.cols());
        n_obs_  = locs.rows();
	n_locs_ = n_obs_;
	n_covs_ = 0;
        eval_basis_at_(locs);   // update \Psi matrix
        update_response_and_weights(y, W);
        return;
    }
    // semi-parametric fit
    // \sum_i w_i * (y_i - x_i^\top * \beta - f(p_i))^2 + \int_D (Lf - u)^2
    template <typename DataLocs, typename WeightMatrix>
        requires(std::is_same_v<DataLocs, matrix_t> || std::is_same_v<DataLocs, binary_t>)
    void analyze_data(const DataLocs& locs, const matrix_t& y, const matrix_t& X, const WeightMatrix& W) {
        fdapde_assert(
          locs.rows() > 0 && y.rows() == locs.rows() && y.cols() == 1 && X.rows() == locs.rows() &&
          W.rows() == locs.rows() && W.rows() == W.cols());
        n_obs_  = locs.rows();
	n_locs_ = n_obs_;
        bool require_woodbury_realloc = n_covs_ != X.cols();
        n_covs_ = X.cols();
        eval_basis_at_(locs);   // update \Psi matrix
        if (require_woodbury_realloc) { U_ = matrix_t::Zero(2 * n_dofs_, n_covs_); }
        if (require_woodbury_realloc) { V_ = matrix_t::Zero(n_covs_, 2 * n_dofs_); }
        update_response_and_weights(y, X, W);
        return;
    }
    // fit from formula
    template <typename GeoFrame, typename WeightMatrix>
    void analyze_data(const std::string& formula, const GeoFrame& gf, const WeightMatrix& W) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        fdapde_assert(gf.n_layers() == 1);
        n_obs_  = gf[0].rows();
	n_locs_ = n_obs_;
        // basis system evaluation
        switch (gf.category(0)[0]) {
        case ltype::point: {
            const auto& spatial_index = geo_index_cast<0, POINT>(gf[0]);	  
            if (spatial_index.points_at_dofs()) {
                Psi_.resize(n_locs_, n_dofs_);
                Psi_.setIdentity();
            } else {
                Psi_ = point_eval_(spatial_index.coordinates());
            }
            D_ = vector_t::Ones(n_locs_).asDiagonal();
            break;
        }
        case ltype::areal: {
            const auto& spatial_index = geo_index_cast<0, POLYGON>(gf[0]);
            const auto& [psi, measure_vect] = areal_eval_(spatial_index.incidence_matrix());
            Psi_ = psi;
            D_ = measure_vect.asDiagonal();
            break;
        }
        }
        // parse formula, extract response vector and design matrix
        Formula formula_(formula);
        std::vector<std::string> covs;
        for (const std::string& token : formula_.rhs()) {
            if (gf.contains(token)) { covs.push_back(token); }
        }
	bool require_woodbury_realloc = n_covs_ != covs.size();
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
        update_response_and_weights(y_, W);   // this updates also design_matrix releated matrices
        return;
    }
  
    // modifiers
    void update_response(const vector_t& y) {
        fdapde_assert(Psi_.rows() > 0 && y.rows() == n_locs_ && y.cols() == 1);
        y_ = y;
        // correct \Psi for missing observations
        auto nan_pattern = na_matrix(y);
	int old_n_obs = n_obs_;
        if (nan_pattern.any()) {
            n_obs_ = n_locs_ - nan_pattern.count();
            B_ = nan_pattern.repeat(1, n_dofs_).select(Psi_, 0);
            y_ = nan_pattern.select(y_, 0);
        }
        if (old_n_obs != n_obs_) { W_ *= (double)old_n_obs / n_obs_; }
        b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * W_ * y;
        return;
    }
    template <typename WeightMatrix> void update_weights(const WeightMatrix& W) {
        fdapde_assert(Psi_.rows() > 0 && W.rows() == n_locs_ && W.rows() == W.cols());
        W_ = W;
	W_ /= n_obs_;
        if (n_covs_ == 0) {
            b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * W_ * y_;
        } else {
            XtWX_ = X_.transpose() * W_ * X_;
            invXtWX_ = XtWX_.partialPivLu();
            invXtWXXtW_ = invXtWX_.solve(X_.transpose() * W_);   // (X^\top * W * X)^{-1} * (X^\top * W)
            // woodbury decomposition matrices
            U_.block(0, 0, n_dofs_, n_covs_) = Psi_.transpose() * D_ * W_ * X_;
            V_.block(0, 0, n_covs_, n_dofs_) = X_.transpose() * W_ * Psi_;
            b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * internals::lmbQ(W_, X_, invXtWX_, y_);
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
            B_ = nan_pattern.repeat(1, n_dofs_).select(Psi_, 0);
            y_ = nan_pattern.select(y_, 0);
        }
        update_weights(W);
        return;
    }

    // main fit entry point
    std::pair<vector_t, vector_t> fit(double lambda) {
        fdapde_assert(lambda > 0 && n_dofs_ > 0 && n_obs_ > 0);
        if (lambda_saved_.value() != lambda || W_changed_) {
            // assemble and factorize system matrix for nonparameteric part
            SparseBlockMatrix<double, 2, 2> A(
              -Psi_.transpose() * D_ * W_ * Psi_, lambda * R1_.transpose(), lambda * R1_, lambda * R0_);
            invA_.compute(A);
	    W_changed_ = false;
        }
        if (lambda_saved_.value() != lambda) {
            // update linear system rhs
            b_.block(n_dofs_, 0, n_dofs_, 1) = lambda * u_;
        }
        lambda_saved_ = lambda;
        vector_t x;
        if (n_covs_ == 0) {
            x = invA_.solve(b_);
            f_ = x.topRows(n_dofs_);
        } else {
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
        return fit(lambda[0]);
    }
    // perform a nonparametric_fit, e.g. discarding possible covariates
    vector_t nonparametric_fit(double lambda) {
        fdapde_assert(lambda > 0 && n_dofs_ > 0 && n_obs_ > 0);
        if (lambda_saved_.value() != lambda) {
            // assemble and factorize system matrix for nonparameteric part
            SparseBlockMatrix<double, 2, 2> A(
              -Psi_.transpose() * D_ * W_ * Psi_, lambda * R1_.transpose(), lambda * R1_, lambda * R0_);
            invA_.compute(A);
        }
        vector_t x;
        if (n_covs_ == 0) {   // equivalent to calling fit(lambda)
            if (lambda_saved_.value() != lambda) { b_.block(n_dofs_, 0, n_dofs_, 1) = lambda * u_; }
            x = invA_.solve(b_);
        } else {
            vector_t b(2 * n_dofs_);
            // assemble nonparametric linear system rhs
            b.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * W_ * y_;
            b.block(n_dofs_, 0, n_dofs_, 1) = lambda * u_;   
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
            Bs_->topRows(n_dofs_) = -Psi_.transpose() * D_ * W_ * (*Us_);
        } else {
            Bs_->topRows(n_dofs_) = -Psi_.transpose() * D_ * internals::lmbQ(W_, X_, invXtWX_, *Us_);
        }
        matrix_t x = n_covs_ == 0 ? invA_.solve(*Bs_) : woodbury_system_solve(invA_, U_, XtWX_, V_, *Bs_);
        double trS = 0;   // monte carlo Tr[S] approximation
        for (int i = 0; i < r; ++i) { trS += Ys_->row(i).dot(x.col(i).head(n_dofs_)); }
        return trS / r;
    }
    template <typename LambdaT>
        requires(internals::is_vector_like_v<LambdaT> || std::is_floating_point_v<LambdaT>)
    double edf(const LambdaT& lambda, int r = 100, int seed = random_seed) {
        double lambda_;
        if constexpr (internals::is_vector_like_v<LambdaT>) {
            fdapde_assert(lambda.size() == n_lambda && lambda[0] > 0);
            lambda_ = lambda[0];
        } else {
            fdapde_assert(lambda > 0);
            lambda_ = lambda;
        }
        if (lambda_saved_.value() != lambda_) {
            SparseBlockMatrix<double, 2, 2> A_(
              -Psi_.transpose() * D_ * W_ * Psi_, lambda_ * R1_.transpose(), lambda_ * R1_, lambda_ * R0_);
            invA_.compute(A_);
            lambda_saved_ = lambda_;
        }
        return edf(r, seed);
    }
    // penalty matrix: \lambda * R1^\top * (R0)^{-1} * R1
    matrix_t P(double lambda) const {
        if (!invR0_.has_value()) { invR0_.compute(R0_); }
        return lambda * R1_.transpose() * invR0_.solve(R1_);
    }
    template <typename LambdaT>
        requires(internals::is_vector_like_v<LambdaT>)
    matrix_t P(const LambdaT& lambda) const {
        fdapde_assert(lambda.size() == n_lambda);
        return P(lambda[0]);
    }
    matrix_t P() const { return P(1.0); }
    // efficient evaluation of f^\top * P * f = g^\top * R0 * g
    double ftPf(double lambda) {
        if (lambda_saved_.value() != lambda || W_changed_) { fit(lambda); }
        return lambda * g_.dot(R0_ * g_);
    }
    template <typename LambdaT>
        requires(internals::is_vector_like_v<LambdaT>)
    double ftPf(const LambdaT& lambda) {
        fdapde_assert(lambda.size() == n_lambda);
        return ftPf(lambda[0]);
    }
    // left multiplication by \Psi
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
    std::optional<double> lambda_saved_ = -1;
    sparse_solver_t invA_;
    matrix_t b_;
    // matrices for hutchinson stochastic estimation of Tr[S]
    std::optional<matrix_t> Ys_;
    std::optional<matrix_t> Bs_;
    std::optional<matrix_t> Us_;

    int n_dofs_ = 0, n_locs_ = 0, n_obs_ = 0, n_covs_ = 0;
    sparse_matrix_t R0_;    // n_dofs x n_dofs matrix [R0]_{ij} = \int_D \psi_i * \psi_j
    sparse_matrix_t R1_;    // n_dofs x n_dofs matrix [R1]_{ij} = \int_D a(\psi_i, \psi_j)
    sparse_matrix_t Psi_;   // n_obs x n_dofs matrix [Psi]_{ij} = \psi_j(p_i)
    vector_t u_;            // n_dofs x 1 vector u_i = \int_D u * \psi_i
    diag_matrix_t D_;       // vector of regions' measures (areal sampling)
    mutable sparse_solver_t invR0_;
    std::optional<sparse_matrix_t> B_;   // \Psi matrix corrected for missing observations
    vector_t f_, beta_, g_;
    // basis system evaluation handles
    std::function<sparse_matrix_t(const matrix_t& locs)> point_eval_;
    std::function<std::pair<sparse_matrix_t, vector_t>(const binary_t& locs)> areal_eval_;

    matrix_t X_;               // n_obs x n_covs design matrix
    vector_t y_;               // n_obs x 1 observation vector
    sparse_matrix_t W_;        // n_obs x n_obs matrix of observation weights
    matrix_t U_, V_;           // (2 * n_dofs) x n_covs matrices [\Psi^\top * D * W * y, 0] and [X^\top * W * \Psi, 0]
    matrix_t XtWX_;            // n_covs x n_covs matrix X^\top * W * X
    dense_solver_t invXtWX_;   // factorization of n_covs x n_covs matrix X^\top * W * X
    matrix_t invXtWXXtW_;      // n_covs x n_obs matrix (X^\top * X)^{-1} * (X^\top W)
    bool W_changed_;
};

}   // namespace internals

// general non-parametrized elliptic solver factory method
template <typename BilinearForm_, typename LinearForm_> struct fe_elliptic_penalty {
    using BilinearForm = std::decay_t<BilinearForm_>;
    using LinearForm = std::decay_t<LinearForm_>;
    using Triangulation = std::tuple<typename BilinearForm::Triangulation>;
    using solver_t = internals::fe_elliptic_solver;
    fdapde_static_assert(
      std::is_same_v<typename BilinearForm::discretization_category FDAPDE_COMMA finite_element_tag>&&
        std::is_same_v<typename LinearForm::discretization_category FDAPDE_COMMA finite_element_tag>,
      FE_ELLIPTIC_PENALTY_IS_FOR_FINITE_ELEMENT_DISCRETIZATIONS_ONLY);

    fe_elliptic_penalty(const BilinearForm_& bilinear_form, const LinearForm_& linear_form) :
        penalty_(std::make_pair(bilinear_form, linear_form)) { }
    const std::tuple<BilinearForm, LinearForm>& get() const { return penalty_; }
   private:
    std::tuple<BilinearForm, LinearForm> penalty_;
};
template <typename BilinearForm, typename LinearForm>
    requires(internals::is_bilinear_form_v<BilinearForm> && internals::is_linear_form_v<LinearForm>)
fe_elliptic_penalty<BilinearForm, LinearForm>
fe_elliptic(const BilinearForm& bilinear_form, const LinearForm& linear_form) {
    return fe_elliptic_penalty(bilinear_form, linear_form);
}
template <typename BilinearForm>
    requires(internals::is_bilinear_form_v<BilinearForm>)
auto fe_elliptic(const BilinearForm& bilinear_form) {   // implicit homogeneous forcing
    using FeSpace = typename BilinearForm::TrialSpace;
    const FeSpace& Vh = bilinear_form.trial_space();
    static constexpr int embed_dim = FeSpace::embed_dim;
    ScalarField<embed_dim, decltype([](const Eigen::Matrix<double, embed_dim, 1>&) { return 0; })> u;
    TestFunction v(Vh);
    auto linear_form = integral(Vh.triangulation())(u * v);
    return fe_elliptic(bilinear_form, linear_form);
}

// catalogue of standard elliptic penalizations
template <typename Functor> struct fe_elliptic_factory {
    using solver_t = internals::fe_elliptic_solver;
    fe_elliptic_factory(const Functor& f) : f_(f) { }
    template <typename Triangulation> auto operator()(const Triangulation& D) const { return f_(D); }
  private:
    Functor f_;
};
// the laplace equation: -\Delta f = 0
auto fe_laplace() {
    return fe_elliptic_factory([]<typename Triangulation>(const Triangulation& D) {
        auto Vh = std::make_shared<FeSpace<Triangulation, FeP<1, 1>>>(D, P1<1>);
        TrialFunction f(Vh);
        TestFunction  v(Vh);
        auto a = integral(D)(dot(grad(f), grad(v)));
        return fe_elliptic(a);
    });
}
// the poisson equation: -\Delta f = u
template <typename Force> auto fe_poisson(Force&& u) {
    return fe_elliptic_factory([u]<typename Triangulation>(const Triangulation& D) {
        auto Vh = std::make_shared<FeSpace<Triangulation, FeP<1, 1>>>(D, P1<1>);
        TrialFunction f(Vh);
        TestFunction  v(Vh);
        auto a = integral(D)(dot(grad(f), grad(v)));
	auto F = integral(D)(u * v);
        return fe_elliptic(a, F);
    });
}
// the general homogeneous diffusion-transport-reaction equation: -div[K + grad(f)] + b \cdot grad(f) + c * f = 0
template <typename Diffusion, typename Transport, typename Reaction>
auto fe_diffusion_transport_reaction(Diffusion&& K, Transport&& b, Reaction&& c) {
    return fe_elliptic_factory([K, b, c]<typename Triangulation>(const Triangulation& D) {
        auto Vh = std::make_shared<FeSpace<Triangulation, FeP<1, 1>>>(D, P1<1>);
        TrialFunction f(Vh);
        TestFunction  v(Vh);
        auto a = integral(D)(dot(K * grad(f), grad(v)) + dot(b, grad(f)) * v + c * f * v);
        return fe_elliptic(a);
    });
}
// the general non-homogeneous diffusion-transport-reaction equation: -div[K + grad(f)] + b \cdot grad(f) + c * f = u
template <typename Diffusion, typename Transport, typename Reaction, typename Force>
auto fe_diffusion_transport_reaction(Diffusion&& K, Transport&& b, Reaction&& c, Force&& u) {
  return fe_elliptic_factory([K, b, c, u]<typename Triangulation>(const Triangulation& D) {
        auto Vh = std::make_shared<FeSpace<Triangulation, FeP<1, 1>>>(D, P1<1>);
        TrialFunction f(Vh);
        TestFunction  v(Vh);
        auto a = integral(D)(dot(K * grad(f), grad(v)) + dot(b, grad(f)) * v + c * f * v);
        auto F = integral(D)(u * v);
        return fe_elliptic(a, F);
    });
}

}   // namespace fdapde

#endif // __FE_ELLIPTIC_SOLVER_H__

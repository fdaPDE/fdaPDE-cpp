#ifndef __FDAPDE_LS_GRAPH_H__
#define __FDAPDE_LS_GRAPH_H__

#include "header_check.h"

namespace fdapde {
namespace internals {

// Spline-only smoother: solves (Psi^T W Psi + lambda * R1) f = Psi^T W y
struct ls_graph {
   private:
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using binary_t = BinaryMatrix<Dynamic, Dynamic>;
    using diag_matrix_t   = Eigen::DiagonalMatrix<double, Dynamic, Dynamic>;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using sparse_solver_t = eigen_sparse_solver_movable_wrap<Eigen::SparseLU<sparse_matrix_t>>;
    using dense_solver_t  = Eigen::PartialPivLU<matrix_t>;

    template <typename DataLocs>
    static constexpr bool is_valid_data_locs_descriptor_v =
      std::is_same_v<DataLocs, matrix_t> || std::is_same_v<DataLocs, binary_t>;
    template <typename Penalty> struct is_valid_penalty {
        static constexpr bool value = requires(Penalty penalty) {
            penalty.bilinear_form();
            penalty.linear_form();
        };
    };
    template <typename Penalty> static constexpr bool is_valid_penalty_v = is_valid_penalty<Penalty>::value;

    // evaluate basis at locations (point or areal)
    template <typename DataLocs>
        requires(is_valid_data_locs_descriptor_v<DataLocs>)
    void eval_basis_at_(const DataLocs& locs) {
        fdapde_assert(n_locs_ == locs.rows());
        if constexpr (std::is_same_v<DataLocs, matrix_t>) {
            Psi_ = point_eval_(locs);
            D_ = vector_t::Ones(n_locs_).asDiagonal();
        } else {
            auto pr = areal_eval_(locs);
            Psi_ = pr.first;
            D_ = pr.second.asDiagonal();
        }
    }
    // optimized basis evaluation at geoframe
    template <typename GeoFrame> void eval_basis_at_(const GeoFrame& gf) {
        switch (gf.category(0)[0]) {
            case ltype::point: {
                const auto& spatial_index = geo_index_cast<0, POINT>(gf[0]);
                Psi_ = point_eval_(spatial_index.coordinates());
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
        return;
    }

   public:
    static constexpr int n_lambda = 1;
    using solver_category = ls_solver;

    ls_graph() noexcept = default;
    // construct from formula + geoframe
    template <typename GeoFrame, typename Penalty, typename WeightMatrix>
        requires(is_valid_penalty_v<Penalty>)
    ls_graph(const std::string& formula, const GeoFrame& gf, Penalty&& penalty, const WeightMatrix& W) : W_(W) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        fdapde_assert(gf.n_layers() == 1);
        n_obs_  = gf[0].rows();
        n_locs_ = n_obs_;

        discretize(penalty);
        analyze_data(formula, gf, W);
    }
    template <typename GeoFrame, typename Penalty>
        requires(is_valid_penalty_v<Penalty>)
    ls_graph(const std::string& formula, const GeoFrame& gf, Penalty&& penalty) :
        ls_graph(formula, gf, penalty, vector_t::Ones(gf[0].rows()).asDiagonal()) { }
    // construct with no data
    template <typename GeoFrame, typename Penalty, typename WeightMatrix>
        requires(is_valid_penalty_v<Penalty>)
    ls_graph(const GeoFrame& gf, Penalty&& penalty, const WeightMatrix& W) : W_(W) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        fdapde_assert(gf.n_layers() == 1);
        n_obs_  = gf[0].rows();
        n_locs_ = n_obs_;

        discretize(penalty);
        eval_basis_at_(gf);
    }
    template <typename GeoFrame, typename Penalty>
        requires(is_valid_penalty_v<Penalty>)
    ls_graph(const GeoFrame& gf, Penalty&& penalty) :
        ls_graph(gf, penalty, vector_t::Ones(gf[0].rows()).asDiagonal()) { }

    // discretization: assemble R1 and basis evaluation handles
    template <typename Penalty> void discretize(Penalty&& penalty) {
        using BilinearForm = typename std::decay_t<Penalty>::BilinearForm;
        using LinearForm = typename std::decay_t<Penalty>::LinearForm;
        fdapde_static_assert(internals::is_valid_penalty_pair_v<BilinearForm FDAPDE_COMMA LinearForm>, INVALID_PENALTY_DESCRIPTION);
        using BsSpace      = typename BilinearForm::TrialSpace;
        // discretization
        const BilinearForm& bilinear_form = penalty.bilinear_form();
        const LinearForm&   linear_form   = penalty.linear_form();
        n_dofs_ = bilinear_form.n_dofs();
        // assemble penalty (stiffness) R1 and forcing u (if provided)
        auto& space = bilinear_form.trial_space();
        TrialFunction u(space);
        TestFunction  v(space);
        R0_ = integral(space.triangulation())(u * v).assemble();
        R1_ = bilinear_form.assemble();
        u_  = linear_form.assemble();

        // store handles to evaluate basis at locations
        point_eval_ = [bs_space = bilinear_form.trial_space()](const matrix_t& locs) -> decltype(auto) {
            return internals::point_basis_eval(bs_space, locs);
        };
        areal_eval_ = [bs_space = bilinear_form.trial_space()](const binary_t& locs) -> decltype(auto) {
            return internals::areal_basis_eval(bs_space, locs);
        };
        // preallocate
        f_.resize(n_dofs_);
        return;
    }

    // non-parametric fit
    // \sum_i w_i * (y_i - f(p_i))^2 + \int_D (Lf - u)^2
    template <typename DataLocs, typename WeightMatrix>
        requires(is_valid_data_locs_descriptor_v<DataLocs>)
    void analyze_data(const DataLocs& locs, const matrix_t& y, const WeightMatrix& W) {
        fdapde_assert(locs.rows() > 0 && y.rows() == locs.rows() && y.cols() == 1);
        n_obs_ = locs.rows();
        n_locs_ = n_obs_;
        n_covs_ = 0;
        eval_basis_at_(locs);
        update_response_and_weights(y, W);
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
        if (require_woodbury_realloc) { U_ = matrix_t::Zero(n_dofs_, n_covs_); }
        if (require_woodbury_realloc) { V_ = matrix_t::Zero(n_covs_, n_dofs_); }
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
        eval_basis_at_(gf);   // update \Psi matrix

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

        if (b_.cols() != y_.cols()) { b_.resize(n_dofs_, y_.cols()); }
        if (n_covs_ != 0) {
            if (require_woodbury_realloc) { U_ = matrix_t::Zero(n_dofs_, n_covs_); }
            if (require_woodbury_realloc) { V_ = matrix_t::Zero(n_covs_, n_dofs_); }
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
        nan_pattern_ = na_matrix(y);
        int old_n_obs = n_obs_;
        if (nan_pattern_.any()) {
            n_obs_ = n_locs_ - nan_pattern_.count();
            B_ = (~nan_pattern_).repeat(1, n_dofs_).select(Psi_, 0);
            y_ = (~nan_pattern_).select(y_, 0);
        }
        if (old_n_obs != n_obs_) { W_ *= (double)old_n_obs / n_obs_; }
        b_ = PsiNA().transpose() * D_ * W_ * y_;
        return;
    }
    template <typename WeightMatrix> void update_weights(const WeightMatrix& W) {
        fdapde_assert(Psi_.rows() > 0 && W.rows() == n_locs_ && W.rows() == W.cols());
        W_ = W;
        W_ /= n_obs_;
        if (n_covs_ == 0) {
            b_ = PsiNA().transpose() * D_ * W_ * y_;
        } else {
            XtWX_ = X_.transpose() * W_ * X_;
            invXtWX_ = XtWX_.partialPivLu();
            invXtWXXtW_ = invXtWX_.solve(X_.transpose() * W_);   // (X^\top * W * X)^{-1} * (X^\top * W)
            // woodbury decomposition matrices
            U_ = PsiNA().transpose() * D_ * W_ * X_;
            V_ = X_.transpose() * W_ * PsiNA();
            b_ = PsiNA().transpose() * D_ * internals::lmbQ(W_, X_, invXtWX_, y_);
        }
        W_changed_ = true;
        return;
    }
    // update response and weights
    template <typename WeightMatrix>
    void update_response_and_weights(const vector_t& y, const WeightMatrix& W) {
        fdapde_assert(Psi_.rows() > 0 && y.rows() == n_locs_);
        y_ = y;
        nan_pattern_ = na_matrix(y_);
        if (nan_pattern_.any()) {
            n_obs_ = n_locs_ - nan_pattern_.count();
            B_ = (~nan_pattern_).repeat(1, n_dofs_).select(Psi_, 0);
            y_ = (~nan_pattern_).select(y_, 0);
        }
        update_weights(W);
        return;
    }
    // fit: solve (Psi^T W Psi + lambda R1) f = Psi^T W y
    std::pair<vector_t, vector_t> fit(double lambda) {
        fdapde_assert(lambda > 0 && n_dofs_ > 0 && n_obs_ > 0);
        if ( lambda_saved_.value() != lambda || W_changed_) {
            // assemble spline system: A = Psi^T W Psi + lambda * R1
            sparse_matrix_t A = Psi_.transpose() * W_ * Psi_ + lambda * R1_;
            invA_.compute(A);
            W_changed_ = false;
        }
        if (lambda_saved_.value() != lambda) {
            // update linear system rhs
            b_ = PsiNA().transpose() * (W_ * y_);
        }
        lambda_saved_ = lambda;
        vector_t x;
        if (n_covs_ == 0) {
            f_ = invA_.solve(b_);
        }else {
            // semi-parametric: include covariates using Woodbury
            //matrix_t U = Psi_.transpose() * W_ * X_;
            //matrix_t Cinv = (X_.transpose() * W_ * X_).partialPivLu().inverse(); // small n_covs x n_covs
            //matrix_t V = X_.transpose() * W_ * Psi_;
            //matrix_t y_vec = b_;
            f_ = woodbury_system_solve(invA_, U_, XtWX_, V_, b_);
            beta_ = invXtWXXtW_ * (y_ - Psi_ * f_);
        }
        return std::make_pair(f_, beta_);;
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
            sparse_matrix_t A = Psi_.transpose() * W_ * Psi_ + lambda * R1_;
            invA_.compute(A);
        }
        vector_t x;
        if (n_covs_ == 0) {   // equivalent to calling fit(lambda)
            if (lambda_saved_.value() != lambda) {
                b_ = PsiNA().transpose() * (W_ * y_);
            }
            x = invA_.solve(b_);
        } else {
            vector_t b(n_dofs_);
            // assemble nonparametric linear system rhs
            b = -PsiNA().transpose() * D_ * W_ * y_;
            x = invA_.solve(b);
        }
        lambda_saved_ = lambda;
        f_ = x;
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
            Bs_ = matrix_t::Zero(n_dofs_, r);   // implicitly enforce homogeneous forcing
        }
        if (n_covs_ == 0) {
            Bs_ = PsiNA().transpose() * D_ * W_ * (*Us_);
        } else {
            Bs_ = PsiNA().transpose() * D_ * internals::lmbQ(W_, X_, invXtWX_, *Us_);
        }
        matrix_t x = n_covs_ == 0 ? invA_.solve(*Bs_) : woodbury_system_solve(invA_, U_, XtWX_, V_, *Bs_);
        double trS = 0;   // monte carlo Tr[S] approximation
        for (int i = 0; i < r; ++i) { trS += Ys_->row(i).dot(x.col(i)); }
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
            sparse_matrix_t A = Psi_.transpose() * W_ * Psi_ + lambda_ * R1_;
            invA_.compute(A);
            lambda_saved_ = lambda_;
        }
        return edf(r, seed);
    }
    // penalty matrix: \lambda * R1
    matrix_t P(double lambda) const {
        return lambda * R1_;
    }
    template <typename LambdaT>
        requires(internals::is_vector_like_v<LambdaT>)
    matrix_t P(const LambdaT& lambda) const {
        fdapde_assert(lambda.size() == n_lambda);
        return P(lambda[0]);
    }
    matrix_t P() const { return P(1.0); }
    // efficient evaluation of f^\top * P * f = ...
    double ftPf(double lambda) {
        if (lambda_saved_.value() != lambda || W_changed_) { fit(lambda); }
        return lambda * f_.dot(R1_ * f_);
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
    int n_obs() const { return n_obs_;}
    const binary_t& nan_pattern() const { return nan_pattern_; }
    const sparse_matrix_t mass() const { return R0_; }
    const sparse_matrix_t& stiff() const { return R1_; }
    const sparse_matrix_t& Psi() const { return Psi_; }
    const sparse_matrix_t& PsiNA() const { return B_.has_value() ? *B_ : Psi_; }
    const vector_t& force() const { return u_; }
    const vector_t& f() const { return f_; }
    const vector_t& beta() const { return beta_; }
    const vector_t& misfit() const { return g_; }
    const matrix_t& design_matrix() const { return X_; }
    const vector_t& response() const { return y_; }
    const sparse_matrix_t& weights() const { return W_; }
    double lambda() const { return *lambda_saved_; }

    const matrix_t& U() const { return U_; }
    const matrix_t& V() const { return V_; }

   protected:
    std::optional<double> lambda_saved_ = -1;
    sparse_solver_t invA_;
    sparse_matrix_t A_;
    matrix_t b_;
    // matrices for Hutchinson stochastic estimation of Tr[S]
    std::optional<matrix_t> Ys_, Bs_, Us_;

    int n_dofs_ = 0, n_locs_ = 0, n_obs_ = 0, n_covs_ = 0;
    sparse_matrix_t R0_;    // n_dofs x n_dofs matrix [R0]_{ij} = \int_D \psi_i * \psi_j
    sparse_matrix_t R1_;    // n_dofs x n_dofs matrix [R1]_{ij} = \int_D a(\psi_i, \psi_j)
    vector_t u_;            // n_dofs x 1 vector u_i = \int_D u * \psi_i (not used for spline smoother?)
    sparse_matrix_t Psi_;   // n_obs x n_dofs
    std::optional<sparse_matrix_t> B_; // \Psi matrix corrected for missing observations
    diag_matrix_t D_; // vector of regions' measures (areal sampling)

    vector_t f_, beta_, g_;

    matrix_t X_;               // n_obs x n_covs design matrix
    vector_t y_;               // n_obs x 1 observation vector
    binary_t nan_pattern_;     // n_obs x 1 indicator vector for NaNs
    sparse_matrix_t W_;        // n_obs x n_obs matrix of observation weights
    matrix_t U_, V_;           // (2 * n_dofs) x n_covs matrices [\Psi^\top * D * W * y, 0] and [X^\top * W * \Psi, 0]
    matrix_t XtWX_;            // n_covs x n_covs matrix X^\top * W * X
    dense_solver_t invXtWX_;   // factorization of n_covs x n_covs matrix X^\top * W * X
    matrix_t invXtWXXtW_;      // n_covs x n_obs matrix (X^\top * X)^{-1} * (X^\top W)
    bool W_changed_;

    // basis eval handles
    std::function<sparse_matrix_t(const matrix_t& locs)> point_eval_;
    std::function<std::pair<sparse_matrix_t, vector_t>(const binary_t& locs)> areal_eval_;
};
} // namespace internals

// solver factory
template <typename BilinearForm_, typename LinearForm_> struct ls_graph {
using solver_t = internals::ls_graph;
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
    ls_graph(const BilinearForm_& bilinear_form, const LinearForm_& linear_form) :
        penalty_(bilinear_form, linear_form) { }
    const penalty_packet& get() const { return penalty_; }
private:
    penalty_packet penalty_;
};
} // namespace fdapde

#endif // __FDAPDE_LS_GRAPH_H__

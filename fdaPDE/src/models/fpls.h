#ifndef __FPLS_H__
#define __FPLS_H__

#include "header_check.h"

namespace fdapde {

enum class fPLSMode { Regression, ModeA, SymmetricBlock };
template <fPLSMode Mode> using fPLSModeTag = std::integral_constant<fPLSMode, Mode>;
inline constexpr fPLSModeTag<fPLSMode::Regression> fPLS_R {};
inline constexpr fPLSModeTag<fPLSMode::ModeA> fPLS_A {};
inline constexpr fPLSModeTag<fPLSMode::SymmetricBlock> fPLS_SB {};

template <typename DirectionSolver, typename LoadingSolver, fPLSMode Mode = fPLSMode::Regression> class fPLS {
   private:
    using direction_solver_t = std::decay_t<DirectionSolver>;
    using loading_solver_t = std::decay_t<LoadingSolver>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    static constexpr int direction_n_lambda = direction_solver_t::n_lambda;
    static constexpr int loading_n_lambda = loading_solver_t::n_lambda;

    struct direction_fit_result {
        vector_t f;
        vector_t v;
        std::vector<double> objective_history;
        int iterations = 0;
        bool monotone = true;
    };
   public:
    fPLS() noexcept = default;

    template <typename GeoFrame, typename DirectionPenalty, typename LoadingPenalty>
    fPLS(
      const std::string& colname, const matrix_t& Y, const GeoFrame& gf, DirectionPenalty&& direction_penalty,
      LoadingPenalty&& loading_penalty, fPLSModeTag<Mode> = {}) {
        direction_solver_.discretize(direction_penalty.get());
        loading_solver_.discretize(loading_penalty.get());
        analyze_data(colname, Y, gf);
    }

    template <typename GeoFrame> void analyze_data(const std::string& colname, const matrix_t& Y, const GeoFrame& gf) {
        fdapde_assert(gf.n_layers() == 1);
        data_ = gf[0].data().template col<double>(colname).as_matrix();
        n_locs_ = data_.rows();
        n_units_ = data_.cols();
        fdapde_assert(Y.rows() == n_units_);
        Y_ = Y;
        auto W = vector_t::Ones(n_locs_).asDiagonal();
        direction_solver_.analyze_data(gf, W);
        loading_solver_.analyze_data(gf, W);
    }

    template <typename DirectionLambda, typename LoadingLambda>
        requires(internals::is_subscriptable<DirectionLambda, int> && internals::is_subscriptable<LoadingLambda, int>)
    void fit(
      int n_comp, const DirectionLambda& direction_lambda, const LoadingLambda& loading_lambda, int max_iter = 20,
      double tol = 1e-6) {
        fdapde_assert(direction_lambda.size() == direction_n_lambda);
        fdapde_assert(loading_lambda.size() == loading_n_lambda);
        n_comp_ = n_comp;

        matrix_t X_h = data_.transpose();
        matrix_t Y_h = Y_;
        W_.resize(direction_solver_.n_dofs(), n_comp_);
        C_.resize(loading_solver_.n_dofs(), n_comp_);
        V_.resize(Y_.cols(), n_comp_);
        D_.resize(Y_.cols(), n_comp_);
        T_.resize(n_units_, n_comp_);
        U_.resize(n_units_, n_comp_);
        sigma_.resize(n_comp_);
        direction_objective_history_.resize(n_comp_);
        direction_iterations_.assign(n_comp_, 0);
        direction_monotone_.assign(n_comp_, true);

        matrix_t M_h = Y_h.transpose() * X_h;
        for (int h = 0; h < n_comp_; ++h) {
            fit_direction_(M_h, direction_lambda, max_iter, tol, h);
            project_(X_h, Y_h, h);

            fit_loadings_(X_h, Y_h, loading_lambda, h);

            deflate_(X_h, Y_h, M_h, h);
        }
        if constexpr (Mode == fPLSMode::Regression) { B_ = coefficient_(n_comp_); }
    }

    void fit(
      int n_comp, const std::vector<double>& direction_lambda_grid, const std::vector<double>& loading_lambda_grid,
      int flag, int max_iter = 20, double tol = 1e-6, int edf_r = 100, int seed = random_seed) {
        fdapde_assert(direction_lambda_grid.size() > 0 && direction_lambda_grid.size() % direction_n_lambda == 0);
        fdapde_assert(loading_lambda_grid.size() > 0 && loading_lambda_grid.size() % loading_n_lambda == 0);
        n_comp_ = n_comp;

        matrix_t X_h = data_.transpose();
        matrix_t Y_h = Y_;
        W_.resize(direction_solver_.n_dofs(), n_comp_);
        C_.resize(loading_solver_.n_dofs(), n_comp_);
        V_.resize(Y_.cols(), n_comp_);
        D_.resize(Y_.cols(), n_comp_);
        T_.resize(n_units_, n_comp_);
        U_.resize(n_units_, n_comp_);
        sigma_.resize(n_comp_);
        direction_objective_history_.resize(n_comp_);
        direction_iterations_.assign(n_comp_, 0);
        direction_monotone_.assign(n_comp_, true);
        direction_lambda_.resize(n_comp_, direction_n_lambda);
        loading_lambda_.resize(n_comp_, loading_n_lambda);

        int calibration = (flag & 0b11110);
        matrix_t M_h = Y_h.transpose() * X_h;
        for (int h = 0; h < n_comp_; ++h) {
            Eigen::Matrix<double, direction_n_lambda, 1> direction_lambda;
            Eigen::Matrix<double, loading_n_lambda, 1> loading_lambda;
            switch (calibration) {
            case 0: {
                fdapde_assert(direction_lambda_grid.size() == direction_n_lambda);
                fdapde_assert(loading_lambda_grid.size() == loading_n_lambda);
                std::copy(direction_lambda_grid.begin(), direction_lambda_grid.end(), direction_lambda.begin());
                std::copy(loading_lambda_grid.begin(), loading_lambda_grid.end(), loading_lambda.begin());
            } break;
            case OptimizeGCV: {
                Eigen::JacobiSVD<matrix_t> svd(M_h, Eigen::ComputeThinU | Eigen::ComputeThinV);
                vector_t f0 = svd.matrixV().col(0);
                auto direction_gcv = [&](auto lambda) {
                    return direction_gcv_(M_h, lambda, f0, max_iter, tol, edf_r, seed);
                };
                auto loading_gcv = [&](auto lambda) {
                    return loading_gcv_(X_h, T_.col(h), lambda, edf_r, seed);
                };
                GridSearch<direction_n_lambda> direction_optimizer;
                direction_lambda = direction_optimizer.optimize(direction_gcv, direction_lambda_grid);

                fit_direction_(M_h, direction_lambda, f0, max_iter, tol, h);
                project_(X_h, Y_h, h);

                if constexpr (Mode == fPLSMode::SymmetricBlock) {
                    std::copy(loading_lambda_grid.begin(), loading_lambda_grid.begin() + loading_n_lambda, loading_lambda.begin());
                }
                if constexpr (Mode == fPLSMode::Regression || Mode == fPLSMode::ModeA) {
                    GridSearch<loading_n_lambda> loading_optimizer;
                    loading_lambda = loading_optimizer.optimize(loading_gcv, loading_lambda_grid);
                }
            } break;
            default: {
                throw std::runtime_error("Unrecognized calibration option.");
            }
            }

            if (calibration != OptimizeGCV) {
                Eigen::JacobiSVD<matrix_t> svd(M_h, Eigen::ComputeThinU | Eigen::ComputeThinV);
                fit_direction_(M_h, direction_lambda, svd.matrixV().col(0), max_iter, tol, h);
                project_(X_h, Y_h, h);
            }

            fit_loadings_(X_h, Y_h, loading_lambda, h);
            for (int j = 0; j < direction_n_lambda; ++j) { direction_lambda_(h, j) = direction_lambda[j]; }
            for (int j = 0; j < loading_n_lambda; ++j) { loading_lambda_(h, j) = loading_lambda[j]; }

            deflate_(X_h, Y_h, M_h, h);
        }
        if constexpr (Mode == fPLSMode::Regression) { B_ = coefficient_(n_comp_); }
    }

    static constexpr fPLSMode mode() { return Mode; }
    const matrix_t& X_space_directions() const { return W_; }
    const matrix_t& Y_space_directions() const { return V_; }
    const matrix_t& X_latent() const { return T_; }
    const matrix_t& X_latent_scores() const { return T_; }
    const matrix_t& Y_latent_scores() const { return U_; }
    const matrix_t& X_loadings() const { return C_; }
    const matrix_t& Y_loadings() const { return D_; }
    matrix_t fitted() const { return fitted(n_comp_); }
    matrix_t fitted(int h) const {
        h = components_(h);
        if constexpr (Mode == fPLSMode::Regression) { return T_.leftCols(h) * D_.leftCols(h).transpose(); }
        if constexpr (Mode == fPLSMode::ModeA || Mode == fPLSMode::SymmetricBlock) {
            return U_.leftCols(h) * D_.leftCols(h).transpose();
        }
    }
    matrix_t reconstructed() const { return reconstructed(n_comp_); }
    matrix_t reconstructed(int h) const {
        h = components_(h);
        return T_.leftCols(h) * (loading_solver_.Psi() * C_.leftCols(h)).transpose();
    }
    const matrix_t& B() const requires(Mode == fPLSMode::Regression) { return B_; }
    matrix_t B(int h) const requires(Mode == fPLSMode::Regression) {
        h = components_(h);
        if (h == n_comp_) return B_;
        return coefficient_(h);
    }
    const matrix_t& Beta() const requires(Mode == fPLSMode::Regression) { return B(); }
    matrix_t Beta(int h) const requires(Mode == fPLSMode::Regression) { return B(h); }
    const matrix_t& direction_lambda() const { return direction_lambda_; }
    const matrix_t& loading_lambda() const { return loading_lambda_; }
    const std::vector<std::vector<double>>& direction_objective_history() const {
        return direction_objective_history_;
    }
    const std::vector<int>& direction_iterations() const { return direction_iterations_; }
    const std::vector<bool>& direction_monotone() const { return direction_monotone_; }

   private:
    int components_(int h) const {
        if (h == 0) h = n_comp_;
        fdapde_assert(h > 0 && h <= n_comp_);
        return h;
    }
    matrix_t coefficient_(int h) const {
        static_assert(Mode == fPLSMode::Regression);
        const auto W_h = W_.leftCols(h);
        const auto C_h = C_.leftCols(h);
        const auto D_h = D_.leftCols(h);
        return W_h * (C_h.transpose() * loading_solver_.Psi().transpose() * loading_solver_.Psi() * W_h)
                       .partialPivLu()
                       .solve(D_h.transpose());
    }
    template <typename Lambda>
        requires(internals::is_subscriptable<Lambda, int>)
    void fit_loadings_(const matrix_t& X_h, const matrix_t& Y_h, const Lambda& loading_lambda, int h) {
        if constexpr (Mode == fPLSMode::SymmetricBlock) {
            C_.col(h) = W_.col(h);
            D_.col(h) = V_.col(h);
            return;
        }
        const double t_norm = T_.col(h).squaredNorm();
        if (!std::isfinite(t_norm) || t_norm <= 0) {
            throw std::runtime_error("fPLS loading update has a non-finite or zero X score");
        }
        const vector_t x_response = X_h.transpose() * T_.col(h) / t_norm;
        if (!x_response.array().isFinite().all()) {
            throw std::runtime_error("fPLS loading update produced a non-finite response");
        }
        loading_solver_.update_response(x_response);
        loading_solver_.fit(loading_lambda);
        if (!loading_solver_.f().array().isFinite().all()) {
            throw std::runtime_error("fPLS loading solver produced a non-finite solution");
        }
        C_.col(h) = loading_solver_.f();
        if constexpr (Mode == fPLSMode::Regression) { D_.col(h) = Y_h.transpose() * T_.col(h) / t_norm; }
        if constexpr (Mode == fPLSMode::ModeA) {
            const double u_norm = U_.col(h).squaredNorm();
            if (!std::isfinite(u_norm) || u_norm <= 0) {
                throw std::runtime_error("fPLS loading update has a non-finite or zero Y score");
            }
            D_.col(h) = Y_h.transpose() * U_.col(h) / u_norm;
        }
        if (!D_.col(h).array().isFinite().all()) { throw std::runtime_error("fPLS response loading is not finite"); }
    }
    void project_(const matrix_t& X_h, const matrix_t& Y_h, int h) {
        T_.col(h) = X_h * direction_solver_.Psi() * W_.col(h);
        U_.col(h) = Y_h * V_.col(h);
    }
    void deflate_(matrix_t& X_h, matrix_t& Y_h, matrix_t& M_h, int h) {
        if constexpr (Mode == fPLSMode::SymmetricBlock) {
            M_h -= sigma_[h] * V_.col(h) * W_.col(h).transpose();
            return;
        }
        X_h -= T_.col(h) * (loading_solver_.Psi() * C_.col(h)).transpose();
        if constexpr (Mode == fPLSMode::Regression) {
            Y_h -= T_.col(h) * D_.col(h).transpose();
        }
        if constexpr (Mode == fPLSMode::ModeA) {
            Y_h -= U_.col(h) * D_.col(h).transpose();
        }
        M_h = Y_h.transpose() * X_h;
    }
    template <typename Lambda>
        requires(internals::is_subscriptable<Lambda, int>)
    void fit_direction_(const matrix_t& M, const Lambda& lambda, int max_iter, double tol, int h) {
        Eigen::JacobiSVD<matrix_t> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
        fit_direction_(M, lambda, svd.matrixV().col(0), max_iter, tol, h);
    }
    template <typename Lambda, typename Init>
        requires(internals::is_subscriptable<Lambda, int>)
    auto solve_direction_(const matrix_t& M, const Lambda& lambda, const Init& f0, int max_iter, double tol) {
        vector_t fn = f0;
        vector_t v(M.rows());
        double Jold = std::numeric_limits<double>::max(), Jnew = 1.0;
        direction_fit_result result;
        result.objective_history.reserve(max_iter);

        for (int i = 0; !almost_equal(Jnew, Jold, tol) && i < max_iter; ++i) {
            v = M * fn;
            const double v_norm = v.norm();
            if (!v.array().isFinite().all() || !std::isfinite(v_norm) || v_norm <= 0) {
                throw std::runtime_error("fPLS direction update is non-finite or numerically singular");
            }
            v /= v_norm;
            const vector_t response = M.transpose() * v;
            if (!response.array().isFinite().all()) {
                throw std::runtime_error("fPLS direction update produced a non-finite response");
            }
            direction_solver_.update_response(response);
            direction_solver_.fit(lambda);
            Jold = Jnew;
            fn = direction_solver_.fn();
            if (!fn.array().isFinite().all()) {
                throw std::runtime_error("fPLS direction solver produced a non-finite solution");
            }
            Jnew = (M - v * fn.transpose()).squaredNorm() + direction_solver_.ftPf(lambda);
            if (!std::isfinite(Jnew)) { throw std::runtime_error("fPLS direction objective is not finite"); }
            result.objective_history.push_back(Jnew);
            result.iterations = i + 1;
            if (i > 0 && (Jnew - Jold) / (1.0 + std::abs(Jold)) > tol) { result.monotone = false; }
        }

        result.f = direction_solver_.f();
        result.v = std::move(v);
        return result;
    }
    template <typename Lambda, typename Init>
        requires(internals::is_subscriptable<Lambda, int>)
    void fit_direction_(const matrix_t& M, const Lambda& lambda, const Init& f0, int max_iter, double tol, int h) {
        auto result = solve_direction_(M, lambda, f0, max_iter, tol);
        const double w_norm = (direction_solver_.Psi() * result.f).norm();
        const double v_norm = result.v.norm();
        if (!std::isfinite(w_norm) || !std::isfinite(v_norm) || w_norm <= 0 || v_norm <= 0) {
            throw std::runtime_error("fPLS direction normalization is non-finite or numerically singular");
        }
        W_.col(h) = result.f / w_norm;
        V_.col(h) = result.v / v_norm;
        sigma_[h] = w_norm * v_norm;
        direction_objective_history_[h] = std::move(result.objective_history);
        direction_iterations_[h] = result.iterations;
        direction_monotone_[h] = result.monotone;
    }
    template <typename Lambda, typename Init>
        requires(internals::is_subscriptable<Lambda, int>)
    double direction_gcv_(
      const matrix_t& M, const Lambda& lambda, const Init& f0, int max_iter, double tol, int edf_r, int seed) {
        const auto result = solve_direction_(M, lambda, f0, max_iter, tol);
        double dor = n_locs_ - direction_solver_.edf(lambda, edf_r, seed);
        return (n_locs_ / std::pow(dor, 2)) *
               ((direction_solver_.Psi() * result.f) - direction_solver_.response()).squaredNorm();
    }
    template <typename Lambda>
        requires(internals::is_subscriptable<Lambda, int>)
    double loading_gcv_(const matrix_t& X, const vector_t& t, const Lambda& lambda, int edf_r, int seed) {
        loading_solver_.update_response(X.transpose() * t / t.squaredNorm());
        loading_solver_.fit(lambda);
        double dor = n_locs_ - loading_solver_.edf(lambda, edf_r, seed);
        return (n_locs_ / std::pow(dor, 2)) *
               (loading_solver_.fn() - loading_solver_.response()).squaredNorm();
    }

    matrix_t data_, Y_;
    direction_solver_t direction_solver_;
    loading_solver_t loading_solver_;
    int n_locs_ = 0, n_units_ = 0, n_comp_ = 0;

    matrix_t W_;   // X directions
    matrix_t V_;   // Y directions
    matrix_t T_;   // X scores
    matrix_t U_;   // Y scores
    matrix_t C_;   // X loadings
    matrix_t D_;   // Y loadings
    matrix_t B_;   // regression operator
    vector_t sigma_;
    matrix_t direction_lambda_;
    matrix_t loading_lambda_;
    std::vector<std::vector<double>> direction_objective_history_;
    std::vector<int> direction_iterations_;
    std::vector<bool> direction_monotone_;
};

template <typename GeoFrame, typename DirectionPenalty, typename LoadingPenalty>
fPLS(
  const std::string& colname, const Eigen::Matrix<double, Dynamic, Dynamic>& Y, const GeoFrame& gf,
  DirectionPenalty&& direction_penalty, LoadingPenalty&& loading_penalty)
  -> fPLS<typename DirectionPenalty::solver_t, typename LoadingPenalty::solver_t, fPLSMode::Regression>;

template <typename GeoFrame, typename DirectionPenalty, typename LoadingPenalty, fPLSMode Mode>
fPLS(
  const std::string& colname, const Eigen::Matrix<double, Dynamic, Dynamic>& Y, const GeoFrame& gf,
  DirectionPenalty&& direction_penalty, LoadingPenalty&& loading_penalty, fPLSModeTag<Mode>)
  -> fPLS<typename DirectionPenalty::solver_t, typename LoadingPenalty::solver_t, Mode>;

}   // namespace fdapde

#endif   // __FPLS_H__

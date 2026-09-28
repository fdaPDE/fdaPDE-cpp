#ifndef __FPLS_H__
#define __FPLS_H__

#include "header_check.h"

namespace fdapde {

/// @brief selects regression, separate block deflation, or cross-covariance deflation
enum class fPLSMode {
    Regression,
    ModeA,
    SymmetricBlock
};
template <fPLSMode Mode> using fPLSModeTag = std::integral_constant<fPLSMode, Mode>;
inline constexpr fPLSModeTag<fPLSMode::Regression> fPLS_R {};
inline constexpr fPLSModeTag<fPLSMode::ModeA> fPLS_A {};
inline constexpr fPLSModeTag<fPLSMode::SymmetricBlock> fPLS_SB {};

/// @brief extracts smooth functional partial least squares components from centered predictor and response data
template <typename DirectionSolver, typename LoadingSolver, fPLSMode Mode = fPLSMode::Regression> class fPLS {
   private:
    using direction_solver_t = std::decay_t<DirectionSolver>;
    using loading_solver_t = std::decay_t<LoadingSolver>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    static constexpr int direction_n_lambda = direction_solver_t::n_lambda;
    static constexpr int loading_n_lambda = loading_solver_t::n_lambda;

    /// @brief stores a direction solution and its objective convergence diagnostics
    struct direction_fit_result {
        vector_t f;
        vector_t v;
        std::vector<double> objective_history;
        int iterations = 0;
        bool monotone = true;
    };
   public:
    /// @brief creates an uninitialized model without discretized penalties or data
    fPLS() noexcept = default;

    /// @brief discretizes direction and loading penalties and binds the centered data
    template <typename GeoFrame, typename DirectionPenalty, typename LoadingPenalty>
    fPLS(
      const std::string& colname,             // name of the centered functional predictor column in the sole data layer
      const matrix_t& Y,                      // centered responses with one row per statistical unit
      const GeoFrame& gf,                     // single-layer predictor data: locations by statistical units
      DirectionPenalty&& direction_penalty,   // variational penalty packet for smoothing predictor directions
      LoadingPenalty&& loading_penalty,       // variational penalty packet for smoothing predictor loadings
      fPLSModeTag<Mode> = {}                  // compile-time choice of the component deflation mode
    ) {
        direction_solver_.discretize(direction_penalty.get());
        loading_solver_.discretize(loading_penalty.get());
        analyze_data(colname, Y, gf);
    }

    /// @brief binds centered predictor and response data to both discretized smoothers
    template <typename GeoFrame>
    void analyze_data(
      const std::string& colname,   // name of the centered functional predictor column in the sole data layer
      const matrix_t& Y,            // centered responses with one row per statistical unit
      const GeoFrame& gf            // single-layer predictor data: locations by statistical units
    ) {
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

    /// @brief extracts components using fixed direction and loading penalties
    template <typename DirectionLambda, typename LoadingLambda>
        requires(internals::is_subscriptable<DirectionLambda, int> && internals::is_subscriptable<LoadingLambda, int>)
    void fit(
      int n_comp,                                // number of components to extract
      const DirectionLambda& direction_lambda,   // fixed direction smoothing parameters shared by all components
      const LoadingLambda& loading_lambda,       // fixed loading smoothing parameters shared by all components
      int max_iter = 1000,                       // maximum number of alternating direction updates per component
      double tol = 1e-8                          // objective tolerance for stopping and relative increase diagnostics
    ) {
        fdapde_assert(direction_lambda.size() == direction_n_lambda);
        fdapde_assert(loading_lambda.size() == loading_n_lambda);
        n_comp_ = n_comp;

        // start from centered blocks and allocate component coefficients and diagnostics
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

        // extract each component from the current cross-product before mode-specific deflation
        matrix_t M_h = Y_h.transpose() * X_h;
        for (int h = 0; h < n_comp_; ++h) {
            fit_direction_(M_h, direction_lambda, max_iter, tol, h);
            project_(X_h, Y_h, h);

            fit_loadings_(X_h, Y_h, loading_lambda, h);

            deflate_(X_h, Y_h, M_h, h);
        }
        if constexpr (Mode == fPLSMode::Regression) { B_ = coefficient_(n_comp_); }
    }

    /// @brief extracts components using fixed penalty schedules or componentwise GCV selection
    /// fixed schedules contain one shared tuple or one tuple per component; GCV grids contain candidate tuples
    void fit(
      int n_comp,                                         // number of components to extract
      const std::vector<double>& direction_lambda_grid,   // candidate tuples or fixed penalty schedule
      const std::vector<double>& loading_lambda_grid,     // candidate tuples or fixed penalty schedule
      int flag,                // calibration bits: zero for fixed penalties or OptimizeGCV for grid search
      int max_iter = 1000,     // maximum number of alternating direction updates per component
      double tol = 1e-8,       // objective tolerance for stopping and relative increase diagnostics
      int edf_r = 100,         // number of random probes for the effective degrees of freedom estimate
      int seed = random_seed   // seed for the effective degrees of freedom estimate
    ) {
        fdapde_assert(direction_lambda_grid.size() > 0 && direction_lambda_grid.size() % direction_n_lambda == 0);
        fdapde_assert(loading_lambda_grid.size() > 0 && loading_lambda_grid.size() % loading_n_lambda == 0);
        n_comp_ = n_comp;

        // start from centered blocks and allocate component coefficients and diagnostics
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
        direction_gcv_values_.assign(n_comp_, {});
        loading_gcv_values_.assign(n_comp_, {});
        direction_lambda_.resize(n_comp_, direction_n_lambda);
        loading_lambda_.resize(n_comp_, loading_n_lambda);

        int calibration = (flag & 0b11110);
        // extract each component from the current cross-product before mode-specific deflation
        matrix_t M_h = Y_h.transpose() * X_h;
        for (int h = 0; h < n_comp_; ++h) {
            Eigen::Matrix<double, direction_n_lambda, 1> direction_lambda;
            Eigen::Matrix<double, loading_n_lambda, 1> loading_lambda;
            switch (calibration) {
            case 0: {
                fdapde_assert(
                  direction_lambda_grid.size() == direction_n_lambda ||
                  direction_lambda_grid.size() == n_comp_ * direction_n_lambda);
                fdapde_assert(
                  loading_lambda_grid.size() == loading_n_lambda ||
                  loading_lambda_grid.size() == n_comp_ * loading_n_lambda);
                const auto direction_begin =
                  direction_lambda_grid.begin() +
                  (direction_lambda_grid.size() == direction_n_lambda ? 0 : h * direction_n_lambda);
                const auto loading_begin = loading_lambda_grid.begin() +
                                           (loading_lambda_grid.size() == loading_n_lambda ? 0 : h * loading_n_lambda);
                std::copy(direction_begin, direction_begin + direction_n_lambda, direction_lambda.begin());
                std::copy(loading_begin, loading_begin + loading_n_lambda, loading_lambda.begin());
            } break;
            case OptimizeGCV: {
                // compare direction penalties from a common initialization before selecting the loading penalty
                Eigen::JacobiSVD<matrix_t> svd(M_h, Eigen::ComputeThinU | Eigen::ComputeThinV);
                vector_t f0 = svd.matrixV().col(0);
                auto direction_gcv = [&](auto lambda) {
                    return direction_gcv_(M_h, lambda, f0, max_iter, tol, edf_r, seed);
                };
                auto loading_gcv = [&](auto lambda) { return loading_gcv_(X_h, T_.col(h), lambda, edf_r, seed); };
                GridSearch<direction_n_lambda> direction_optimizer;
                direction_lambda = direction_optimizer.optimize(direction_gcv, direction_lambda_grid);
                direction_gcv_values_[h] = direction_optimizer.values();

                fit_direction_(M_h, direction_lambda, f0, max_iter, tol, h);
                project_(X_h, Y_h, h);

                if constexpr (Mode == fPLSMode::SymmetricBlock) {
                    std::copy(
                      loading_lambda_grid.begin(), loading_lambda_grid.begin() + loading_n_lambda,
                      loading_lambda.begin());
                }
                if constexpr (Mode == fPLSMode::Regression || Mode == fPLSMode::ModeA) {
                    GridSearch<loading_n_lambda> loading_optimizer;
                    loading_lambda = loading_optimizer.optimize(loading_gcv, loading_lambda_grid);
                    loading_gcv_values_[h] = loading_optimizer.values();
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

    /// @brief returns the compile-time component deflation mode
    static constexpr fPLSMode mode() { return Mode; }
    /// @brief returns predictor direction coefficients with one column per component
    const matrix_t& X_space_directions() const { return W_; }
    /// @brief returns response directions with one column per component
    const matrix_t& Y_space_directions() const { return V_; }
    /// @brief returns predictor scores with statistical units in rows and components in columns
    const matrix_t& X_latent() const { return T_; }
    /// @brief returns predictor scores through the explicit score accessor
    const matrix_t& X_latent_scores() const { return T_; }
    /// @brief returns response scores with statistical units in rows and components in columns
    const matrix_t& Y_latent_scores() const { return U_; }
    /// @brief returns predictor loading coefficients with one column per component
    const matrix_t& X_loadings() const { return C_; }
    /// @brief returns response loadings with one column per component
    const matrix_t& Y_loadings() const { return D_; }
    /// @brief returns per-component direction GCV curves from the last grid or schedule fit; empty without GCV
    const std::vector<std::vector<double>>& direction_gcv_values() const { return direction_gcv_values_; }
    /// @brief returns per-component loading GCV curves; empty without GCV or in symmetric mode
    const std::vector<std::vector<double>>& loading_gcv_values() const { return loading_gcv_values_; }
    /// @brief reconstructs centered responses using all fitted components
    matrix_t fitted() const { return fitted(n_comp_); }
    /// @brief reconstructs responses from a prefix of the fitted components
    matrix_t fitted(
      int h   // component count; zero selects the full fit
    ) const {
        h = components_(h);
        if constexpr (Mode == fPLSMode::Regression) { return T_.leftCols(h) * D_.leftCols(h).transpose(); }
        if constexpr (Mode == fPLSMode::ModeA || Mode == fPLSMode::SymmetricBlock) {
            return U_.leftCols(h) * D_.leftCols(h).transpose();
        }
    }
    /// @brief reconstructs centered predictors at observation locations using all fitted components
    matrix_t reconstructed() const { return reconstructed(n_comp_); }
    /// @brief reconstructs predictors at observation locations from a component prefix
    matrix_t reconstructed(
      int h   // component count; zero selects the full fit
    ) const {
        h = components_(h);
        return T_.leftCols(h) * (loading_solver_.Psi() * C_.leftCols(h)).transpose();
    }
    /// @brief returns the cached regression operator mapping predictors in the direction basis to responses
    const matrix_t& B() const
        requires(Mode == fPLSMode::Regression)
    {
        return B_;
    }
    /// @brief returns the regression operator for a prefix of the fitted components
    // clang-format off
    matrix_t B(
      int h   // component count; zero selects the full fit
    ) const
        requires(Mode == fPLSMode::Regression)
    {
        // clang-format on
        h = components_(h);
        if (h == n_comp_) return B_;
        return coefficient_(h);
    }
    /// @brief returns the cached full regression operator through the B accessor
    const matrix_t& Beta() const
        requires(Mode == fPLSMode::Regression)
    {
        return B();
    }
    /// @brief returns the component-prefix regression operator through the B accessor
    matrix_t Beta(
      int h   // component count; zero selects the full fit
    ) const
        requires(Mode == fPLSMode::Regression)
    {
        return B(h);
    }
    /// @brief returns direction penalties recorded by the grid or schedule overload, one row per component
    const matrix_t& direction_lambda() const { return direction_lambda_; }
    /// @brief returns loading penalties recorded by the grid or schedule overload, one row per component
    const matrix_t& loading_lambda() const { return loading_lambda_; }
    /// @brief returns penalized direction objectives after each update, grouped by component
    const std::vector<std::vector<double>>& direction_objective_history() const { return direction_objective_history_; }
    /// @brief returns the number of completed direction updates for each component
    const std::vector<int>& direction_iterations() const { return direction_iterations_; }
    /// @brief reports whether each direction objective avoided increases beyond the relative tolerance
    const std::vector<bool>& direction_monotone() const { return direction_monotone_; }
   private:
    /// @brief resolves zero to the full component count and checks the requested prefix
    int components_(
      int h   // component count; zero selects the full fit
    ) const {
        if (h == 0) h = n_comp_;
        fdapde_assert(h > 0 && h <= n_comp_);
        return h;
    }
    /// @brief converts sequential score regression to coefficients on the original predictors
    matrix_t coefficient_(
      int h   // number of leading fitted components
    ) const {
        static_assert(Mode == fPLSMode::Regression);
        const auto W_h = W_.leftCols(h);
        const auto C_h = C_.leftCols(h);
        const auto D_h = D_.leftCols(h);
        // sequential deflation couples each score only to preceding scores, with unit diagonal
        const matrix_t coupling = C_h.transpose() * loading_solver_.Psi().transpose() * direction_solver_.Psi() * W_h;
        return W_h * coupling.template triangularView<Eigen::UnitUpper>().solve(D_h.transpose());
    }
    /// @brief stores direction loadings in symmetric mode or fits mode-specific block loadings
    template <typename Lambda>
        requires(internals::is_subscriptable<Lambda, int>)
    void fit_loadings_(
      const matrix_t& X_h,            // current predictor residuals with statistical units in rows
      const matrix_t& Y_h,            // current response residuals with statistical units in rows
      const Lambda& loading_lambda,   // loading smoothing parameters for the current component
      int h                           // zero-based component index
    ) {
        if constexpr (Mode == fPLSMode::SymmetricBlock) {
            C_.col(h) = W_.col(h);
            D_.col(h) = V_.col(h);
            return;
        }
        // smooth predictor loadings on predictor scores; the mode selects the response score block
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
    /// @brief projects the current predictor and response blocks onto their component directions
    void project_(
      const matrix_t& X_h,   // current predictor residuals with statistical units in rows
      const matrix_t& Y_h,   // current response residuals with statistical units in rows
      int h                  // zero-based component index
    ) {
        T_.col(h) = X_h * direction_solver_.Psi() * W_.col(h);
        U_.col(h) = Y_h * V_.col(h);
    }
    /// @brief removes the component from the blocks or directly from the cross-product in symmetric mode
    void deflate_(
      matrix_t& X_h,   // current predictor residuals with statistical units in rows
      matrix_t& Y_h,   // current response residuals with statistical units in rows
      matrix_t& M_h,   // current response-predictor cross-product, updated in place
      int h            // zero-based component index
    ) {
        if constexpr (Mode == fPLSMode::SymmetricBlock) {
            M_h -= sigma_[h] * V_.col(h) * W_.col(h).transpose();
            return;
        }
        X_h -= T_.col(h) * (loading_solver_.Psi() * C_.col(h)).transpose();
        if constexpr (Mode == fPLSMode::Regression) { Y_h -= T_.col(h) * D_.col(h).transpose(); }
        if constexpr (Mode == fPLSMode::ModeA) { Y_h -= U_.col(h) * D_.col(h).transpose(); }
        M_h = Y_h.transpose() * X_h;
    }
    /// @brief initializes a direction from the leading right singular vector and fits the component
    template <typename Lambda>
        requires(internals::is_subscriptable<Lambda, int>)
    void fit_direction_(
      const matrix_t& M,      // response-predictor cross-product used for direction estimation
      const Lambda& lambda,   // smoothing parameters for the current solver
      int max_iter,           // maximum number of alternating direction updates per component
      double tol,             // objective tolerance for stopping and relative increase diagnostics
      int h                   // zero-based component index
    ) {
        Eigen::JacobiSVD<matrix_t> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
        fit_direction_(M, lambda, svd.matrixV().col(0), max_iter, tol, h);
    }
    /// @brief alternates response direction and predictor smoothing until convergence or the iteration limit
    template <typename Lambda, typename Init>
        requires(internals::is_subscriptable<Lambda, int>)
    auto solve_direction_(
      const matrix_t& M,      // response-predictor cross-product used for direction estimation
      const Lambda& lambda,   // smoothing parameters for the current solver
      const Init& f0,         // initial predictor direction evaluated at the observation locations
      int max_iter,           // maximum number of alternating direction updates per component
      double tol              // objective tolerance for stopping and relative increase diagnostics
    ) {
        vector_t fn = f0;
        vector_t v(M.rows());
        double Jold = std::numeric_limits<double>::max(), Jnew = 1.0;
        direction_fit_result result;
        result.objective_history.reserve(max_iter);

        // alternate a unit response direction with a penalized predictor fit and track objective changes
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
    /// @brief fits and normalizes a direction pair and records its component diagnostics
    template <typename Lambda, typename Init>
        requires(internals::is_subscriptable<Lambda, int>)
    void fit_direction_(
      const matrix_t& M,      // response-predictor cross-product used for direction estimation
      const Lambda& lambda,   // smoothing parameters for the current solver
      const Init& f0,         // initial predictor direction evaluated at the observation locations
      int max_iter,           // maximum number of alternating direction updates per component
      double tol,             // objective tolerance for stopping and relative increase diagnostics
      int h                   // zero-based component index
    ) {
        auto result = solve_direction_(M, lambda, f0, max_iter, tol);
        const double w_norm = (direction_solver_.Psi() * result.f).norm();
        const double v_norm = result.v.norm();
        if (!std::isfinite(w_norm) || !std::isfinite(v_norm) || w_norm <= 0 || v_norm <= 0) {
            throw std::runtime_error("fPLS direction normalization is non-finite or numerically singular");
        }
        // normalize in observation space while retaining coefficient-space predictor directions
        W_.col(h) = result.f / w_norm;
        V_.col(h) = result.v / v_norm;
        sigma_[h] = w_norm * v_norm;
        direction_objective_history_[h] = std::move(result.objective_history);
        direction_iterations_[h] = result.iterations;
        direction_monotone_[h] = result.monotone;
    }
    /// @brief evaluates direction GCV after fitting from the supplied initial direction
    template <typename Lambda, typename Init>
        requires(internals::is_subscriptable<Lambda, int>)
    double direction_gcv_(
      const matrix_t& M,      // response-predictor cross-product used for direction estimation
      const Lambda& lambda,   // smoothing parameters for the current solver
      const Init& f0,         // initial predictor direction evaluated at the observation locations
      int max_iter,           // maximum number of alternating direction updates per component
      double tol,             // objective tolerance for stopping and relative increase diagnostics
      int edf_r,              // number of random probes for the effective degrees of freedom estimate
      int seed                // seed for the effective degrees of freedom estimate
    ) {
        const auto result = solve_direction_(M, lambda, f0, max_iter, tol);
        double dor = n_locs_ - direction_solver_.edf(lambda, edf_r, seed);
        return (n_locs_ / std::pow(dor, 2)) *
               ((direction_solver_.Psi() * result.f) - direction_solver_.response()).squaredNorm();
    }
    /// @brief evaluates loading GCV for the predictor residuals and current scores
    template <typename Lambda>
        requires(internals::is_subscriptable<Lambda, int>)
    double loading_gcv_(
      const matrix_t& X,      // predictor residuals with statistical units in rows
      const vector_t& t,      // predictor scores for the current component
      const Lambda& lambda,   // smoothing parameters for the current solver
      int edf_r,              // number of random probes for the effective degrees of freedom estimate
      int seed                // seed for the effective degrees of freedom estimate
    ) {
        loading_solver_.update_response(X.transpose() * t / t.squaredNorm());
        loading_solver_.fit(lambda);
        double dor = n_locs_ - loading_solver_.edf(lambda, edf_r, seed);
        return (n_locs_ / std::pow(dor, 2)) * (loading_solver_.fn() - loading_solver_.response()).squaredNorm();
    }

    matrix_t data_, Y_;
    direction_solver_t direction_solver_;
    loading_solver_t loading_solver_;
    int n_locs_ = 0, n_units_ = 0, n_comp_ = 0;

    matrix_t W_;   // predictor direction coefficients
    matrix_t V_;   // response directions
    matrix_t T_;   // predictor scores
    matrix_t U_;   // response scores
    matrix_t C_;   // predictor loading coefficients
    matrix_t D_;   // response loadings
    matrix_t B_;   // regression operator
    vector_t sigma_;
    matrix_t direction_lambda_;
    matrix_t loading_lambda_;
    std::vector<std::vector<double>> direction_objective_history_;
    std::vector<std::vector<double>> direction_gcv_values_;
    std::vector<std::vector<double>> loading_gcv_values_;
    std::vector<int> direction_iterations_;
    std::vector<bool> direction_monotone_;
};

/// @brief deduces solver types with regression as the default deflation mode
template <typename GeoFrame, typename DirectionPenalty, typename LoadingPenalty>
fPLS(
  const std::string& colname,   // name of the centered functional predictor column in the sole data layer
  const Eigen::Matrix<double, Dynamic, Dynamic>& Y,   // centered responses with one row per statistical unit
  const GeoFrame& gf,                                 // single-layer predictor data: locations by statistical units
  DirectionPenalty&& direction_penalty,               // variational penalty packet for smoothing predictor directions
  LoadingPenalty&& loading_penalty                    // variational penalty packet for smoothing predictor loadings
  ) -> fPLS<typename DirectionPenalty::solver_t, typename LoadingPenalty::solver_t, fPLSMode::Regression>;

/// @brief deduces solver types and the explicit deflation mode
template <typename GeoFrame, typename DirectionPenalty, typename LoadingPenalty, fPLSMode Mode>
fPLS(
  const std::string& colname,   // name of the centered functional predictor column in the sole data layer
  const Eigen::Matrix<double, Dynamic, Dynamic>& Y,   // centered responses with one row per statistical unit
  const GeoFrame& gf,                                 // single-layer predictor data: locations by statistical units
  DirectionPenalty&& direction_penalty,               // variational penalty packet for smoothing predictor directions
  LoadingPenalty&& loading_penalty,                   // variational penalty packet for smoothing predictor loadings
  fPLSModeTag<Mode>                                   // compile-time choice of the component deflation mode
  ) -> fPLS<typename DirectionPenalty::solver_t, typename LoadingPenalty::solver_t, Mode>;

}   // namespace fdapde

#endif   // __FPLS_H__

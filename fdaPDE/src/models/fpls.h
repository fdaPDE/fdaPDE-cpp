#ifndef __FPLS_H__
#define __FPLS_H__

#include "header_check.h"

namespace fdapde {

template <typename DirectionSolver, typename LoadingSolver> class fPLS {
   private:
    using direction_solver_t = std::decay_t<DirectionSolver>;
    using loading_solver_t = std::decay_t<LoadingSolver>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    static constexpr int direction_n_lambda = direction_solver_t::n_lambda;
    static constexpr int loading_n_lambda = loading_solver_t::n_lambda;
   public:
    fPLS() noexcept = default;

    template <typename GeoFrame, typename DirectionPenalty, typename LoadingPenalty>
    fPLS(
      const std::string& colname, const matrix_t& Y, const GeoFrame& gf, DirectionPenalty&& direction_penalty,
      LoadingPenalty&& loading_penalty) {
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

        for (int h = 0; h < n_comp_; ++h) {
            fit_direction_(Y_h.transpose() * X_h, direction_lambda, max_iter, tol, h);
            T_.col(h) = X_h * direction_solver_.Psi() * W_.col(h);

            loading_solver_.update_response(X_h.transpose() * T_.col(h) / T_.col(h).squaredNorm());
            loading_solver_.fit(loading_lambda);
            C_.col(h) = loading_solver_.f();
            D_.col(h) = Y_h.transpose() * T_.col(h) / T_.col(h).squaredNorm();

            X_h -= T_.col(h) * (loading_solver_.Psi() * C_.col(h)).transpose();
            Y_h -= T_.col(h) * D_.col(h).transpose();
        }
        B_ = W_ * (C_.transpose() * loading_solver_.Psi().transpose() * loading_solver_.Psi() * W_)
                   .partialPivLu()
                   .solve(D_.transpose());
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
        direction_lambda_.resize(n_comp_, direction_n_lambda);
        loading_lambda_.resize(n_comp_, loading_n_lambda);

        int calibration = (flag & 0b11110);
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
                matrix_t M = Y_h.transpose() * X_h;
                Eigen::JacobiSVD<matrix_t> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
                vector_t f0 = svd.matrixV().col(0);
                auto direction_gcv = [&](auto lambda) {
                    return direction_gcv_(M, lambda, f0, max_iter, tol, edf_r, seed);
                };
                auto loading_gcv = [&](auto lambda) {
                    return loading_gcv_(X_h, T_.col(h), lambda, edf_r, seed);
                };
                GridSearch<direction_n_lambda> direction_optimizer;
                direction_lambda = direction_optimizer.optimize(direction_gcv, direction_lambda_grid);

                fit_direction_(M, direction_lambda, f0, max_iter, tol, h);
                T_.col(h) = X_h * direction_solver_.Psi() * W_.col(h);

                GridSearch<loading_n_lambda> loading_optimizer;
                loading_lambda = loading_optimizer.optimize(loading_gcv, loading_lambda_grid);
            } break;
            default: {
                throw std::runtime_error("Unrecognized calibration option.");
            }
            }

            if (calibration != OptimizeGCV) {
                Eigen::JacobiSVD<matrix_t> svd(Y_h.transpose() * X_h, Eigen::ComputeThinU | Eigen::ComputeThinV);
                fit_direction_(Y_h.transpose() * X_h, direction_lambda, svd.matrixV().col(0), max_iter, tol, h);
                T_.col(h) = X_h * direction_solver_.Psi() * W_.col(h);
            }

            loading_solver_.update_response(X_h.transpose() * T_.col(h) / T_.col(h).squaredNorm());
            loading_solver_.fit(loading_lambda);
            C_.col(h) = loading_solver_.f();
            D_.col(h) = Y_h.transpose() * T_.col(h) / T_.col(h).squaredNorm();
            for (int j = 0; j < direction_n_lambda; ++j) { direction_lambda_(h, j) = direction_lambda[j]; }
            for (int j = 0; j < loading_n_lambda; ++j) { loading_lambda_(h, j) = loading_lambda[j]; }

            X_h -= T_.col(h) * (loading_solver_.Psi() * C_.col(h)).transpose();
            Y_h -= T_.col(h) * D_.col(h).transpose();
        }
        B_ = W_ * (C_.transpose() * loading_solver_.Psi().transpose() * loading_solver_.Psi() * W_)
                   .partialPivLu()
                   .solve(D_.transpose());
    }

    const matrix_t& X_space_directions() const { return W_; }
    const matrix_t& Y_space_directions() const { return V_; }
    const matrix_t& X_latent() const { return T_; }
    const matrix_t& X_loadings() const { return C_; }
    const matrix_t& Y_loadings() const { return D_; }
    matrix_t fitted() const { return T_ * D_.transpose(); }
    matrix_t reconstructed() const { return T_ * (loading_solver_.Psi() * C_).transpose(); }
    const matrix_t& B() const { return B_; }
    const matrix_t& direction_lambda() const { return direction_lambda_; }
    const matrix_t& loading_lambda() const { return loading_lambda_; }

   private:
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

        for (int i = 0; !almost_equal(Jnew, Jold, tol) && i < max_iter; ++i) {
            v = M * fn;
            v.normalize();
            direction_solver_.update_response(M.transpose() * v);
            direction_solver_.fit(lambda);
            Jold = Jnew;
            fn = direction_solver_.fn();
            Jnew = (M - v * fn.transpose()).squaredNorm() + direction_solver_.ftPf(lambda);
        }

        return std::make_pair(direction_solver_.f(), v);
    }
    template <typename Lambda, typename Init>
        requires(internals::is_subscriptable<Lambda, int>)
    void fit_direction_(const matrix_t& M, const Lambda& lambda, const Init& f0, int max_iter, double tol, int h) {
        const auto& [f, v] = solve_direction_(M, lambda, f0, max_iter, tol);
        const double norm = std::sqrt(direction_solver_.f().dot(direction_solver_.mass() * direction_solver_.f()));
        W_.col(h) = f / norm;
        V_.col(h) = v;
    }
    template <typename Lambda, typename Init>
        requires(internals::is_subscriptable<Lambda, int>)
    double direction_gcv_(
      const matrix_t& M, const Lambda& lambda, const Init& f0, int max_iter, double tol, int edf_r, int seed) {
        const auto& [f, v] = solve_direction_(M, lambda, f0, max_iter, tol);
        double dor = n_locs_ - direction_solver_.edf(lambda, edf_r, seed);
        return (n_locs_ / std::pow(dor, 2)) *
               ((direction_solver_.Psi() * f) - direction_solver_.response()).squaredNorm();
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
    matrix_t C_;   // X loadings
    matrix_t D_;   // Y loadings
    matrix_t B_;   // regression operator
    matrix_t direction_lambda_;
    matrix_t loading_lambda_;
};

template <typename GeoFrame, typename DirectionPenalty, typename LoadingPenalty>
fPLS(
  const std::string& colname, const Eigen::Matrix<double, Dynamic, Dynamic>& Y, const GeoFrame& gf,
  DirectionPenalty&& direction_penalty, LoadingPenalty&& loading_penalty)
  -> fPLS<typename DirectionPenalty::solver_t, typename LoadingPenalty::solver_t>;

}   // namespace fdapde

#endif   // __FPLS_H__

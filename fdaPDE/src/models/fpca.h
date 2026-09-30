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

#ifndef __FPCA_H__
#define __FPCA_H__

#include "header_check.h"

namespace fdapde {

[[maybe_unused]] constexpr int ComputeRandSVD = 0x1;
[[maybe_unused]] constexpr int ComputeXactSVD = 0x0;

// bits 1 - 5 reserved to calibration strategies
[[maybe_unused]] constexpr int NoCalibration = 0x0;
[[maybe_unused]] constexpr int OptimizeGCV = 0x1 << 1;
[[maybe_unused]] constexpr int OptimizeMSRE = 0x2 << 1;

namespace internals {

/// @brief estimates sequential smooth principal components by alternating scores and loading smoothing
template <typename VariationalSolver> class fpca_power_iteration_impl {
   private:
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;

    /// @brief retains a converged component and its objective diagnostics
    struct fit_result {
        vector_t f;
        vector_t s;
        std::vector<double> objective_history;
        int iterations = 0;
        bool monotone = true;
    };
   public:
    using smoother_t = std::decay_t<VariationalSolver>;
    static constexpr int n_lambda = smoother_t::n_lambda;

    /// @brief creates an unbound power iteration solver
    fpca_power_iteration_impl() noexcept = default;
    /// @brief binds the smoother with default iteration controls
    fpca_power_iteration_impl(VariationalSolver& smoother) noexcept :
        smoother_(std::addressof(smoother)), n_dofs_(smoother.n_dofs()) { }
    /// @brief binds the smoother and explicit objective stopping controls
    fpca_power_iteration_impl(VariationalSolver& smoother, int max_iter, double tol) noexcept :
        smoother_(std::addressof(smoother)), n_dofs_(smoother.n_dofs()), max_iter_(max_iter), tol_(tol) { }

    /// @brief estimates sequential components with fixed penalties or explicit-probe GCV calibration
    template <typename DataT>
    auto fit(
      const DataT& data, int rank, const std::vector<double>& lambda_grid, int flag, int edf_r = 100,
      int seed = random_seed) {
        fdapde_assert(lambda_grid.size() > 0 && lambda_grid.size() % n_lambda == 0);
        edf_r_ = edf_r;
        seed_ = seed;
        edf_map_.clear();
        gcv_values_.clear();
        gcv_edf_.clear();
        selected_indices_.clear();
        matrix_t X = data.transpose();
        n_locs_ = X.cols(), n_units_ = X.rows();
        // first guess of PCs set to a multivariate PCA (SVD)
        matrix_t V;
        if (flag & ComputeRandSVD) {
            RSI<matrix_t> svd(X, rank);
            V = std::move(svd.matrixV());
        } else {
            Eigen::JacobiSVD<matrix_t> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
            V = std::move(svd.matrixV());
        }
        // allocate memory
        f_.resize(n_dofs_, rank);
        s_.resize(n_units_, rank);
        f_norm_.resize(rank);
        lambda_.resize(rank, n_lambda);
        objective_history_.resize(rank);
        iterations_.assign(rank, 0);
        monotone_.assign(rank, true);

        int calibration = (flag & 0b11110);   // detect calibration strategy
        if (calibration == OptimizeGCV) {
            gcv_values_.resize(rank);
            gcv_edf_.resize(rank);
            selected_indices_.resize(rank);
        }
        for (int i = 0; i < rank; ++i) {
            // select optimal smoothing level for i-th component
            Eigen::Matrix<double, n_lambda, 1> opt_lambda;
            switch (calibration) {
            case 0: {   // no calibration
                fdapde_assert(lambda_grid.size() == n_lambda);
                std::copy(lambda_grid.begin(), lambda_grid.end(), opt_lambda.begin());
            } break;
            case OptimizeGCV: {
                auto gcv_functor = [&](auto lambda) {
                    const double value = gcv_(X, lambda, V.col(i));
                    gcv_edf_[i].push_back(last_edf_);
                    return value;
                };
                GridSearch<n_lambda> optimizer;
                opt_lambda = optimizer.optimize(gcv_functor, lambda_grid);
                gcv_values_[i] = optimizer.values();
                selected_indices_[i] =
                  std::distance(gcv_values_[i].begin(), std::min_element(gcv_values_[i].begin(), gcv_values_[i].end()));
            } break;
            case OptimizeMSRE: {
            } break;
            default: {
                throw std::runtime_error("Unrecognized calibration option.");
            }
            }
            // fit with optimal lambda
            auto result = solve_(X, opt_lambda, V.col(i));
            for (int j = 0; j < n_lambda; ++j) { lambda_(i, j) = opt_lambda[j]; }
            // store results
            f_norm_[i] = std::sqrt(result.f.dot(smoother_->mass() * result.f));
            if (!std::isfinite(f_norm_[i]) || f_norm_[i] <= 0) {
                throw std::runtime_error("fPCA loading has zero or non-finite functional norm");
            }
            f_.col(i) = result.f / f_norm_[i];
            s_.col(i) = result.s * f_norm_[i];
            objective_history_[i] = std::move(result.objective_history);
            iterations_[i] = result.iterations;
            monotone_[i] = result.monotone;
            // deflate
            X = X - s_.col(i) * (smoother_->Psi() * f_.col(i)).transpose();
        }
        return std::tie(f_, s_);
    }
    /// @brief returns normalized subject scores
    const matrix_t& scores() const { return s_; }
    /// @brief returns functional L2-normalized loading coefficients
    const matrix_t& loading() const { return f_; }
    /// @brief returns loading norms transferred to the scores
    const std::vector<double>& loadings_norm() const { return f_norm_; }
    /// @brief returns the fitted penalty tuple for each component
    const matrix_t& lambda() const { return lambda_; }
    /// @brief returns the bound smoother
    const smoother_t* smoother() const { return smoother_; }
    /// @brief returns objective traces for the retained component fits
    const std::vector<std::vector<double>>& objective_history() const { return objective_history_; }
    /// @brief returns completed power updates per component
    const std::vector<int>& iterations() const { return iterations_; }
    /// @brief reports whether each objective trace avoided relative increases beyond tolerance
    const std::vector<bool>& monotone() const { return monotone_; }
    /// @brief returns actual GCV evaluations per component in candidate order
    const std::vector<std::vector<double>>& gcv_values() const { return gcv_values_; }
    /// @brief returns the EDF estimates used by each recorded GCV evaluation
    const std::vector<std::vector<double>>& gcv_edf() const { return gcv_edf_; }
    /// @brief returns the zero-based first minimizing candidate for each component
    const std::vector<int>& selected_indices() const { return selected_indices_; }
   private:
    /// @brief alternates a unit score with a smoothed loading while monitoring the penalized objective
    template <typename LambdaT, typename InitT>
        requires(internals::is_subscriptable<LambdaT, int>)
    auto solve_(const matrix_t& X, const LambdaT& lambda, const InitT& f0) {
        // initialization
        vector_t fn = f0;
        vector_t s(n_units_);
        double Jold = std::numeric_limits<double>::max(), Jnew = 1.0;
        int n_iter = 0;
        fit_result result;
        result.objective_history.reserve(max_iter_);
        while (!almost_equal(Jnew, Jold, tol_) && n_iter < max_iter_) {
            // project the current loading and normalize its subject score
            s = X * fn;
            const double score_norm = s.norm();
            if (!s.allFinite() || !std::isfinite(score_norm) || score_norm <= 0) {
                throw std::runtime_error("fPCA score update is zero or non-finite");
            }
            s /= score_norm;
            // smooth the conditional loading response for the current unit score
            smoother_->update_response(X.transpose() * s);
            smoother_->fit(lambda);
            // prepare for next iteration
            n_iter++;
            fn = smoother_->Psi() * smoother_->f();
            if (!fn.allFinite() || !smoother_->f().allFinite()) {
                throw std::runtime_error("fPCA smoother produced a non-finite loading");
            }
            Jold = Jnew;
            Jnew = (X - s * fn.transpose()).squaredNorm() + smoother_->ftPf(lambda);
            if (!std::isfinite(Jnew)) throw std::runtime_error("fPCA objective is non-finite");
            result.objective_history.push_back(Jnew);
            result.iterations = n_iter;
            if (!std::isfinite(Jnew) || (n_iter > 1 && (Jnew - Jold) / (1.0 + std::abs(Jold)) > tol_)) {
                result.monotone = false;
            }
        }
        result.f = smoother_->f();
        result.s = std::move(s);
        return result;
    }
    /// @brief evaluates a candidate and retains its single cached EDF estimate
    template <typename LambdaT, typename InitT>
        requires(internals::is_subscriptable<LambdaT, int>)
    double gcv_(const matrix_t& X, const LambdaT lambda, const InitT& f0) {
        const auto result = solve_(X, lambda, f0);
        // evaluate GCV index at convergence
        std::array<double, n_lambda> lambda_vec;
        for (int i = 0; i < n_lambda; ++i) { lambda_vec[i] = lambda[i]; }
        if (edf_map_.find(lambda_vec) == edf_map_.end()) {   // cache Tr[S]
            edf_map_[lambda_vec] = smoother_->edf(edf_r_, seed_);
        }
        last_edf_ = edf_map_.at(lambda_vec);
        const double dor = n_locs_ - last_edf_;
        if (!std::isfinite(last_edf_) || !std::isfinite(dor) || dor <= 0) {
            throw std::runtime_error("fPCA GCV has no positive residual degrees of freedom");
        }
        const double value =
          (n_locs_ / std::pow(dor, 2)) * ((smoother_->Psi() * result.f) - smoother_->response()).squaredNorm();
        if (!std::isfinite(value)) throw std::runtime_error("fPCA GCV evaluation is non-finite");
        return value;
    }
    std::unordered_map<std::array<double, n_lambda>, double, internals::std_array_hash<double, n_lambda>> edf_map_;
    std::vector<std::vector<double>> gcv_values_, gcv_edf_;
    std::vector<int> selected_indices_;
    int edf_r_ = 100, seed_ = random_seed;
    double last_edf_ = 0;
    int n_locs_ = 0, n_units_ = 0, n_dofs_ = 0;
    smoother_t* smoother_;         // smoothing variational solver
    matrix_t f_;                   // PCs expansion coefficient vector
    matrix_t s_;                   // PCs scores
    std::vector<double> f_norm_;   // L^2 norm of estimated PCs
    matrix_t lambda_;              // selected PCs smoothing level
    std::vector<std::vector<double>> objective_history_;
    std::vector<int> iterations_;
    std::vector<bool> monotone_;

    // power iteration algorithm parameters
    double tol_ = 1e-6;
    int max_iter_ = 20;
};

template <typename VariationalSolver> class fpca_subspace_iteration_impl {
   private:
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using svd_t    = Eigen::JacobiSVD<matrix_t>;
   public:
    using smoother_t = std::decay_t<VariationalSolver>;
    static constexpr int n_lambda = smoother_t::n_lambda;
    
    fpca_subspace_iteration_impl() noexcept = default;
    fpca_subspace_iteration_impl(VariationalSolver& smoother) noexcept :
        smoother_(std::addressof(smoother)), n_dofs_(smoother.n_dofs()) { }
    fpca_subspace_iteration_impl(VariationalSolver& smoother, int max_iter, double tol) noexcept :
        smoother_(std::addressof(smoother)), n_dofs_(smoother.n_dofs()), max_iter_(max_iter), tol_(tol) { }

    /// @brief jointly estimates the requested components with fixed or GCV-selected penalties
    template <typename DataT> auto fit(const DataT& data, int rank, const std::vector<double>& lambda_grid, int flag) {
        fdapde_assert(lambda_grid.size() > 0 && lambda_grid.size() % n_lambda == 0);
        matrix_t X = data.transpose();
        n_locs_ = X.cols(), n_units_ = X.rows();
        // first guess of PCs set to a multivariate PCA (SVD)
        matrix_t V;
        if (flag & ComputeRandSVD) {
            RSI<matrix_t> svd(X, rank);
            V = std::move(svd.matrixV());
        } else {
            Eigen::JacobiSVD<matrix_t> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
            V = svd.matrixV().leftCols(rank);
        }
        // allocate memory
        f_.resize(n_dofs_, rank);
        s_.resize(n_units_, rank);
        f_norm_.resize(rank);
        lambda_.resize(rank, n_lambda);

        int calibration = (flag & 0b11110);   // detect calibration strategy
        Eigen::Matrix<double, n_lambda, 1> opt_lambda;
        switch (calibration) {
        case 0: {   // no calibration
            fdapde_assert(lambda_grid.size() == n_lambda);
            std::copy(lambda_grid.begin(), lambda_grid.end(), opt_lambda.begin());
        } break;
        case OptimizeGCV: {
            auto gcv_functor = [&](auto lambda) { return gcv_(X, rank, lambda, V); };
            GridSearch<n_lambda> optimizer;
            auto opt_ = optimizer.optimize(gcv_functor, lambda_grid);
            for (int i = 0; i < n_lambda; ++i) { opt_lambda[i] = opt_[i]; }
        } break;
        case OptimizeMSRE: {
        } break;
        default: {
            throw std::runtime_error("Unrecognized calibration option.");
        }
        }
        // fit with optimal lambda
        const auto& [F, S] = solve_(X, rank, opt_lambda, V);
	// store results
        for (int i = 0; i < rank; ++i) {
            for (int j = 0; j < n_lambda; ++j) { lambda_(i, j) = opt_lambda[j]; }
            f_norm_[i] = std::sqrt(F.col(i).dot(smoother_->mass() * F.col(i)));   // L^2 norm
            f_.col(i) = F.col(i) / f_norm_[i];
	    s_.col(i) = S.col(i) * f_norm_[i];
        }
        return std::tie(f_, s_);
    }
    // observers
    const matrix_t& scores() const { return s_; }
    const matrix_t& loading() const { return f_; }
    const std::vector<double>& loadings_norm() const { return f_norm_; }
    const matrix_t& lambda() const { return lambda_; }
    const smoother_t* smoother() const { return smoother_; }
  private:
    // finds matrices S, F minimizing \norm{X - S * F^\top}_F^2 + \sum_{i=1}^rank P_{\lambda_i}(f_i)
    template <typename LambdaT, typename InitT>
        requires(internals::is_subscriptable<LambdaT, int>)
    auto solve_(const matrix_t& X, int rank, const LambdaT& lambda, const InitT& F0) {
        // initialization
        matrix_t Fn = F0;
	matrix_t F(n_dofs_, rank);
        matrix_t S(n_units_, rank);
        double Jold = std::numeric_limits<double>::max(), Jnew = 1.0;
        int n_iter = 0;
        while (!almost_equal(Jnew, Jold, tol_) && n_iter < max_iter_) {
            // solve the orthogonal procrustes problem
            // S = \argmin \| X - S * F^\top \|_F^2 subject to S^\top * S = I
            S = X * Fn;
	    svd_t svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
            S = svd.matrixU();
            // f_j = \argmin_f \sum_i (y_i - f_j(p_i))^2 + \int_D (\Delta f_j)^2, with y = X^\top * S_j,
	    // j = 1, ..., rank
            double pen = 0;
            for (int j = 0; j < rank; ++j) {
                smoother_->update_response(X.transpose() * S.col(j));
                smoother_->fit(lambda);
		F .col(j) = smoother_->f();
		Fn.col(j) = smoother_->Psi() * smoother_->f();
		pen = pen + smoother_->ftPf(lambda);
            }
            // prepare for next iteration
            n_iter++;
            Jold = Jnew;
            Jnew = (X - S * Fn.transpose()).squaredNorm() + pen;
        }
        return std::make_pair(F, S);
    }
    // finds vectors s, f minimizing \norm{X - s * f^\top}_F^2 + P_{\lambda}(f) and returns the GCV index
    template <typename LambdaT>
        requires(internals::is_subscriptable<LambdaT, int>)
    double gcv_(const matrix_t& X, int rank, const LambdaT lambda, const matrix_t F0) {
        const auto& [F, S] = solve_(X, rank, lambda, F0);
        // evaluate GCV index at convergence
        int dor = n_locs_ - smoother_->edf(lambda);
        return (n_locs_ / std::pow(dor, 2)) * (X.transpose() * S - (smoother_->Psi() * F)).squaredNorm();
    }
  
    int n_locs_ = 0, n_units_ = 0, n_dofs_ = 0;
    smoother_t* smoother_;         // smoothing variational solver
    matrix_t f_;                   // PCs expansion coefficient vector
    matrix_t s_;                   // PCs scores
    std::vector<double> f_norm_;   // L^2 norm of estimated PCs
    matrix_t lambda_;              // selected PCs smoothing level

    // subspace iteration algorithm parameters
    double tol_ = 1e-6;
    int max_iter_ = 20;
};

// direct fPCA
template <typename VariationalSolver> class fpca_direct_impl {
   private:
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
   public:
    using smoother_t = std::decay_t<VariationalSolver>;
    static constexpr int n_lambda = smoother_t::n_lambda;
  
    fpca_direct_impl() noexcept = default;
    fpca_direct_impl(VariationalSolver& smoother) noexcept :
        smoother_(std::addressof(smoother)), n_dofs_(smoother.n_dofs()) { }

    template <typename DataT> auto fit(const DataT& data, int rank, const std::vector<double>& lambda_grid, int flag) {
        fdapde_assert(lambda_grid.size() > 0 && lambda_grid.size() % n_lambda == 0);
        matrix_t X = data.transpose();
        n_locs_ = X.cols(), n_units_ = X.rows();
        // allocate memory
        f_.resize(n_dofs_, rank);
        s_.resize(n_units_, rank);
        f_norm_.resize(rank);
        lambda_.resize(rank, n_lambda);
	
        int calibration = (flag & 0b11110);   // detect calibration strategy
        Eigen::Matrix<double, n_lambda, 1> opt_lambda;
        switch (calibration) {
        case 0: {   // no calibration
            fdapde_assert(lambda_grid.size() == n_lambda);
            std::copy(lambda_grid.begin(), lambda_grid.end(), opt_lambda.begin());
        } break;
        case OptimizeGCV: {
            auto gcv_functor = [&](auto lambda) { return gcv_(X, rank, lambda, flag); };
            GridSearch<n_lambda> optimizer;
            auto opt_ = optimizer.optimize(gcv_functor, lambda_grid);
            for (int i = 0; i < n_lambda; ++i) { opt_lambda[i] = opt_[i]; }
        } break;
        case OptimizeMSRE: {
        } break;
        default: {
            throw std::runtime_error("Unrecognized calibration option.");
        }
        }
        // fit with optimal lambda
        const auto& [f, s] = solve_(X, rank, opt_lambda, flag);
        // store results
        for (int i = 0; i < rank; ++i) {
            for (int j = 0; j < n_lambda; ++j) { lambda_(i, j) = opt_lambda[j]; }
            f_norm_[i] = std::sqrt(f.col(i).dot(smoother_->mass() * f.col(i)));   // L^2 norm
            f_.col(i) = f.col(i) / f_norm_[i];
	    s_.col(i) = s.col(i) * f_norm_[i];
        }
        return std::tie(f_, s_);
    }
    // observers
    const matrix_t& scores() const { return s_; }
    const matrix_t& loading() const { return f_; }
    const std::vector<double>& loadings_norm() const { return f_norm_; }
    const matrix_t& lambda() const { return lambda_; }
    const smoother_t* smoother() const { return smoother_; }
   private:
    // finds vectors s, f minimizing \norm{X - s * f^\top}_F^2 + P_{\lambda}(f)
    template <typename LambdaT>
        requires(internals::is_subscriptable<LambdaT, int>)
    auto solve_(const matrix_t& X, int rank, const LambdaT& lambda, int flag) {
        for (int i = 0; i < lambda.size(); ++i) { fdapde_assert(lambda[i] > 0); }
        matrix_t C = smoother_->Psi().transpose() * smoother_->Psi() + smoother_->P(lambda);
        // given the cholesky decomposition of C as C = D * D^\top, compute D^{-1}
        Eigen::LLT<matrix_t> chol(C);
        invD_ = chol.matrixL().solve(matrix_t::Identity(n_dofs_, n_dofs_));
        // compute SVD of X * \Psi * (D^{-1})^\top
        matrix_t V, s;
	vector_t singularValues;
        if (flag & ComputeRandSVD) {
            RSI<matrix_t> svd(X * smoother_->Psi() * invD_.transpose(), rank);
            V = std::move(svd.matrixV());
            singularValues = std::move(svd.singularValues());
	    s = svd.matrixU().leftCols(rank);
        } else {
            Eigen::JacobiSVD<matrix_t> svd(
              X * smoother_->Psi() * invD_.transpose(), Eigen::ComputeThinU | Eigen::ComputeThinV);
            V = std::move(svd.matrixV());
            singularValues = std::move(svd.singularValues());
	    s = svd.matrixU().leftCols(rank);
        }
        matrix_t f = (singularValues.head(rank).asDiagonal() * V.leftCols(rank).transpose() * invD_).transpose();
	return std::make_pair(f, s);
    }
    template <typename LambdaT>
        requires(internals::is_subscriptable<LambdaT, int>)
    double gcv_(const matrix_t& X, int rank, const LambdaT lambda, int flag) {
        const auto& [F, S] = solve_(X, rank, lambda, flag);
        // evaluate GCV index at convergence (note that Tr[S] = \|D^(-1)\|_F^2)
        int dor = n_locs_ - invD_.squaredNorm();
        return (n_locs_ / std::pow(dor, 2)) * (X.transpose() * S - (smoother_->Psi() * F)).squaredNorm();
    }
    int n_locs_ = 0, n_units_ = 0, n_dofs_ = 0;
    matrix_t invD_;                // inverse of the cholesky factor of \Psi^\top * \Psi + P_{\lambda}
    smoother_t* smoother_;         // smoothing variational solver
    matrix_t f_;                   // PCs expansion coefficient vector
    matrix_t s_;                   // PCs scores
    std::vector<double> f_norm_;   // L^2 norm of estimated PCs
    matrix_t lambda_;              // selected PCs smoothing level
};
  
// class for handling nan
template <typename fPCASolver> class fpca_na_impl {
   private:
    using fpca_t = std::decay_t<fPCASolver>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using binary_t = BinaryMatrix<Dynamic, Dynamic>;
   public:
    using smoother_t = typename fPCASolver::smoother_t;
    static constexpr int n_lambda = smoother_t::n_lambda;

    fpca_na_impl() noexcept = default;
    fpca_na_impl(fPCASolver& fpca) noexcept :
        fpca_(std::addressof(fpca)), smoother_(fpca.smoother()), n_dofs_(smoother_->n_dofs()) { }

    template <typename DataT>
    auto fit(const DataT& data, int rank, const std::vector<double>& lambda_grid, int flag) {
        matrix_t X = data.transpose();         // create temporary of mapped data
        binary_t nan_pattern = na_matrix(X);   // compute missingness pattern

        int n_locs_ = X.cols(), n_units_ = X.rows();
        // initialization
        f_.resize(n_dofs_, rank);
        s_.resize(n_units_, rank);
        f_norm_.resize(rank);

        int calibration = (flag & 0b11110);   // detect calibration strategy
        std::vector<double> opt_lambda(n_lambda);
        switch (calibration) {
        case 0: {   // no calibration
            fdapde_assert(lambda_grid.size() == n_lambda);
            std::copy(lambda_grid.begin(), lambda_grid.end(), opt_lambda.begin());
        } break;
        case OptimizeMSRE: {
        } break;
        default: {
            throw std::runtime_error("Unrecognized calibration option.");
        }
        }
        // fit with optimal lambda
        const auto& [F, S] = solve_(X, nan_pattern, rank, opt_lambda, flag);
        // store results
        for (int i = 0; i < rank; ++i) {
            f_norm_[i] = std::sqrt(F.col(i).dot(smoother_->mass() * F.col(i)));   // L^2 norm
            f_.col(i) = F.col(i) / f_norm_[i];
	    s_.col(i) = S.col(i) * f_norm_[i];
        }
        return std::tie(f_, s_);
    }
    // observers
    const matrix_t& scores() const { return s_; }
    const matrix_t& loading() const { return f_; }
    const std::vector<double>& loadings_norm() const { return f_norm_; }
   private:
    template <typename LambdaT>
        requires(internals::is_subscriptable<LambdaT, int>)
    auto solve_(const matrix_t& X, const binary_t& nan, int rank, const LambdaT& lambda, int flag) {
        for (int i = 0; i < lambda.size(); ++i) { fdapde_assert(lambda[i] > 0); }
	matrix_t F(n_dofs_ , rank);
	matrix_t S(n_units_, rank);
        matrix_t U  = matrix_t::Zero(n_units_, n_dofs_);
	matrix_t Un = matrix_t::Zero(n_units_, n_dofs_);
	// repeat for increasing rank
        for (int k = 1; k <= rank; ++k) {
            int n_iter = 0;
            double Jold = std::numeric_limits<double>::max(), Jnew = 1.0;
            while (!almost_equal(Jnew, Jold, tol_) && n_iter < max_iter_) {
                // imputation update
                matrix_t Xn = (~nan).select(X, Un);
                Xn.rowwise() -= Xn.colwise().mean();   // re-center
                // rank-k fPCA on imputed data
                const auto& [f, s] = fpca_->fit(X, k, lambda, flag);
                U = s * f.transpose();   // reconstruction update
                // prepare for next iteration
                n_iter++;
                Jold = Jnew;
                Un = U * smoother_->Psi().transpose();
                Jnew = ((~nan).select(X - Un, 0)).squaredNorm() + (U * smoother_->P(lambda) * U.transpose()).trace();
                if (almost_equal(Jnew, Jold, tol_) || n_iter == max_iter_) {
                    F.leftCols(k) = f;
                    S.leftCols(k) = s;
                }
            }
        }
	return std::make_pair(F, S);
    }
    int n_locs_ = 0, n_units_ = 0, n_dofs_ = 0;
    fpca_t* fpca_;
    const smoother_t* smoother_;
    matrix_t f_;                   // PCs expansion coefficient vector
    matrix_t s_;                   // PCs scores
    std::vector<double> f_norm_;   // L^2 norm of estimated PCs

    // MM scheme parameters
    double tol_ = 1e-6;
    int max_iter_ = 100;
};

}   // namespace internals

class fpca_power_solver {
    template <typename Smoother> using impl_t = internals::fpca_power_iteration_impl<Smoother>;
   public:
    fpca_power_solver() noexcept : max_iter_(20), tol_(1e-6) { }
    fpca_power_solver(int max_iter, double tol) noexcept : max_iter_(max_iter), tol_(tol) { }
    template <typename Solver> [[nodiscard]] auto get(Solver&& solver) const {
        return impl_t<Solver>(solver, max_iter_, tol_);
    }
   private:
    int max_iter_;
    double tol_;
};
class fpca_subspace_solver {
    template <typename Smoother> using impl_t = internals::fpca_subspace_iteration_impl<Smoother>;
   public:
    fpca_subspace_solver() noexcept : max_iter_(20), tol_(1e-6) { }
    fpca_subspace_solver(int max_iter, double tol) noexcept : max_iter_(max_iter), tol_(tol) { }
    template <typename Solver> [[nodiscard]] auto get(Solver&& solver) const {
        return impl_t<Solver>(solver, max_iter_, tol_);
    }
   private:
    int max_iter_;
    double tol_;
};
class fpca_direct_solver {
    template <typename Smoother> using impl_t = internals::fpca_direct_impl<Smoother>;
   public:
    fpca_direct_solver() noexcept = default;
    template <typename Smoother> [[nodiscard]] auto get(Smoother&& solver) const { return impl_t<Smoother>(solver); }
};

/// @brief estimates smooth principal components through a variational smoother and a solver policy
template <typename VariationalSolver> class fPCA {
   private:
    using smoother_t = std::decay_t<VariationalSolver>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using binary_t = BinaryMatrix<Dynamic, Dynamic>;
    static constexpr int n_lambda = smoother_t::n_lambda;
   public:
    /// @brief creates an uninitialized principal component model
    fPCA() noexcept = default;
    /// @brief discretizes the supplied penalty and binds functional observations
    template <typename GeoFrame, typename Penalty>
    fPCA(const std::string& colname, const GeoFrame& gf, Penalty&& penalty) noexcept : smoother_(), data_() {
        discretize(penalty.get());
        analyze_data(colname, gf);
    }
    /// @brief discretizes the variational penalty and records its coefficient-space size
    template <typename... Args> void discretize(Args&&... args) {
        smoother_.discretize(std::forward<Args>(args)...);
        n_dofs_ = smoother_.n_dofs();
        return;
    }
    /// @brief binds locations-by-subjects observations and configures equal observation weights
    template <typename GeoFrame> void analyze_data(const std::string& colname, const GeoFrame& gf) {
        fdapde_assert(gf.n_layers() == 1);
        data_ = gf[0].data().template col<double>(colname).as_matrix();
        n_locs_ = data_.rows();
        n_units_ = data_.cols();
        smoother_.analyze_data(gf, vector_t::Ones(gf[0].rows()).asDiagonal());
        // detect if data_ has at least one missing value
        has_nan_ = false;
        for (int i = 0; i < n_locs_; ++i) {
            for (int j = 0; j < n_units_; ++j) {
                if (std::isnan(data_(i, j))) {
                    has_nan_ = true;
                    break;
                }
            }
        }
        return;
    }

    /// @brief fits components, retaining actual power-solver GCV curves and explicitly controlled EDF estimates
    /// explicit EDF controls require complete data and the power solver
    template <typename LambdaT, typename Policy = fpca_power_solver>
        requires(internals::is_vector_like_v<LambdaT>)
    auto fit(
      int rank, const LambdaT& lambda_grid, int flag = ComputeRandSVD, Policy policy = Policy(), int edf_r = 100,
      int seed = random_seed) {
        if (
          rank < 1 || rank > std::min(n_locs_, n_units_) || lambda_grid.size() == 0 ||
          lambda_grid.size() % n_lambda != 0 || edf_r < 1) {
            throw std::invalid_argument("invalid fPCA component count, lambda grid or EDF probe count");
        }
        for (int i = 0; i < lambda_grid.size(); ++i) {
            if (!std::isfinite(lambda_grid[i]) || lambda_grid[i] <= 0) {
                throw std::invalid_argument("fPCA smoothing parameters must be finite and positive");
            }
        }
        lambda_grid_.resize(0, n_lambda);
        gcv_values_.clear();
        gcv_edf_.clear();
        selected_indices_.clear();
        auto solver_ = policy.get(smoother_);   // instantiate solver implementation
        f_.resize(n_dofs_, rank);
        s_.resize(n_units_, rank);
        f_norm_.resize(rank);
        lambda_.resize(n_lambda, rank);
        // dispatch to processing logic
        if (has_nan_) {
            // the imputation solver cannot apply explicit stochastic EDF controls
            if (edf_r != 100 || seed != random_seed) {
                throw std::invalid_argument("explicit fPCA EDF controls require complete observations");
            }
            // default to OptimMSRE calibration, if no calibration provided
            if (lambda_grid.size() > n_lambda && (flag & 0b11110) == 0) { flag = flag | OptimizeMSRE; }
            internals::fpca_na_impl mm_scheme(solver_);
            const auto& [f, s] = mm_scheme.fit(data_, rank, lambda_grid, flag);
            f_ = std::move(f);
            s_ = std::move(s);
            f_norm_ = solver_.loadings_norm();
        } else {
            // default to OptimGCV calibration, if no calibration provided
            if (lambda_grid.size() > n_lambda && (flag & 0b11110) == 0) { flag = flag | OptimizeGCV; }

            auto fit_solver = [&]() {
                if constexpr (requires { solver_.fit(data_, rank, lambda_grid, flag, edf_r, seed); }) {
                    return solver_.fit(data_, rank, lambda_grid, flag, edf_r, seed);
                } else {
                    if (edf_r != 100 || seed != random_seed) {
                        throw std::invalid_argument("explicit fPCA EDF controls require the power solver");
                    }
                    return solver_.fit(data_, rank, lambda_grid, flag);
                }
            };
            const auto& [f, s] = fit_solver();
            f_ = std::move(f);
            s_ = std::move(s);
            f_norm_ = solver_.loadings_norm();
            if constexpr (requires {
                              solver_.gcv_values();
                              solver_.gcv_edf();
                              solver_.selected_indices();
                          }) {
                gcv_values_ = solver_.gcv_values();
                gcv_edf_ = solver_.gcv_edf();
                selected_indices_ = solver_.selected_indices();
                if (!gcv_values_.empty()) {
                    lambda_grid_.resize(lambda_grid.size() / n_lambda, n_lambda);
                    for (int i = 0; i < lambda_grid_.rows(); ++i) {
                        for (int j = 0; j < n_lambda; ++j) { lambda_grid_(i, j) = lambda_grid[i * n_lambda + j]; }
                    }
                }
            }
        }
        lambda_ = solver_.lambda();
        if constexpr (requires {
                          solver_.objective_history();
                          solver_.iterations();
                          solver_.monotone();
                      }) {
            objective_history_ = solver_.objective_history();
            iterations_ = solver_.iterations();
            monotone_ = solver_.monotone();
        } else {
            objective_history_.assign(rank, {});
            iterations_.assign(rank, 0);
            monotone_.assign(rank, true);
        }
        return std::tie(f_, s_);
    }
    /// @brief returns fitted subject scores
    const matrix_t& S() const { return s_; }
    /// @brief returns functional loading coefficients
    const matrix_t& F() const { return f_; }
    /// @brief evaluates loadings at the bound observation locations
    matrix_t Fn() const { return smoother_.Psi() * f_; }
    /// @brief returns loading norms transferred to subject scores
    const std::vector<double>& loadings_norm() const { return f_norm_; }
    /// @brief returns the selected penalty tuple per component
    const matrix_t& lambda() const { return lambda_; }
    /// @brief returns actual power-solver GCV evaluations in candidate order; empty without a power GCV fit
    const std::vector<std::vector<double>>& gcv_values() const { return gcv_values_; }
    /// @brief returns the EDF estimates used by the recorded power-solver GCV evaluations
    const std::vector<std::vector<double>>& gcv_edf() const { return gcv_edf_; }
    /// @brief returns candidate penalty tuples in evaluation order; empty without a power GCV fit
    const matrix_t& lambda_grid() const { return lambda_grid_; }
    /// @brief returns the zero-based first minimizing candidate per component
    const std::vector<int>& selected_indices() const { return selected_indices_; }
    /// @brief returns the objective trace for each retained component
    const std::vector<std::vector<double>>& objective_history() const { return objective_history_; }
    /// @brief returns the completed power updates per component, or zero for other policies
    const std::vector<int>& iterations() const { return iterations_; }
    /// @brief reports whether each retained power objective avoided relative increases beyond tolerance
    const std::vector<bool>& monotone() const { return monotone_; }
   private:
    matrix_t data_;         // mapped geoframe data
    smoother_t smoother_;   // variational solver used in the smoothing step
    bool has_nan_ = false;

    int n_locs_ = 0, n_units_ = 0, n_dofs_ = 0;
    matrix_t f_;                   // PCs expansion coefficient vector
    matrix_t s_;                   // PCs scores
    std::vector<double> f_norm_;   // L^2 norm of estimated components
    matrix_t lambda_;              // selected level of smoothing for each component
    matrix_t lambda_grid_;
    std::vector<std::vector<double>> gcv_values_, gcv_edf_;
    std::vector<int> selected_indices_;
    std::vector<std::vector<double>> objective_history_;
    std::vector<int> iterations_;
    std::vector<bool> monotone_;
};

// deduction guide
template <typename GeoFrame, typename Penalty>
fPCA(const std::string& colname, const GeoFrame& gf, Penalty&& solver) -> fPCA<typename Penalty::solver_t>;

}   // namespace fdapde

#endif   // __FPCA_H__

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
[[maybe_unused]] constexpr int OptimizeGCV  = 0x1 << 1;
[[maybe_unused]] constexpr int OptimizeMSRE = 0x2 << 1;
  
namespace internals {

// power iteration based fPCA
// finds vectors s, f minimizing \norm{X - s * f^\top}_F^2 + P_{\lambda}(f)
template <typename VariationalSolver> class fpca_power_iteration_impl {
   private:
    using smoother_t = std::decay_t<VariationalSolver>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    static constexpr int n_lambda = smoother_t::n_lambda;
   public:
    fpca_power_iteration_impl() noexcept = default;
    fpca_power_iteration_impl(VariationalSolver& solver) noexcept : smoother_(std::addressof(solver)) { }
    fpca_power_iteration_impl(VariationalSolver& solver, int max_iter, double tol) noexcept :
        smoother_(std::addressof(solver)), max_iter_(max_iter), tol_(tol) { }

    template <typename DataT> auto fit(const DataT& data, int rank, const std::vector<double>& lambda_grid, int flag) {
        fdapde_assert(lambda_grid.size() > 0 && lambda_grid.size() % n_lambda == 0);
        matrix_t X = data.transpose();
        n_locs_ = X.cols();
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
        f_.resize(smoother_->n_dofs(), rank);
        s_.resize(X.rows(), rank);
        f_norm_.resize(rank);
        lambda_.resize(rank, n_lambda);

	int calibration = (flag & 0b11110); // detect calibration strategy
        for (int i = 0; i < rank; ++i) {
	    // select optimal smoothing level for i-th component
            std::array<double, n_lambda> opt_lambda;
            switch (calibration) {
            case 0: {   // no calibration
                fdapde_assert(lambda_grid.size() == n_lambda);
		std::copy(lambda_grid.begin(), lambda_grid.end(), opt_lambda.begin());
            } break;
            case OptimizeGCV: {
                auto gcv_functor = [&](auto lambda) { return gcv_(X, lambda, V.col(i)); };
                // GridOptimizer<n_lambda> optimizer;
		GridOptimizer<n_lambda> optimizer;
                opt_lambda = optimizer.optimize(gcv_functor, lambda_grid);
            } break;
            case OptimizeMSRE: {
            } break;
            }
            // fit with optimal lambda
            const auto& [f, s] = solve_(X, opt_lambda, V.col(i));
            for (int j = 0; j < n_lambda; ++j) { lambda_(i, j) = opt_lambda[j]; }
            // store results
            f_norm_[i] = std::sqrt(f.dot(smoother_->mass() * f));
            f_.col(i) = f / f_norm_[i];
            s_.col(i) = s;
            // deflate
            X = X - s_.col(i) * (smoother_->Psi() * f_.col(i)).transpose() * f_norm_[i];
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
    template <typename LambdaT, typename InitT>
        requires(internals::is_subscriptable<LambdaT, int>)
    auto solve_(const matrix_t& X, const LambdaT& lambda, const InitT& f0) {
        // initialization
        vector_t fn = f0;
        vector_t s(X.rows());
        double Jold = std::numeric_limits<double>::max(), Jnew = 1.0;
        int n_iter = 0;
        while (!almost_equal(Jnew, Jold, tol_) && n_iter < max_iter_) {
            // s = X * fn / \norm(X * fn)
            s = X * fn;
            s = s / s.norm();
            // f = \argmin_f \sum_i (y_i - f(p_i))^2 + \int_D (\Delta f)^2, with y = X^\top * s
            smoother_->update_response(X.transpose() * s);
            smoother_->fit(lambda);
            // prepare for next iteration
            n_iter++;
            fn = smoother_->Psi() * smoother_->f();
            Jold = Jnew;
            Jnew = (X - s * fn.transpose()).squaredNorm() + smoother_->ftPf(lambda);
        }
        return std::make_pair(smoother_->f(), s);
    }
    // finds vectors s, f minimizing \norm{X - s * f^\top}_F^2 + P_{\lambda}(f) and returns the GCV index
    template <typename LambdaT, typename InitT>
        requires(internals::is_subscriptable<LambdaT, int>)
    double gcv_(const matrix_t& X, const LambdaT lambda, const InitT& f0) {
        const auto& [f, s] = solve_(X, lambda, f0);
        // evaluate GCV index at convergence
        if (edf_map_.find(lambda) == edf_map_.end()) {   // cache Tr[S]
            edf_map_[lambda] = smoother_->edf();
        }
        int dor = n_locs_ - edf_map_.at(lambda);
        return (n_locs_ / std::pow(dor, 2)) * ((smoother_->Psi() * f) - smoother_->response()).squaredNorm();
    }
    std::unordered_map<std::array<double, n_lambda>, double, internals::std_array_hash<double, n_lambda>> edf_map_;
    int n_locs_ = 0;

    smoother_t* smoother_;         // smoothing variational solver
    matrix_t f_;                   // PCs expansion coefficient vector
    matrix_t s_;                   // PCs scores
    std::vector<double> f_norm_;   // L^2 norm of estimated PCs
    matrix_t lambda_;              // selected PCs smoothing level
  
    // power iteration algorithm parameters
    double tol_ = 1e-6;
    int max_iter_ = 20;
};

  template <typename VariationalSolver> class fpca_subspace_iteration_impl {
   private:
    using smoother_t = std::decay_t<VariationalSolver>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using svd_t    = Eigen::JacobiSVD<matrix_t>;
    static constexpr int n_lambda = smoother_t::n_lambda;
   public:
    fpca_subspace_iteration_impl() noexcept = default;
    fpca_subspace_iteration_impl(VariationalSolver& solver) noexcept : smoother_(std::addressof(solver)) { }

    template <typename DataT> auto fit(const DataT& data, int rank, const std::vector<double>& lambda_grid, int flag) {
        fdapde_assert(lambda_grid.size() > 0 && lambda_grid.size() % n_lambda == 0);
        matrix_t X = data.transpose();
        n_locs_ = X.cols();
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
        f_.resize(smoother_->n_dofs(), rank);
        s_.resize(X.rows(), rank);
        f_norm_.resize(rank);
        lambda_.resize(rank, n_lambda);

        int calibration = (flag & 0b11110);   // detect calibration strategy	
        std::array<double, n_lambda> opt_lambda;
        switch (calibration) {
        case 0: {   // no calibration
            fdapde_assert(lambda_grid.size() == n_lambda);
            std::copy(lambda_grid.begin(), lambda_grid.end(), opt_lambda.begin());
        } break;
        case OptimizeGCV: {
            auto gcv_functor = [&](auto lambda) { return gcv_(X, rank, lambda, V); };
            // GridOptimizer<n_lambda> optimizer;
            GridOptimizer<n_lambda> optimizer;
            opt_lambda = optimizer.optimize(gcv_functor, lambda_grid);
        } break;
        case OptimizeMSRE: {
        } break;
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
	matrix_t F(smoother_->n_dofs(), rank);
        matrix_t S(X.rows(), rank);
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
    int n_locs_ = 0;

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
    using smoother_t = std::decay_t<VariationalSolver>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    static constexpr int n_lambda = smoother_t::n_lambda;
   public:
    fpca_direct_impl() noexcept = default;
    fpca_direct_impl(VariationalSolver& solver) noexcept : smoother_(std::addressof(solver)) { }

    template <typename DataT> auto fit(const DataT& data, int rank, const std::vector<double>& lambda_grid, int flag) {
        fdapde_assert(lambda_grid.size() > 0 && lambda_grid.size() % n_lambda == 0);
        matrix_t X = data.transpose();
        n_locs_ = X.cols();	
	// allocate memory
        f_.resize(smoother_->n_dofs(), rank);
        s_.resize(X.rows(), rank);
        f_norm_.resize(rank);
        lambda_.resize(rank, n_lambda);
	
        int calibration = (flag & 0b11110);   // detect calibration strategy
        std::array<double, n_lambda> opt_lambda;
        switch (calibration) {
        case 0: {   // no calibration
            fdapde_assert(lambda_grid.size() == n_lambda);
            std::copy(lambda_grid.begin(), lambda_grid.end(), opt_lambda.begin());
        } break;
        case OptimizeGCV: {
            auto gcv_functor = [&](auto lambda) { return gcv_(X, rank, lambda, flag); };
            // GridOptimizer<n_lambda> optimizer;
            GridOptimizer<n_lambda> optimizer;
            opt_lambda = optimizer.optimize(gcv_functor, lambda_grid);
        } break;
        case OptimizeMSRE: {
        } break;
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
        invD_ = chol.matrixL().solve(matrix_t::Identity(smoother_->n_dofs(), smoother_->n_dofs()));
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
        // evaluate GCV index at convergence
	int dor = n_locs_ - smoother_->edf();	
        return (n_locs_ / std::pow(dor, 2)) * (X.transpose() * S - (smoother_->Psi() * F)).squaredNorm();
    }

    int n_locs_ = 0;   // number of observations
    matrix_t invD_;   // inverse of the cholesky factor of \Psi^\top * \Psi + P_{\lambda}

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
    using smoother_t = typename fPCASolver::smoother_t;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using binary_t = BinaryMatrix<Dynamic, Dynamic>;
    static constexpr int n_lambda = smoother_t::n_lambda;
   public:
    fpca_na_impl() noexcept = default;
    fpca_na_impl(fPCASolver& fpca) noexcept : fpca_(std::addressof(fpca)), smoother_(fpca.smoother()) { }

    template <typename DataT, typename MissingT, typename... LambdaT>
        requires((std::is_convertible_v<LambdaT, double> && ...) && (sizeof...(LambdaT) == n_lambda))
    auto fit(const DataT& data, const binary_t& nan, int rank, int flag, LambdaT... lambda) {
        matrix_t X = data.transpose();   // create temporary of mapped data
        int n_units_ = X.rows();
        int n_dofs_ = smoother_->n_dofs();

        // initialization
        f_.resize(n_dofs_, rank);
        s_.resize(n_units_, rank);
        f_norm_.resize(rank);
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
                const auto& [f, s] = fpca_->fit(X, k, flag, lambda...);
                U = s * f.transpose();   // reconstruction update
                // prepare for next iteration
                n_iter++;
                Jold = Jnew;
                Un = U * smoother_->Psi().transpose();
                Jnew = ((~nan).select(X - Un, 0)).squaredNorm() + (U * smoother_->P(lambda...) * U.transpose()).trace();
            }
        }
        // store results
        for (int i = 0; i < rank; ++i) {
            f_norm_[i] = std::sqrt(f_.col(i).dot(smoother_->mass() * f_.col(i)));   // L^2 norm
            f_.col(i) = f_.col(i) / f_norm_[i];
	    s_.col(i) = s_.col(i) * f_norm_[i];
        }
        return std::tie(f_, s_);
    }
    // observers
    const matrix_t& scores() const { return s_; }
    const matrix_t& loading() const { return f_; }
    const std::vector<double>& loadings_norm() const { return f_norm_; }
   private:
    fpca_t* fpca_;
    smoother_t* smoother_;

    matrix_t f_;                   // PCs expansion coefficient vector
    matrix_t s_;                   // PCs scores
    std::vector<double> f_norm_;   // L^2 norm of estimated PCs

    // MM scheme parameters
    double tol_ = 1e-6;
    int max_iter_ = 100;  
};

}   // namespace internals

template <typename VariationalSolver> class fPCA {
   private:
    using smoother_t = std::decay_t<VariationalSolver>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using data_t   = Eigen::Map<const Eigen::Matrix<double, Dynamic, Dynamic, Eigen::ColMajor>>;
    static constexpr int n_lambda = smoother_t::n_lambda;
   public:
    fPCA() noexcept = default;
    template <typename GeoFrame, typename Penalty>
    fPCA(const std::string& colname, const GeoFrame& gf, Penalty&& penalty) noexcept :
        smoother_(), data_(gf[0].data().template col<double>(colname).as_matrix()) {
        fdapde_assert(gf.n_layers() == 1);
        n_locs_ = data_.rows();
	n_units_ = data_.cols();
        if constexpr (requires(Penalty p) { p.get(); }) {
            smoother_ = smoother_t(gf, penalty.get());
        } else {
            smoother_ = smoother_t(gf, penalty(gf.template triangulation<0>()).get());
        }
	n_dofs_ = smoother_.n_dofs();
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
    }

    template <typename LambdaT, typename Policy = internals::fpca_direct_impl<smoother_t>>
        requires(internals::is_subscriptable<LambdaT, int>)
    auto fit(int rank, const LambdaT& lambda, int flag, Policy solver = Policy()) {
        Policy solver_(smoother_);
        f_.resize(n_dofs_ , rank);
        s_.resize(n_units_, rank);
        f_norm_.resize(rank);
        lambda_.resize(n_lambda, rank);
	// dispatch to processing logic
        if (has_nan_) {
        } else {
            const auto& [f, s] = solver_.fit(data_, rank, lambda, flag);
            f_ = std::move(f);
            s_ = std::move(s);
            f_norm_ = solver_.loadings_norm();
        }
        return std::tie(f_, s_);
    }
    // template <typename... LambdaT>
    //     requires((std::is_convertible_v<LambdaT, double> && ...) && (sizeof...(LambdaT) == n_lambda))
    // std::pair<const matrix_t&, const matrix_t&> fit(int rank, LambdaT... lambda) {
    //     return fit(rank, ComputeRandSVD, lambda...);
    // }
    // observers
    const matrix_t& scores() const { return s_; }
    const matrix_t& loading() const { return f_; }
    const std::vector<double>& loadings_norm() const { return f_norm_; }
    const matrix_t& lambda() const { return lambda_; }
   private:
    data_t data_;           // mapped geoframe data
    smoother_t smoother_;   // variational solver used in the smoothing step
    bool has_nan_;

    int n_locs_ = 0, n_units_ = 0, n_dofs_ = 0;
    matrix_t f_;                   // PCs expansion coefficient vector
    matrix_t s_;                   // PCs scores
    std::vector<double> f_norm_;   // L^2 norm of estimated components
    matrix_t lambda_;              // selected level of smoothing for each component
};

// deduction guide
template <typename GeoFrame, typename Penalty>
fPCA(const std::string& colname, const GeoFrame& gf, Penalty&& solver) -> fPCA<typename Penalty::solver_t>;

}   // namespace fdapde

#endif   // __FPCA_H__

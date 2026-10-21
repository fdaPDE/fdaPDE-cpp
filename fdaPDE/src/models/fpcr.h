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

#ifndef __FPCR_H__
#define __FPCR_H__

#include "header_check.h"

namespace fdapde {

/// @brief regresses centered scalar or multivariate responses on sequential smooth principal component scores
/// predictor data is bound through GeoFrame with locations in rows and statistical units in columns
/// each component prefix refits OLS; beta coefficients map the original centered predictors after sequential deflation
/// supported bindings are single-parameter elliptic finite-element and spline smoothers
/// calibration uses the native fPCA power solver with fixed penalties or explicitly controlled stochastic EDF estimates
///
/// @code
/// fPCR model("X", Y, data, fe_ls_elliptic(penalty, load));
/// model.fit(3, std::vector<double> {0.01, 0.1, 1.0}, ComputeXactSVD | OptimizeGCV, 100, 1e-8, 100, 42);
/// auto fitted = model.fitted(2);
/// auto beta = model.Beta(2);
/// @endcode
template <typename VariationalSolver> class fPCR {
   private:
    using smoother_t = std::decay_t<VariationalSolver>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    static constexpr int n_lambda = smoother_t::n_lambda;
    static_assert(n_lambda == 1, "fPCR requires a single-parameter elliptic finite-element or spline smoother");

    fPCA<smoother_t> fpca_;
    matrix_t X_, Y_;
    matrix_t scores_, directions_, loading_directions_, projection_scales_, sampled_directions_;
    std::vector<double> projection_errors_;
    std::vector<matrix_t> score_coefficients_, beta_;
    int n_comp_ = 0;
   public:
    /// @brief creates an uninitialized model without discretized penalties or data
    fPCR() noexcept = default;

    /// @brief discretizes the principal component penalty and binds centered predictors and responses
    template <typename GeoFrame, typename Penalty>
    fPCR(const std::string& colname, const matrix_t& Y, const GeoFrame& gf, Penalty&& penalty) {
        fpca_.discretize(penalty.get());
        analyze_data(colname, Y, gf);
    }

    /// @brief discretizes the principal component smoothing penalty
    template <typename... Args> void discretize(Args&&... args) { fpca_.discretize(std::forward<Args>(args)...); }

    /// @brief binds finite centered predictors and all response columns to the discretized model
    template <typename GeoFrame> void analyze_data(const std::string& colname, const matrix_t& Y, const GeoFrame& gf) {
        if (gf.n_layers() != 1) throw std::invalid_argument("fPCR requires a single functional data layer");
        matrix_t X = gf[0].data().template col<double>(colname).as_matrix().transpose();
        if (
          X.rows() == 0 || X.cols() == 0 || Y.rows() != X.rows() || Y.cols() == 0 || !X.allFinite() || !Y.allFinite()) {
            throw std::invalid_argument("fPCR requires finite predictor and response matrices with matching units");
        }
        fpca_.analyze_data(colname, gf);
        X_ = std::move(X);
        Y_ = Y;
        clear_();
    }

    /// @brief fits native smooth principal components, sequential score projections and OLS for every prefix
    /// fixed calibration requires one shared penalty tuple; GCV candidates are consecutive tuples of solver width
    void fit(
      int n_comp,                               // number of sequential principal components
      const std::vector<double>& lambda_grid,   // shared fixed penalty tuple or consecutive GCV candidate tuples
      int flag = ComputeXactSVD,                // SVD initialization and NoCalibration or OptimizeGCV option bits
      int max_iter = 20,                        // maximum power updates per candidate or retained component
      double tol = 1e-6,                        // relative objective stopping tolerance
      int edf_r = 100,                          // number of Hutchinson probes for GCV effective degrees of freedom
      int seed = random_seed                    // seed shared across the fit's cached EDF estimates
    ) {
        clear_();
        const int calibration = flag & 0b11110;
        if (
          n_comp < 1 || n_comp > std::min(X_.rows(), X_.cols()) || max_iter < 1 || !std::isfinite(tol) || tol <= 0 ||
          edf_r < 1 || lambda_grid.empty() || lambda_grid.size() % n_lambda != 0 ||
          (calibration != NoCalibration && calibration != OptimizeGCV) ||
          (calibration == NoCalibration && lambda_grid.size() != n_lambda)) {
            throw std::invalid_argument("invalid fPCR component count, solver controls or calibration grid");
        }

        // reuse the native power solver's component fits and calibration diagnostics
        fpca_.fit(n_comp, lambda_grid, flag, fpca_power_solver(max_iter, tol), edf_r, seed);
        const matrix_t sampled_loadings = fpca_.Fn();
        scores_.resize(X_.rows(), n_comp);
        directions_.resize(fpca_.F().rows(), n_comp);
        loading_directions_.resize(fpca_.F().rows(), n_comp);
        sampled_directions_.resize(X_.cols(), n_comp);
        projection_scales_.resize(n_comp, 1);
        projection_errors_.resize(n_comp);
        score_coefficients_.resize(n_comp);
        beta_.resize(n_comp);
        matrix_t residual = X_;

        // project final loadings on successive residuals and undo previous deflations in the original-X score map
        for (int h = 0; h < n_comp; ++h) {
            const matrix_t raw_score = residual * sampled_loadings.col(h);
            const double projected_norm = raw_score.norm();
            const double fitted_norm = fpca_.S().col(h).norm();
            if (
              !std::isfinite(projected_norm) || projected_norm <= 0 || !std::isfinite(fitted_norm) ||
              fitted_norm <= 0) {
                throw std::runtime_error("fPCR cannot project a zero or non-finite component");
            }
            const double scale = fitted_norm / projected_norm;
            projection_scales_(h, 0) = scale;
            loading_directions_.col(h) = scale * fpca_.F().col(h);
            directions_.col(h) = loading_directions_.col(h);
            const matrix_t coupling = sampled_loadings.leftCols(h).transpose() * (scale * sampled_loadings.col(h));
            directions_.col(h) -= directions_.leftCols(h) * coupling;
            sampled_directions_.col(h) = scale * sampled_loadings.col(h);
            sampled_directions_.col(h) -= sampled_directions_.leftCols(h) * coupling;
            scores_.col(h) = scale * raw_score;
            projection_errors_[h] = (scores_.col(h) - fpca_.S().col(h)).cwiseAbs().maxCoeff();
            residual.noalias() -= scores_.col(h) * sampled_loadings.col(h).transpose();
        }
        if (!scores_.allFinite() || !directions_.allFinite()) {
            throw std::runtime_error("fPCR score projection or direction is non-finite");
        }

        // refit OLS on each prefix because sequential smooth principal component scores need not be orthogonal
        for (int h = 1; h <= n_comp; ++h) {
            score_coefficients_[h - 1] = scores_.leftCols(h).colPivHouseholderQr().solve(Y_);
            beta_[h - 1] = directions_.leftCols(h) * score_coefficients_[h - 1];
            if (!score_coefficients_[h - 1].allFinite() || !beta_[h - 1].allFinite()) {
                throw std::runtime_error("fPCR score regression produced non-finite coefficients");
            }
        }
        n_comp_ = n_comp;
    }

    /// @brief returns functional L2-normalized principal component loading coefficients
    const matrix_t& X_loadings() const { return fpca_.F(); }
    /// @brief returns principal component loading coefficients through the native fPCA accessor
    const matrix_t& F() const { return fpca_.F(); }
    /// @brief evaluates loading coefficients at the training observation locations
    matrix_t Fn() const { return fpca_.Fn(); }
    /// @brief returns native power-iteration scores before final-loading sequential projection
    const matrix_t& S() const { return fpca_.S(); }
    /// @brief returns sequential final-loading scores with subjects in rows and components in columns
    const matrix_t& X_latent_scores() const { return scores_; }
    /// @brief returns sequential scores through the fPLS-compatible latent accessor
    const matrix_t& X_latent() const { return scores_; }
    /// @brief returns coefficient-space directions mapping the original predictors to sequential scores
    const matrix_t& X_space_directions() const { return directions_; }
    /// @brief returns scaled final loadings used for score projection on each component's residual
    const matrix_t& loading_directions() const { return loading_directions_; }
    /// @brief returns the score rescaling factor per component as a single-column matrix
    const matrix_t& projection_scales() const { return projection_scales_; }
    /// @brief returns maximum absolute differences between projected and native fitted scores per component
    const std::vector<double>& projection_errors() const { return projection_errors_; }
    /// @brief returns the separately fitted OLS coefficients for a prefix, with one column per response
    const matrix_t& score_coefficients(int h = 0) const { return score_coefficients_[components_(h) - 1]; }
    /// @brief returns response loadings of the full score regression, with response variables in rows
    matrix_t Y_loadings() const { return score_coefficients().transpose(); }
    /// @brief returns original-X functional regression coefficients for a prefix, with one column per response
    const matrix_t& B(int h = 0) const { return beta_[components_(h) - 1]; }
    /// @brief returns original-X regression coefficients through the fPLS-compatible accessor
    const matrix_t& Beta(int h = 0) const { return B(h); }
    /// @brief returns training response predictions for a prefix using its own OLS fit
    matrix_t fitted(int h = 0) const {
        h = components_(h);
        return scores_.leftCols(h) * score_coefficients_[h - 1];
    }
    /// @brief projects finite centered predictors at the training locations onto sequential scores for a prefix
    matrix_t transform(const matrix_t& X, int h = 0) const {
        h = components_(h);
        if (X.cols() != X_.cols() || !X.allFinite()) {
            throw std::invalid_argument("fPCR projection requires finite predictors at the training locations");
        }
        return X * sampled_directions_.leftCols(h);
    }
    /// @brief predicts every response column from centered new predictors using a separately fitted prefix
    matrix_t predict(const matrix_t& X, int h = 0) const {
        h = components_(h);
        return transform(X, h) * score_coefficients_[h - 1];
    }
    /// @brief reconstructs training predictors at the observation locations from a component prefix
    matrix_t reconstructed(int h = 0) const {
        h = components_(h);
        return scores_.leftCols(h) * fpca_.Fn().leftCols(h).transpose();
    }
    /// @brief returns the functional norms transferred to native principal component scores
    const std::vector<double>& loadings_norm() const { return fpca_.loadings_norm(); }
    /// @brief returns the selected penalty tuple per component
    const matrix_t& lambda() const { return fpca_.lambda(); }
    /// @brief returns actual GCV evaluations per component in candidate order; empty after fixed calibration
    const std::vector<std::vector<double>>& gcv_values() const { return fpca_.gcv_values(); }
    /// @brief returns the EDF estimates used by the actual GCV evaluations; empty after fixed calibration
    const std::vector<std::vector<double>>& gcv_edf() const { return fpca_.gcv_edf(); }
    /// @brief returns candidate penalty tuples in evaluation order; empty after fixed calibration
    const matrix_t& lambda_grid() const { return fpca_.lambda_grid(); }
    /// @brief returns the zero-based first minimizing candidate per component; empty after fixed calibration
    const std::vector<int>& selected_indices() const { return fpca_.selected_indices(); }
    /// @brief returns retained principal component objective traces
    const std::vector<std::vector<double>>& objective_history() const { return fpca_.objective_history(); }
    /// @brief returns completed power updates per retained principal component
    const std::vector<int>& iterations() const { return fpca_.iterations(); }
    /// @brief reports whether each principal component objective avoided relative increases beyond tolerance
    const std::vector<bool>& monotone() const { return fpca_.monotone(); }
   private:
    /// @brief clears the regression state before rebinding observations or starting a new fit
    void clear_() {
        n_comp_ = 0;
        scores_.resize(0, 0);
        directions_.resize(0, 0);
        loading_directions_.resize(0, 0);
        sampled_directions_.resize(0, 0);
        projection_scales_.resize(0, 0);
        projection_errors_.clear();
        score_coefficients_.clear();
        beta_.clear();
    }

    /// @brief resolves zero to all fitted components and validates a requested prefix
    int components_(int h) const {
        if (h == 0) h = n_comp_;
        if (h < 1 || h > n_comp_) throw std::out_of_range("fPCR component prefix is outside the fitted model");
        return h;
    }
};

/// @brief deduces the principal component smoother from its variational penalty packet
template <typename GeoFrame, typename Penalty>
fPCR(
  const std::string& colname, const Eigen::Matrix<double, Dynamic, Dynamic>& Y, const GeoFrame& gf, Penalty&& penalty)
  -> fPCR<typename std::decay_t<Penalty>::solver_t>;

}   // namespace fdapde

#endif   // __FPCR_H__

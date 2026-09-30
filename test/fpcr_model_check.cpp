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

#include <fdaPDE/models.h>

#include <cassert>
#include <iostream>

using namespace fdapde;

/// @brief counts native smoother fits and EDF estimates without altering their numerical behavior
struct counted_smoother : internals::fe_ls_elliptic {
    static inline int fits = 0, edfs = 0;
    /// @brief delegates the fit and counts one actual solver invocation
    template <typename Lambda> void fit(const Lambda& lambda) {
        ++fits;
        internals::fe_ls_elliptic::fit(lambda);
    }
    /// @brief delegates explicit-probe EDF estimation and counts the actual trace invocation
    double edf(int probes, int seed) {
        ++edfs;
        return internals::fe_ls_elliptic::edf(probes, seed);
    }
};

/// @brief reports whether a public API operation rejects invalid input
template <typename Callable> bool rejects(Callable&& callable) {
    try {
        callable();
    } catch (const std::exception&) { return true; }
    return false;
}

/// @brief checks public fPCR regression, calibration provenance, repeated fits and invalid inputs on a small FEM domain
int main() {
    auto domain = Triangulation<2, 2>::Rectangle(0, 1, 0, 1, 4, 4);
    FeSpace space(domain, P1<1>);
    TrialFunction f(space);
    TestFunction v(space);
    const auto penalty = integral(domain)(dot(grad(f), grad(v)));
    ZeroField<2> forcing;
    const auto load = integral(domain)(forcing * v);
    const auto smoother_penalty = fe_ls_elliptic(penalty, load);
    Eigen::MatrixXd X(9, domain.n_nodes());
    for (int i = 0; i < X.rows(); ++i) {
        for (int j = 0; j < X.cols(); ++j) {
            X(i, j) = std::sin(0.2 * j + 0.5 * i) + std::cos(0.3 * j - 0.4 * i) + 0.2 * std::sin(0.7 * j + 0.8 * i);
        }
    }
    X.rowwise() -= X.colwise().mean().eval();
    Eigen::MatrixXd Y(X.rows(), 2);
    Y.col(0) = X.col(2) + 0.7 * X.col(5);
    Y.col(1) = X.col(6) - 0.4 * X.col(1);
    GeoFrame data(domain);
    auto& layer = data.insert_scalar_layer<POINT>("observations", MESH_NODES);
    layer.load_blk("X", X.transpose());
    internals::fe_ls_elliptic sampler;
    sampler.discretize(smoother_penalty.get());
    sampler.analyze_data(data, Eigen::VectorXd::Ones(X.cols()).asDiagonal());

    // finite-element beta and unseen predictions must reproduce independently projected scores for every prefix
    fPCR model("X", Y, data, smoother_penalty);
    model.fit(3, std::vector<double> {0.05}, ComputeXactSVD, 100, 1e-9, 7, 42);
    Eigen::MatrixXd unseen = 0.6 * X.topRows(4);
    unseen.array() += 0.02;
    Eigen::MatrixXd residual = unseen;
    Eigen::MatrixXd unseen_scores(unseen.rows(), 3);
    for (int h = 1; h <= 3; ++h) {
        unseen_scores.col(h - 1) = residual * model.Fn().col(h - 1) * model.projection_scales()(h - 1, 0);
        residual -= unseen_scores.col(h - 1) * model.Fn().col(h - 1).transpose();
        const Eigen::MatrixXd scores = model.X_latent_scores().leftCols(h);
        const Eigen::MatrixXd coefficients = scores.colPivHouseholderQr().solve(Y);
        // each prefix coefficient must equal an independent OLS fit on all response columns
        assert((model.score_coefficients(h) - coefficients).norm() < 1e-10);
        // the original-X coefficient map must equal the sequential training scores
        assert((X * sampler.Psi() * model.X_space_directions().leftCols(h) - scores).norm() < 1e-9);
        // beta must reproduce both response columns of prefix OLS on the original predictors
        assert((X * sampler.Psi() * model.Beta(h) - model.fitted(h)).norm() < 1e-9);
        // unseen score projection must match explicit residual deflation using fitted training scales
        assert((model.transform(unseen, h) - unseen_scores.leftCols(h)).norm() < 1e-9);
        // unseen predictions must equal the independently projected prefix scores times their own OLS coefficients
        assert((model.predict(unseen, h) - unseen_scores.leftCols(h) * coefficients).norm() < 1e-9);
        // reconstructed predictors must equal the sum of the first h rank-one projected components
        assert((model.reconstructed(h) - scores * model.Fn().leftCols(h).transpose()).norm() < 1e-9);
    }
    // scalar response fitting must reproduce its column in the multivariate regression
    fPCR scalar("X", Eigen::MatrixXd(Y.leftCols(1)), data, smoother_penalty);
    scalar.fit(3, std::vector<double> {0.05}, ComputeXactSVD, 100, 1e-9, 7, 42);
    assert((scalar.fitted() - model.fitted().leftCols(1)).norm() < 1e-9);
    // nonorthogonal score prefixes require an OLS refit rather than truncating the full regression
    assert((model.score_coefficients(1) - model.score_coefficients().topRows(1)).norm() > 1e-4);

    // reversing predictor signs must preserve regressions when prediction data is reversed consistently
    GeoFrame reversed_data(domain);
    auto& reversed_layer = reversed_data.insert_scalar_layer<POINT>("observations", MESH_NODES);
    reversed_layer.load_blk("X", -X.transpose());
    fPCR reversed("X", Y, reversed_data, smoother_penalty);
    reversed.fit(3, std::vector<double> {0.05}, ComputeXactSVD, 100, 1e-9, 7, 42);
    // the score signs induced by reversed predictors must cancel in every prefix's OLS prediction
    for (int h = 1; h <= 3; ++h) assert((reversed.predict(-unseen, h) - model.predict(unseen, h)).norm() < 1e-9);

    // instrumentation proves that calibration records real evaluations without diagnostic refits or repeated EDF
    fPCR<counted_smoother> calibrated;
    calibrated.discretize(smoother_penalty.get());
    calibrated.analyze_data("X", Y, data);
    const std::vector<double> grid {0.01, 0.1, 0.01};
    counted_smoother::fits = counted_smoother::edfs = 0;
    calibrated.fit(2, grid, ComputeXactSVD | OptimizeGCV, 1, 1e-9, 7, 0);
    // one update per candidate plus the retained fit gives H*(G+1) actual smoother fits
    assert(counted_smoother::fits == 8);
    // a duplicated lambda and subsequent components must reuse the same two EDF estimates
    assert(counted_smoother::edfs == 2);
    const int fits = counted_smoother::fits, edfs = counted_smoother::edfs;
    for (int h = 0; h < 2; ++h) {
        // GCV and EDF curves must retain every evaluation, including repeated candidate tuples
        assert(calibrated.gcv_values()[h].size() == grid.size() && calibrated.gcv_edf()[h].size() == grid.size());
        const auto& values = calibrated.gcv_values()[h];
        const int minimum = std::distance(values.begin(), std::min_element(values.begin(), values.end()));
        // selected indices must identify the first actual minimizing candidate used by GridSearch
        assert(calibrated.selected_indices()[h] == minimum && calibrated.lambda()(h, 0) == grid[minimum]);
        // duplicate grid positions must retain identical cached EDF and deterministic candidate GCV values
        assert(calibrated.gcv_edf()[h][0] == calibrated.gcv_edf()[h][2] && values[0] == values[2]);
    }
    // diagnostic accessors must not perform additional smoothing fits or EDF estimates
    assert(counted_smoother::fits == fits && counted_smoother::edfs == edfs);
    const auto seed_zero_edf = calibrated.gcv_edf();
    calibrated.fit(2, grid, ComputeXactSVD | OptimizeGCV, 1, 1e-9, 7, 42);
    fPCR reference("X", Y, data, smoother_penalty);
    reference.fit(2, grid, ComputeXactSVD | OptimizeGCV, 1, 1e-9, 7, 42);
    // changed-seed refits must match fresh models with that seed and refresh the EDF sample
    assert(calibrated.gcv_edf() == reference.gcv_edf() && calibrated.gcv_edf() != seed_zero_edf);
    // identical explicit seeds must reproduce all GCV evaluations and fitted multivariate responses
    assert(
      calibrated.gcv_values() == reference.gcv_values() && (calibrated.fitted() - reference.fitted()).norm() < 1e-10);
    calibrated.fit(1, std::vector<double> {0.05}, ComputeXactSVD, 2, 1e-9, 7, 42);
    // fixed refits must clear candidate grids, GCV/EDF curves and selection indices from the previous fit
    assert(
      calibrated.gcv_values().empty() && calibrated.gcv_edf().empty() && calibrated.lambda_grid().rows() == 0 &&
      calibrated.selected_indices().empty());
    // reduced-rank refits must replace component and objective diagnostic dimensions
    assert(calibrated.X_latent_scores().cols() == 1 && calibrated.objective_history().size() == 1);

    // permanent public validation must reject invalid prefixes, dimensions, non-finite inputs and solver controls
    assert(rejects([&] { model.fitted(4); }));
    // prediction locations must agree with the fitted sampling design
    assert(rejects([&] { model.predict(Eigen::MatrixXd::Zero(2, X.cols() - 1)); }));
    // non-finite new predictor observations must be rejected before projection
    assert(rejects([&] { model.predict(Eigen::MatrixXd::Constant(2, X.cols(), NAN)); }));
    // every response column must have finite observations and matching subject count
    assert(rejects([&] { fPCR invalid("X", Eigen::MatrixXd::Constant(Y.rows(), 2, NAN), data, smoother_penalty); }));
    // smoothing penalties must be finite and strictly positive
    assert(rejects([&] { model.fit(1, std::vector<double> {-0.1}); }));
    // zero predictors must fail the score update rather than producing singular latent components
    GeoFrame zero_data(domain);
    auto& zero_layer = zero_data.insert_scalar_layer<POINT>("observations", MESH_NODES);
    zero_layer.load_blk("X", Eigen::MatrixXd::Zero(X.cols(), X.rows()));
    fPCR zero("X", Y, zero_data, smoother_penalty);
    assert(rejects([&] { zero.fit(1, std::vector<double> {0.1}); }));
    // EDF probe counts must be positive even when supplied to fixed calibration
    assert(rejects([&] { model.fit(1, std::vector<double> {0.1}, ComputeXactSVD, 2, 1e-9, 0, 42); }));
    std::cout << "fPCR model checks passed\n";
}

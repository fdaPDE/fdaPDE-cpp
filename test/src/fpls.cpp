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

using namespace fdapde;
using fdapde::test::almost_equal;

namespace {

/// @brief counts solver evaluations independently for direction and loading calibration
template <typename Solver, int Block> struct counted_fpls_solver : Solver {
    inline static int edf_calls = 0, fit_calls = 0;
    /// @brief delegates fitting while recording each solve requested by the model
    template <typename... Args> decltype(auto) fit(Args&&... args) {
        ++fit_calls;
        return Solver::fit(std::forward<Args>(args)...);
    }
    /// @brief delegates EDF estimation while recording each trace evaluation
    template <typename... Args> double edf(Args&&... args) {
        ++edf_calls;
        return Solver::edf(std::forward<Args>(args)...);
    }
};

/// @brief decorates an existing penalty packet with calibration evaluation counters
template <typename Penalty, int Block> struct counted_fpls_penalty {
    using solver_t = counted_fpls_solver<typename Penalty::solver_t, Block>;
    const Penalty& penalty;
    /// @brief forwards the existing discretization packet without duplicating its assembly
    decltype(auto) get() const { return penalty.get(); }
};

// check saved calibration values against independent EDF and residual oracles, with no observer evaluations
TEST(fpls, calibration_snapshot) {
    auto D = Triangulation<2, 2>::Rectangle(0, 1, 0, 1, 4, 4);
    GeoFrame data(D);
    auto& layer = data.insert_scalar_layer<POINT>("l1", MESH_NODES);
    Eigen::MatrixXd X(8, D.n_nodes()), Y(8, 2);
    for (int i = 0; i < X.rows(); ++i) {
        for (int j = 0; j < X.cols(); ++j) {
            X(i, j) = std::sin((i + 1) * D.nodes()(j, 0)) + std::cos((i + 2) * D.nodes()(j, 1));
        }
        Y(i, 0) = std::sin(i + 1);
        Y(i, 1) = std::cos(i + 1);
    }
    layer.load_blk("X", X.transpose());
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    const auto penalty = fe_ls_elliptic(a, F);
    counted_fpls_penalty<decltype(penalty), 0> direction {penalty};
    counted_fpls_penalty<decltype(penalty), 1> loading {penalty};
    using direction_t = decltype(direction)::solver_t;
    using loading_t = decltype(loading)::solver_t;
    fPLS<direction_t, loading_t> model("X", Y, data, direction, loading);
    const std::vector<double> grid {0.001, 0.1, 10.0};
    for (int seed : {0, 42, 0}) {
        direction_t::edf_calls = loading_t::edf_calls = loading_t::fit_calls = 0;
        model.fit(2, grid, grid, OptimizeGCV, 15, 1e-8, 12, seed);
        // every component and candidate must trigger exactly one direction EDF evaluation
        EXPECT_EQ(direction_t::edf_calls, 6);
        // every component and candidate must trigger exactly one loading EDF evaluation
        EXPECT_EQ(loading_t::edf_calls, 6);
        // loading solves consist only of candidate evaluation and the selected fit per component
        EXPECT_EQ(loading_t::fit_calls, 8);
        Eigen::MatrixXd residual = X;
        for (int h = 0; h < 2; ++h) {
            const auto& curve = model.loading_gcv_values()[h];
            const int selected = model.loading_selected_indices()[h];
            // stored selection must use the first minimum in the actual evaluated curve
            EXPECT_EQ(selected, std::min_element(curve.begin(), curve.end()) - curve.begin());
            // selected lambda must correspond to the saved candidate row
            EXPECT_DOUBLE_EQ(model.loading_lambda()(h, 0), model.loading_lambda_grid()(selected, 0));
            auto reference = penalty.get();
            internals::fe_ls_elliptic smoother;
            smoother.discretize(reference);
            smoother.analyze_data(data, Eigen::VectorXd::Ones(X.cols()).asDiagonal());
            const Eigen::VectorXd response =
              residual.transpose() * model.X_latent().col(h) / model.X_latent().col(h).squaredNorm();
            smoother.update_response(response);
            for (int j = 0; j < int(grid.size()); ++j) {
                smoother.fit(grid[j]);
                const double edf = smoother.edf(grid[j], 12, seed);
                const double gcv = X.cols() * (smoother.fn() - response).squaredNorm() / std::pow(X.cols() - edf, 2);
                // saved loading EDF must equal an independent same-seed solver evaluation
                EXPECT_DOUBLE_EQ(model.loading_gcv_edfs()[h][j], edf);
                // saved direction EDF uses the same smoother operator and explicit probe sample
                EXPECT_DOUBLE_EQ(model.direction_gcv_edfs()[h][j], edf);
                // saved loading GCV must use the actual component residual and the saved EDF
                EXPECT_NEAR(curve[j], gcv, 1e-12);
            }
            residual -= model.X_latent().col(h) * model.X_loadings().col(h).transpose();
        }
        // observing all diagnostics must leave the model's EDF evaluation count unchanged
        EXPECT_EQ(direction_t::edf_calls + loading_t::edf_calls, 12);
    }
    const Eigen::Matrix<double, 1, 1> fixed = Eigen::Matrix<double, 1, 1>::Constant(0.1);
    for (bool schedule : {false, true}) {
        if (schedule) model.fit(1, grid, grid, OptimizeGCV, 15, 1e-8, 12, 42);
        if (schedule)
            model.fit(1, std::vector<double> {0.1}, std::vector<double> {0.1}, NoCalibration, 15, 1e-8);
        else
            model.fit(1, fixed, fixed, 15, 1e-8);
        // fixed fits must clear both EDF curves and candidate grids from preceding calibration
        EXPECT_TRUE(
          model.direction_gcv_edfs().empty() && model.loading_gcv_edfs().empty() &&
          model.direction_lambda_grid().rows() == 0 && model.loading_lambda_grid().rows() == 0);
        // fixed fits must clear both candidate selections from preceding calibration
        EXPECT_TRUE(model.direction_selected_indices().empty() && model.loading_selected_indices().empty());
        // fixed overloads must refresh selected penalty matrices to the new component count and values
        EXPECT_TRUE(
          model.direction_lambda().rows() == 1 && model.loading_lambda().rows() == 1 &&
          model.direction_lambda()(0, 0) == 0.1 && model.loading_lambda()(0, 0) == 0.1);
    }
}

/// @brief estimates the predictor mean with a fixed smoothing penalty
Eigen::RowVectorXd smooth_mean(
  const Eigen::Matrix<double, Dynamic, Dynamic>& X,   // predictors with statistical units in rows
  fdapde::internals::fe_ls_elliptic& smoother,        // discretized smoother bound to the observation locations
  double lambda                                       // fixed mean smoothing parameter
) {
    smoother.update_response(X.transpose() * Eigen::VectorXd::Ones(X.rows()) / X.rows());
    smoother.fit(lambda);
    return smoother.fn().transpose();
}

/// @brief estimates the predictor mean with a GCV-selected smoothing penalty
Eigen::RowVectorXd smooth_mean(
  const Eigen::Matrix<double, Dynamic, Dynamic>& X,   // predictors with statistical units in rows
  fdapde::internals::fe_ls_elliptic& smoother,        // discretized smoother bound to the observation locations
  const std::vector<double>& lambda_grid,             // candidate mean smoothing parameters
  int edf_r,                                          // number of random probes for effective degrees of freedom
  int seed                                            // seed for the effective degrees of freedom estimate
) {
    auto gcv = [&](auto lambda) {
        smoother.update_response(X.transpose() * Eigen::VectorXd::Ones(X.rows()) / X.rows());
        smoother.fit(lambda);
        double dor = X.cols() - smoother.edf(lambda, edf_r, seed);
        return (X.cols() / std::pow(dor, 2)) * (smoother.fn() - smoother.response()).squaredNorm();
    };
    GridSearch<1> optimizer;
    Eigen::Matrix<double, 1, 1> lambda = optimizer.optimize(gcv, lambda_grid);
    smoother.update_response(X.transpose() * Eigen::VectorXd::Ones(X.rows()) / X.rows());
    smoother.fit(lambda);
    return smoother.fn().transpose();
}

/// @brief checks coefficient predictions against the independent sequential score reconstruction
template <typename Model>
void check_coefficient_predictions(
  const Model& model,             // fitted regression model
  const Eigen::MatrixXd& X_Psi,   // original centered predictors multiplied by the direction evaluation matrix
  int n_comp                      // number of fitted component prefixes to check
) {
    for (int h = 1; h <= n_comp; ++h) {
        SCOPED_TRACE(h);
        // each prefix must predict the same response as sequential score regression
        EXPECT_TRUE((X_Psi * model.Beta(h)).isApprox(model.fitted(h), 1e-10));
        // both coefficient accessors must expose the same operator
        EXPECT_TRUE(model.B(h).isApprox(model.Beta(h), 1e-12));
    }
    // the cached full fit must agree with predictions from the cached coefficients
    EXPECT_TRUE((X_Psi * model.Beta()).isApprox(model.fitted(), 1e-10));
    // zero selects the full fit in the explicit component overload
    EXPECT_TRUE(model.Beta(0).isApprox(model.Beta(), 1e-12));
    // the default B accessor must expose the cached full coefficient matrix
    EXPECT_TRUE(model.B().isApprox(model.Beta(), 1e-12));
}

/// @brief checks fixed-penalty fits against reference predictions and verifies component penalty schedules
void check_fpls_case(
  const std::string& data_path,   // directory containing input data and reference fit outputs
  double lambda                   // smoothing parameter for centering, directions, and loadings
) {
    std::string mesh_path = "../data/mesh/unit_square_60/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);

    Eigen::Matrix<double, Dynamic, Dynamic> X = read_csv<double>(data_path + "X.csv").as_matrix();
    Eigen::Matrix<double, Dynamic, Dynamic> Y = read_csv<double>(data_path + "Y.csv").as_matrix();

    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);

    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);

    fdapde::internals::fe_ls_elliptic x_centering;
    x_centering.discretize(fe_ls_elliptic(a, F).get());
    x_centering.analyze_data(data, Eigen::VectorXd::Ones(data[0].rows()).asDiagonal());
    Eigen::RowVectorXd X_mean = smooth_mean(X, x_centering, lambda);
    Eigen::Matrix<double, Dynamic, Dynamic> X_centered = X.rowwise() - X_mean;
    Eigen::RowVectorXd Y_mean = Y.colwise().mean();
    Eigen::Matrix<double, Dynamic, Dynamic> Y_centered = Y.rowwise() - Y_mean;

    l1.load_blk("X", X_centered.transpose());

    fPLS m("X", Y_centered, data, fe_ls_elliptic(a, F), fe_ls_elliptic(a, F));
    Eigen::Matrix<double, 1, 1> lambda_vec;
    lambda_vec << lambda;
    m.fit(3, lambda_vec, lambda_vec, 20, 1e-2);

    // compare every coefficient prefix with the independent score-based predictions
    check_coefficient_predictions(m, X_centered * x_centering.Psi(), 3);

    // compare uncentered response predictions with the stored reference using the test tolerance
    EXPECT_TRUE(almost_equal<double>(m.fitted().rowwise() + Y_mean, data_path + "Y_hat.csv"));
    // compare uncentered predictor reconstructions with the stored reference using the test tolerance
    EXPECT_TRUE(almost_equal<double>(m.reconstructed().rowwise() + X_mean, data_path + "X_hat.csv"));
    // compare uncentered response predictions with the stored reference using the test tolerance
    EXPECT_TRUE(almost_equal<double>(m.fitted(3).rowwise() + Y_mean, data_path + "Y_hat.csv"));
    // compare uncentered predictor reconstructions with the stored reference using the test tolerance
    EXPECT_TRUE(almost_equal<double>(m.reconstructed(3).rowwise() + X_mean, data_path + "X_hat.csv"));
    // match output rows to the number of input statistical units
    EXPECT_EQ(m.X_latent_scores().rows(), Y.rows());
    // require one score column for each of the three requested components
    EXPECT_EQ(m.X_latent_scores().cols(), 3);
    // match output rows to the number of input statistical units
    EXPECT_EQ(m.Y_latent_scores().rows(), Y.rows());
    // require one score column for each of the three requested components
    EXPECT_EQ(m.Y_latent_scores().cols(), 3);
    // match output columns to the number of response variables
    EXPECT_EQ(m.fitted(1).cols(), Y.cols());
    // match output rows to the number of input statistical units
    EXPECT_EQ(m.reconstructed(1).rows(), Y.rows());
    // match output columns to the number of response variables
    EXPECT_EQ(m.B(1).cols(), Y.cols());
    // match the iteration diagnostic count to the requested component count
    EXPECT_EQ(m.direction_iterations().size(), 3);
    // match the objective history count to the requested component count
    EXPECT_EQ(m.direction_objective_history().size(), 3);
    // match the monotonicity diagnostic count to the requested component count
    EXPECT_EQ(m.direction_monotone().size(), 3);
    for (int h = 0; h < 3; ++h) {
        // require at least one completed direction update for each component
        EXPECT_GT(m.direction_iterations()[h], 0);
        // match the number of recorded objectives to the completed direction updates
        EXPECT_EQ(m.direction_objective_history()[h].size(), m.direction_iterations()[h]);
        // require the fit diagnostic to report no objective increase beyond tolerance
        EXPECT_TRUE(m.direction_monotone()[h]);
    }

    const std::vector<double> direction_schedule {lambda, lambda / 2, lambda / 4};
    const std::vector<double> loading_schedule {lambda / 4, lambda / 2, lambda};
    m.fit(3, direction_schedule, loading_schedule, NoCalibration, 20, 1e-2);
    for (int h = 0; h < 3; ++h) {
        // each component receives its own direction penalty from the input schedule
        EXPECT_DOUBLE_EQ(m.direction_lambda()(h, 0), direction_schedule[h]);
        // the loading schedule is indexed independently from the direction schedule
        EXPECT_DOUBLE_EQ(m.loading_lambda()(h, 0), loading_schedule[h]);
    }
    // component-specific penalties preserve the coefficient prediction identity
    check_coefficient_predictions(m, X_centered * x_centering.Psi(), 3);
}

/// @brief checks GCV-selected fits against reference predictions and convergence diagnostics
void check_fpls_gcv_case(
  const std::string& data_path   // directory containing input data and reference GCV-fit outputs
) {
    std::string mesh_path = "../data/mesh/unit_square_60/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);

    Eigen::Matrix<double, Dynamic, Dynamic> X = read_csv<double>(data_path + "X.csv").as_matrix();
    Eigen::Matrix<double, Dynamic, Dynamic> Y = read_csv<double>(data_path + "Y.csv").as_matrix();

    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);

    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);

    std::vector<double> lambda_grid(5);
    for (int i = 0; i < 5; ++i) { lambda_grid[i] = std::pow(10, -4.0 + i); }
    int seed = 476813;

    fdapde::internals::fe_ls_elliptic x_centering;
    x_centering.discretize(fe_ls_elliptic(a, F).get());
    x_centering.analyze_data(data, Eigen::VectorXd::Ones(data[0].rows()).asDiagonal());
    Eigen::RowVectorXd X_mean = smooth_mean(X, x_centering, lambda_grid, 1000, seed);
    Eigen::Matrix<double, Dynamic, Dynamic> X_centered = X.rowwise() - X_mean;
    Eigen::RowVectorXd Y_mean = Y.colwise().mean();
    Eigen::Matrix<double, Dynamic, Dynamic> Y_centered = Y.rowwise() - Y_mean;

    l1.load_blk("X", X_centered.transpose());

    fPLS m("X", Y_centered, data, fe_ls_elliptic(a, F), fe_ls_elliptic(a, F));
    m.fit(3, lambda_grid, lambda_grid, OptimizeGCV, 20, 1e-2, 1000, seed);

    // compare every coefficient prefix with the independent score-based predictions
    check_coefficient_predictions(m, X_centered * x_centering.Psi(), 3);

    // compare uncentered response predictions with the stored reference using the test tolerance
    EXPECT_TRUE(almost_equal<double>(m.fitted().rowwise() + Y_mean, data_path + "Y_hat.csv"));
    // compare uncentered predictor reconstructions with the stored reference using the test tolerance
    EXPECT_TRUE(almost_equal<double>(m.reconstructed().rowwise() + X_mean, data_path + "X_hat.csv"));
    // compare uncentered response predictions with the stored reference using the test tolerance
    EXPECT_TRUE(almost_equal<double>(m.fitted(3).rowwise() + Y_mean, data_path + "Y_hat.csv"));
    // compare uncentered predictor reconstructions with the stored reference using the test tolerance
    EXPECT_TRUE(almost_equal<double>(m.reconstructed(3).rowwise() + X_mean, data_path + "X_hat.csv"));
    // match output rows to the number of input statistical units
    EXPECT_EQ(m.X_latent_scores().rows(), Y.rows());
    // require one score column for each of the three requested components
    EXPECT_EQ(m.X_latent_scores().cols(), 3);
    // match output rows to the number of input statistical units
    EXPECT_EQ(m.Y_latent_scores().rows(), Y.rows());
    // require one score column for each of the three requested components
    EXPECT_EQ(m.Y_latent_scores().cols(), 3);
    // match output columns to the number of response variables
    EXPECT_EQ(m.fitted(1).cols(), Y.cols());
    // match output rows to the number of input statistical units
    EXPECT_EQ(m.reconstructed(1).rows(), Y.rows());
    // match output columns to the number of response variables
    EXPECT_EQ(m.B(1).cols(), Y.cols());
    // match the iteration diagnostic count to the requested component count
    EXPECT_EQ(m.direction_iterations().size(), 3);
    // match the objective history count to the requested component count
    EXPECT_EQ(m.direction_objective_history().size(), 3);
    // match the monotonicity diagnostic count to the requested component count
    EXPECT_EQ(m.direction_monotone().size(), 3);
    for (int h = 0; h < 3; ++h) {
        // require at least one completed direction update for each component
        EXPECT_GT(m.direction_iterations()[h], 0);
        // match the number of recorded objectives to the completed direction updates
        EXPECT_EQ(m.direction_objective_history()[h].size(), m.direction_iterations()[h]);
        // require the fit diagnostic to report no objective increase beyond tolerance
        EXPECT_TRUE(m.direction_monotone()[h]);
    }
}

/// @brief checks that mode A and symmetric block fits produce finite outputs with the expected dimensions
void check_restored_modes_smoke(
  const std::string& data_path,   // directory containing predictor and response data
  double lambda                   // fixed smoothing parameter for centering and both component solvers
) {
    std::string mesh_path = "../data/mesh/unit_square_60/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);

    Eigen::Matrix<double, Dynamic, Dynamic> X = read_csv<double>(data_path + "X.csv").as_matrix();
    Eigen::Matrix<double, Dynamic, Dynamic> Y = read_csv<double>(data_path + "Y.csv").as_matrix();

    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);

    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);

    fdapde::internals::fe_ls_elliptic x_centering;
    x_centering.discretize(fe_ls_elliptic(a, F).get());
    x_centering.analyze_data(data, Eigen::VectorXd::Ones(data[0].rows()).asDiagonal());
    Eigen::RowVectorXd X_mean = smooth_mean(X, x_centering, lambda);
    Eigen::Matrix<double, Dynamic, Dynamic> X_centered = X.rowwise() - X_mean;
    Eigen::RowVectorXd Y_mean = Y.colwise().mean();
    Eigen::Matrix<double, Dynamic, Dynamic> Y_centered = Y.rowwise() - Y_mean;

    l1.load_blk("X", X_centered.transpose());

    Eigen::Matrix<double, 1, 1> lambda_vec;
    lambda_vec << lambda;

    fPLS mode_a("X", Y_centered, data, fe_ls_elliptic(a, F), fe_ls_elliptic(a, F), fPLS_A);
    mode_a.fit(3, lambda_vec, lambda_vec, 20, 1e-2);
    // verify that the constructor tag selects the expected deflation mode
    EXPECT_EQ(mode_a.mode(), fPLSMode::ModeA);
    // match output rows to the number of input statistical units
    EXPECT_EQ(mode_a.fitted().rows(), Y.rows());
    // match output columns to the number of response variables
    EXPECT_EQ(mode_a.fitted().cols(), Y.cols());
    // match output rows to the number of input statistical units
    EXPECT_EQ(mode_a.reconstructed().rows(), X.rows());
    // require one score column for each of the three requested components
    EXPECT_EQ(mode_a.Y_latent_scores().cols(), 3);
    // check every entry of the response predictions for NaN or infinity
    EXPECT_TRUE(mode_a.fitted().array().isFinite().all());
    // check every entry of the predictor reconstructions for NaN or infinity
    EXPECT_TRUE(mode_a.reconstructed().array().isFinite().all());

    fPLS mode_sb("X", Y_centered, data, fe_ls_elliptic(a, F), fe_ls_elliptic(a, F), fPLS_SB);
    const std::vector<double> loading_grid {};
    mode_sb.fit(3, lambda_vec, loading_grid, 20, 1e-2);
    // verify that the constructor tag selects the expected deflation mode
    EXPECT_EQ(mode_sb.mode(), fPLSMode::SymmetricBlock);
    // match output rows to the number of input statistical units
    EXPECT_EQ(mode_sb.fitted().rows(), Y.rows());
    // match output columns to the number of response variables
    EXPECT_EQ(mode_sb.fitted().cols(), Y.cols());
    // match output rows to the number of input statistical units
    EXPECT_EQ(mode_sb.reconstructed().rows(), X.rows());
    // require one score column for each of the three requested components
    EXPECT_EQ(mode_sb.Y_latent_scores().cols(), 3);
    // check every entry of the response predictions for NaN or infinity
    EXPECT_TRUE(mode_sb.fitted().array().isFinite().all());
    // check every entry of the predictor reconstructions for NaN or infinity
    EXPECT_TRUE(mode_sb.reconstructed().array().isFinite().all());

    const Eigen::MatrixXd fixed_fitted = mode_sb.fitted();
    const std::vector<double> direction_grid {lambda};
    for (int calibration : {NoCalibration, OptimizeGCV}) {
        mode_sb.fit(3, direction_grid, loading_grid, calibration, 20, 1e-2, 10, 476813);
        // symmetric block fits must not report unused loading penalties for either calibration strategy
        EXPECT_EQ(mode_sb.loading_lambda().size(), 0);
        // reusing directions as loadings requires no loading GCV curve
        EXPECT_TRUE(mode_sb.loading_gcv_values().empty());
        // the common loading update must copy the predictor directions exactly
        EXPECT_TRUE(mode_sb.X_loadings().isApprox(mode_sb.X_space_directions(), 1e-12));
        // the common loading update must copy the response directions exactly
        EXPECT_TRUE(mode_sb.Y_loadings().isApprox(mode_sb.Y_space_directions(), 1e-12));
        // one direction candidate reproduces the fixed fit without any loading penalty
        EXPECT_TRUE(mode_sb.fitted().isApprox(fixed_fitted, 1e-12));
    }
}

/// @brief checks spline regression outputs, coefficient equivalence, and GCV curve reset on fixed-penalty refits
void check_spline_smoke() {
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 1, 21);
    GeoFrame data(T);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);

    const int n_units = 8;
    const int n_locs = data[0].rows();
    Eigen::Matrix<double, Dynamic, Dynamic> X(n_units, n_locs);
    Eigen::Matrix<double, Dynamic, Dynamic> Y(n_units, 2);
    for (int i = 0; i < n_units; ++i) {
        for (int j = 0; j < n_locs; ++j) {
            double x = T.nodes()(j, 0);
            X(i, j) = std::sin((i + 1) * x) + 0.2 * std::cos((j + 1) * x);
        }
        Y(i, 0) = 0.3 * X.row(i).mean() + i * 0.05;
        Y(i, 1) = X(i, n_locs / 2) - X(i, 0);
    }

    Eigen::RowVectorXd X_mean = X.colwise().mean();
    Eigen::Matrix<double, Dynamic, Dynamic> X_centered = X.rowwise() - X_mean;
    Eigen::RowVectorXd Y_mean = Y.colwise().mean();
    Eigen::Matrix<double, Dynamic, Dynamic> Y_centered = Y.rowwise() - Y_mean;
    l1.load_blk("X", X_centered.transpose());

    BsSpace Bh(T, 3);
    TrialFunction f(Bh);
    TestFunction v(Bh);
    auto a = integral(T)(dxx(f) * dxx(v));
    ZeroField<1> u;
    auto F = integral(T)(u * v);

    fPLS m("X", Y_centered, data, bs_ls_elliptic(a, F), bs_ls_elliptic(a, F));
    Eigen::Matrix<double, 1, 1> lambda_vec;
    lambda_vec << 1e-3;
    m.fit(2, lambda_vec, lambda_vec, 20, 1e-6);

    // compare every coefficient prefix with the independent score-based predictions
    check_coefficient_predictions(m, X_centered * internals::point_basis_eval(Bh, T.nodes()), 2);

    // match output rows to the number of input statistical units
    EXPECT_EQ(m.fitted().rows(), Y.rows());
    // match output columns to the number of response variables
    EXPECT_EQ(m.fitted().cols(), Y.cols());
    // check every entry of the response predictions for NaN or infinity
    EXPECT_TRUE(m.fitted().array().isFinite().all());
    // check every entry of the predictor reconstructions for NaN or infinity
    EXPECT_TRUE(m.reconstructed().array().isFinite().all());
    // check every entry of the regression coefficients for NaN or infinity
    EXPECT_TRUE(m.Beta().array().isFinite().all());
    // match the iteration diagnostic count to the requested component count
    EXPECT_EQ(m.direction_iterations().size(), 2);
    // match the objective history count to the requested component count
    EXPECT_EQ(m.direction_objective_history().size(), 2);
    // match the monotonicity diagnostic count to the requested component count
    EXPECT_EQ(m.direction_monotone().size(), 2);

    const std::vector<double> lambda_grid {1e-3};
    for (bool use_schedule : {false, true}) {
        m.fit(2, lambda_grid, lambda_grid, OptimizeGCV, 20, 1e-6, 10, 476813);
        // both components must retain a direction curve after GCV calibration
        ASSERT_EQ(m.direction_gcv_values().size(), 2);
        // both components must retain a loading curve after GCV calibration
        ASSERT_EQ(m.loading_gcv_values().size(), 2);
        for (int h = 0; h < 2; ++h) {
            // a single candidate produces exactly one direction GCV value per component
            EXPECT_EQ(m.direction_gcv_values()[h].size(), 1);
            // a single candidate produces exactly one loading GCV value per component
            EXPECT_EQ(m.loading_gcv_values()[h].size(), 1);
        }
        if (use_schedule) {
            m.fit(2, lambda_grid, lambda_grid, NoCalibration, 20, 1e-6);
        } else {
            m.fit(2, lambda_vec, lambda_vec, 20, 1e-6);
        }
        // either fixed-penalty overload must discard direction curves from the preceding GCV fit
        EXPECT_TRUE(m.direction_gcv_values().empty());
        // either fixed-penalty overload must discard loading curves from the preceding GCV fit
        EXPECT_TRUE(m.loading_gcv_values().empty());
    }
}

/// @brief checks that a zero predictor block fails with an exception during direction estimation
void check_singular_direction_failure() {
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 1, 11);
    GeoFrame data(T);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);

    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(6, data[0].rows());
    Eigen::MatrixXd Y = Eigen::MatrixXd::Random(6, 1);
    l1.load_blk("X", X.transpose());

    BsSpace Bh(T, 3);
    TrialFunction f(Bh);
    TestFunction v(Bh);
    auto a = integral(T)(dxx(f) * dxx(v));
    ZeroField<1> u;
    auto F = integral(T)(u * v);

    fPLS m("X", Y, data, bs_ls_elliptic(a, F), bs_ls_elliptic(a, F));
    Eigen::Matrix<double, 1, 1> lambda;
    lambda << 1e-9;
    // require direction estimation on zero predictors to throw a runtime error
    EXPECT_THROW(m.fit(1, lambda, lambda), std::runtime_error);
}

}   // namespace

// check fixed-penalty predictions and component-specific penalty schedules on the reference data
TEST(fpls, test_01) { check_fpls_case("../data/models/fpls/2D_test1/", 10.0); }

// check GCV-selected predictions against the reference data and verify component diagnostics
TEST(fpls, test_02) { check_fpls_gcv_case("../data/models/fpls/2D_test2/"); }

// check mode A and symmetric block selection, output dimensions, and finite reconstructions
TEST(fpls, restored_modes_smoke) { check_restored_modes_smoke("../data/models/fpls/2D_test1/", 10.0); }

// check spline regression and ensure both fixed-penalty overloads clear curves from preceding GCV fits
TEST(fpls, spline_smoke) { check_spline_smoke(); }

// check that zero predictors trigger a runtime error instead of non-finite directions
TEST(fpls, singular_direction_failure) { check_singular_direction_failure(); }

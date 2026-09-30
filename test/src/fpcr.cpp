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
using fdapde::test::check_coefficient_predictions;
using fdapde::test::smooth_mean;

namespace {

// reference outputs use the same MatrixMarket format and tolerance as the fPLS regression cases
/// @brief checks frozen fPCR predictions, reconstructions, coefficients and component diagnostics
template <typename Model>
void check_fpcr_outputs(
  const Model& model,                  // fitted functional principal component regression model
  const Eigen::MatrixXd& X_centered,   // centered training predictors before component deflation
  const Eigen::MatrixXd& X_Psi,        // original centered predictors evaluated in the functional basis
  const Eigen::RowVectorXd& X_mean,    // smooth predictor mean restored in the reference reconstruction
  const Eigen::RowVectorXd& Y_mean,    // response means restored in the reference predictions
  const std::string& reference_path    // directory containing the frozen fPCR regression outputs
) {
    // every prefix must map original predictors to the same response as its separately fitted score regression
    check_coefficient_predictions(model, X_Psi, 3);

    // compare uncentered response predictions with the frozen current-implementation reference
    EXPECT_TRUE(almost_equal<double>(model.fitted().rowwise() + Y_mean, reference_path + "Y_hat.csv"));
    // compare uncentered predictor reconstructions with the frozen current-implementation reference
    EXPECT_TRUE(almost_equal<double>(model.reconstructed().rowwise() + X_mean, reference_path + "X_hat.csv"));
    // the explicit full prefix must reproduce the same frozen response predictions
    EXPECT_TRUE(almost_equal<double>(model.fitted(3).rowwise() + Y_mean, reference_path + "Y_hat.csv"));
    // the explicit full prefix must reproduce the same frozen predictor reconstructions
    EXPECT_TRUE(almost_equal<double>(model.reconstructed(3).rowwise() + X_mean, reference_path + "X_hat.csv"));
    // original-X regression coefficients must match the frozen coefficient matrix for both responses
    EXPECT_TRUE(almost_equal<double>(model.Beta(), reference_path + "B_hat.csv"));
    // predicting the training observations must agree with the cached separately fitted score regression
    EXPECT_TRUE(model.predict(X_centered).isApprox(model.fitted(), 1e-10));
    // score rows must match the number of input statistical units
    EXPECT_EQ(model.X_latent_scores().rows(), X_centered.rows());
    // scores must retain the three requested principal components
    EXPECT_EQ(model.X_latent_scores().cols(), 3);
    // response loadings must retain both response columns of the multivariate regression
    EXPECT_EQ(model.Y_loadings().rows(), Y_mean.size());
    // response loadings must retain one column for each fitted component
    EXPECT_EQ(model.Y_loadings().cols(), 3);
    // a response prediction prefix must retain every input response variable
    EXPECT_EQ(model.fitted(1).cols(), Y_mean.size());
    // predictor reconstruction rows must match the number of input statistical units
    EXPECT_EQ(model.reconstructed(1).rows(), X_centered.rows());
    // a coefficient prefix must retain every input response variable
    EXPECT_EQ(model.B(1).cols(), Y_mean.size());
    // iteration diagnostics must have one entry for each of the three fitted components
    ASSERT_EQ(model.iterations().size(), 3);
    // objective diagnostics must have one trace for each of the three fitted components
    ASSERT_EQ(model.objective_history().size(), 3);
    // monotonicity diagnostics must have one flag for each of the three fitted components
    ASSERT_EQ(model.monotone().size(), 3);
    for (int h = 0; h < 3; ++h) {
        SCOPED_TRACE(h);
        // each retained principal component must complete at least one power update
        EXPECT_GT(model.iterations()[h], 0);
        // each recorded objective must correspond to one completed power update
        EXPECT_EQ(model.objective_history()[h].size(), model.iterations()[h]);
        // every retained objective trace must avoid increases beyond the fit tolerance
        EXPECT_TRUE(model.monotone()[h]);
    }
}

/// @brief checks fixed-penalty fPCR on the same centered data and variational penalty as the fixed fPLS case
void check_fpcr_case(
  const std::string& data_path,        // directory containing the existing fPLS predictor and response inputs
  const std::string& reference_path,   // directory containing the frozen fixed-penalty fPCR outputs
  double lambda                        // shared smoothing parameter for centering and principal components
) {
    std::string mesh_path = "../data/mesh/unit_square_60/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);

    Eigen::MatrixXd X = read_csv<double>(data_path + "X.csv").as_matrix();
    Eigen::MatrixXd Y = read_csv<double>(data_path + "Y.csv").as_matrix();

    GeoFrame data(D);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);

    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);

    // apply the same smooth-X and ordinary-Y centering used by the fixed fPLS reference case
    internals::fe_ls_elliptic x_centering;
    x_centering.discretize(fe_ls_elliptic(a, F).get());
    x_centering.analyze_data(data, Eigen::VectorXd::Ones(data[0].rows()).asDiagonal());
    Eigen::RowVectorXd X_mean = smooth_mean(X, x_centering, lambda);
    Eigen::MatrixXd X_centered = X.rowwise() - X_mean;
    Eigen::RowVectorXd Y_mean = Y.colwise().mean();
    Eigen::MatrixXd Y_centered = Y.rowwise() - Y_mean;

    l1.load_blk("X", X_centered.transpose());

    fPCR m("X", Y_centered, data, fe_ls_elliptic(a, F));
    m.fit(
      3,                              // number of retained principal components
      std::vector<double> {lambda},   // fixed smoothing parameter shared by all components
      ComputeXactSVD,                 // exact SVD initialization without calibration
      20,                             // maximum power updates per component, as in the fPLS case
      1e-2                            // objective stopping tolerance, as in the fPLS case
    );

    // check the frozen fixed-penalty outputs and coefficient prediction identities
    check_fpcr_outputs(m, X_centered, X_centered * x_centering.Psi(), X_mean, Y_mean, reference_path);
    for (int h = 0; h < 3; ++h) {
        // fixed calibration must record the supplied shared penalty for every component
        EXPECT_DOUBLE_EQ(m.lambda()(h, 0), lambda);
    }
    // fixed calibration must not report candidate grids or GCV/EDF evaluations
    EXPECT_TRUE(m.lambda_grid().rows() == 0 && m.gcv_values().empty() && m.gcv_edf().empty());
}

/// @brief checks GCV-selected fPCR on the same centered data, grid and EDF controls as the GCV fPLS case
void check_fpcr_gcv_case(
  const std::string& data_path,       // directory containing the existing fPLS predictor and response inputs
  const std::string& reference_path   // directory containing the frozen GCV-selected fPCR outputs
) {
    std::string mesh_path = "../data/mesh/unit_square_60/";
    Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);

    Eigen::MatrixXd X = read_csv<double>(data_path + "X.csv").as_matrix();
    Eigen::MatrixXd Y = read_csv<double>(data_path + "Y.csv").as_matrix();

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
    const int seed = 476813;

    // select the smooth-X mean with the same GCV grid and stochastic EDF controls as the fPLS case
    internals::fe_ls_elliptic x_centering;
    x_centering.discretize(fe_ls_elliptic(a, F).get());
    x_centering.analyze_data(data, Eigen::VectorXd::Ones(data[0].rows()).asDiagonal());
    Eigen::RowVectorXd X_mean = smooth_mean(X, x_centering, lambda_grid, 1000, seed);
    Eigen::MatrixXd X_centered = X.rowwise() - X_mean;
    Eigen::RowVectorXd Y_mean = Y.colwise().mean();
    Eigen::MatrixXd Y_centered = Y.rowwise() - Y_mean;

    l1.load_blk("X", X_centered.transpose());

    fPCR m("X", Y_centered, data, fe_ls_elliptic(a, F));
    m.fit(
      3,                              // number of retained principal components
      lambda_grid,                    // the same five candidate lambdas as the GCV fPLS case
      ComputeXactSVD | OptimizeGCV,   // exact SVD initialization and componentwise GCV selection
      20,                             // maximum power updates per candidate and retained component
      1e-2,                           // objective stopping tolerance, as in the fPLS case
      1000,                           // Hutchinson probes, as in the GCV fPLS case
      seed                            // reproducible EDF seed, as in the GCV fPLS case
    );

    // check the frozen GCV-selected outputs and coefficient prediction identities
    check_fpcr_outputs(m, X_centered, X_centered * x_centering.Psi(), X_mean, Y_mean, reference_path);
    // every retained component must expose its actual calibration curve
    ASSERT_EQ(m.gcv_values().size(), 3);
    // every retained component must expose the EDF estimates used by its calibration curve
    ASSERT_EQ(m.gcv_edf().size(), 3);
    // every retained component must expose the selected candidate position
    ASSERT_EQ(m.selected_indices().size(), 3);
    for (int h = 0; h < 3; ++h) {
        SCOPED_TRACE(h);
        // each GCV curve must contain the five evaluated candidates
        ASSERT_EQ(m.gcv_values()[h].size(), lambda_grid.size());
        // each EDF curve must retain the same candidate order as the GCV curve
        EXPECT_EQ(m.gcv_edf()[h].size(), lambda_grid.size());
        const auto& curve = m.gcv_values()[h];
        const int selected = std::min_element(curve.begin(), curve.end()) - curve.begin();
        // the saved selection must identify the first actual minimizing candidate
        EXPECT_EQ(m.selected_indices()[h], selected);
        // the component's retained penalty must equal the selected grid candidate
        EXPECT_DOUBLE_EQ(m.lambda()(h, 0), lambda_grid[selected]);
    }
}

}   // namespace

// check fixed-penalty fPCR outputs against frozen references using the same multivariate inputs as fpls.test_01
TEST(fpcr, test_01) { check_fpcr_case("../data/models/fpls/2D_test1/", "../data/models/fpcr/2D_test1/", 10.0); }

// check GCV-selected fPCR outputs against frozen references using the same multivariate inputs as fpls.test_02
TEST(fpcr, test_02) { check_fpcr_gcv_case("../data/models/fpls/2D_test2/", "../data/models/fpcr/2D_test2/"); }

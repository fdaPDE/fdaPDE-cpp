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

#include <fdaPDE/models.h>   // fdaPDE
#include <gtest/gtest.h>     // testing framework

#include <unsupported/Eigen/SparseExtra>

namespace fdapde {
namespace test {

  [[maybe_unused]] constexpr double testing_double_tolerance = 1e-7;

  // floating point comparison utilities
  template <typename Scalar>
  bool almost_equal(
    const Eigen::Matrix<Scalar, Dynamic, Dynamic>& op1,
    const Eigen::Matrix<Scalar, Dynamic, Dynamic>& op2,
    double epsilon) {
    return (op1 - op2).template lpNorm<Eigen::Infinity>() < epsilon ||
           (op1 - op2).template lpNorm<Eigen::Infinity>() <
             (std::max(op1.template lpNorm<Eigen::Infinity>(), op2.template lpNorm<Eigen::Infinity>()) * epsilon);
  }

  template <typename Scalar>
  bool almost_equal(
    const Eigen::Matrix<Scalar, Dynamic, Dynamic>& op1,
    const Eigen::Matrix<Scalar, Dynamic, Dynamic>& op2) {
    return almost_equal(op1, op2, testing_double_tolerance);
  }

  template <typename Scalar>
  bool almost_equal(const Eigen::SparseMatrix<Scalar>& op1, std::string op2) {
    Eigen::SparseMatrix<Scalar> mem_buff;
    Eigen::loadMarket(mem_buff, op2);
    return almost_equal(op1, mem_buff);
  }

  template <typename Scalar>
  bool almost_equal(const Eigen::Matrix<Scalar, Dynamic, Dynamic>& op1, std::string op2) {
    Eigen::SparseMatrix<Scalar> mem_buff;
    Eigen::loadMarket(mem_buff, op2);
    return almost_equal(op1, Eigen::Matrix<Scalar, Dynamic, Dynamic>(mem_buff));
  }

  template <typename Scalar>
  bool almost_equal(const std::vector<Scalar>& op1, std::string op2) {
    Eigen::SparseMatrix<Scalar> mem_buff;
    Eigen::Matrix<double, Dynamic, Dynamic> values(op1.size(), 1);

    Eigen::loadMarket(mem_buff, op2);
    for (int i = 0; i < op1.size(); ++i) {
      values(i, 0) = op1[i];
    }
    return almost_equal(values, Eigen::Matrix<Scalar, Dynamic, Dynamic>(mem_buff));
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

}   // namespace test
}   // namespace fdapde


// #include "src/sr.cpp"
#include "src/fpls.cpp"
#include "src/fpcr.cpp"
#include "src/fpca.cpp"
#include "src/edf.cpp"
#include "src/formula.cpp"
// #include "src/gsr.cpp"
// #include "src/qsr.cpp"
// #include "src/de.cpp"

int main(int argc, char **argv){
  // start testing
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();

  return 0;
}

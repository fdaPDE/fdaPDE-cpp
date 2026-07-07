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

#ifndef __FDAPDE_RGCCA_TYPES_H__
#define __FDAPDE_RGCCA_TYPES_H__

#include "../header_check.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <functional>
#include <limits>
#include <mutex>
#include <numeric>
#include <optional>
#include <ostream>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "validation.h"

namespace fdapde {
namespace rgcca {

using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using BoolMatrix = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;
using IndexVector = Eigen::Vector<int, Eigen::Dynamic>;
using SparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;

// Public RGCCA options shared by blocks, model fitting, and bootstrap
enum class InitStrategy { None, SVD, Uniform, WarmStart };
enum class DesignMode { Empty, Custom, FullyConnected };
enum class LambdaSelection { Manual, Automatic };
enum class Mode { CorMax, Regularized, CovMax };
enum class Deflation { None, Scores };
enum class WeightSignConstraint { None, NonNegative };
enum class ResamplingStrategy { Ordinary, Stationary };
enum class InactiveBlockSignalAction { None, KeptInactive, Deactivated, KeptActive, Reactivated };

inline const char* to_string(InitStrategy x) {
    switch (x) {
    case InitStrategy::None:      return "None";
    case InitStrategy::SVD:       return "SVD";
    case InitStrategy::Uniform:   return "Uniform";
    case InitStrategy::WarmStart: return "WarmStart";
    }
    return "Unknown";
}

inline const char* to_string(LambdaSelection x) {
    switch (x) {
    case LambdaSelection::Manual:    return "Manual";
    case LambdaSelection::Automatic: return "Automatic";
    }
    return "Unknown";
}

inline const char* to_string(Mode x) {
    switch (x) {
    case Mode::CorMax:      return "CorMax";
    case Mode::Regularized: return "Regularized";
    case Mode::CovMax:      return "CovMax";
    }
    return "Unknown";
}

inline const char* to_string(WeightSignConstraint x) {
    switch (x) {
    case WeightSignConstraint::None:        return "None";
    case WeightSignConstraint::NonNegative: return "NonNegative";
    }
    return "Unknown";
}

inline const char* to_string(Deflation x) {
    switch (x) {
    case Deflation::None:   return "None";
    case Deflation::Scores: return "Scores";
    }
    return "Unknown";
}

inline const char* to_string(ResamplingStrategy x) {
    switch (x) {
    case ResamplingStrategy::Ordinary:   return "Ordinary";
    case ResamplingStrategy::Stationary: return "Stationary";
    }
    return "Unknown";
}

inline const char* to_string(InactiveBlockSignalAction x) {
    switch (x) {
    case InactiveBlockSignalAction::None:         return "none";
    case InactiveBlockSignalAction::KeptInactive: return "kept inactive";
    case InactiveBlockSignalAction::Deactivated:  return "deactivated";
    case InactiveBlockSignalAction::KeptActive:   return "kept active";
    case InactiveBlockSignalAction::Reactivated:  return "reactivated";
    }
    return "unknown";
}

struct Scheme {
    std::function<double(double)> g;   // g(t)
    std::function<double(double)> w;   // w(t)
    double phi = 1.0;
    const char* name = "custom";
    bool sign_invariant = false;

    static Scheme Horst() {
        return {[](double t) { return t; }, [](double) { return 1.0; }, 1.0, "Horst", false};
    }
    static Scheme Centroid() {
        return {
            [](double t) { return std::abs(t); }, [](double t) { return t >= 0 ? 1.0 : -1.0; }, 1.0, "Centroid", true};
    }
    static Scheme Factorial() {
        return {[](double t) { return t * t; }, [](double t) { return t; }, 2.0, "Factorial", true};
    }
};

struct Options {
    int max_iter;
    double tol;
    bool cache_covariances;
    bool bias;
    InitStrategy init_strategy;
    LambdaSelection lambda_selection_weights;
    LambdaSelection lambda_selection_components;
    bool component_significance;
    bool block_importance;
    bool inactive_block_signal_test;
    bool block_deactivation;
    bool connection_deactivation;
    Mode mode;
    WeightSignConstraint weight_sign_constraint;
    Deflation deflation_mode;
    Scheme scheme;

    explicit Options(
      const int max_iter_ = 1000, const double tol_ = 1e-8, const bool bias_ = true,
      const InitStrategy init_strategy_ = InitStrategy::SVD, const Mode mode_ = Mode::CovMax,
      const WeightSignConstraint weight_sign_constraint_ = WeightSignConstraint::None,
      const LambdaSelection lambda_selection_weights_ = LambdaSelection::Manual,
      const LambdaSelection lambda_selection_components_ = LambdaSelection::Automatic,
      const bool component_significance_ = false,
      const Deflation deflation_mode_ = Deflation::Scores, const Scheme& scheme_ = Scheme::Factorial(),
      const bool cache_ = true,
      const bool block_deactivation_ = false, const bool connection_deactivation_ = false,
      const bool block_importance_ = false,
      const bool inactive_block_signal_test_ = false) :
        max_iter(max_iter_),
        tol(tol_),
        cache_covariances(cache_),
        bias(bias_),
        init_strategy(init_strategy_),
        lambda_selection_weights(lambda_selection_weights_),
        lambda_selection_components(lambda_selection_components_),
        component_significance(component_significance_),
        block_importance(block_importance_),
        inactive_block_signal_test(inactive_block_signal_test_),
        block_deactivation(block_deactivation_),
        connection_deactivation(connection_deactivation_),
        mode(mode_),
        weight_sign_constraint(weight_sign_constraint_),
        deflation_mode(deflation_mode_),
        scheme(scheme_) { }
};

struct BootstrapConfig {
    unsigned seed = 12345;
    int max_threads = 12;

    int B_min = 500;
    int B_max = 1000;
    int check_every = 50;
    int check_every_block_deactivation = 1;
    int check_every_connection_deactivation = 100;
    int fit_max_iter = -1; // negative means use Options::max_iter

    bool adaptive = true;
    double adaptive_tol = 1e-3;
    int stable_checks_required = 3;

    double active_block_tol = 1e-8;

    double active_connection_sign_stability = 0.95;
    double active_connection_min_abs_corr = 0.05;

    double ci_level = 0.95;
    int patience = 1;

    ResamplingStrategy resampling_strategy = ResamplingStrategy::Ordinary;
    double stationary_block_length = 10.0;

    int component_significance_resamples = 100;
    double component_significance_alpha = 0.05;
    int block_importance_resamples = 100;
    double block_importance_alpha = 0.05;

    int inactive_block_signal_resamples = 100;
    double inactive_block_signal_alpha = 0.05;
};

std::ostream& operator<<(std::ostream& os, const Options& opt);
std::ostream& operator<<(std::ostream& os, const BootstrapConfig& config);

struct Result {
    using Matrix = fdapde::rgcca::Matrix;
    using BoolMatrix = fdapde::rgcca::BoolMatrix;

    int h = 0;
    int n_blocks = 0;
    std::vector<double> obj_history;
    bool monotone = true;
    bool cancelled = false;
    int iters = 0;
    BoolMatrix C;
    Matrix covariance_matrix;
    Matrix correlation_matrix;
    std::vector<double> tau_values;
    std::vector<double> lambda_components_values;
    std::vector<double> lambda_weights_values;
    std::vector<bool> active_blocks;
    double rho_tot = std::numeric_limits<double>::quiet_NaN();
    double rho_tot_p_value = std::numeric_limits<double>::quiet_NaN();
    int rho_tot_bootstrap_count = 0;
    bool component_significant = true;
    std::vector<double> block_importance;
    std::vector<double> block_importance_p_values;
    std::vector<bool> block_importance_significant;
    int block_importance_bootstrap_count = 0;
    std::vector<InactiveBlockSignalAction> inactive_block_signal_actions;

    explicit Result(const int n_blocks_) : n_blocks(n_blocks_), C(n_blocks_, n_blocks_), covariance_matrix(n_blocks_, n_blocks_),
    tau_values(n_blocks_), lambda_components_values(n_blocks_), lambda_weights_values(n_blocks_), active_blocks(n_blocks_),
    block_importance(n_blocks_, std::numeric_limits<double>::quiet_NaN()),
    block_importance_p_values(n_blocks_, std::numeric_limits<double>::quiet_NaN()),
    block_importance_significant(n_blocks_, false),
    inactive_block_signal_actions(n_blocks_, InactiveBlockSignalAction::None) {}
};

struct BootstrapResult {
    using Matrix = fdapde::rgcca::Matrix;
    using Vector = fdapde::rgcca::Vector;

    int h = 0;
    int B = 0;

    std::vector<double> lambda_grid;
    std::vector<double> criterion;

    double lambda_opt = std::numeric_limits<double>::quiet_NaN();
    int lambda_opt_index = -1;
    double ci_level = std::numeric_limits<double>::quiet_NaN();

    std::vector<std::string> block_names;
    std::vector<bool> active_blocks;

    // [lambda][block] -> vector/matrix
    std::vector<std::vector<Vector>> w_fit_by_lambda;
    std::vector<std::vector<Matrix>> w_boot_by_lambda;
    std::vector<std::vector<Vector>> w_min_by_lambda;
    std::vector<int> B_used_by_lambda;

    // [lambda] -> matrix
    std::vector<Matrix> corr_boot_by_lambda;
    std::vector<Matrix> corr_min_by_lambda;

    // [lambda] -> n_blocks x n_blocks confidence intervals
    std::vector<Matrix> corr_ci_low_by_lambda;
    std::vector<Matrix> corr_ci_high_by_lambda;

    BootstrapResult() = default;

    BootstrapResult(
        const int h_,
        const int B_,
        const std::vector<double>& lambda_grid_,
        const std::vector<std::string>& block_names_,
        const std::vector<int>& block_dims_,
        const double ci_level_
    ) :
        h(h_),
        B(B_),
        lambda_grid(lambda_grid_),
        criterion(lambda_grid_.size(), -std::numeric_limits<double>::infinity()),
        ci_level(ci_level_),
        block_names(block_names_)
    {
        const std::size_t n_lambda = lambda_grid.size();
        const std::size_t n_blocks = block_dims_.size();

        w_fit_by_lambda.resize(n_lambda);
        w_boot_by_lambda.resize(n_lambda);
        w_min_by_lambda.resize(n_lambda);
        B_used_by_lambda.resize(n_lambda);

        active_blocks.resize(n_blocks);

        corr_min_by_lambda.resize(n_lambda);
        corr_boot_by_lambda.resize(n_lambda);
        corr_ci_low_by_lambda.resize(n_lambda);
        corr_ci_high_by_lambda.resize(n_lambda);

        for (std::size_t i = 0; i < n_lambda; ++i) {
            corr_min_by_lambda[i].setZero(n_blocks, n_blocks);
            corr_ci_low_by_lambda[i].setZero(n_blocks, n_blocks);
            corr_ci_high_by_lambda[i].setZero(n_blocks, n_blocks);
        }
    }
};

std::ostream& operator<<(std::ostream& os, const Result& r);
std::ostream& operator<<(std::ostream& os, const std::vector<Result>& results);

} // namespace rgcca
} // namespace fdapde

#endif // __FDAPDE_RGCCA_TYPES_H__

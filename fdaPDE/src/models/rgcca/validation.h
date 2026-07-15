// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>

#ifndef __FDAPDE_RGCCA_VALIDATION_UTILS_H__
#define __FDAPDE_RGCCA_VALIDATION_UTILS_H__

namespace fdapde {
namespace rgcca {
namespace internals {

// validates a positive finite regularization parameter
inline void validate_positive_regularization_lambda(const double lambda, const char* name) {
    if (!(lambda > 0.0) || !std::isfinite(lambda)) {
        throw std::invalid_argument(std::string("RGCCA: ") + name + " must be finite and positive");
    }
}

} // namespace internals
} // namespace rgcca
} // namespace fdapde

#endif // __FDAPDE_RGCCA_VALIDATION_UTILS_H__

#ifdef __FDAPDE_RGCCA_DEFINE_MODEL_VALIDATION__
#ifndef __FDAPDE_RGCCA_MODEL_VALIDATION_H__
#define __FDAPDE_RGCCA_MODEL_VALIDATION_H__

namespace fdapde {

// validates the public fit entry point and enabled optional workflows
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::validate_fit_() const {
    if (n_blocks() < 2) throw std::runtime_error("RGCCA: need ≥ 2 blocks");

    // bootstrap-backed workflows share the same sampling restriction
    const bool run_model_selection = bootstrap_model_selection_requested_();
    if (run_model_selection || opt_.component_significance || opt_.block_importance)
        validate_bootstrap_support_();
    if (run_model_selection)
        validate_bootstrap_config_();
    if (opt_.component_significance)
        validate_component_significance_config_();
    if (opt_.block_importance)
        validate_block_importance_config_();
    if (opt_.lambda_selection_weights == rgcca::LambdaSelection::Automatic)
        validate_lambda_grid_weights_();
}

// rejects bootstrap workflows for sampling strategies that cannot resample rows
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::validate_bootstrap_support_() const {
    if constexpr (std::same_as<SamplingStrategy, rgcca::TimeDependentSampling>) {
        throw std::runtime_error(
            "RGCCA: bootstrap is not supported "
            "for rgcca::TimeDependentSampling"
        );
    }
}

// validates adaptive bootstrap and model-selection configuration
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::validate_bootstrap_config_() const {
    // execution and sample budget
    if (bootstrap_config_.max_threads <= 0)
        throw std::invalid_argument("RGCCA: bootstrap max_threads must be positive");
    if (bootstrap_config_.B_min <= 0)
        throw std::invalid_argument("RGCCA: bootstrap B_min must be positive");
    if (bootstrap_config_.B_max <= 0)
        throw std::invalid_argument("RGCCA: bootstrap B_max must be positive");
    if (bootstrap_config_.B_min > bootstrap_config_.B_max)
        throw std::invalid_argument("RGCCA: bootstrap B_max must be greater than B_min");
    if (bootstrap_config_.check_every <= 0)
        throw std::invalid_argument("RGCCA: bootstrap check_every must be positive");
    if (bootstrap_config_.check_every_block_deactivation <= 0)
        throw std::invalid_argument("RGCCA: bootstrap check_every_block_deactivation must be positive");
    if (bootstrap_config_.check_every_connection_deactivation <= 0)
        throw std::invalid_argument("RGCCA: bootstrap check_every_connection_deactivation must be positive");
    if (bootstrap_config_.fit_max_iter == 0 || bootstrap_config_.fit_max_iter < -1)
        throw std::invalid_argument("RGCCA: bootstrap fit_max_iter must be positive or -1");

    // adaptive stopping and deactivation thresholds
    if (bootstrap_config_.stable_checks_required <= 0)
        throw std::invalid_argument("RGCCA: bootstrap stable_checks_required must be positive");
    if (!(bootstrap_config_.adaptive_tol >= 0.0) || !std::isfinite(bootstrap_config_.adaptive_tol))
        throw std::invalid_argument("RGCCA: bootstrap adaptive_tol must be finite and nonnegative");
    if (!(bootstrap_config_.active_block_tol >= 0.0) || !std::isfinite(bootstrap_config_.active_block_tol))
        throw std::invalid_argument("RGCCA: bootstrap active_block_tol must be finite and nonnegative");
    if (
        !(bootstrap_config_.active_connection_sign_stability >= 0.0) ||
        bootstrap_config_.active_connection_sign_stability > 1.0 ||
        !std::isfinite(bootstrap_config_.active_connection_sign_stability)
    ) {
        throw std::invalid_argument(
            "RGCCA: bootstrap active_connection_sign_stability must be finite and in [0, 1]"
        );
    }
    if (
        !(bootstrap_config_.active_connection_min_abs_corr >= 0.0) ||
        !std::isfinite(bootstrap_config_.active_connection_min_abs_corr)
    ) {
        throw std::invalid_argument(
            "RGCCA: bootstrap active_connection_min_abs_corr must be finite and nonnegative"
        );
    }
    // confidence intervals and lambda-level stopping
    if (
        !(bootstrap_config_.ci_level > 0.0) ||
        bootstrap_config_.ci_level >= 1.0 ||
        !std::isfinite(bootstrap_config_.ci_level)
    ) {
        throw std::invalid_argument("RGCCA: bootstrap ci_level must be finite and in (0, 1)");
    }
    if (bootstrap_config_.patience <= 0)
        throw std::invalid_argument("RGCCA: bootstrap patience must be positive");
    if (
        bootstrap_config_.resampling_strategy == rgcca::ResamplingStrategy::Stationary &&
        (!(bootstrap_config_.stationary_block_length > 0.0) ||
         !std::isfinite(bootstrap_config_.stationary_block_length))
    ) {
        throw std::invalid_argument("RGCCA: bootstrap stationary_block_length must be finite and positive");
    }
}

// validates the permutation bootstrap used for component significance
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::validate_component_significance_config_() const {
    if (bootstrap_config_.max_threads <= 0)
        throw std::invalid_argument("RGCCA: component significance max_threads must be positive");
    if (bootstrap_config_.component_significance_resamples <= 0)
        throw std::invalid_argument("RGCCA: component significance resamples must be positive");
    if (bootstrap_config_.fit_max_iter == 0 || bootstrap_config_.fit_max_iter < -1)
        throw std::invalid_argument("RGCCA: component significance fit_max_iter must be positive or -1");
    if (
        !(bootstrap_config_.component_significance_alpha > 0.0) ||
        bootstrap_config_.component_significance_alpha >= 1.0 ||
        !std::isfinite(bootstrap_config_.component_significance_alpha)
    ) {
        throw std::invalid_argument(
            "RGCCA: component significance alpha must be finite and in (0, 1)"
        );
    }
    if (
        bootstrap_config_.resampling_strategy == rgcca::ResamplingStrategy::Stationary &&
        (!(bootstrap_config_.stationary_block_length > 0.0) ||
         !std::isfinite(bootstrap_config_.stationary_block_length))
    ) {
        throw std::invalid_argument("RGCCA: component significance stationary_block_length must be finite and positive");
    }
}

// validates the bootstrap used for block importance
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::validate_block_importance_config_() const {
    if (bootstrap_config_.max_threads <= 0)
        throw std::invalid_argument("RGCCA: block importance max_threads must be positive");
    if (bootstrap_config_.block_importance_resamples <= 0)
        throw std::invalid_argument("RGCCA: block importance resamples must be positive");
    if (
        !(bootstrap_config_.block_importance_alpha > 0.0) ||
        bootstrap_config_.block_importance_alpha >= 1.0 ||
        !std::isfinite(bootstrap_config_.block_importance_alpha)
    ) {
        throw std::invalid_argument(
            "RGCCA: block importance alpha must be finite and in (0, 1)"
        );
    }
    if (
        bootstrap_config_.resampling_strategy == rgcca::ResamplingStrategy::Stationary &&
        (!(bootstrap_config_.stationary_block_length > 0.0) ||
         !std::isfinite(bootstrap_config_.stationary_block_length))
    ) {
        throw std::invalid_argument("RGCCA: block importance stationary_block_length must be finite and positive");
    }
}

// validates a block index against the model block list
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::validate_index_(const int j) const {
    if (j < 0 || j >= static_cast<int>(blocks_.size()))
        throw std::out_of_range("block index");
}

// validates an explicit per-component weight-lambda grid
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::validate_lambda_grid_weights_(
    const std::vector<std::vector<double>>& lambda_grid
) const {
    if (static_cast<int>(lambda_grid.size()) != n_comp_)
        throw std::invalid_argument("RGCCA: weight lambda grid must have size 1 or n_comp");

    for (const auto& grid : lambda_grid) {
        if (grid.empty()) {
            throw std::invalid_argument("RGCCA: weight lambda grid contains an empty component grid");
        }

        for (std::size_t i = 0; i < grid.size(); ++i) {
            const double lambda = grid[i];

            if (!(lambda > 0.0) || !std::isfinite(lambda)) {
                throw std::invalid_argument("RGCCA: weight lambda grid values must be finite and positive");
            }

            if (i > 0 && grid[i] < grid[i - 1]) {
                throw std::invalid_argument("RGCCA: weight lambda grid must be sorted in nondecreasing order");
            }
        }
    }
}

// validates the stored weight-lambda grid before automatic selection
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::validate_lambda_grid_weights_() const {
    if (static_cast<int>(lambda_grid_weights_.size()) != n_comp_) {
        throw std::runtime_error(
            "RGCCA: automatic weight lambda selection requires set_lambda_grid_weights(...) "
            "with one grid or one grid per component"
        );
    }

    validate_lambda_grid_weights_(lambda_grid_weights_);
}

// validates bootstrap weight-ci indices, confidence level and evaluation matrix
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::validate_bootstrap_weights_ci_(
    const typename RGCCA<SamplingStrategy>::BootstrapResult& boot_results,
    const int lambda_i,
    const int block_j,
    const rgcca::SparseMatrix& Psi
) const {
    if (lambda_i < 0 || lambda_i >= static_cast<int>(boot_results.lambda_grid.size()))
        throw std::out_of_range("RGCCA: bootstrap lambda index");
    if (block_j < 0 || block_j >= static_cast<int>(boot_results.block_names.size()))
        throw std::out_of_range("RGCCA: bootstrap block index");
    if (
        !(boot_results.ci_level > 0.0) ||
        boot_results.ci_level >= 1.0 ||
        !std::isfinite(boot_results.ci_level)
    ) {
        throw std::logic_error("RGCCA: bootstrap CI level is unavailable");
    }

    const rgcca::Matrix& w_boot =
        boot_results.w_boot_by_lambda[lambda_i][block_j];
    if (Psi.rows() <= 0 || Psi.cols() != w_boot.rows()) {
        throw std::invalid_argument(
            "RGCCA: Psi must have one column per bootstrap weight coefficient"
        );
    }
}

} // namespace fdapde

#endif // __FDAPDE_RGCCA_MODEL_VALIDATION_H__
#endif // __FDAPDE_RGCCA_DEFINE_MODEL_VALIDATION__

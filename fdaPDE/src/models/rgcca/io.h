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


#ifndef __FDAPDE_RGCCA_IO_H__
#define __FDAPDE_RGCCA_IO_H__

#include "model.h"

namespace fdapde {
namespace rgcca {

inline const char* bool_text_(const bool value) { return value ? "true" : "false"; }

inline std::ostream& operator<<(std::ostream& os, const Options& opt) {
    os << "RGCCA::Options {\n"
       << "  max_iter                    = " << opt.max_iter << '\n'
       << "  tol                         = " << opt.tol << '\n'
       << "  cache_covariances           = " << bool_text_(opt.cache_covariances) << '\n'
       << "  bias                        = " << bool_text_(opt.bias) << '\n'
       << "  init_strategy               = " << to_string(opt.init_strategy) << '\n'
       << "  lambda_selection_weights    = " << to_string(opt.lambda_selection_weights) << '\n'
       << "  lambda_selection_components = " << to_string(opt.lambda_selection_components) << '\n'
       << "  component_significance      = " << bool_text_(opt.component_significance) << '\n'
       << "  block_deactivation          = " << bool_text_(opt.block_deactivation) << '\n'
       << "  connection_deactivation     = " << bool_text_(opt.connection_deactivation) << '\n'
       << "  mode                        = " << to_string(opt.mode) << '\n'
       << "  weight_sign_constraint      = " << to_string(opt.weight_sign_constraint) << '\n'
       << "  deflation_mode              = " << to_string(opt.deflation_mode) << '\n'
       << "  scheme                      = " << opt.scheme.name << '\n'
       << "}";
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const BootstrapConfig& config) {
    os << "RGCCA::BootstrapConfig {\n"
       << "  seed                               = " << config.seed << '\n'
       << "  max_threads                        = " << config.max_threads << '\n'
       << "  B_min                              = " << config.B_min << '\n'
       << "  B_max                              = " << config.B_max << '\n'
       << "  check_every                        = " << config.check_every << '\n'
       << "  fit_max_iter                       = " << config.fit_max_iter << '\n'
       << "  adaptive                           = " << bool_text_(config.adaptive) << '\n'
       << "  adaptive_tol                       = " << config.adaptive_tol << '\n'
       << "  stable_checks_required             = " << config.stable_checks_required << '\n'
       << "  active_block_tol                   = " << config.active_block_tol << '\n'
       << "  active_connection_sign_stability   = " << config.active_connection_sign_stability << '\n'
       << "  active_connection_min_abs_corr     = " << config.active_connection_min_abs_corr << '\n'
       << "  aggressive_connection_deactivation = " << bool_text_(config.aggressive_connection_deactivation) << '\n'
       << "  min_boots_before_connection_deactivation = " << config.min_boots_before_connection_deactivation << '\n'
       << "  ci_level                           = " << config.ci_level << '\n'
       << "  patience                           = " << config.patience << '\n'
       << "  resampling_strategy                = " << to_string(config.resampling_strategy) << '\n'
       << "  stationary_block_length            = " << config.stationary_block_length << '\n'
       << "  component_significance_resamples   = " << config.component_significance_resamples << '\n'
       << "  component_significance_alpha       = " << config.component_significance_alpha << '\n'
       << "}";
    return os;
}

// pretty printer for a single Result
inline std::ostream& operator<<(std::ostream& os, const Result& r) {
    const bool minimal = false;
    if (!minimal) {
        os << "shrinkage parameters used : " << std::endl;
        for (size_t i = 0; i < r.tau_values.size(); ++i) {
            os << "- Block " << i+1  << ": tau = "<< r.tau_values[i] << "\n";
        }
        os << std::endl;
        os << "active blocks :\n";
        for (size_t i = 0; i < r.active_blocks.size(); ++i) {
            os << "- Block " << i+1  << ": " << (r.active_blocks[i] ? "active    " : "non-active" ) << "\n";
        }
        os << std::endl;
        if (r.C.rows() <= 20) {
            os << "Updated design matrix:\n";
            os << r.C << std::endl;
        } else {
            os << "Updated design matrix: skipped (" << r.C.rows() << " blocks)\n";
        }
        os << std::endl;

        os << "regularization parameters used : " << std::endl;
        os << std::scientific;
        for (size_t i = 0; i < r.tau_values.size(); ++i) {
            os << "- Block " << i+1  << ": lambda_c = "<< r.lambda_components_values[i]
               << ", lambda_l = "<< r.lambda_weights_values[i] << "\n";
        }
        os << std::fixed;
        os << std::endl;
    }
    os << "n_iters: " << r.iters << "\n";
    os << "monotone: " << (r.monotone ? "yes" : "no") << "\n";
    os << std::endl;
    os << "objective :\n";
    double prev_obj = r.obj_history[0];
    for (size_t i = 1; i < r.obj_history.size(); ++i) {
        const double obj = r.obj_history[i];
        os << "- iter " << std::setw(3) << (i)
           << "   |   fit = " << std::setw(12) << std::setprecision(8) << std::fixed << obj
           << "   |   overall diff = " << std::setw(7) << (obj - prev_obj) << "\n";
        os << std::fixed << std::setprecision(8);
        prev_obj = obj;
    }
    os << std::endl;
    if (!minimal) {
        os << "covariance matrix :\n";
        os << std::fixed << std::setprecision(2);
        os << r.covariance_matrix << std::endl;
        os << std::fixed << std::setprecision(8);
        os << "\ncorrelation matrix :\n";
        os << std::fixed << std::setprecision(2);
        os << r.correlation_matrix << std::endl;
        os << std::fixed << std::setprecision(8);
    }
    os << std::endl;
    if (std::isfinite(r.rho_tot_p_value)) {
        os << "significance:\n";
        os << "- rho_tot: " << std::fixed << r.rho_tot << "\n";
        os << "- p-value: " << std::fixed << r.rho_tot_p_value
           << " (" << r.rho_tot_bootstrap_count << " resamples)\n";
        os << "- signif.: " << (r.component_significant ? "yes" : "no") << "\n";
    }

    return os;
}

// pretty printer for a vector of Result (components)
inline std::ostream& operator<<(std::ostream& os, const std::vector<Result>& results) {
    for (size_t h = 0; h < results.size(); ++h) {
        os << "\n";
        os << "========================================\n";
        os << "Component " << (h + 1) << "\n";
        os << "----------------------------------------\n";
        os << results[h]; // delegate to the single-result printer
    }
    os << "\n";
    return os;
}

} // namespace rgcca
} // namespace fdapde

#endif // __FDAPDE_RGCCA_IO_H__

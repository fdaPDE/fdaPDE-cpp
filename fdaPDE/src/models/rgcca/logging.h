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

#ifndef __FDAPDE_RGCCA_LOGGING_H__
#define __FDAPDE_RGCCA_LOGGING_H__

#include "../header_check.h"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <ios>
#include <iostream>
#include <ostream>
#include <string>
#include <string_view>

namespace fdapde {

#ifdef FDAPDE_ENABLE_COUT

inline std::ostream& cout = std::cout;

#else

class null_ostream {
   public:
    using ostream_manipulator = std::ostream& (*)(std::ostream&);
    using ios_manipulator = std::ios& (*)(std::ios&);
    using ios_base_manipulator = std::ios_base& (*)(std::ios_base&);

    template <typename T>
    constexpr const null_ostream& operator<<(T&&) const noexcept {
        return *this;
    }

    constexpr const null_ostream& operator<<(ostream_manipulator) const noexcept {
        return *this;
    }
    constexpr const null_ostream& operator<<(ios_manipulator) const noexcept {
        return *this;
    }
    constexpr const null_ostream& operator<<(ios_base_manipulator) const noexcept {
        return *this;
    }
};

inline constexpr null_ostream cout {};

#endif

}   // namespace fdapde

#endif   // __FDAPDE_RGCCA_LOGGING_H__

#ifdef __FDAPDE_RGCCA_DEFINE_MODEL_LOGGING__
#ifndef __FDAPDE_RGCCA_MODEL_LOGGING_H__
#define __FDAPDE_RGCCA_MODEL_LOGGING_H__

namespace fdapde {

inline void log_header_(const std::string_view title) {
    fdapde::cout << '\n';
    for (std::size_t i = 0; i < title.size(); ++i) fdapde::cout << '=';
    fdapde::cout << '\n' << title << '\n';
    for (std::size_t i = 0; i < title.size(); ++i) fdapde::cout << '=';
    fdapde::cout << "\n\n";
}

// starts a timed log step and returns its start time
template <typename SamplingStrategy>
auto RGCCA<SamplingStrategy>::log_step_start_(std::string_view label) const
    -> std::chrono::high_resolution_clock::time_point {
    fdapde::cout << label << " --> ";
    return std::chrono::high_resolution_clock::now();
}

// closes a timed log step
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::log_step_end_(
    const std::chrono::high_resolution_clock::time_point start
) const {
    const auto end = std::chrono::high_resolution_clock::now();
    const double elapsed_sec = std::chrono::duration<double>(end - start).count();
    fdapde::cout << "<-- " << std::fixed << std::setprecision(3) << elapsed_sec << std::defaultfloat << "s\n";
}

// logs the active weight regularization value
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::log_weight_lambda_(const double lambda) const {
    fdapde::cout << "lambda=";
    if (std::isfinite(lambda))
        fdapde::cout << lambda;
    else
        fdapde::cout << "none";
}

// logs the fit entry header and active configuration
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::log_fit_header_(const bool run_model_selection) const {
    log_header_("RGCCA fit");
    fdapde::cout << opt_ << '\n';
    if (run_model_selection || opt_.component_significance || opt_.block_importance || inactive_block_signal_test_requested_())
        fdapde::cout << bootstrap_config_ << '\n';
}

// logs the bootstrap model-selection section header
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::log_bootstrap_model_selection_header_() const {
    log_header_(std::string("Bootstrap model selection for component ") + std::to_string(h_ + 1));
}

// logs the component-significance bootstrap section header
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::log_bootstrap_component_significance_header_() const {
    log_header_(std::string("Bootstrap component significance for component ") + std::to_string(h_ + 1));
}

// logs the block-importance bootstrap section header
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::log_bootstrap_block_importance_header_() const {
    log_header_(std::string("Bootstrap block importance for component ") + std::to_string(h_ + 1));
}

// logs the inactive-block signal test section header
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::log_inactive_block_signal_test_header_() const {
    log_header_(std::string("Inactive-block signal test for component ") + std::to_string(h_ + 1));
}

// logs the lambda candidate currently being evaluated
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::log_bootstrap_lambda_candidate_(
    const bool lambda_selected,
    const double lambda
) const {
    if (lambda_selected) {
        fdapde::cout << "- lambda = " << lambda << '\n';
        return;
    }

    fdapde::cout << "- fixed weight regularization ";
    log_weight_lambda_(lambda);
    fdapde::cout << '\n';
}

// logs the component-significance test result
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::log_component_significance_(
    const typename RGCCA<SamplingStrategy>::ComponentSignificanceResult& significance
) const {
    fdapde::cout << "\nSignificance:\n"
                 << "  - rho_tot = " << significance.rho_tot << '\n'
                 << "  - p-value = " << significance.p_value << '\n'
                 << "  - significant = " << std::boolalpha << significance.significant << std::noboolalpha << "\n\n";
}

// logs the block-importance test result
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::log_block_importance_(
    const typename RGCCA<SamplingStrategy>::BlockImportanceResult& importance
) const {
    fdapde::cout << "\nBlock importance:\n";
    for (int j = 0; j < static_cast<int>(importance.rho.size()); ++j) {
        fdapde::cout << "  - block[" << j << "]"
                     << ": rho = " << importance.rho[j]
                     << ", p-value = " << importance.p_value[j]
                     << ", significant = " << std::boolalpha << importance.significant[j] << std::noboolalpha
                     << '\n';
    }
    fdapde::cout << '\n';
}

// logs one inactive-block gate decision before fitting the next component
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::log_inactive_block_signal_gate_(
    const int block,
    const InactiveBlockSignalAction action,
    const double statistic,
    const double p_value
) const {
    fdapde::cout << "Inactive-block signal gate for component " << h_ + 1
                 << ": block[" << block << "] " << rgcca::to_string(action);
    if (std::isfinite(statistic))
        fdapde::cout << ", stat=" << std::setprecision(4) << statistic;
    if (std::isfinite(p_value))
        fdapde::cout << ", p=" << std::setprecision(4) << p_value;
    fdapde::cout << std::defaultfloat << '\n';
}

// logs the final bootstrap summary for one lambda candidate
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::log_bootstrap_lambda_summary_(
    const typename RGCCA<SamplingStrategy>::AdaptiveBootstrapState& state,
    const typename RGCCA<SamplingStrategy>::BootstrapTimingSummary& timing_summary,
    const double elapsed_sec,
    const int n_active_blocks,
    const int n_active_connections,
    const double lambda
) const {
    // sample accounting
    const int B_discarded =
        state.B_total -
        state.B_design -
        state.B_stale -
        state.B_cancelled -
        state.B_done;

    // timing and iteration recap
    fdapde::cout << "  Bootstrap used: total=" << state.B_total
                 << ", design=" << state.B_design
                 << ", stale=" << state.B_stale
                 << ", cancelled=" << state.B_cancelled
                 << ", good=" << state.B_done
                 << ", discarded=" << B_discarded << '\n';
    fdapde::cout << "  Time: total=" << std::fixed << std::setprecision(3) << elapsed_sec
                 << "s, avg_fit=" << timing_summary.avg_fit_time()
                 << " ± " << timing_summary.sd_fit_time()
                 << "s, efficiency=" << std::setprecision(1)
                 << 100.0 * timing_summary.efficiency() << "%\n";
    fdapde::cout << "  Iters: avg=" << std::setprecision(1) << timing_summary.avg_iters()
                 << ", max=" << timing_summary.max_iters
                 << ", capped=" << timing_summary.capped_fits
                 << "/" << timing_summary.n_fits << '\n';

    // design and criterion recap
    fdapde::cout << "  Design: active_blocks=" << n_active_blocks
                 << ", active_connections=" << n_active_connections
                 << ", crit=" << std::setprecision(3) << state.crit
                 << std::defaultfloat << ", ";
    log_weight_lambda_(lambda);
    fdapde::cout << std::defaultfloat << "\n\n";
}

// logs the selected bootstrap model and resulting active design
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::log_bootstrap_selection_result_(
    const typename RGCCA<SamplingStrategy>::BootstrapResult& boot_results,
    const rgcca::BoolMatrix& C_best,
    const bool lambda_selected
) const {
    // lambda selection result
    if (lambda_selected) {
        fdapde::cout << "\nOptimal lambda: " << boot_results.lambda_opt << '\n';
    } else {
        fdapde::cout << "\nNo weight lambda selection requested\n";
    }

    // block deactivation result
    fdapde::cout << "\nBlock deactivation:\n";
    for (int j = 0; j < static_cast<int>(boot_results.active_blocks.size()); ++j) {
        const double nrm = boot_results.w_min_by_lambda[boot_results.lambda_opt_index][j].norm();
        fdapde::cout << "- block[" << std::setw(2) << j << "]: "
                     << "||w_min|| = " << std::fixed << std::setprecision(4) << nrm
                     << std::defaultfloat
                     << ", active = " << (boot_results.active_blocks[j] ? "yes" : "no")
                     << '\n';
    }

    // design matrix recap
    const int J = static_cast<int>(C_best.rows());
    if (J <= 20) {
        fdapde::cout << "\nUpdated design matrix:\n";
        fdapde::cout << C_best << "\n\n";
    } else {
        fdapde::cout << "\nUpdated design matrix: skipped (" << J << " blocks)\n\n";
    }
}

// logs lambda-grid early stopping
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::log_bootstrap_early_stop_(const int patience) const {
    fdapde::cout << "  Early stop: no improvement for "
                 << patience
                 << " consecutive lambdas\n\n";
}

// logs streaming bootstrap progress
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::log_bootstrap_progress_(
    const typename RGCCA<SamplingStrategy>::AdaptiveBootstrapState& state,
    const int n_active_blocks,
    const int n_active_connections,
    const typename RGCCA<SamplingStrategy>::AdaptiveStopInfo& stop_info,
    const double elapsed_since_last_log
) const {
    fdapde::cout << "  " << std::left << std::setw(13) << "progress" << std::right
                 << " total=" << std::setw(6) << state.B_total
                 << " design=" << std::setw(6) << state.B_design
                 << " stale=" << std::setw(4) << state.B_stale
                 << " cancelled=" << std::setw(6) << state.B_cancelled
                 << " good=" << std::setw(6) << state.B_done
                 << " | dt=" << std::fixed << std::setprecision(3) << std::setw(8)
                 << elapsed_since_last_log << "s"
                 << std::defaultfloat
                 << " | ab=" << std::setw(4) << n_active_blocks
                 << ", ac=" << std::setw(6) << n_active_connections
                 << " | crit=" << std::setprecision(3) << state.crit
                 << ", rel=";
    if (std::isfinite(stop_info.rel_change)) {
        fdapde::cout << std::setprecision(3) << stop_info.rel_change;
    } else {
        fdapde::cout << "  -  ";
    }
    if (stop_info.stop)
        fdapde::cout << " | stop=" << stop_info.reason;
    fdapde::cout << std::defaultfloat << '\n';
}

// logs a bootstrap design reset after block or connection deactivation
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::log_bootstrap_design_reset_(
    const typename RGCCA<SamplingStrategy>::AdaptiveBootstrapState& state,
    const int n_active_blocks,
    const int n_active_connections,
    const double elapsed_since_last_log
) const {
    fdapde::cout << "  " << std::left << std::setw(13) << "design reset" << std::right
                 << " total=" << std::setw(6) << state.B_total
                 << " design=" << std::setw(6) << state.B_design
                 << " stale=" << std::setw(4) << state.B_stale
                 << " cancelled=" << std::setw(6) << state.B_cancelled
                 << " good=" << std::setw(6) << state.B_done
                 << " | dt=" << std::fixed << std::setprecision(3) << std::setw(8)
                 << elapsed_since_last_log << "s"
                 << std::defaultfloat
                 << " | ab=" << std::setw(4) << n_active_blocks
                 << ", ac=" << std::setw(6) << n_active_connections
                 << std::defaultfloat
                 << " | epoch=" << state.design_epoch
                 << '\n';
}

} // namespace fdapde

#endif // __FDAPDE_RGCCA_MODEL_LOGGING_H__
#endif // __FDAPDE_RGCCA_DEFINE_MODEL_LOGGING__

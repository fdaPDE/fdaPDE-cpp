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

#ifndef __FDAPDE_RGCCA_BOOTSTRAP_UTILS_H__
#define __FDAPDE_RGCCA_BOOTSTRAP_UTILS_H__

namespace fdapde {

// returns bootstrap weight confidence intervals for a selected component, lambda and block
template <typename SamplingStrategy>
auto RGCCA<SamplingStrategy>::bootstrap_weights_ci(
    const int h,
    const int lambda_i,
    const int block_j,
    const rgcca::SparseMatrix& Psi
) const -> std::pair<rgcca::Vector, rgcca::Vector> {
    return bootstrap_weights_ci_(h, lambda_i, block_j, Psi);
}

// returns bootstrap weight confidence intervals at the selected lambda for a component and block
template <typename SamplingStrategy>
auto RGCCA<SamplingStrategy>::bootstrap_weights_ci(
    const int h,
    const int block_j,
    const rgcca::SparseMatrix& Psi
) const -> std::pair<rgcca::Vector, rgcca::Vector> {
    const auto& boot_results = bootstrap_selection_result_(h);
    return bootstrap_weights_ci_(h, bootstrap_lambda_opt_index_(boot_results), block_j, Psi);
}

// evaluates bootstrap weight confidence intervals on locations for an explicit lambda index
template <typename SamplingStrategy>
template <typename DataLocs>
requires(!std::same_as<std::decay_t<DataLocs>, rgcca::SparseMatrix>)
auto RGCCA<SamplingStrategy>::bootstrap_weights_ci(
    const int h,
    const int lambda_i,
    const int block_j,
    const DataLocs& locs
) const -> std::pair<rgcca::Vector, rgcca::Vector> {
    validate_index_(block_j);
    const rgcca::SparseMatrix Psi = blocks_[block_j]->Psi_at(locs);
    return bootstrap_weights_ci(h, lambda_i, block_j, Psi);
}

// evaluates bootstrap weight confidence intervals on locations at the selected lambda
template <typename SamplingStrategy>
template <typename DataLocs>
requires(!std::same_as<std::decay_t<DataLocs>, rgcca::SparseMatrix>)
auto RGCCA<SamplingStrategy>::bootstrap_weights_ci(
    const int h,
    const int block_j,
    const DataLocs& locs
) const -> std::pair<rgcca::Vector, rgcca::Vector> {
    validate_index_(block_j);
    const rgcca::SparseMatrix Psi = blocks_[block_j]->Psi_at(locs);
    return bootstrap_weights_ci(h, block_j, Psi);
}

// finds stored bootstrap model-selection results for one component
template <typename SamplingStrategy>
auto RGCCA<SamplingStrategy>::bootstrap_selection_result_(const int h) const
    -> const typename RGCCA<SamplingStrategy>::BootstrapResult& {
    if (bootstrap_selection_results_.empty()) {
        throw std::logic_error(
            "RGCCA: bootstrap weight CIs require automatic weight lambda selection results"
        );
    }

    for (const auto& result : bootstrap_selection_results_) {
        if (result.h == h) return result;
    }

    throw std::out_of_range("RGCCA: bootstrap component index");
}

// resolves the selected lambda index from bootstrap model-selection results
template <typename SamplingStrategy>
int RGCCA<SamplingStrategy>::bootstrap_lambda_opt_index_(
    const typename RGCCA<SamplingStrategy>::BootstrapResult& boot_results
) const {
    if (
        boot_results.lambda_opt_index >= 0 &&
        boot_results.lambda_opt_index < static_cast<int>(boot_results.lambda_grid.size())
    ) {
        return boot_results.lambda_opt_index;
    }

    for (int i = 0; i < static_cast<int>(boot_results.lambda_grid.size()); ++i) {
        if (boot_results.lambda_grid[i] == boot_results.lambda_opt) return i;
    }

    throw std::logic_error("RGCCA: bootstrap optimal lambda index is unavailable");
}

// computes pointwise weight confidence intervals from stored bootstrap weights
template <typename SamplingStrategy>
auto RGCCA<SamplingStrategy>::bootstrap_weights_ci_(
    const int h,
    const int lambda_i,
    const int block_j,
    const rgcca::SparseMatrix& Psi
) const -> std::pair<rgcca::Vector, rgcca::Vector> {
    const auto& boot_results = bootstrap_selection_result_(h);
    validate_bootstrap_weights_ci_(boot_results, lambda_i, block_j, Psi);
    const rgcca::Matrix& w_boot =
        boot_results.w_boot_by_lambda[lambda_i][block_j];

    const double alpha_low = (1.0 - boot_results.ci_level) / 2.0;
    const double alpha_high = 1.0 - alpha_low;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const int B_eff = std::min(
        boot_results.B_used_by_lambda[lambda_i],
        static_cast<int>(w_boot.cols())
    );

    rgcca::Vector ci_low(Psi.rows());
    rgcca::Vector ci_high(Psi.rows());

    if (B_eff <= 0) {
        ci_low.setConstant(nan);
        ci_high.setConstant(nan);
        return {ci_low, ci_high};
    }

    const rgcca::Matrix w_eval = Psi * w_boot.leftCols(B_eff);

    for (int r = 0; r < w_eval.rows(); ++r) {
        std::vector<double> values;
        values.reserve(B_eff);

        for (int b = 0; b < B_eff; ++b)
            values.push_back(w_eval(r, b));

        ci_low[r] = internals::empirical_quantile(values, alpha_low);
        ci_high[r] = internals::empirical_quantile(values, alpha_high);
    }

    return {ci_low, ci_high};
}

// picks the bootstrap fit iteration cap, falling back to the main fit cap
template <typename SamplingStrategy>
int RGCCA<SamplingStrategy>::bootstrap_fit_max_iter_() const {
    return bootstrap_config_.fit_max_iter > 0 ? bootstrap_config_.fit_max_iter : opt_.max_iter;
}

// configures and returns the number of bootstrap worker threads
template <typename SamplingStrategy>
int RGCCA<SamplingStrategy>::bootstrap_n_threads_() const {
    const int max_threads = bootstrap_config_.max_threads;
    if (max_threads <= 0) {
        throw std::invalid_argument("RGCCA: bootstrap max_threads must be positive");
    }

    const int available_threads = std::max(1, static_cast<int>(fdapde::available_concurrency()));
    const int requested_threads = std::min(max_threads, available_threads);
    parallel_set_num_threads(requested_threads);
    const int actual_threads = static_cast<int>(internals::threaded_executor::instance().size());

    if (actual_threads > max_threads) {
        throw std::runtime_error(
            "RGCCA: bootstrap allows at most " + std::to_string(max_threads) +
            " threads, but the executor is already initialized with " +
            std::to_string(actual_threads) + " threads"
        );
    }

    return actual_threads;
}

// fits one bootstrap resample and returns its weights, correlations and timing
template <typename SamplingStrategy>
auto RGCCA<SamplingStrategy>::fit_bootstrap_sample_(
    typename RGCCA<SamplingStrategy>::BootstrapBlocks& boot_blocks,
    const int b,
    const unsigned seed,
    const rgcca::BoolMatrix& C_active,
    const std::vector<rgcca::Vector>& w_fit,
    const std::function<bool()>& cancelled
) -> typename RGCCA<SamplingStrategy>::BootstrapSampleResult {
    typename RGCCA<SamplingStrategy>::BootstrapSampleResult out;
    const auto active_blocks = active_blocks_from_C_(C_active);
    const int fit_max_iter = bootstrap_fit_max_iter_();

    if (cancelled && cancelled()) {
        out.cancelled = true;
        return out;
    }

    copy_weights_snapshot_(boot_blocks.refs, w_fit);
    set_row_index_all_(boot_blocks.refs, bootstrap_index_(n_, seed + static_cast<unsigned>(b)));

    init_comp_(boot_blocks.refs, rgcca::InitStrategy::WarmStart, true, &active_blocks);
    if (cancelled && cancelled()) {
        out.cancelled = true;
        clear_row_index_all_(boot_blocks.refs);
        return out;
    }

    const auto fit_start = std::chrono::high_resolution_clock::now();
    out.fit_started = true;
    const typename RGCCA<SamplingStrategy>::Result fit_result =
        fit_component_(boot_blocks.refs, C_active, true, fit_max_iter, cancelled);
    const auto fit_end = std::chrono::high_resolution_clock::now();

    out.fit_time = std::chrono::duration<double>(fit_end - fit_start).count();
    out.fit_iters = fit_result.iters;
    out.capped = fit_result.iters >= fit_max_iter;
    out.cancelled = fit_result.cancelled;
    if (out.cancelled) {
        clear_row_index_all_(boot_blocks.refs);
        return out;
    }

    out.w = snapshot_weights_(boot_blocks.refs);
    for (int j = 0; j < static_cast<int>(out.w.size()); ++j) {
        if (out.w[j].dot(w_fit[j]) < 0.0) out.w[j] *= -1.0;
    }
    correlation_matrix_raw_data_(boot_blocks.refs, out.w, out.corr);
    clear_row_index_all_(boot_blocks.refs);

    return out;
}

// runs adaptive bootstrap resampling for one lambda value
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::run_bootstrap_stream_(
    int lambda_i,
    typename RGCCA<SamplingStrategy>::AdaptiveBootstrapState& bootstrap_state,
    std::vector<typename RGCCA<SamplingStrategy>::BootstrapBlocks>& thread_boot_worker,
    rgcca::BoolMatrix& C_active,
    const std::vector<rgcca::Vector>& w_fit,
    std::vector<rgcca::Vector>& w_min,
    typename RGCCA<SamplingStrategy>::BootstrapResult& boot_results,
    const typename RGCCA<SamplingStrategy>::BlockRefList& blocks,
    typename RGCCA<SamplingStrategy>::BootstrapTimingSummary& timing_summary
) {
    const int J = static_cast<int>(w_fit.size());
    const int n_threads = bootstrap_state.n_threads;
    const unsigned seed = static_cast<unsigned>(bootstrap_state.seed);
    std::atomic<int> next_boot {0};
    std::atomic<bool> stop {false};
    std::mutex merge_mutex;

    const auto parallel_start = std::chrono::high_resolution_clock::now();
    auto last_log_time = parallel_start;
    auto elapsed_since_last_log = [&]() {
        const auto now = std::chrono::high_resolution_clock::now();
        const double elapsed = std::chrono::duration<double>(now - last_log_time).count();
        last_log_time = now;
        return elapsed;
    };

    parallel_for(0, n_threads, 1, [&](int) {
        while (!stop.load(std::memory_order_acquire)) {
            const int b = next_boot.fetch_add(1, std::memory_order_relaxed);

            // take a consistent design snapshot for this sample
            rgcca::BoolMatrix C_snapshot;
            int epoch = 0;
            {
                std::lock_guard<std::mutex> lock(merge_mutex);
                if (bootstrap_state.stop || bootstrap_state.B_done >= bootstrap_state.B_max) {
                    stop.store(true, std::memory_order_release);
                    break;
                }
                C_snapshot = C_active;
                epoch = bootstrap_state.design_epoch;
            }

            // fit outside the merge lock; cancellation handles design resets
            const int tid = this_thread_id();
            auto cancelled = [&]() {
                return stop.load(std::memory_order_acquire) ||
                    epoch != bootstrap_state.design_epoch_signal.load(std::memory_order_acquire);
            };
            auto sample = fit_bootstrap_sample_(
                thread_boot_worker[tid], b, seed, C_snapshot, w_fit, cancelled
            );

            // merge sample accounting and discard stale work
            std::lock_guard<std::mutex> lock(merge_mutex);
            ++bootstrap_state.B_total;
            timing_summary.add_sample(sample.fit_time, sample.fit_iters, sample.capped, sample.fit_started);

            if (sample.cancelled) {
                ++bootstrap_state.B_cancelled;
                continue;
            }

            if (bootstrap_state.stop || epoch != bootstrap_state.design_epoch) {
                if (epoch != bootstrap_state.design_epoch)
                    ++bootstrap_state.B_stale;
                continue;
            }

            if (bootstrap_state.B_done >= bootstrap_state.B_max) {
                bootstrap_state.stop = true;
                stop.store(true, std::memory_order_release);
                continue;
            }

            // store accepted weights and correlations
            const int b_good = bootstrap_state.B_done;
            for (int j = 0; j < J; ++j) {
                boot_results.w_boot_by_lambda[lambda_i][j].col(b_good) = sample.w[j];
                update_w_min_(w_min[j], w_fit[j], sample.w[j]);
            }
            boot_results.corr_boot_by_lambda[lambda_i].col(b_good) =
                Eigen::Map<const rgcca::Vector>(sample.corr.data(), J * J);

            for (int j = 0; j < J; ++j) {
                for (int k = j + 1; k < J; ++k) {
                    const double c = sample.corr(j, k);
                    if (c > 0.0) {
                        ++bootstrap_state.corr_pos_count(j, k);
                        ++bootstrap_state.corr_pos_count(k, j);
                    } else if (c < 0.0) {
                        ++bootstrap_state.corr_neg_count(j, k);
                        ++bootstrap_state.corr_neg_count(k, j);
                    }
                }
            }
            ++bootstrap_state.B_done;

            // block deactivation invalidates previously accepted samples
            if (opt_.block_deactivation) {
                const rgcca::BoolMatrix C_before = C_active;
                threshold_inactive_blocks_(w_min, C_active);
                if (!same_design_(C_before, C_active)) {
                    reset_good_bootstrap_(bootstrap_state, w_min, w_fit, C_active);
                    const int n_active_blocks = count_active_blocks_(C_active);
                    const int n_active_connections = count_active_connections_(C_active);
                    log_bootstrap_design_reset_(
                        bootstrap_state,
                        n_active_blocks,
                        n_active_connections,
                        elapsed_since_last_log()
                    );
                    if (n_active_connections == 0) {
                        bootstrap_state.crit = 0.0;
                        bootstrap_state.stop = true;
                        stop.store(true, std::memory_order_release);
                    }
                    continue;
                }
            }

            // expensive checks run only at configured checkpoints
            const bool force_check = bootstrap_state.B_done >= bootstrap_state.B_max;
            const bool due_check =
                force_check ||
                bootstrap_state.B_done - bootstrap_state.last_check_B_done >= bootstrap_state.check_every;
            if (!due_check) continue;
            bootstrap_state.last_check_B_done = bootstrap_state.B_done;

            // connection deactivation also restarts the accepted-sample stream
            if (connection_deactivation_ready_(bootstrap_state)) {
                const rgcca::BoolMatrix C_before = C_active;
                threshold_inactive_connections_(lambda_i, bootstrap_state, boot_results, C_active);
                if (!same_design_(C_before, C_active)) {
                    reset_good_bootstrap_(bootstrap_state, w_min, w_fit, C_active);
                    const int n_active_blocks = count_active_blocks_(C_active);
                    const int n_active_connections = count_active_connections_(C_active);
                    log_bootstrap_design_reset_(
                        bootstrap_state,
                        n_active_blocks,
                        n_active_connections,
                        elapsed_since_last_log()
                    );
                    if (n_active_connections == 0) {
                        bootstrap_state.crit = 0.0;
                        bootstrap_state.stop = true;
                        stop.store(true, std::memory_order_release);
                    }
                    continue;
                }
            }

            // criterion and adaptive stopping use the current minimum-weight envelope
            bootstrap_state.crit = criterion_score_with_weights_(blocks, w_min, C_);
            const typename RGCCA<SamplingStrategy>::AdaptiveStopInfo stop_info =
                bootstrap_config_.adaptive ?
                    adaptive_stop_(bootstrap_state, bootstrap_config_) :
                    typename RGCCA<SamplingStrategy>::AdaptiveStopInfo{};
            log_bootstrap_progress_(
                bootstrap_state,
                count_active_blocks_(C_active),
                count_active_connections_(C_active),
                stop_info,
                elapsed_since_last_log()
            );

            if (stop_info.stop || bootstrap_state.B_done >= bootstrap_state.B_max) {
                bootstrap_state.stop = true;
                stop.store(true, std::memory_order_release);
            }
        }
    });

    const auto parallel_end = std::chrono::high_resolution_clock::now();
    const double parallel_time = std::chrono::duration<double>(parallel_end - parallel_start).count();
    timing_summary.set_parallel_capacity(parallel_time, n_threads);
}

// compares two bootstrap design matrices entry by entry
template <typename SamplingStrategy>
bool RGCCA<SamplingStrategy>::same_design_(
    const rgcca::BoolMatrix& lhs,
    const rgcca::BoolMatrix& rhs
) const {
    if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) return false;
    for (int i = 0; i < lhs.rows(); ++i)
        for (int j = 0; j < lhs.cols(); ++j)
            if (lhs(i, j) != rhs(i, j)) return false;
    return true;
}

// discards accepted samples after a design change and keeps inactive blocks zeroed
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::reset_good_bootstrap_(
    typename RGCCA<SamplingStrategy>::AdaptiveBootstrapState& state,
    std::vector<rgcca::Vector>& w_min,
    const std::vector<rgcca::Vector>& w_fit,
    const rgcca::BoolMatrix& C_active
) const {
    state.B_design += state.B_done;
    state.reset_good();
    ++state.design_epoch;
    state.design_epoch_signal.store(state.design_epoch, std::memory_order_release);
    w_min = w_fit; // !!!!! Not sure about this

    const auto active_blocks = active_blocks_from_C_(C_active);
    for (int j = 0; j < static_cast<int>(w_min.size()); ++j) {
        if (!active_blocks[j])
            w_min[j].setZero();
    }
}

// checks whether connection deactivation can run for the current bootstrap state
template <typename SamplingStrategy>
bool RGCCA<SamplingStrategy>::connection_deactivation_ready_(
    const typename RGCCA<SamplingStrategy>::AdaptiveBootstrapState& state
) const {
    if (!opt_.connection_deactivation) return false;
    if (state.B_done < bootstrap_config_.min_boots_before_connection_deactivation) return false;
    return bootstrap_config_.aggressive_connection_deactivation || state.B_done >= state.B_min;
}

// allocates bootstrap result matrices for one lambda when first needed
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::ensure_bootstrap_lambda_storage_(
    typename RGCCA<SamplingStrategy>::BootstrapResult& boot_results,
    const int lambda_i,
    const std::vector<int>& block_dims,
    const int J
) const {
    if (boot_results.corr_boot_by_lambda[lambda_i].size() != 0) return;

    boot_results.w_boot_by_lambda[lambda_i].resize(J);
    for (int j = 0; j < J; ++j)
        boot_results.w_boot_by_lambda[lambda_i][j].setZero(block_dims[j], boot_results.B);
    boot_results.corr_boot_by_lambda[lambda_i].setZero(J * J, boot_results.B);
}

// decides whether adaptive bootstrap stopping has stabilized
template <typename SamplingStrategy>
auto RGCCA<SamplingStrategy>::adaptive_stop_(
    typename RGCCA<SamplingStrategy>::AdaptiveBootstrapState& state,
    const typename RGCCA<SamplingStrategy>::BootstrapConfig& config
) const -> typename RGCCA<SamplingStrategy>::AdaptiveStopInfo {
    typename RGCCA<SamplingStrategy>::AdaptiveStopInfo out;
    if (state.crit == 0.0) {
        out.stop = true;
        out.reason = "crit=0";
        return out;
    }

    if (!std::isfinite(state.crit_prev_check)) {
        state.crit_prev_check = state.crit;
        return out;
    }

    const double rel_change =
        std::abs(state.crit_prev_check - state.crit) /
        (std::abs(state.crit_prev_check) + 1e-12);
    out.rel_change = rel_change;

    if (rel_change < config.adaptive_tol) {
        ++state.stable_checks;
    } else {
        state.stable_checks = 0;
    }

    state.crit_prev_check = state.crit;

    if (state.stable_checks >= config.stable_checks_required && state.B_done >= state.B_min) {
        out.stop = true;
        out.reason = "stable";
    }

    return out;
}

// updates the best lambda candidate and applies lambda-level early stopping
template <typename SamplingStrategy>
bool RGCCA<SamplingStrategy>::early_stop_lambda_(
    typename RGCCA<SamplingStrategy>::AdaptiveBootstrapState& state,
    int lambda_i,
    const typename RGCCA<SamplingStrategy>::BootstrapConfig& config
) const {
    if (state.crit > state.best_criterion) {
        state.best_criterion = state.crit;
        state.best_i = lambda_i;
        state.no_improve = 0;
    } else {
        ++state.no_improve;
    }

    if (state.no_improve >= config.patience) {
        log_bootstrap_early_stop_(config.patience);
        return true;
    }

    return false;
}

// deactivates blocks whose minimum bootstrap weights are below tolerance
template <typename SamplingStrategy>
int RGCCA<SamplingStrategy>::threshold_inactive_blocks_(
    std::vector<rgcca::Vector>& w_min,
    rgcca::BoolMatrix& C_active
) const {
    const int J = static_cast<int>(w_min.size());

    for (int j = 0; j < J; ++j) {
        const double nrm = w_min[j].norm();

        if (nrm < bootstrap_config_.active_block_tol) {
            w_min[j].setZero();
            C_active.row(j).setConstant(false);
            C_active.col(j).setConstant(false);
        }
    }

    int n_active_blocks = 0;
    int last_active_block = -1;
    auto active_blocks = active_blocks_from_C_(C_active);
    for (int j = 0; j < J; ++j) {
        if (active_blocks[j]) {
            ++n_active_blocks;
            last_active_block = j;
        }
    }

    if (n_active_blocks == 1) {
        w_min[last_active_block].setZero();
        C_active.row(last_active_block).setConstant(false);
        C_active.col(last_active_block).setConstant(false);
        n_active_blocks = 0;
    }

    return n_active_blocks;
}

// deactivates unstable or weak bootstrap connections
template <typename SamplingStrategy>
int RGCCA<SamplingStrategy>::threshold_inactive_connections_(
    int lambda_i,
    typename RGCCA<SamplingStrategy>::AdaptiveBootstrapState& state,
    const typename RGCCA<SamplingStrategy>::BootstrapResult& boot_results,
    rgcca::BoolMatrix& C_active,
    const int B_eff_override
) {
    const int J = static_cast<int>(C_active.rows());
    const int B_eff = B_eff_override > 0 ? B_eff_override : state.B_done;
    if (B_eff <= 0)
        return count_active_connections_(C_active);

    for (int j = 0; j < J; ++j) {
        for (int k = j + 1; k < J; ++k) {
            if (!C_active(j, k)) continue;

            std::vector<double> abs_corr;
            abs_corr.reserve(B_eff);

            const int row_jk = j + k * J;

            for (int b = 0; b < B_eff; ++b) {
                const double corr_jk = boot_results.corr_boot_by_lambda[lambda_i](row_jk, b);
                abs_corr.push_back(std::abs(corr_jk));
            }

            const int n_pos = state.corr_pos_count(j, k);
            const int n_neg = state.corr_neg_count(j, k);
            const double sign_stability = B_eff > 0 ?
                static_cast<double>(std::max(n_pos, n_neg)) / static_cast<double>(B_eff) :
                0.0;

            const double med_abs_corr = internals::median(abs_corr);

            const bool active =
                sign_stability >= bootstrap_config_.active_connection_sign_stability &&
                med_abs_corr >= bootstrap_config_.active_connection_min_abs_corr;

            if (!active) {
                C_active(j, k) = false;
                C_active(k, j) = false;
            }
        }
    }

    deactivate_isolated_blocks_(C_active);

    return count_active_connections_(C_active);
}

// restores the full design while preserving already inactive blocks
template <typename SamplingStrategy>
auto RGCCA<SamplingStrategy>::reset_connections_keep_inactive_blocks_(
    const rgcca::BoolMatrix& C_full,
    const rgcca::BoolMatrix& C_current
) const -> rgcca::BoolMatrix {
    rgcca::BoolMatrix C_reset = C_full;

    const auto active_blocks = active_blocks_from_C_(C_current);
    const int J = static_cast<int>(C_reset.rows());

    for (int j = 0; j < J; ++j) {
        if (!active_blocks[j]) {
            C_reset.row(j).setConstant(false);
            C_reset.col(j).setConstant(false);
        }
    }

    for (int j = 0; j < J; ++j)
        C_reset(j, j) = false;

    return C_reset;
}

// removes blocks that have no remaining active connections
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::deactivate_isolated_blocks_(
    rgcca::BoolMatrix& C_active
) const {
    const int J = static_cast<int>(C_active.rows());

    for (int j = 0; j < J; ++j) {
        bool active = false;

        for (int k = 0; k < J; ++k) {
            if (C_active(j, k)) {
                active = true;
                break;
            }
        }

        if (!active) {
            C_active.row(j).setConstant(false);
            C_active.col(j).setConstant(false);
        }
    }
}

// counts active upper-triangular design connections
template <typename SamplingStrategy>
int RGCCA<SamplingStrategy>::count_active_connections_(
    const rgcca::BoolMatrix& C_active
) const {
    const int J = static_cast<int>(C_active.rows());

    int n_active_connections = 0;

    for (int j = 0; j < J; ++j) {
        for (int k = j + 1; k < J; ++k) {
            if (C_active(j, k))
                ++n_active_connections;
        }
    }

    return n_active_connections;
}

// counts blocks with at least one active connection
template <typename SamplingStrategy>
int RGCCA<SamplingStrategy>::count_active_blocks_(
    const rgcca::BoolMatrix& C_active
) const {
    const auto active_blocks = active_blocks_from_C_(C_active);
    return static_cast<int>(std::count(active_blocks.begin(), active_blocks.end(), true));
}

// returns a block activity mask from a design matrix
template <typename SamplingStrategy>
std::vector<bool> RGCCA<SamplingStrategy>::active_blocks_from_C_(
    const rgcca::BoolMatrix& C_active
) const {
    const int J = static_cast<int>(C_active.rows());

    std::vector<bool> active_blocks(J, false);

    for (int j = 0; j < J; ++j) {
        for (int k = 0; k < J; ++k) {
            if (C_active(j, k)) {
                active_blocks[j] = true;
                break;
            }
        }
    }

    return active_blocks;
}

// computes a Fisher z confidence interval for one bootstrap correlation row
template <typename SamplingStrategy>
std::pair<double, double> RGCCA<SamplingStrategy>::fisher_z_corr_ci_(
    const rgcca::Matrix& corr_boot,
    const int row,
    const int B_eff,
    const double alpha_low,
    const double alpha_high
) const {
    constexpr double eps = 1e-12;

    std::vector<double> z_values;
    z_values.reserve(B_eff);

    for (int b = 0; b < B_eff; ++b) {
        const double corr = corr_boot(row, b);
        if (!std::isfinite(corr)) continue;

        const double corr_clamped = std::clamp(corr, -1.0 + eps, 1.0 - eps);
        z_values.push_back(std::atanh(corr_clamped));
    }

    const double z_low = internals::empirical_quantile(z_values, alpha_low);
    const double z_high = internals::empirical_quantile(z_values, alpha_high);

    if (!std::isfinite(z_low) || !std::isfinite(z_high)) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return {nan, nan};
    }

    return {std::tanh(z_low), std::tanh(z_high)};
}

// updates the elementwise minimum bootstrap weight envelope for one block
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::update_w_min_(
    rgcca::Vector& w_min_j,
    const rgcca::Vector& w_fit_j,
    const rgcca::Vector& w_bj
) const {
    for (int r = 0; r < w_min_j.size(); ++r) {
        if (w_fit_j[r] > 0.0) {
            w_min_j[r] = std::max(0.0, std::min(w_min_j[r], w_bj[r]));
        } else if (w_fit_j[r] < 0.0) {
            w_min_j[r] = std::min(0.0, std::max(w_min_j[r], w_bj[r]));
        } else {
            w_min_j[r] = 0.0;
        }
    }
}

// shrinks bootstrap storage to the number of accepted samples per lambda
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::resize_bootstrap_results_(
    typename RGCCA<SamplingStrategy>::BootstrapResult& boot_results,
    int n_lambdas,
    int J,
    const std::vector<int>& block_dims
) const {
    for (int i = 0; i < n_lambdas; ++i) {
        const int B_eff = boot_results.B_used_by_lambda[i];

        if (boot_results.corr_boot_by_lambda[i].size() == 0) {
            boot_results.w_boot_by_lambda[i].resize(J);
            for (int j = 0; j < J; ++j)
                boot_results.w_boot_by_lambda[i][j].setZero(block_dims[j], 0);
            boot_results.corr_boot_by_lambda[i].setZero(J * J, 0);
            continue;
        }

        for (int j = 0; j < J; ++j) {
            boot_results.w_boot_by_lambda[i][j].conservativeResize(
                Eigen::NoChange, B_eff
            );
        }

        boot_results.corr_boot_by_lambda[i].conservativeResize(
            Eigen::NoChange, B_eff
        );
    }
}

// computes bootstrap correlation confidence intervals for every lambda and block pair
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::compute_bootstrap_corr_cis_(
    typename RGCCA<SamplingStrategy>::BootstrapResult& boot_results,
    int J
) const {
    const double alpha_low = (1.0 - boot_results.ci_level) / 2.0;
    const double alpha_high = 1.0 - alpha_low;
    const double nan = std::numeric_limits<double>::quiet_NaN();

    for (std::size_t i = 0; i < boot_results.lambda_grid.size(); ++i) {
        const int B_eff = boot_results.B_used_by_lambda[i];

        if (B_eff <= 0) {
            boot_results.corr_ci_low_by_lambda[i].setConstant(J, J, nan);
            boot_results.corr_ci_high_by_lambda[i].setConstant(J, J, nan);

            continue;
        }

        boot_results.corr_ci_low_by_lambda[i].setZero(J, J);
        boot_results.corr_ci_high_by_lambda[i].setZero(J, J);

        for (int j = 0; j < J; ++j) {
            boot_results.corr_ci_low_by_lambda[i](j, j) = 1.0;
            boot_results.corr_ci_high_by_lambda[i](j, j) = 1.0;

            for (int k = j + 1; k < J; ++k) {
                const int row_jk = j + k * J;
                const auto [ci_low, ci_high] = fisher_z_corr_ci_(
                    boot_results.corr_boot_by_lambda[i],
                    row_jk,
                    B_eff,
                    alpha_low,
                    alpha_high
                );

                boot_results.corr_ci_low_by_lambda[i](j, k) = ci_low;
                boot_results.corr_ci_low_by_lambda[i](k, j) = ci_low;
                boot_results.corr_ci_high_by_lambda[i](j, k) = ci_high;
                boot_results.corr_ci_high_by_lambda[i](k, j) = ci_high;
            }
        }
    }
}

// applies the same row index to all block views
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::set_row_index_all_(
    const typename RGCCA<SamplingStrategy>::BlockRefList& blocks,
    const rgcca::IndexVector& idx
) {
    for (auto* b : blocks) b->set_row_index(idx);
}

// applies independent permutation row indices to all block views
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::set_permuted_row_index_all_(
    const typename RGCCA<SamplingStrategy>::BlockRefList& blocks,
    const unsigned seed
) {
    for (int j = 0; j < static_cast<int>(blocks.size()); ++j) {
        std::mt19937_64 rng(seed + static_cast<unsigned>(104729 * (j + 1)));
        blocks[j]->set_row_index(permutation_indices_(blocks[j]->n_raw(), rng));
    }
}

// clears row-index views on all blocks
template <typename SamplingStrategy>
void RGCCA<SamplingStrategy>::clear_row_index_all_(
    const typename RGCCA<SamplingStrategy>::BlockRefList& blocks
) {
    for (auto* b : blocks) b->clear_row_index();
}

// builds one bootstrap row index vector from a seed
template <typename SamplingStrategy>
auto RGCCA<SamplingStrategy>::bootstrap_index_(
    int n,
    unsigned seed
) const -> rgcca::IndexVector {
    std::mt19937_64 rng(seed);
    return bootstrap_indices_(n, rng);
}

// dispatches to the configured bootstrap resampling scheme
template <typename SamplingStrategy>
auto RGCCA<SamplingStrategy>::bootstrap_indices_(
    int n,
    std::mt19937_64& rng
) const -> rgcca::IndexVector {
    switch (bootstrap_config_.resampling_strategy) {
        case rgcca::ResamplingStrategy::Ordinary:
            return ordinary_bootstrap_indices_(n, rng);
        case rgcca::ResamplingStrategy::Stationary:
            return stationary_bootstrap_indices_(n, bootstrap_config_.stationary_block_length, rng);
    }

    throw std::logic_error("unsupported resampling strategy");
}

// samples ordinary bootstrap indices with replacement
template <typename SamplingStrategy>
auto RGCCA<SamplingStrategy>::ordinary_bootstrap_indices_(
    const int n,
    std::mt19937_64& rng
) const -> rgcca::IndexVector {
    if (n <= 0)
        throw std::invalid_argument("n must be positive");

    std::uniform_int_distribution<int> U(0, n - 1);

    rgcca::IndexVector idx(n);
    for (int i = 0; i < n; ++i)
        idx(i) = U(rng);

    return idx;
}

// samples a random permutation of row indices
template <typename SamplingStrategy>
auto RGCCA<SamplingStrategy>::permutation_indices_(
    const int n,
    std::mt19937_64& rng
) const -> rgcca::IndexVector {
    if (n <= 0)
        throw std::invalid_argument("n must be positive");

    rgcca::IndexVector idx(n);
    std::iota(idx.data(), idx.data() + idx.size(), 0);
    std::shuffle(idx.data(), idx.data() + idx.size(), rng);

    return idx;
}

// samples stationary bootstrap indices with geometric block lengths
template <typename SamplingStrategy>
auto RGCCA<SamplingStrategy>::stationary_bootstrap_indices_(
    const int n,
    const double mean_block_length,
    std::mt19937_64& rng
) const -> rgcca::IndexVector {
    if (n <= 0)
        throw std::invalid_argument("n must be positive");

    if (!(mean_block_length > 0.0) || !std::isfinite(mean_block_length))
        throw std::invalid_argument("stationary block length must be positive");

    const double p = std::clamp(1.0 / mean_block_length, 0.0, 1.0);

    std::uniform_int_distribution<int> U_index(0, n - 1);
    std::bernoulli_distribution start_new_block(p);

    rgcca::IndexVector idx(n);

    int current = U_index(rng);
    idx(0) = current;

    for (int i = 1; i < n; ++i) {
        if (start_new_block(rng)) {
            current = U_index(rng);
        } else {
            current = (current + 1) % n;
        }

        idx(i) = current;
    }

    return idx;
}

} // namespace fdapde

#endif // __FDAPDE_RGCCA_BOOTSTRAP_UTILS_H__

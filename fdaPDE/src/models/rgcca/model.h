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


#ifndef __FDAPDE_RGCCA_MODEL_H__
#define __FDAPDE_RGCCA_MODEL_H__

#include "blocks.h"
#include "logging.h"
#include "statistics.h"
#include "fdaPDE/execution.h"

#include <algorithm>
#include <condition_variable>
#include <map>
#include <numeric>

namespace fdapde {

// rgcca owns block wiring, component fitting, bootstrap selection, and result storage
template <typename SamplingStrategy>
class RGCCA {
public:
    using Block = rgcca::internals::BaseBlock<SamplingStrategy>;
    using BlockPtr = std::unique_ptr<Block>;
    using BlockRefList = std::vector<Block*>;
    using BlockOwnerList = std::vector<std::unique_ptr<Block>>;
    using Matrix = rgcca::Matrix;
    using BoolMatrix = rgcca::BoolMatrix;
    using SparseMatrix = rgcca::SparseMatrix;
    using Vector = rgcca::Vector;
    using IndexVector = rgcca::IndexVector;
    using IndependentSampling = rgcca::IndependentSampling;
    using TimeDependentSampling = rgcca::TimeDependentSampling;
    using InitStrategy = rgcca::InitStrategy;
    using DesignMode = rgcca::DesignMode;
    using LambdaSelection = rgcca::LambdaSelection;
    using Mode = rgcca::Mode;
    using Deflation = rgcca::Deflation;
    using WeightSignConstraint = rgcca::WeightSignConstraint;
    using ResamplingStrategy = rgcca::ResamplingStrategy;
    using SignificanceStatus = rgcca::SignificanceStatus;
    using ComponentStatus = rgcca::ComponentStatus;
    using Scheme = rgcca::Scheme;
    using Options = rgcca::Options;
    using BootstrapConfig = rgcca::BootstrapConfig;
    using Result = rgcca::Result;
    using BootstrapResult = rgcca::BootstrapResult;
    using SamplingDomain = std::conditional_t<std::same_as<SamplingStrategy, TimeDependentSampling>, Triangulation<1, 1>, internals::empty_t>;
    using ComponentCallback = std::function<void(RGCCA&, const Result&)>;

    // constructors
    template <typename S = SamplingStrategy>
    requires std::same_as<S, IndependentSampling>
    explicit RGCCA(const int n, const Options& opt = Options(), const int n_comp = 1) : n_(n), opt_(opt), n_comp_(n_comp) {}
    template <typename S = SamplingStrategy>
    requires std::same_as<S, TimeDependentSampling>
    explicit RGCCA(const int n, const Triangulation<1, 1>& T, const Options& opt = Options(), const int n_comp = 1) : n_(n), T_(T), opt_(opt), n_comp_(n_comp) {}

    // blocks management
    int add_block(BlockPtr b);
    template<typename S = SamplingStrategy>
    requires std::same_as<S, IndependentSampling>
    int add_multivariate_block(std::string block_name, Matrix&& X);
    template<typename S = SamplingStrategy>
    requires std::same_as<S, TimeDependentSampling>
    int add_multivariate_block(std::string block_name, const Vector& times, Matrix&& X);
    template <typename GeoFrame, typename WeightsPenaltyType>
    requires std::same_as<SamplingStrategy, IndependentSampling>
    int add_functional_block(std::string block_name, const GeoFrame& gf, Matrix&& X, WeightsPenaltyType&& weights_penalty);
    template <typename GeoFrame, typename WeightsPenaltyType>
    requires std::same_as<SamplingStrategy, TimeDependentSampling>
    int add_functional_block(std::string block_name, const Vector& times, const GeoFrame& gf, Matrix&& X, WeightsPenaltyType&& weights_penalty);

    // design matrix
    void connect(int j, int k, bool on = true) {
        ensure_design_initialized_();

        validate_index_(j);
        validate_index_(k);

        if (j == k) return;

        C_(j, k) = on;
        C_(k, j) = on;

        design_mode_ = DesignMode::Custom;
    }

    // weights and components regularization utilities
    void set_lambda_weights_all(const double lambda) const {
        rgcca::internals::validate_positive_regularization_lambda(lambda, "weight lambda");
        for (auto& b : blocks_) b->set_lambda_weights(lambda);
    }
    void set_lambda_components_all(const double lambda) const {
        rgcca::internals::validate_positive_regularization_lambda(lambda, "component lambda");
        for (auto& b : blocks_) b->set_lambda_components(lambda);
    }
    void set_lambda_grid_weights(const std::vector<double>& lambda_grid) {
        std::vector<std::vector<double>> normalized(n_comp(), lambda_grid);
        validate_lambda_grid_weights_(normalized);
        lambda_grid_weights_ = std::move(normalized);
    }
    void set_lambda_grid_weights(const std::vector<std::vector<double>>& lambda_grid) {
        std::vector<std::vector<double>> normalized;
        if (static_cast<int>(lambda_grid.size()) == 1) {
            normalized.assign(n_comp(), lambda_grid.front());
        } else {
            normalized = lambda_grid;
        }
        validate_lambda_grid_weights_(normalized);
        lambda_grid_weights_ = std::move(normalized);
    }

    // setters
    void set_n_comp(const int n_comp) {
        auto blocks = main_blocks_();
        set_n_comp_(blocks, n_comp);
    }
    void set_bootstrap_config(const BootstrapConfig bootstrap_config) {
        bootstrap_config_ = bootstrap_config;
        if (!bootstrap_config_.adaptive) {
            bootstrap_config_.B_max = bootstrap_config_.B_min;
            bootstrap_config_.check_every = bootstrap_config_.B_min;
            bootstrap_config_.check_every_block_deactivation = bootstrap_config_.B_min;
            bootstrap_config_.check_every_connection_deactivation = bootstrap_config_.B_min;
        }
    }

    // initialization
    void init() {

        // init design matrix
        ensure_design_initialized_();
        if (design_mode_ == DesignMode::Empty) {
            set_fully_connected_design_();
            design_mode_ = DesignMode::FullyConnected;
        }

        // init components evaluation matrix (only in the TimeDependentSampling scenario)
        if constexpr (std::same_as<SamplingStrategy, TimeDependentSampling>) compute_Psi_T_();

        // set the model as initialized
        initialized_ = true;
    }

    // fit
    std::vector<Result> fit(ComponentCallback component_callback = {}) {

        // initialization
        if (!initialized_) init();

        // validation and logging
        validate_fit_();
        const bool run_model_selection = bootstrap_model_selection_requested_();
        log_fit_header_(run_model_selection);

        // room for results
        std::vector<Result> results;
        results.reserve(n_comp());
        bootstrap_selection_results_.clear();
        bootstrap_selection_results_.reserve(n_comp());
        const auto initial_block_variance = block_variance_trace_(main_blocks_());
        auto current_block_variance = initial_block_variance;

        // components loop
        n_comp_effective_ = 0;
        n_comp_attempted_ = 0;

        for (int hh = 0; hh < n_comp(); ++hh) {
            set_h_(hh);

            BoolMatrix C_active = C_;

            // bootstrap model selection
            if (run_model_selection && count_active_connections_(C_active) > 0) {
                auto selection = bootstrap_model_selection_(C_active);
                C_active = std::move(selection.C_active);
                if (selection.lambda_selected)
                    set_lambda_weights_all(selection.lambda);
            }

            // final fit
            auto step_start = log_step_start_("Final component fit");
            auto blocks = main_blocks_();
            const auto active_blocks = active_blocks_from_C_(C_active);
            init_comp_(blocks, InitStrategy::None, true, &active_blocks);
            Result component_result = fit_component_(blocks, C_active);
            log_step_end_(step_start);
            const auto block_variance_before = current_block_variance;
            const auto observed_correlation = rho_tot_maxvar_(blocks, C_active);
            component_result.rho_tot = observed_correlation.normalized;
            component_result.rho_tot_raw = observed_correlation.raw;
            component_result.inner_ave = observed_correlation.inner_ave;

            // structural stop: no active design left
            if (count_active_connections_(C_active) == 0) {
                component_result.status = ComponentStatus::RejectedInactiveDesign;
                annotate_component_significance_(component_result, inactive_component_significance_());
                annotate_explained_variance_(
                    component_result, initial_block_variance, block_variance_before, block_variance_before
                );
                results.push_back(std::move(component_result));
                finish_component_(results.back(), component_callback);
                break;
            }

            // bootstrap component significance
            if (opt_.component_significance) {
                const auto significance = bootstrap_test_component_significance_(C_active, observed_correlation);
                annotate_component_significance_(component_result, significance);
                if (significance.status != SignificanceStatus::Significant) {
                    component_result.status =
                        significance.status == SignificanceStatus::NotSignificant ?
                            ComponentStatus::RejectedNotSignificant :
                            ComponentStatus::RejectedSignificanceUnavailable;
                    annotate_explained_variance_(
                        component_result, initial_block_variance, block_variance_before, block_variance_before
                    );
                    results.push_back(std::move(component_result));
                    finish_component_(results.back(), component_callback);
                    break;
                }
            }

            component_result.status = ComponentStatus::Retained;

            // bootstrap block importance
            if (opt_.block_importance) {
                const auto importance = bootstrap_test_block_importance_(C_active);
                annotate_block_importance_(component_result, importance);
            }

            // store results
            results.push_back(std::move(component_result));

            // deflation
            step_start = log_step_start_("Deflate blocks");
            deflate_all_();
            log_step_end_(step_start);
            current_block_variance = block_variance_trace_(blocks);
            annotate_explained_variance_(
                results.back(), initial_block_variance, block_variance_before, current_block_variance
            );

            // component post-processing
            finish_component_(results.back(), component_callback);
        }

        return results;
    }

    // getters
    void get_tau(const BlockRefList& blocks, std::vector<double>& tau_values) const {
        for (int j = 0; j < n_blocks(); ++j)
            tau_values[j] = blocks[j]->tau();
    }
    void get_lambdas(const BlockRefList& blocks, std::vector<double> & lambda_components_values, std::vector<double> & lambda_weights_values) const {
        for (int j = 0; j < n_blocks(); ++j) {
            lambda_components_values[j] = blocks[j]->lambda_components();
            lambda_weights_values[j] = blocks[j]->lambda_weights();
        }
    }

    // observers
    [[nodiscard]] int n() const { return n_; }
    [[nodiscard]] int n_comp() const { return n_comp_; }
    [[nodiscard]] int n_comp_effective() const { return n_comp_effective_; }
    [[nodiscard]] int n_comp_attempted() const { return n_comp_attempted_; }
    [[nodiscard]] int n_blocks() const { return J_; }
    [[nodiscard]] const Options& options() const { return opt_; }
    [[nodiscard]] const Scheme& scheme() const { return opt_.scheme; }
    [[nodiscard]] const std::vector<BlockPtr>& blocks() const { return blocks_; }
    [[nodiscard]] const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& C() const { return C_; }
    template <typename S = SamplingStrategy>
    requires std::same_as<S, TimeDependentSampling>
    [[nodiscard]] const SparseMatrix& Psi_T() const { return Psi_T_; };

    // bootstrap results
    [[nodiscard]] const std::vector<BootstrapResult>& bootstrap_selection_results() const { return bootstrap_selection_results_; }
    void clear_bootstrap_selection_results() {
        bootstrap_selection_results_.clear();
        bootstrap_selection_results_.shrink_to_fit();
    }

    // bootstrap CI
    [[nodiscard]] std::pair<Vector, Vector> bootstrap_weights_ci(const int h, const int lambda_i, const int block_j, const SparseMatrix& Psi) const;
    [[nodiscard]] std::pair<Vector, Vector> bootstrap_weights_ci(const int h, const int block_j, const SparseMatrix& Psi) const;
    template <typename DataLocs>
    requires(!std::same_as<std::decay_t<DataLocs>, SparseMatrix>)
    [[nodiscard]] std::pair<Vector, Vector> bootstrap_weights_ci(const int h, const int lambda_i, const int block_j, const DataLocs& locs) const;
    template <typename DataLocs>
    requires(!std::same_as<std::decay_t<DataLocs>, SparseMatrix>)
    [[nodiscard]] std::pair<Vector, Vector> bootstrap_weights_ci(const int h, const int block_j, const DataLocs& locs) const;

private:
    // private state types
    struct FitWorkspace {
        Matrix Cov;
        Eigen::ArrayXXi dirty;
        Vector means;

        explicit FitWorkspace(int n_blocks) {
            Cov.setZero(n_blocks, n_blocks);
            dirty.setOnes(n_blocks, n_blocks);
            means.setZero(n_blocks);
            for (int j = 0; j < n_blocks; ++j) {
                Cov(j, j) = 1.0;
                dirty(j, j) = 0;
            }
        }
    };

    // bootstrap worker types
    struct BootstrapBlocks {
        BlockOwnerList owners;
        BlockRefList refs;
    };

    // bootstrap state types
    struct AdaptiveBootstrapState {

        explicit AdaptiveBootstrapState(
            const BootstrapConfig bootstrap_config,
            const int n_threads_,
            const int h,
            const int n_blocks_
        ) : n_threads(n_threads_), n_blocks(n_blocks_) {
            seed = bootstrap_config.seed + static_cast<unsigned>(h);
            B_min = bootstrap_config.B_min;
            B_max = bootstrap_config.adaptive ? bootstrap_config.B_max : B_min;
            check_every = bootstrap_config.adaptive ? bootstrap_config.check_every : B_min;
            check_every_block_deactivation =
                bootstrap_config.adaptive ? bootstrap_config.check_every_block_deactivation : B_min;
            check_every_connection_deactivation =
                bootstrap_config.adaptive ? bootstrap_config.check_every_connection_deactivation : B_min;
            corr_pos_count.setZero(n_blocks, n_blocks);
            corr_neg_count.setZero(n_blocks, n_blocks);
        }

        void reset() {
            B_done = 0;
            B_total = 0;
            B_design = 0;
            B_stale = 0;
            B_cancelled = 0;
            B_final_capped = 0;
            design_epoch = 0;
            last_check_B_done = 0;
            last_block_deactivation_check_B_done = 0;
            last_connection_deactivation_check_B_done = 0;
            stop = false;
            reset_good();
        }
        void reset_good() {
            B_done = 0;
            stable_checks = 0;
            crit_prev_check = std::numeric_limits<double>::infinity();
            crit = std::numeric_limits<double>::quiet_NaN();
            last_check_B_done = 0;
            last_block_deactivation_check_B_done = 0;
            last_connection_deactivation_check_B_done = 0;
            corr_pos_count.setZero(n_blocks, n_blocks);
            corr_neg_count.setZero(n_blocks, n_blocks);
        }

        // config
        int n_threads;
        int n_blocks;
        int seed;
        int B_min;
        int B_max;
        int check_every;
        int check_every_block_deactivation;
        int check_every_connection_deactivation;

        // state
        int B_done = 0;
        int B_total = 0;
        int B_design = 0;
        int B_stale = 0;
        int B_cancelled = 0;
        int B_final_capped = 0;
        int design_epoch = 0;
        int last_check_B_done = 0;
        int last_block_deactivation_check_B_done = 0;
        int last_connection_deactivation_check_B_done = 0;
        bool stop = false;

        Eigen::MatrixXi corr_pos_count;
        Eigen::MatrixXi corr_neg_count;

        // adaptive checks
        int stable_checks = 0;
        double crit_prev_check = std::numeric_limits<double>::infinity();
        double crit = std::numeric_limits<double>::quiet_NaN();

        // early stop
        double best_criterion = -std::numeric_limits<double>::infinity();
        int best_i = -1;
        int no_improve = 0;
    };

    struct ModelSelectionResult {
        bool lambda_selected = false;
        double lambda = std::numeric_limits<double>::quiet_NaN();
        BoolMatrix C_active;
    };
    struct AdaptiveStopInfo {
        bool stop = false;
        double rel_change = std::numeric_limits<double>::quiet_NaN();
        const char* reason = "";
    };

    // bootstrap result types
    struct BootstrapTimingSummary {
        double fit_time = 0.0;
        double fit_time_sq = 0.0;
        double worker_time = 0.0;
        double claim_wait_time = 0.0;
        double merge_wait_time = 0.0;
        double parallel_capacity = 0.0;
        double iters = 0.0;
        int capped_fits = 0;
        int max_iters = 0;
        int n_fits = 0;
        ::fdapde::internals::NonNegativeWeightSolveStats nn_stats;

        void add_sample(
            const double fit_time_sec,
            const int fit_iters,
            const bool capped,
            const bool fit_started,
            const ::fdapde::internals::NonNegativeWeightSolveStats& sample_nn_stats
        ) {
            if (!fit_started) return;
            fit_time += fit_time_sec;
            fit_time_sq += fit_time_sec * fit_time_sec;
            iters += static_cast<double>(fit_iters);
            capped_fits += capped ? 1 : 0;
            max_iters = std::max(max_iters, fit_iters);
            nn_stats += sample_nn_stats;
            ++n_fits;
        }
        void set_parallel_capacity(const double wall_time_sec, const int n_threads) {
            parallel_capacity = wall_time_sec * static_cast<double>(n_threads);
        }
        void add_worker_time(const double worker_time_sec) {
            worker_time += worker_time_sec;
        }
        void set_scheduler_wait_times(const double claim_wait_sec, const double merge_wait_sec) {
            claim_wait_time = claim_wait_sec;
            merge_wait_time = merge_wait_sec;
        }
        [[nodiscard]] double efficiency() const {
            return parallel_capacity > 0.0 ? fit_time / parallel_capacity : 0.0;
        }
        [[nodiscard]] double fraction_of_capacity(const double time_sec) const {
            return parallel_capacity > 0.0 ? time_sec / parallel_capacity : 0.0;
        }
        [[nodiscard]] double worker_overhead_fraction() const {
            return fraction_of_capacity(std::max(0.0, worker_time - fit_time));
        }
        [[nodiscard]] double coordinator_fraction() const {
            return std::max(
                0.0,
                1.0 - fraction_of_capacity(worker_time + claim_wait_time + merge_wait_time)
            );
        }
        [[nodiscard]] double avg_fit_time() const {
            return n_fits > 0 ? fit_time / static_cast<double>(n_fits) : 0.0;
        }
        [[nodiscard]] double sd_fit_time() const {
            if (n_fits < 2) return 0.0;
            const double mean = avg_fit_time();
            const double var = (fit_time_sq - static_cast<double>(n_fits) * mean * mean) /
                static_cast<double>(n_fits - 1);
            return std::sqrt(std::max(0.0, var));
        }
        [[nodiscard]] double avg_iters() const {
            return n_fits > 0 ? iters / static_cast<double>(n_fits) : 0.0;
        }
    };
    struct MaxvarCorrelationResult {
        double raw = std::numeric_limits<double>::quiet_NaN();
        double normalized = std::numeric_limits<double>::quiet_NaN();
        double inner_ave = std::numeric_limits<double>::quiet_NaN();
    };
    struct ComponentSignificanceResult {
        double rho_tot = std::numeric_limits<double>::quiet_NaN();
        double rho_tot_raw = std::numeric_limits<double>::quiet_NaN();
        double p_value = std::numeric_limits<double>::quiet_NaN();
        double null_mean = std::numeric_limits<double>::quiet_NaN();
        double null_q95 = std::numeric_limits<double>::quiet_NaN();
        double null_max = std::numeric_limits<double>::quiet_NaN();
        int B = 0;
        int null_valid_count = 0;
        SignificanceStatus status = SignificanceStatus::NotTested;
    };
    struct BlockImportanceResult {
        std::vector<double> rho;
        std::vector<double> p_value;
        std::vector<bool> significant;
        int B = 0;
    };
    struct BootstrapSampleResult {
        std::vector<Vector> w;
        Matrix corr;
        double fit_time = 0.0;
        int fit_iters = 0;
        bool capped = false;
        bool cancelled = false;
        bool fit_started = false;
        ::fdapde::internals::NonNegativeWeightSolveStats nn_stats;
    };

    // component initialization
    void init_comp_(
        const BlockRefList& blocks,
        InitStrategy init_strategy = InitStrategy::None,
        const bool update_regularization = true,
        const std::vector<bool>* active_blocks = nullptr
    ) {
        if (update_regularization) {
            if (opt_.mode == Mode::Regularized) set_tau_auto_all_(blocks);
            if (opt_.lambda_selection_components == LambdaSelection::Automatic) set_lambda_components_auto_all_(blocks);
        }

        if (active_blocks != nullptr && active_blocks->size() != blocks.size())
            throw std::logic_error("init_comp_: active block mask size mismatch");

        if (init_strategy == InitStrategy::None) init_strategy = opt_.init_strategy;

        // blocks loop
        for (int j = 0; j < n_blocks(); ++j) {
            auto* b = blocks[j];
            b->set_h(h_);
            if (active_blocks != nullptr && !(*active_blocks)[j])
                continue;

            if (init_strategy == InitStrategy::WarmStart) {
                b->refresh_component();
            } else {
                b->init_weight_uniform();
                switch (init_strategy) {
                    case InitStrategy::Uniform: {
                        const auto info = b->uniform_init();
                        b->compute(info.nu);
                        break;
                    }
                    case InitStrategy::SVD: {
                        const auto info = b->svd_init();
                        b->compute(info.nu);
                        break;
                    }
                    default: throw std::logic_error("unsupported initialization strategy");
                }
            }
        }
    }
    void init_comp_() {
        auto blocks = main_blocks_();
        init_comp_(blocks);
    }

    // components fit
    Result fit_component_(
        const BlockRefList& blocks,
        const BoolMatrix& C_active,
        const bool update_component_lambdas = true,
        const int max_iter_override = -1,
        const std::function<bool()>& cancelled = {},
        const bool collect_summary = true
    ) {
        const int max_iter = max_iter_override > 0 ? max_iter_override : opt_.max_iter;
        FitWorkspace ws(n_blocks());

        // room for results
        Result res(n_blocks());
        res.h = h_;
        res.obj_history.reserve(max_iter);

        // design update according to current active blocks
        res.C = C_active;
        res.active_blocks = active_blocks_from_C_(C_active);
        for (int j = 0; j < n_blocks(); ++j) {
            if (!res.active_blocks[j]){
                blocks[j]->weights().col(h_).setZero();
                blocks[j]->components().col(h_).setZero();
            }
        }

        // initialization
        std::vector<Vector> eta_cache = eta_(blocks);
        for (int j = 0; j < n_blocks(); ++j)
            ws.means[j] = eta_cache[j].mean();
        res.obj_history.push_back(objective_(ws, res.C, eta_cache));
        auto w_prev = snapshot_weights_(blocks);
        if (update_component_lambdas && opt_.lambda_selection_components == LambdaSelection::Automatic)
            set_lambda_components_auto_all_(blocks);

        // main loop
        for (int s = 0; s < max_iter; ++s) {
            for (int l = 0; l < n_blocks(); ++l) {

                // allow callers to stop long fits between block updates
                if (cancelled && cancelled()) {
                    res.cancelled = true;
                    return res;
                }

                // skip deactivated blocks
                if (!res.active_blocks[l]) continue;

                // inner-component assembler
                Vector nu_l = Vector::Zero(blocks[l]->n());
                const Vector& eta_l = eta_cache[l];
                for (int k = 0; k < n_blocks(); ++k) {
                    if (!res.C(l, k)) continue;
                    const Vector& eta_k = eta_cache[k];
                    const double cov_lk = cov_value_(ws, l, k, eta_l, eta_k);
                    const double w_lk = opt_.scheme.w(cov_lk);
                    if constexpr (std::same_as<SamplingStrategy, IndependentSampling>) {
                        nu_l.noalias() += w_lk * eta_k;
                    } else {
                        nu_l.noalias() += w_lk * eta_(*blocks[k], *blocks[l]);
                    }
                }

                // block update
                blocks[l]->compute(nu_l);
                eta_cache[l] = eta_(*blocks[l]);
                ws.means[l] = eta_cache[l].mean();
                mark_cov_rowcol_dirty_(ws, l);
            }

            // update metrics
            const double f_obj = objective_(ws, res.C, eta_cache);
            const double obj_prev = res.obj_history.back();
            res.obj_history.push_back(f_obj);
            res.iters = s + 1;

            // check monotonicity
            const double rel_obj_drop = (obj_prev - f_obj) / (1.0 + std::abs(obj_prev));
            if (rel_obj_drop > opt_.tol)
                res.monotone = false;

            // stopping criteria
            const double delta_obj = std::abs(f_obj - obj_prev);
            const double delta_w = weights_variation_(blocks, w_prev);
            if (delta_obj < opt_.tol || delta_w < opt_.tol)
                break;

            // save weights snapshot
            w_prev = snapshot_weights_(blocks);
        }

        // allow callers to stop long fits before cov and cor matrices computation
        if (cancelled && cancelled()) {
            res.cancelled = true;
            return res;
        }

        if (!collect_summary) return res;

        // save results
        covariance_matrix_(blocks, res.covariance_matrix);
        correlation_matrix_(blocks, res.correlation_matrix);
        get_tau(blocks, res.tau_values);
        get_lambdas(blocks, res.lambda_components_values, res.lambda_weights_values);

        return res;
    }
    Result fit_component_(const BlockRefList& blocks) {
        return fit_component_(blocks, C_);
    }
    Result fit_component_(const BoolMatrix& C_active) {
        auto blocks = main_blocks_();
        return fit_component_(blocks, C_active);
    }
    Result fit_component_() {
        auto blocks = main_blocks_();
        return fit_component_(blocks, C_);
    }

    // bootstrap model selection
    ModelSelectionResult bootstrap_model_selection_(const BoolMatrix& C_initial) {
        log_bootstrap_model_selection_header_();

        // lambda selection
        const bool select_lambda = weight_lambda_selection_requested_();
        const std::vector<double> lambda_grid = model_selection_lambda_grid_();

        // set the number of threads
        const int n_threads = bootstrap_n_threads_();

        // original blocks
        auto blocks = main_blocks_();
        BoolMatrix C_active = C_initial;
        BoolMatrix C_best = C_initial;

        // init bootstrap
        auto step_start = log_step_start_("Init bootstrap");
        AdaptiveBootstrapState bootstrap_state(bootstrap_config_, n_threads, h_, n_blocks());
        const auto block_dims = block_dims_(blocks);
        BootstrapResult boot_results(
            h_, bootstrap_state.B_max, lambda_grid,
            block_names_(blocks), block_dims, bootstrap_config_.ci_level
        );
        log_step_end_(step_start);

        // preliminary fit
        step_start = log_step_start_("Preliminary fit");
        if (select_lambda) set_lambda_weights_all(lambda_grid.back());
        const auto preliminary_active_blocks = active_blocks_from_C_(C_active);
        init_comp_(blocks, InitStrategy::None, true, &preliminary_active_blocks);
        fit_component_(blocks, C_active);
        auto preliminary_w_fit = snapshot_weights_(blocks);
        log_step_end_(step_start);

        // clone blocks for bootstrap workers
        step_start = log_step_start_("Clone worker blocks");
        std::vector<BootstrapBlocks> thread_boot_worker(n_threads);
        for (int t = 0; t < n_threads; ++t) {
            thread_boot_worker[t] = clone_blocks_();
        }
        log_step_end_(step_start);

        // model selection loop
        int n_lambda = static_cast<int>(lambda_grid.size());
        for (int lambda_i = n_lambda - 1; lambda_i >= 0; --lambda_i) {
            ensure_bootstrap_lambda_storage_(boot_results, lambda_i, block_dims, n_blocks());
            const bool reuse_preliminary_fit = lambda_i == n_lambda - 1;

            // set the current lambda (if needed)
            if (select_lambda) {
                const double lambda = lambda_grid[lambda_i];
                log_bootstrap_lambda_candidate_(select_lambda, lambda);
                if (!reuse_preliminary_fit) {
                    set_lambda_weights_all(lambda);
                    for (auto& worker : thread_boot_worker)
                        set_lambda_weights_all_(worker.refs, lambda);
                }
            } else {
                log_bootstrap_lambda_candidate_(select_lambda, fixed_weight_lambda_(blocks));
            }

            // init warm-start at lambda
            std::vector<Vector> w_fit;
            if (reuse_preliminary_fit) {
                step_start = log_step_start_("  Warm-start fit (reuse preliminary)");
                w_fit = preliminary_w_fit;
            } else {
                step_start = log_step_start_("  Warm-start fit");
                const auto active_blocks = active_blocks_from_C_(C_active);
                init_comp_(blocks, InitStrategy::WarmStart, true, &active_blocks);
                fit_component_(blocks, C_active);
                w_fit = snapshot_weights_(blocks);
            }
            auto w_min = w_fit;
            if (opt_.block_deactivation)
                threshold_inactive_blocks_(w_min, C_active);
            log_step_end_(step_start);

            // start bootstrap timer
            auto start = std::chrono::high_resolution_clock::now();

            // streaming bootstrap
            bootstrap_state.reset();
            BootstrapTimingSummary bootstrap_timing_summary;
            run_bootstrap_stream_(
                lambda_i,
                bootstrap_state,
                thread_boot_worker,
                C_active,
                w_fit,
                w_min,
                boot_results,
                blocks,
                bootstrap_timing_summary
            );
            BoolMatrix C_lambda = C_active;
            int n_active_connections = count_active_connections_(C_lambda);
            const int n_active_blocks = count_active_blocks_(C_lambda);

            // end bootstrap timer
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            const double elapsed_sec = duration.count() / 1000.0;

            // save results
            boot_results.w_fit_by_lambda[lambda_i] = w_fit;
            boot_results.w_min_by_lambda[lambda_i] = w_min;
            step_start = log_step_start_("  Final lambda correlation");
            correlation_matrix_(blocks, w_min, boot_results.corr_min_by_lambda[lambda_i]);
            boot_results.corr_min_by_lambda[lambda_i].array() *=
                (C_lambda.cast<double>() + Matrix::Identity(n_blocks(), n_blocks())).array();
            log_step_end_(step_start);
            if (!std::isfinite(bootstrap_state.crit))
                bootstrap_state.crit = criterion_score_with_weights_(blocks, w_min, C_initial);
            boot_results.criterion[lambda_i] = bootstrap_state.crit;
            boot_results.B_used_by_lambda[lambda_i] = bootstrap_state.B_done;
            boot_results.B_total_by_lambda[lambda_i] = bootstrap_state.B_total;
            boot_results.B_design_by_lambda[lambda_i] = bootstrap_state.B_design;
            boot_results.B_stale_by_lambda[lambda_i] = bootstrap_state.B_stale;
            boot_results.B_cancelled_by_lambda[lambda_i] = bootstrap_state.B_cancelled;
            boot_results.B_final_capped_by_lambda[lambda_i] = bootstrap_state.B_final_capped;
            boot_results.design_epochs_by_lambda[lambda_i] = bootstrap_state.design_epoch + 1;

            // Release unused B_max columns before allocating the next lambda.
            resize_bootstrap_lambda_results_(boot_results, lambda_i, n_blocks(), block_dims);

            // store current design if better
            if (bootstrap_state.crit > bootstrap_state.best_criterion) {
                C_best = C_lambda;
            }

            // log bootstrap summary
            log_bootstrap_lambda_summary_(
                bootstrap_state,
                bootstrap_timing_summary,
                elapsed_sec,
                n_active_blocks,
                n_active_connections,
                select_lambda ? lambda_grid[lambda_i] : fixed_weight_lambda_(blocks)
            );

            // early stop
            if (early_stop_lambda_(bootstrap_state, lambda_i, bootstrap_config_)) {
                break;
            }

            // reset connections, but keep fully deactivated blocks off
            C_active = reset_connections_keep_inactive_blocks_(C_initial, C_active);

        }

        // compute correlation matrices CI
        step_start = log_step_start_("Compute bootstrap correlation CIs");
        compute_bootstrap_corr_cis_(boot_results, n_blocks());
        log_step_end_(step_start);

        // save optimal design
        if (bootstrap_state.best_i < 0)
            throw std::runtime_error("No model candidate was evaluated during bootstrap selection");
        boot_results.lambda_opt_index = bootstrap_state.best_i;
        if (select_lambda)
            boot_results.lambda_opt = lambda_grid[bootstrap_state.best_i];
        boot_results.active_blocks = active_blocks_from_C_(C_best);
        log_bootstrap_selection_result_(boot_results, C_best, select_lambda);

        // return model selection results
        bootstrap_selection_results_.push_back(std::move(boot_results));
        ModelSelectionResult out;
        out.lambda_selected = select_lambda;
        out.lambda = bootstrap_selection_results_.back().lambda_opt;
        out.C_active = C_best;
        return out;
    }

    // bootstrap component significance
    ComponentSignificanceResult bootstrap_test_component_significance_(
        const BoolMatrix& C_active,
        const MaxvarCorrelationResult& observed
    ) {
        log_bootstrap_component_significance_header_();
        ComponentSignificanceResult out;

        auto blocks = main_blocks_();
        out.rho_tot = observed.normalized;
        out.rho_tot_raw = observed.raw;

        // A degenerate observed statistic cannot support a significance claim.
        if (count_active_connections_(C_active) == 0 || !std::isfinite(out.rho_tot)) {
            out.B = 0;
            log_component_significance_(out);
            return out;
        }

        const int B = bootstrap_config_.component_significance_resamples;
        const int n_threads = bootstrap_n_threads_();
        out.B = B;
        const auto w_fit = snapshot_weights_(blocks);

        // each worker owns a mutable clone with its own row permutation
        auto step_start = log_step_start_("Clone significance worker blocks");
        std::vector<BootstrapBlocks> thread_boot_worker(n_threads);
        for (int t = 0; t < n_threads; ++t)
            thread_boot_worker[t] = clone_blocks_();
        log_step_end_(step_start);

        std::vector<double> rho_null(B, std::numeric_limits<double>::quiet_NaN());
        const unsigned seed = bootstrap_config_.seed + static_cast<unsigned>(1000003 * (h_ + 1));
        const auto active_blocks = active_blocks_from_C_(C_active);

        // null: break cross-block row alignment within each block, then refit
        step_start = log_step_start_("Significance bootstrap");
        parallel_for(0, B, 1, [&](int b) {
            const int tid = this_thread_id();
            auto& boot_blocks = thread_boot_worker[tid];

            for (std::size_t j = 0; j < boot_blocks.refs.size(); ++j)
                boot_blocks.refs[j]->reset_bootstrap_fit_state(w_fit[j]);
            copy_weights_snapshot_(boot_blocks.refs, w_fit);
            set_component_significance_row_index_all_(
                boot_blocks.refs,
                seed + static_cast<unsigned>(7919 * (b + 1))
            );

            init_comp_(boot_blocks.refs, InitStrategy::WarmStart, false, &active_blocks);
            fit_component_(boot_blocks.refs, C_active, false, bootstrap_fit_max_iter_());

            // count null statistics at least as extreme as the observed one
            rho_null[b] = rho_tot_maxvar_(boot_blocks.refs, C_active).normalized;

            clear_row_index_all_(boot_blocks.refs);
        });
        log_step_end_(step_start);

        std::vector<double> valid_null;
        valid_null.reserve(B);
        double null_sum = 0.0;
        int ge_count = 0;
        for (const double rho_star : rho_null) {
            if (!std::isfinite(rho_star)) continue;
            valid_null.push_back(rho_star);
            null_sum += rho_star;
            if (rho_star >= out.rho_tot) ++ge_count;
        }

        out.null_valid_count = static_cast<int>(valid_null.size());
        if (valid_null.empty()) {
            log_component_significance_(out);
            return out;
        }

        out.null_mean = null_sum / static_cast<double>(valid_null.size());
        out.null_max = *std::max_element(valid_null.begin(), valid_null.end());
        out.null_q95 = internals::empirical_quantile(valid_null, 0.95);
        out.p_value = static_cast<double>(ge_count + 1) /
            static_cast<double>(out.null_valid_count + 1);
        out.status = out.p_value <= bootstrap_config_.component_significance_alpha ?
            SignificanceStatus::Significant : SignificanceStatus::NotSignificant;

        log_component_significance_(out);

        return out;
    }

    // bootstrap block importance for the fitted component
    BlockImportanceResult bootstrap_test_block_importance_(const BoolMatrix& C_active) {
        log_bootstrap_block_importance_header_();

        auto blocks = main_blocks_();
        const int B = bootstrap_config_.block_importance_resamples;
        BlockImportanceResult out;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        out.rho.assign(n_blocks(), nan);
        out.p_value.assign(n_blocks(), nan);
        out.significant.assign(n_blocks(), false);
        out.B = B;

        const std::vector<Vector> eta = eta_(blocks);
        Vector v;
        if (!block_importance_weights_(eta, C_active, v)) {
            log_block_importance_(out);
            return out;
        }

        std::vector<Vector> targets(n_blocks());
        std::vector<int> testable(n_blocks(), 0);
        for (int j = 0; j < n_blocks(); ++j) {
            if (v[j] == 0.0) continue;
            targets[j] = block_importance_target_(eta, C_active, v, j);
            const double target_var = cov_(targets[j], targets[j]);
            if (!(target_var > 0.0) || !std::isfinite(target_var)) continue;
            testable[j] = 1;
            out.rho[j] = block_importance_rho_(eta[j], targets[j]);
        }

        std::vector<int> significant(n_blocks(), 0);
        parallel_for(0, n_blocks(), 1, [&](int j) {
            if (!testable[j]) return;

            auto boot_blocks = clone_blocks_();
            auto* block = boot_blocks.refs[j];
            std::mt19937_64 rng(
                bootstrap_config_.seed +
                static_cast<unsigned>(1000003 * (h_ + 1)) +
                static_cast<unsigned>(9176 * (j + 1))
            );

            int ge_count = 0;
            for (int b = 0; b < B; ++b) {
                block->set_row_index(single_block_null_indices_(block->n_raw(), rng));
                block->compute(targets[j]);
                const double rho_star = block_importance_rho_(eta_(*block), targets[j]);
                if (std::isfinite(rho_star) && rho_star >= out.rho[j])
                    ++ge_count;
            }
            block->clear_row_index();

            out.p_value[j] = static_cast<double>(ge_count + 1) / static_cast<double>(B + 1);
            significant[j] = out.p_value[j] <= bootstrap_config_.block_importance_alpha ? 1 : 0;
        });

        for (int j = 0; j < n_blocks(); ++j)
            out.significant[j] = significant[j] != 0;

        log_block_importance_(out);
        return out;
    }

    // initialization utils
    void ensure_design_initialized_() {
        if (C_.rows() == n_blocks() && C_.cols() == n_blocks()) return;

        C_.resize(n_blocks(), n_blocks());
        C_.setConstant(false);
    }
    void clear_design_() {
        C_.resize(n_blocks(), n_blocks());
        C_.setConstant(false);
    }
    void set_fully_connected_design_() {
        clear_design_();
        for (int j = 0; j < n_blocks(); ++j)
            for (int k = 0; k < n_blocks(); ++k)
                C_(j, k) = (j != k);
    }
    template <typename S = SamplingStrategy>
    requires std::same_as<S, TimeDependentSampling>
    void add_times_(const Vector& t) {
        times_.reserve(times_.size() + static_cast<size_t>(t.size()));
        times_.insert(times_.end(), t.data(), t.data() + t.size());
    }
    void compute_Psi_T_() {
        std::ranges::sort(times_);
        times_.erase(std::ranges::unique(times_).begin(), times_.end());
        Eigen::VectorXd times_eig = Eigen::Map<Eigen::VectorXd>(times_.data(), times_.size());
        SamplingStrategy::compute_Psi(T_, Matrix{times_eig}, Psi_T_);
    }

    // weight snapshot helpers
    void copy_weights_snapshot_(const BlockRefList& blocks, const std::vector<Vector>& weights) const {
        if (blocks.size() != weights.size())
            throw std::logic_error("copy_weights_snapshot_: size mismatch");

        for (std::size_t j = 0; j < blocks.size(); ++j) {
            if (blocks[j]->n_dofs_weights() != weights[j].size())
                throw std::logic_error("copy_weights_snapshot_: incompatible weight size");

            blocks[j]->weights().col(h_) = weights[j];
        }
    }
    std::vector<Vector> snapshot_weights_(const BlockRefList& blocks) const {
        std::vector<Vector> out;
        out.reserve(blocks.size());
        for (auto* b : blocks) out.push_back(b->weights().col(h_));
        return out;
    }

    // fit helpers
    std::vector<double> block_variance_trace_(const BlockRefList& blocks) const {
        std::vector<double> out;
        out.reserve(blocks.size());
        for (const auto* block : blocks) out.push_back(block->variance_trace());
        return out;
    }
    void annotate_explained_variance_(
        Result& result,
        const std::vector<double>& initial,
        const std::vector<double>& before,
        const std::vector<double>& after
    ) const {
        if (initial.size() != blocks_.size() || before.size() != blocks_.size() || after.size() != blocks_.size())
            throw std::logic_error("annotate_explained_variance_: block count mismatch");

        result.block_variance_initial = initial;
        result.block_variance_before = before;
        result.block_variance_after = after;
        for (std::size_t j = 0; j < initial.size(); ++j) {
            if (!(initial[j] > 0.0) || !std::isfinite(initial[j])) continue;
            result.block_variance_explained[j] = (before[j] - after[j]) / initial[j];
            result.block_variance_explained_cumulative[j] = (initial[j] - after[j]) / initial[j];
        }
    }
    void deflate_all_() const {
        for (auto& b : blocks_) b->deflate(opt_.deflation_mode);
    }
    void compute_weights_star_(const int h) {
        for (auto& b : blocks_) b->compute_weights_star(h);
    }
    void finish_component_(const Result& result, const ComponentCallback& component_callback) {
        ++n_comp_attempted_;
        if (result.retained()) {
            ++n_comp_effective_;
            auto step_start = log_step_start_("Compute weights_star");
            compute_weights_star_(result.h);
            log_step_end_(step_start);
        }

        if (component_callback) {
            auto step_start = log_step_start_("Component callback");
            component_callback(*this, result);
            log_step_end_(step_start);
        }
    }
    double weights_variation_(const BlockRefList& blocks, const std::vector<Vector>& w_prev) const {
        double acc = 0.0;

        for (int j = 0; j < n_blocks(); ++j) {
            const auto wj = blocks[j]->weights().col(h_);
            acc += (wj - w_prev[j]).squaredNorm();
        }

        return acc;
    }

    // bootstrap request and logging helpers
    bool bootstrap_model_selection_requested_() const {
        return opt_.block_deactivation ||
            opt_.connection_deactivation ||
            opt_.lambda_selection_weights == LambdaSelection::Automatic;
    }
    bool weight_lambda_selection_requested_() const {
        return opt_.lambda_selection_weights == LambdaSelection::Automatic;
    }
    std::vector<double> model_selection_lambda_grid_() const {
        if (weight_lambda_selection_requested_())
            return lambda_grid_weights_[h_];

        return {std::numeric_limits<double>::quiet_NaN()};
    }
    double fixed_weight_lambda_(const BlockRefList& blocks) const {
        for (const auto* block : blocks) {
            const double lambda = block->lambda_weights();
            if (std::isfinite(lambda)) return lambda;
        }
        return std::numeric_limits<double>::quiet_NaN();
    }

    // logging helpers
    std::chrono::high_resolution_clock::time_point log_step_start_(std::string_view label) const;
    void log_step_end_(const std::chrono::high_resolution_clock::time_point start) const;
    void log_weight_lambda_(double lambda) const;
    void log_fit_header_(bool run_model_selection) const;
    void log_bootstrap_model_selection_header_() const;
    void log_bootstrap_component_significance_header_() const;
    void log_bootstrap_block_importance_header_() const;
    void log_bootstrap_lambda_candidate_(bool lambda_selected, double lambda) const;
    void log_component_significance_(const ComponentSignificanceResult& significance) const;
    void log_block_importance_(const BlockImportanceResult& importance) const;
    void log_bootstrap_lambda_summary_(
        const AdaptiveBootstrapState& state,
        const BootstrapTimingSummary& timing_summary,
        double elapsed_sec,
        int n_active_blocks,
        int n_active_connections,
        double lambda
    ) const;
    void log_bootstrap_selection_result_(
        const BootstrapResult& boot_results,
        const BoolMatrix& C_best,
        bool lambda_selected
    ) const;
    void log_bootstrap_early_stop_(int patience) const;

    // bootstrap CI helpers
    const BootstrapResult& bootstrap_selection_result_(const int h) const;
    int bootstrap_lambda_opt_index_(const BootstrapResult& boot_results) const;
    std::pair<Vector, Vector> bootstrap_weights_ci_(const int h, const int lambda_i, const int block_j, const SparseMatrix& Psi) const;

    // bootstrap sample helpers
    BootstrapSampleResult fit_bootstrap_sample_(
        BootstrapBlocks& boot_blocks,
        const int b,
        const unsigned seed,
        const BoolMatrix& C_active,
        const std::vector<Vector>& w_fit,
        const std::function<bool()>& cancelled = {}
    );
    void run_bootstrap_stream_(
        int lambda_i,
        AdaptiveBootstrapState& bootstrap_state,
        std::vector<BootstrapBlocks>& thread_boot_worker,
        BoolMatrix& C_active,
        const std::vector<Vector>& w_fit,
        std::vector<Vector>& w_min,
        BootstrapResult& boot_results,
        const BlockRefList& blocks,
        BootstrapTimingSummary& timing_summary
    );

    // bootstrap state helpers
    bool same_design_(const BoolMatrix& lhs, const BoolMatrix& rhs) const;
    void reset_good_bootstrap_(
        AdaptiveBootstrapState& state,
        std::vector<Vector>& w_min,
        const std::vector<Vector>& w_fit,
        const BoolMatrix& C_active
    ) const;
    // bootstrap logging helpers
    void log_bootstrap_progress_(
        const AdaptiveBootstrapState& state,
        const int n_active_blocks,
        const int n_active_connections,
        const AdaptiveStopInfo& stop_info,
        const double elapsed_since_last_log
    ) const;
    void log_bootstrap_design_reset_(
        const AdaptiveBootstrapState& state,
        const int n_active_blocks,
        const int n_active_connections,
        const double elapsed_since_last_log
    ) const;

    // block helpers
    BlockRefList main_blocks_() const;
    BootstrapBlocks clone_blocks_() const;
    std::vector<std::string> block_names_(const BlockRefList& blocks) const;
    std::vector<int> block_dims_(const BlockRefList& blocks) const;

    // bootstrap storage and stopping helpers
    void ensure_bootstrap_lambda_storage_(
        BootstrapResult& boot_results,
        const int lambda_i,
        const std::vector<int>& block_dims,
        const int n_blocks
    ) const;
    AdaptiveStopInfo adaptive_stop_(AdaptiveBootstrapState& state, const BootstrapConfig& config) const;
    bool early_stop_lambda_(AdaptiveBootstrapState& state, int lambda_i, const BootstrapConfig& config) const;
    double criterion_score_with_weights_(const BlockRefList& blocks, const std::vector<Vector>& weights, const BoolMatrix& C);

    // bootstrap configuration helpers
    int bootstrap_fit_max_iter_() const;
    int bootstrap_n_threads_() const;

    // bootstrap design helpers
    int threshold_inactive_blocks_(std::vector<Vector>& w_min, BoolMatrix& C_active) const;
    int threshold_inactive_connections_(
        int lambda_i,
        AdaptiveBootstrapState& state,
        const BootstrapResult& boot_results,
        BoolMatrix& C_active,
        const int B_eff_override = -1
    );
    BoolMatrix reset_connections_keep_inactive_blocks_(const BoolMatrix& C_full, const BoolMatrix& C_current) const;
    void deactivate_isolated_blocks_(BoolMatrix& C_active) const;
    int count_active_connections_(const BoolMatrix& C_active) const;
    int count_active_blocks_(const BoolMatrix& C_active) const;
    std::vector<bool> active_blocks_from_C_(const BoolMatrix& C_active) const;

    // bootstrap result helpers
    std::pair<double, double> fisher_z_corr_ci_(
        const Matrix& corr_boot,
        const int row,
        const int B_eff,
        const double alpha_low,
        const double alpha_high
    ) const;
    void update_w_min_(Vector& w_min_j, const Vector& w_fit_j, const Vector& w_bj) const;
    void resize_bootstrap_lambda_results_(
        BootstrapResult& boot_results,
        int lambda_i,
        int n_blocks,
        const std::vector<int>& block_dims
    ) const;
    void compute_bootstrap_corr_cis_(BootstrapResult& boot_results, int n_blocks) const;

    // bootstrap row-index helpers
    void set_row_index_all_(const BlockRefList& blocks, const IndexVector& idx);
    void set_component_significance_row_index_all_(const BlockRefList& blocks, const unsigned seed);
    void clear_row_index_all_(const BlockRefList& blocks);

    // bootstrap resampling helpers
    IndexVector bootstrap_index_(int n, unsigned seed) const;
    IndexVector bootstrap_indices_(int n, std::mt19937_64& rng) const;
    IndexVector ordinary_bootstrap_indices_(const int n, std::mt19937_64& rng) const;
    IndexVector permutation_indices_(const int n, std::mt19937_64& rng) const;
    IndexVector stationary_bootstrap_indices_(
        const int n,
        const double mean_block_length,
        std::mt19937_64& rng
    ) const;
    IndexVector single_block_null_indices_(const int n, std::mt19937_64& rng) const;

    // bootstrap component significance helpers
    void annotate_component_significance_(Result& result, const ComponentSignificanceResult& significance) const;
    ComponentSignificanceResult inactive_component_significance_() const;
    MaxvarCorrelationResult rho_tot_maxvar_(const BlockRefList& blocks, const BoolMatrix& C) const;

    // bootstrap block-importance helpers
    void annotate_block_importance_(Result& result, const BlockImportanceResult& importance) const;
    bool block_importance_weights_(const std::vector<Vector>& eta, const BoolMatrix& C_active, Vector& v) const;
    Vector block_importance_target_(
        const std::vector<Vector>& eta,
        const BoolMatrix& C_active,
        const Vector& v,
        const int j
    ) const;
    double block_importance_rho_(const Vector& z, const Vector& s) const;

    // validation helpers
    void validate_fit_() const;
    void validate_bootstrap_support_() const;
    void validate_stationary_resampling_config_(const char* workflow) const;
    void validate_bootstrap_config_() const;
    void validate_component_significance_config_() const;
    void validate_block_importance_config_() const;
    void validate_index_(int j) const;
    void validate_lambda_grid_weights_(const std::vector<std::vector<double>>& lambda_grid) const;
    void validate_lambda_grid_weights_() const;
    void validate_bootstrap_weights_ci_(const BootstrapResult& boot_results, const int lambda_i, const int block_j, const SparseMatrix& Psi) const;

    // components utils
    void set_h_(const BlockRefList& blocks, const int h) {
        if (h < 0 || h >= n_comp_) throw std::out_of_range("component index");
        h_ = h;
        for (auto* b : blocks) b->set_h(h_);
    }
    void set_h_(const int h) {
        auto blocks = main_blocks_();
        set_h_(blocks, h);
    }
    void set_n_comp_(const BlockRefList& blocks, const int n_comp) {
        if (n_comp <= 0) throw std::invalid_argument("n_comp must be > 0");

        n_comp_ = n_comp;
        for (auto* b : blocks)
            b->set_n_comp(n_comp);

        if (h_ >= n_comp)
            set_h_(blocks, n_comp - 1);
    }

    // private setters
    void set_tau_auto_all_(const BlockRefList& blocks) const {
        for (auto* b : blocks)
            b->select_tau_auto();
    }
    void set_lambda_components_auto_all_(const BlockRefList& blocks) const {
        if constexpr (std::same_as<SamplingStrategy, IndependentSampling>) return;
        for (auto* b : blocks)
            b->set_lambda_components(-1);
    }
    void set_lambda_weights_all_(const BlockRefList& blocks, double lambda) const {
        for (auto* b : blocks)
            b->set_lambda_weights(lambda);
    }

    // eta
    Vector eta_(Block& b) const {
        // using the RGCCA own Psi_T (or components_m for independent)
        if constexpr (std::same_as<SamplingStrategy, TimeDependentSampling>) {
            return Psi_T_ * b.components().col(h_);
        } else {
            return b.components().col(h_);
        }
    }
    Vector eta_(Block& b, const Block& ref) const {
        // using reference block's Psi_T
        if constexpr (std::same_as<SamplingStrategy, TimeDependentSampling>) {
            return ref.Psi_T() * b.components().col(h_);
        } else {
            return b.components().col(h_);
        }
    }
    std::vector<Vector> eta_(const BlockRefList& blocks) const {
        std::vector<Vector> out;
        out.reserve(blocks.size());

        for (auto* b : blocks)
            out.push_back(eta_(*b));

        return out;
    }
    std::vector<Vector> eta_with_weights_(const BlockRefList& blocks, const std::vector<Vector>& weights) const {
        if (static_cast<int>(weights.size()) != n_blocks())
            throw std::logic_error("eta_with_weights_: size mismatch");

        std::vector<Vector> out;
        out.reserve(blocks.size());

        for (int j = 0; j < n_blocks(); ++j) {
            if (weights[j].size() != blocks[j]->n_dofs_weights())
                throw std::logic_error("eta_with_weights_: incompatible weight size");

            out.push_back(blocks[j]->data() * blocks[j]->Psi_D() * weights[j]);
        }

        return out;
    }
    std::vector<Vector> eta_with_weights_for_evaluation_(const BlockRefList& blocks, const std::vector<Vector>& weights) {
        if (static_cast<int>(weights.size()) != n_blocks())
            throw std::logic_error("eta_with_weights_for_evaluation_: size mismatch");

        std::vector<Vector> out;
        out.reserve(blocks.size());

        for (int j = 0; j < n_blocks(); ++j) {
            if (weights[j].size() != blocks[j]->n_dofs_weights())
                throw std::logic_error("eta_with_weights_for_evaluation_: incompatible weight size");

            // Evaluation normalizes with M, not Omega, so the penalty does not affect the reported component
            Vector eta = blocks[j]->normalized_component_for_evaluation(weights[j], h_);
            if constexpr (std::same_as<SamplingStrategy, TimeDependentSampling>) {
                eta = Psi_T_ * eta;
            }
            out.push_back(std::move(eta));
        }

        return out;
    }

    // covariance and correlation
    double cov_(const Vector& u, const Vector& v) const {
        const double den = opt_.bias ? u.size() : std::max<int>(1, u.size() - 1);
        return (u.dot(v) - static_cast<double>(u.size()) * u.mean() * v.mean()) / den;
        // return u.dot(v) / den;
    }
    double corr_(const Vector& u, const Vector& v) const {
        const double var_u = cov_(u, u);
        const double var_v = cov_(v, v);
        if (!(var_u > 0.0) || !(var_v > 0.0) || !std::isfinite(var_u) || !std::isfinite(var_v)) return 0.0;
        const double corr = cov_(u, v) / std::sqrt(var_u * var_v);
        return std::isfinite(corr) ? corr : 0.0;
    }
    double cov_value_(FitWorkspace& ws, int l, int k, const Vector& eta_l, const Vector& eta_k) const {
        if (!opt_.cache_covariances) return cov_(eta_l, eta_k);

        // compute or reuse cov(l,k); when computed, store and mark clean (both (l,k) and (k,l))
        if (!ws.dirty(l, k)) return ws.Cov(l, k);
        const double den = opt_.bias ? eta_l.size() : std::max<int>(1, eta_l.size() - 1);
        const double c = (
            eta_l.dot(eta_k) - static_cast<double>(eta_l.size()) * ws.means[l] * ws.means[k]
        ) / den;
        ws.Cov(l, k) = ws.Cov(k, l) = c;
        ws.dirty(l, k) = ws.dirty(k, l) = 0;
        return c;
    }
    void mark_cov_rowcol_dirty_(FitWorkspace& ws, int l) const {
        if (!opt_.cache_covariances) return;
        for (int k = 0; k < ws.Cov.rows(); ++k) {
            ws.dirty(l, k) = 1;
            ws.dirty(k, l) = 1;
        }
        ws.dirty(l, l) = 0;
        ws.Cov(l, l) = 1.0;
    }
    void covariance_matrix_(const std::vector<Vector>& eta, Matrix& Cov) const {
        Cov.setZero(n_blocks(), n_blocks());

        for (int j = 0; j < n_blocks(); ++j) {
            for (int k = j; k < n_blocks(); ++k) {
                const double c = cov_(eta[j], eta[k]);
                Cov(j, k) = c;
                Cov(k, j) = c;
            }
        }
    }
    void covariance_matrix_(const BlockRefList& blocks, Matrix& Cov) const {
        covariance_matrix_(eta_(blocks), Cov);
    }
    void correlation_matrix_(const BlockRefList& blocks, Matrix& Corr) const {
        correlation_matrix_(eta_(blocks), Corr);
    }
    void correlation_matrix_(const std::vector<Vector>& eta, Matrix& Corr) const {
        Corr.setIdentity(n_blocks(), n_blocks());

        std::vector<double> vars(n_blocks());
        for (int j = 0; j < n_blocks(); ++j) {
            vars[j] = cov_(eta[j], eta[j]);
        }

        for (int j = 0; j < n_blocks(); ++j) {
            for (int k = j + 1; k < n_blocks(); ++k) {
                double corr_jk = 0.0;

                if (vars[j] > 0.0 && vars[k] > 0.0) {
                    corr_jk = cov_(eta[j], eta[k]) / std::sqrt(vars[j] * vars[k]);
                }

                Corr(j, k) = corr_jk;
                Corr(k, j) = corr_jk;
            }
        }
    }
    void correlation_matrix_(const BlockRefList& blocks, const std::vector<Vector>& weights, Matrix& Corr) const {
        correlation_matrix_(eta_with_weights_(blocks, weights), Corr);
    }
    void correlation_matrix_raw_data_(
        const BlockRefList& blocks,
        const std::vector<Vector>& weights,
        Matrix& Corr
    ) const {
        if (static_cast<int>(weights.size()) != n_blocks())
            throw std::logic_error("correlation_matrix_raw_data_: size mismatch");

        Corr.setIdentity(n_blocks(), n_blocks());
        std::vector<Vector> eta(n_blocks());
        std::vector<double> vars(n_blocks(), 0.0);

        for (int j = 0; j < n_blocks(); ++j) {
            if (weights[j].size() != blocks[j]->n_dofs_weights())
                throw std::logic_error("correlation_matrix_raw_data_: incompatible weight size");
            if (weights[j].squaredNorm() == 0.0)
                continue;
            if (!blocks[j]->raw_score_cache_for_weight(weights[j], eta[j]))
                eta[j] = blocks[j]->raw_data() * blocks[j]->Psi_D() * weights[j];
            vars[j] = cov_(eta[j], eta[j]);
        }

        for (int j = 0; j < n_blocks(); ++j) {
            for (int k = j + 1; k < n_blocks(); ++k) {
                double corr_jk = 0.0;
                if (vars[j] > 0.0 && vars[k] > 0.0)
                    corr_jk = cov_(eta[j], eta[k]) / std::sqrt(vars[j] * vars[k]);
                Corr(j, k) = Corr(k, j) = corr_jk;
            }
        }
    }

    // optimization criteria
    double objective_(const BlockRefList& blocks, FitWorkspace& ws, const BoolMatrix& C) const {
        auto eta = eta_(blocks);
        for (int j = 0; j < n_blocks(); ++j)
            ws.means[j] = eta[j].mean();
        return objective_(ws, C, eta);
    }
    double objective_(FitWorkspace& ws, const BoolMatrix& C, const std::vector<Vector>& eta) const {
        double f = 0.0;
        for (int j = 0; j < n_blocks(); ++j) {
            for (int k = j; k < n_blocks(); ++k) {
                if (C(j, k)) {
                    const double cov_jk = cov_value_(ws, j, k, eta[j], eta[k]);
                    const double mult = j == k ? 1.0 : 2.0;
                    f += mult * opt_.scheme.g(cov_jk);
                }
            }
        }
        return f;
    }
private:
    Options opt_;
    DesignMode design_mode_ {DesignMode::Empty};

    int J_ {0}; // number of blocks
    int n_ {0}; // global number of observations

    std::vector<double> times_;
    SamplingDomain T_; // only used by TimeDependentSampling
    SparseMatrix Psi_T_; // only used by TimeDependentSampling

    int h_ {0}; // current component index
    int n_comp_{1};
    int n_comp_effective_{0};
    int n_comp_attempted_{0};

    std::vector<std::unique_ptr<Matrix>> data_blocks_;
    std::vector<BlockPtr> blocks_;

    BootstrapConfig bootstrap_config_;
    std::vector<BootstrapResult> bootstrap_selection_results_;
    std::vector<std::vector<double>> lambda_grid_weights_;

    bool initialized_ {false};
    BoolMatrix C_;
};

} // namespace fdapde

#define __FDAPDE_RGCCA_DEFINE_MODEL_BLOCKS__
#include "blocks.h"
#undef __FDAPDE_RGCCA_DEFINE_MODEL_BLOCKS__
#define __FDAPDE_RGCCA_DEFINE_MODEL_LOGGING__
#include "logging.h"
#undef __FDAPDE_RGCCA_DEFINE_MODEL_LOGGING__
#include "bootstrap_utils.h"
#define __FDAPDE_RGCCA_DEFINE_MODEL_VALIDATION__
#include "validation.h"
#undef __FDAPDE_RGCCA_DEFINE_MODEL_VALIDATION__

#endif // __FDAPDE_RGCCA_MODEL_H__

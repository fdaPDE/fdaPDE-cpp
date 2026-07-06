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


#ifndef __FDAPDE_RGCCA_MODEL_H__
#define __FDAPDE_RGCCA_MODEL_H__

#include "blocks.h"
#include "logging.h"
#include "statistics.h"
#include "fdaPDE/execution.h"

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

    // initialization
    void init() {
        ensure_design_initialized_();

        if (design_mode_ == DesignMode::Empty) {
            set_fully_connected_design_();
            design_mode_ = DesignMode::FullyConnected;
        }

        if constexpr (std::same_as<SamplingStrategy, TimeDependentSampling>) compute_Psi_T_();

        initialized_ = true;
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
        }
    }

    // fit
    std::vector<Result> fit(ComponentCallback on_component = {}) {
        if (!initialized_) init();

        // validation and logging
        validate_fit_();
        const int J = n_blocks();
        const bool run_model_selection = bootstrap_model_selection_requested_();
        log_fit_header_(run_model_selection);

        // room for results
        std::vector<Result> results;
        results.reserve(n_comp());
        bootstrap_selection_results_.clear();
        bootstrap_selection_results_.reserve(n_comp());

        // components loop
        bool completed_components = true;
        std::vector<bool> previous_active_blocks;
        std::vector<bool> carried_inactive_blocks(J, false);
        for (int hh = 0; hh < n_comp(); ++hh) {
            set_h_(hh);

            // bootstrap model selection
            BoolMatrix C_active = C_;
            if (inactive_block_signal_test_requested_() && !previous_active_blocks.empty()) {
                auto blocks = main_blocks_();
                apply_inactive_block_signal_test_(blocks, previous_active_blocks, carried_inactive_blocks, C_active);
            }
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

            // structural stop: no active design left, independent of significance testing
            if (count_active_connections_(C_active) == 0) {
                annotate_component_significance_(component_result, inactive_component_significance_());
                results.push_back(std::move(component_result));
                run_component_callback_(on_component, results.back());
                append_inactive_components_(results, hh + 1, J);
                completed_components = false;
                break;
            }

            // bootstrap component significance
            if (opt_.component_significance) {
                const auto significance = bootstrap_test_component_significance_(C_active);
                annotate_component_significance_(component_result, significance);

                if (!significance.significant) {
                    results.push_back(std::move(component_result));
                    run_component_callback_(on_component, results.back());
                    append_inactive_components_(results, hh + 1, J);
                    completed_components = false;
                    break;
                }
            }

            // store results
            results.push_back(std::move(component_result));

            // deflation
            step_start = log_step_start_("Deflate blocks");
            deflate_all_();
            log_step_end_(step_start);

            // component callback
            step_start = log_step_start_("Component callback");
            run_component_callback_(on_component, results.back());
            log_step_end_(step_start);

            previous_active_blocks = results.back().active_blocks;
        }

        // post-processing weights
        if (!on_component || !completed_components) {
            auto step_start = log_step_start_("Compute weights_star");
            compute_weights_star_();
            log_step_end_(step_start);
        }

        return results;
    }

    // getters
    void get_tau(const BlockRefList& blocks, std::vector<double>& tau_values) const {
        const int J = n_blocks();
        for (std::size_t j = 0; j < J; ++j)
            tau_values[j] = blocks[j]->tau();
    }
    void get_lambdas(const BlockRefList& blocks, std::vector<double> & lambda_components_values, std::vector<double> & lambda_weights_values) const {
        const int J = n_blocks();
        for (std::size_t j = 0; j < J; ++j) {
            lambda_components_values[j] = blocks[j]->lambda_components();
            lambda_weights_values[j] = blocks[j]->lambda_weights();
        }
    }

    // observers
    [[nodiscard]] int n() const { return n_; }
    [[nodiscard]] int n_comp() const { return n_comp_; }
    [[nodiscard]] int n_blocks() const { return J_; }
    [[nodiscard]] const Options& options() const { return opt_; }
    [[nodiscard]] const Scheme& scheme() const { return opt_.scheme; }
    [[nodiscard]] const std::vector<BlockPtr>& blocks() const { return blocks_; }
    [[nodiscard]] const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& C() const { return C_; }
    template <typename S = SamplingStrategy>
    requires std::same_as<S, TimeDependentSampling>
    [[nodiscard]] const SparseMatrix& Psi_T() const { return Psi_T_; };
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

        explicit FitWorkspace(int J) {
            Cov.setZero(J, J);
            dirty.setOnes(J, J);
            for (int j = 0; j < J; ++j) {
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
            const int J_
        ) : n_threads(n_threads_), J(J_) {
            seed = bootstrap_config.seed + static_cast<unsigned>(h);
            B_min = bootstrap_config.B_min;
            B_max = bootstrap_config.adaptive ? bootstrap_config.B_max : B_min;
            check_every = bootstrap_config.adaptive ? bootstrap_config.check_every : B_min;
            corr_pos_count.setZero(J, J);
            corr_neg_count.setZero(J, J);
        }

        void reset() {
            B_done = 0;
            B_total = 0;
            B_design = 0;
            B_stale = 0;
            B_cancelled = 0;
            design_epoch = 0;
            design_epoch_signal.store(0, std::memory_order_release);
            last_check_B_done = 0;
            stop = false;
            reset_good();
        }
        void reset_good() {
            B_done = 0;
            stable_checks = 0;
            crit_prev_check = std::numeric_limits<double>::infinity();
            crit = std::numeric_limits<double>::quiet_NaN();
            last_check_B_done = 0;
            corr_pos_count.setZero(J, J);
            corr_neg_count.setZero(J, J);
        }

        // config
        int n_threads;
        int J;
        int seed;
        int B_min;
        int B_max;
        int check_every;

        // state
        int B_done = 0;
        int B_total = 0;
        int B_design = 0;
        int B_stale = 0;
        int B_cancelled = 0;
        int design_epoch = 0;
        std::atomic<int> design_epoch_signal {0};
        int last_check_B_done = 0;
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
        double parallel_capacity = 0.0;
        double iters = 0.0;
        int capped_fits = 0;
        int max_iters = 0;
        int n_fits = 0;

        void add_sample(const double fit_time_sec, const int fit_iters, const bool capped, const bool fit_started) {
            if (!fit_started) return;
            fit_time += fit_time_sec;
            fit_time_sq += fit_time_sec * fit_time_sec;
            iters += static_cast<double>(fit_iters);
            capped_fits += capped ? 1 : 0;
            max_iters = std::max(max_iters, fit_iters);
            ++n_fits;
        }
        void set_parallel_capacity(const double wall_time_sec, const int n_threads) {
            parallel_capacity = wall_time_sec * static_cast<double>(n_threads);
        }
        [[nodiscard]] double efficiency() const {
            return parallel_capacity > 0.0 ? fit_time / parallel_capacity : 0.0;
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
    struct ComponentSignificanceResult {
        double rho_tot = std::numeric_limits<double>::quiet_NaN();
        double p_value = std::numeric_limits<double>::quiet_NaN();
        int B = 0;
        bool significant = true;
    };
    struct BootstrapSampleResult {
        std::vector<Vector> w;
        Matrix corr;
        double fit_time = 0.0;
        int fit_iters = 0;
        bool capped = false;
        bool cancelled = false;
        bool fit_started = false;
    };

    // components initialization
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
        const std::function<bool()>& cancelled = {}
    ) {
        const int J = n_blocks();
        const int max_iter = max_iter_override > 0 ? max_iter_override : opt_.max_iter;
        FitWorkspace ws(J);

        // room for results
        Result res(J);
        res.h = h_;
        res.obj_history.reserve(max_iter);

        // design update according to current active blocks
        res.C = C_active;
        res.active_blocks = active_blocks_from_C_(C_active);
        for (int j = 0; j < J; ++j) {
            if (!res.active_blocks[j]){
                blocks[j]->weights().col(h_).setZero();
                blocks[j]->components().col(h_).setZero();
            }
        }

        // initialization
        std::vector<Vector> eta_cache = eta_(blocks);
        res.obj_history.push_back(objective_(ws, res.C, eta_cache));
        auto w_prev = snapshot_weights_(blocks);
        if (update_component_lambdas && opt_.lambda_selection_components == LambdaSelection::Automatic)
            set_lambda_components_auto_all_(blocks);

        // main loop
        for (int s = 0; s < max_iter; ++s) {
            if (cancelled && cancelled()) {
                res.cancelled = true;
                return res;
            }
            for (int l = 0; l < J; ++l) {
                if (cancelled && cancelled()) {
                    res.cancelled = true;
                    return res;
                }

                // skip deactivated blocks
                if (!res.active_blocks[l]) continue;

                // inner-component assembler
                Vector nu_l = Vector::Zero(blocks[l]->n());
                const Vector& eta_l = eta_cache[l];
                for (int k = 0; k < J; ++k) {
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

        if (cancelled && cancelled()) {
            res.cancelled = true;
            return res;
        }

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

    void apply_inactive_block_signal_test_(
        const BlockRefList& blocks,
        const std::vector<bool>& previous_active_blocks,
        std::vector<bool>& carried_inactive_blocks,
        BoolMatrix& C_active
    ) {
        const int J = static_cast<int>(blocks.size());
        if (static_cast<int>(previous_active_blocks.size()) != J)
            throw std::logic_error("inactive block signal test: active block size mismatch");
        if (static_cast<int>(carried_inactive_blocks.size()) != J)
            throw std::logic_error("inactive block signal test: carried block size mismatch");

        for (int j = 0; j < J; ++j) {
            if (previous_active_blocks[j])
                carried_inactive_blocks[j] = false;
            if (carried_inactive_blocks[j]) {
                C_active.row(j).setConstant(false);
                C_active.col(j).setConstant(false);
            }
        }

        std::vector<Vector> active_scores(J);
        for (int k = 0; k < J; ++k) {
            if (previous_active_blocks[k])
                active_scores[k] = blocks[k]->svd_init().nu;
        }

        for (int j = 0; j < J; ++j) {
            if (previous_active_blocks[j]) continue;
            if (carried_inactive_blocks[j]) continue;
            const bool has_active_neighbor = inactive_block_has_active_neighbor_(previous_active_blocks, C_active, j);
            if (!inactive_block_has_residual_signal_(blocks, active_scores, previous_active_blocks, C_active, j)) {
                C_active.row(j).setConstant(false);
                C_active.col(j).setConstant(false);
                if (has_active_neighbor)
                    carried_inactive_blocks[j] = true;
            }
        }
    }

    bool inactive_block_has_active_neighbor_(
        const std::vector<bool>& active_blocks,
        const BoolMatrix& C_active,
        const int candidate
    ) const {
        const int J = static_cast<int>(active_blocks.size());
        for (int k = 0; k < J; ++k)
            if (active_blocks[k] && C_active(candidate, k))
                return true;
        return false;
    }

    bool inactive_block_has_residual_signal_(
        const BlockRefList& blocks,
        const std::vector<Vector>& active_scores,
        const std::vector<bool>& active_blocks,
        const BoolMatrix& C_active,
        const int candidate
    ) const {
        std::vector<int> active_neighbors;
        const int J = static_cast<int>(active_blocks.size());
        for (int k = 0; k < J; ++k) {
            if (active_blocks[k] && C_active(candidate, k))
                active_neighbors.push_back(k);
        }
        if (active_neighbors.empty()) return false;

        const Vector candidate_score = blocks[candidate]->svd_init().nu;
        const double observed = residual_signal_stat_(candidate_score, active_scores, active_neighbors);
        if (!(observed > 0.0) || !std::isfinite(observed)) return false;

        const int B = bootstrap_config_.inactive_block_signal_resamples;
        int ge_count = 0;
        std::mt19937_64 rng(
            bootstrap_config_.seed +
            static_cast<unsigned>(1000003 * (h_ + 1)) +
            static_cast<unsigned>(9176 * (candidate + 1))
        );
        for (int b = 0; b < B; ++b) {
            const IndexVector idx = inactive_block_signal_null_indices_(candidate_score.size(), rng);
            Vector null_score(idx.size());
            for (int i = 0; i < idx.size(); ++i)
                null_score[i] = candidate_score[idx[i]];

            const double null_stat = residual_signal_stat_(null_score, active_scores, active_neighbors);
            if (std::isfinite(null_stat) && null_stat >= observed)
                ++ge_count;
        }

        const double p_value = static_cast<double>(ge_count + 1) / static_cast<double>(B + 1);
        return p_value <= bootstrap_config_.inactive_block_signal_alpha;
    }

    double residual_signal_stat_(
        const Vector& candidate_score,
        const std::vector<Vector>& active_scores,
        const std::vector<int>& active_neighbors
    ) const {
        const double candidate_var = cov_(candidate_score, candidate_score);
        if (!(candidate_var > 0.0) || !std::isfinite(candidate_var)) return 0.0;

        double stat = 0.0;
        for (const int k : active_neighbors) {
            const double active_var = cov_(active_scores[k], active_scores[k]);
            if (!(active_var > 0.0) || !std::isfinite(active_var)) continue;
            const double corr = cov_(candidate_score, active_scores[k]) / std::sqrt(candidate_var * active_var);
            if (std::isfinite(corr)) stat = std::max(stat, std::abs(corr));
        }
        return stat;
    }

    IndexVector inactive_block_signal_null_indices_(const int n, std::mt19937_64& rng) const {
        if (bootstrap_config_.resampling_strategy == ResamplingStrategy::Stationary)
            return stationary_bootstrap_indices_(n, bootstrap_config_.stationary_block_length, rng);
        return permutation_indices_(n, rng);
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
        const int J = static_cast<int>(blocks.size());
        BoolMatrix C_active = C_initial;
        BoolMatrix C_best = C_initial;

        // init bootstrap
        auto step_start = log_step_start_("Init bootstrap");
        AdaptiveBootstrapState bootstrap_state(bootstrap_config_, n_threads, h_, J);
        const auto block_dims = block_dims_(blocks);
        BootstrapResult boot_results(
            h_, bootstrap_state.B_max, lambda_grid,
            block_names_(blocks), block_dims, bootstrap_config_.ci_level
        );
        log_step_end_(step_start);

        // preliminary fit
        step_start = log_step_start_("Preliminary fit");
        if (select_lambda)
            set_lambda_weights_all(lambda_grid.back());
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
        auto fixed_weight_lambda = [&]() {
            for (const auto& block : blocks) {
                const double lambda = block->lambda_weights();
                if (std::isfinite(lambda)) return lambda;
            }
            return std::numeric_limits<double>::quiet_NaN();
        };

        // model selection loop
        int n_lambda = static_cast<int>(lambda_grid.size());
        for (int lambda_i = n_lambda - 1; lambda_i >= 0; --lambda_i) {
            ensure_bootstrap_lambda_storage_(boot_results, lambda_i, block_dims, J);
            const bool reuse_preliminary_fit = lambda_i == n_lambda - 1;

            if (select_lambda) {
                // current lambda
                const double lambda = lambda_grid[lambda_i];
                log_bootstrap_lambda_candidate_(select_lambda, lambda);
                if (!reuse_preliminary_fit) {
                    set_lambda_weights_all(lambda);
                    for (auto& worker : thread_boot_worker)
                        set_lambda_weights_all_(worker.refs, lambda);
                }
            } else {
                log_bootstrap_lambda_candidate_(select_lambda, fixed_weight_lambda());
            }

            // init warm start at lambda
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
            boot_results.corr_min_by_lambda[lambda_i].array() *= (C_lambda.cast<double>() + Matrix::Identity(J, J)).array();
            log_step_end_(step_start);
            if (!std::isfinite(bootstrap_state.crit))
                bootstrap_state.crit = criterion_score_with_weights_(blocks, w_min, C_initial);
            boot_results.criterion[lambda_i] = bootstrap_state.crit;
            boot_results.B_used_by_lambda[lambda_i] = bootstrap_state.B_done;

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
                select_lambda ? lambda_grid[lambda_i] : fixed_weight_lambda()
            );

            // early stop
            if (early_stop_lambda_(bootstrap_state, lambda_i, bootstrap_config_)) {
                break;
            }

            // reset connections, but keep fully deactivated blocks off
            C_active = reset_connections_keep_inactive_blocks_(C_initial, C_active);

        }

        // resize allocated space for bootstrap results
        step_start = log_step_start_("Resize bootstrap results");
        resize_bootstrap_results_(boot_results, static_cast<int>(lambda_grid.size()), J, block_dims);
        log_step_end_(step_start);

        // compute correlation matrices CI
        step_start = log_step_start_("Compute bootstrap correlation CIs");
        compute_bootstrap_corr_cis_(boot_results, J);
        log_step_end_(step_start);

        // save optimal design
        if (bootstrap_state.best_i < 0)
            throw std::runtime_error("No model candidate was evaluated during bootstrap selection");
        boot_results.lambda_opt_index = bootstrap_state.best_i;
        if (select_lambda)
            boot_results.lambda_opt = lambda_grid[bootstrap_state.best_i];
        boot_results.active_blocks = active_blocks_from_C_(C_best);
        log_bootstrap_selection_result_(boot_results, C_best, select_lambda);

        bootstrap_selection_results_.push_back(std::move(boot_results));
        ModelSelectionResult out;
        out.lambda_selected = select_lambda;
        out.lambda = bootstrap_selection_results_.back().lambda_opt;
        out.C_active = C_best;
        return out;
    }

    // bootstrap component significance
    ComponentSignificanceResult bootstrap_test_component_significance_(const BoolMatrix& C_active) {
        log_bootstrap_component_significance_header_();
        ComponentSignificanceResult out;

        // observed statistic on the fitted component
        auto blocks = main_blocks_();
        const int J = static_cast<int>(blocks.size());
        out.rho_tot = rho_tot_(blocks, C_active);

        // inactive or degenerate components are declared non-significant
        if (count_active_connections_(C_active) == 0 || !std::isfinite(out.rho_tot)) {
            out.B = 0;
            out.p_value = 1.0;
            out.significant = false;
            return out;
        }

        const int B = bootstrap_config_.component_significance_resamples;
        const int n_threads = bootstrap_n_threads_();
        out.B = B;

        // each worker owns a mutable clone with its own row permutation
        auto step_start = log_step_start_("Clone significance worker blocks");
        std::vector<BootstrapBlocks> thread_boot_worker(n_threads);
        for (int t = 0; t < n_threads; ++t)
            thread_boot_worker[t] = clone_blocks_();
        log_step_end_(step_start);

        std::vector<int> thread_ge_count(n_threads, 0);
        const unsigned seed = bootstrap_config_.seed + static_cast<unsigned>(1000003 * (h_ + 1));
        const auto active_blocks = active_blocks_from_C_(C_active);

        // permutation null: break row alignment within each block, then refit
        step_start = log_step_start_("Significance bootstrap");
        parallel_for(0, B, 1, [&](int b) {
            const int tid = this_thread_id();
            auto& boot_blocks = thread_boot_worker[tid];

            set_permuted_row_index_all_(
                boot_blocks.refs,
                seed + static_cast<unsigned>(7919 * (b + 1))
            );

            init_comp_(boot_blocks.refs, InitStrategy::WarmStart, false, &active_blocks);
            fit_component_(boot_blocks.refs, C_active, false, bootstrap_fit_max_iter_());

            // count null statistics at least as extreme as the observed one
            const double rho_star = rho_tot_(boot_blocks.refs, C_active);
            if (std::isfinite(rho_star) && rho_star >= out.rho_tot)
                ++thread_ge_count[tid];

            clear_row_index_all_(boot_blocks.refs);
        });
        log_step_end_(step_start);

        // corrected empirical upper-tail p-value
        const int ge_count = std::accumulate(thread_ge_count.begin(), thread_ge_count.end(), 0);
        out.p_value = static_cast<double>(ge_count + 1) / static_cast<double>(B + 1);
        out.significant = out.p_value <= bootstrap_config_.component_significance_alpha;

        log_component_significance_(out);

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
    BoolMatrix inactive_design_(const int J) const {
        BoolMatrix C_inactive(J, J);
        C_inactive.setConstant(false);
        return C_inactive;
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
    void deflate_all_() const {
        for (auto& b : blocks_) b->deflate(opt_.deflation_mode);
    }
    void compute_weights_star_() {
        for (auto& b : blocks_) b->compute_weights_star();
    }
    void compute_weights_star_(const int h) {
        for (auto& b : blocks_) b->compute_weights_star(h);
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
    bool inactive_block_signal_test_requested_() const {
        return bootstrap_config_.inactive_block_signal_test;
    }
    bool weight_lambda_selection_requested_() const {
        return opt_.lambda_selection_weights == LambdaSelection::Automatic;
    }
    std::vector<double> model_selection_lambda_grid_() const {
        if (weight_lambda_selection_requested_())
            return lambda_grid_weights_[h_];

        return {std::numeric_limits<double>::quiet_NaN()};
    }

    // component callback runner
    void run_component_callback_(const ComponentCallback& on_component, const Result& result) {
        if (!on_component) return;
        compute_weights_star_(result.h);
        on_component(*this, result);
    }

    // logging helpers
    std::chrono::high_resolution_clock::time_point log_step_start_(std::string_view label) const;
    void log_step_end_(const std::chrono::high_resolution_clock::time_point start) const;
    void log_weight_lambda_(double lambda) const;
    void log_fit_header_(bool run_model_selection) const;
    void log_bootstrap_model_selection_header_() const;
    void log_bootstrap_component_significance_header_() const;
    void log_bootstrap_lambda_candidate_(bool lambda_selected, double lambda) const;
    void log_component_significance_(const ComponentSignificanceResult& significance) const;
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
    bool connection_deactivation_ready_(const AdaptiveBootstrapState& state) const;

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
        const int J
    ) const;
    AdaptiveStopInfo adaptive_stop_(AdaptiveBootstrapState& state, const BootstrapConfig& config) const;
    bool early_stop_lambda_(AdaptiveBootstrapState& state, int lambda_i, const BootstrapConfig& config) const;

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
    void resize_bootstrap_results_(
        BootstrapResult& boot_results,
        int n_lambdas,
        int J,
        const std::vector<int>& block_dims
    ) const;
    void compute_bootstrap_corr_cis_(BootstrapResult& boot_results, int J) const;

    // bootstrap row-index helpers
    void set_row_index_all_(const BlockRefList& blocks, const IndexVector& idx);
    void set_permuted_row_index_all_(const BlockRefList& blocks, const unsigned seed);
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

    // bootstrap component significance helpers
    void annotate_component_significance_(Result& result, const ComponentSignificanceResult& significance) const {
        result.rho_tot = significance.rho_tot;
        result.rho_tot_p_value = significance.p_value;
        result.rho_tot_bootstrap_count = significance.B;
        result.component_significant = significance.significant;
    }
    ComponentSignificanceResult inactive_component_significance_() const {
        ComponentSignificanceResult out;
        out.rho_tot = 0.0;
        out.p_value = 1.0;
        out.B = 0;
        out.significant = false;
        return out;
    }
    void append_inactive_components_(std::vector<Result>& results, const int from_h, const int J) {
        const BoolMatrix C_inactive = inactive_design_(J);
        const ComponentSignificanceResult significance = inactive_component_significance_();

        for (int hh = from_h; hh < n_comp(); ++hh) {
            set_h_(hh);
            const auto step_start = log_step_start_("Inactive component fit");
            Result inactive_result = fit_component_(C_inactive);
            log_step_end_(step_start);
            annotate_component_significance_(inactive_result, significance);
            results.push_back(std::move(inactive_result));
        }
    }

    // validation helpers
    void validate_fit_() const;
    void validate_bootstrap_support_() const;
    void validate_bootstrap_config_() const;
    void validate_component_significance_config_() const;
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
        const int J = static_cast<int>(blocks.size());

        if (static_cast<int>(weights.size()) != J)
            throw std::logic_error("eta_with_weights_: size mismatch");

        std::vector<Vector> out;
        out.reserve(blocks.size());

        for (int j = 0; j < J; ++j) {
            if (weights[j].size() != blocks[j]->n_dofs_weights())
                throw std::logic_error("eta_with_weights_: incompatible weight size");

            out.push_back(blocks[j]->data() * blocks[j]->Psi_D() * weights[j]);
        }

        return out;
    }
    std::vector<Vector> eta_with_weights_for_evaluation_(const BlockRefList& blocks, const std::vector<Vector>& weights) {
        const int J = static_cast<int>(blocks.size());

        if (static_cast<int>(weights.size()) != J)
            throw std::logic_error("eta_with_weights_for_evaluation_: size mismatch");

        std::vector<Vector> out;
        out.reserve(blocks.size());

        for (int j = 0; j < J; ++j) {
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
    double cov_value_(FitWorkspace& ws, int l, int k, const Vector& eta_l, const Vector& eta_k) const {
        if (!opt_.cache_covariances) return cov_(eta_l, eta_k);

        // compute or reuse cov(l,k); when computed, store and mark clean (both (l,k) and (k,l))
        if (!ws.dirty(l, k)) return ws.Cov(l, k);
        const double c = cov_(eta_l, eta_k);
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
        const int J = static_cast<int>(eta.size());
        Cov.setZero(J, J);

        for (int j = 0; j < J; ++j) {
            for (int k = j; k < J; ++k) {
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
        const int J = static_cast<int>(eta.size());
        Corr.setIdentity(J, J);

        std::vector<double> vars(J);
        for (int j = 0; j < J; ++j) {
            vars[j] = cov_(eta[j], eta[j]);
        }

        for (int j = 0; j < J; ++j) {
            for (int k = j + 1; k < J; ++k) {
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
        const int J = static_cast<int>(blocks.size());
        if (static_cast<int>(weights.size()) != J)
            throw std::logic_error("correlation_matrix_raw_data_: size mismatch");

        Corr.setIdentity(J, J);
        std::vector<Vector> eta(J);
        std::vector<double> vars(J, 0.0);

        for (int j = 0; j < J; ++j) {
            if (weights[j].size() != blocks[j]->n_dofs_weights())
                throw std::logic_error("correlation_matrix_raw_data_: incompatible weight size");
            if (weights[j].squaredNorm() == 0.0)
                continue;
            if (!blocks[j]->raw_score_cache_for_weight(weights[j], eta[j]))
                eta[j] = blocks[j]->raw_data() * blocks[j]->Psi_D() * weights[j];
            vars[j] = cov_(eta[j], eta[j]);
        }

        for (int j = 0; j < J; ++j) {
            for (int k = j + 1; k < J; ++k) {
                double corr_jk = 0.0;
                if (vars[j] > 0.0 && vars[k] > 0.0)
                    corr_jk = cov_(eta[j], eta[k]) / std::sqrt(vars[j] * vars[k]);
                Corr(j, k) = Corr(k, j) = corr_jk;
            }
        }
    }

    // optimization criteria
    double objective_(const BlockRefList& blocks, FitWorkspace& ws, const BoolMatrix& C) const {
        return objective_(ws, C, eta_(blocks));
    }
    double objective_(FitWorkspace& ws, const BoolMatrix& C, const std::vector<Vector>& eta) const {
        const int J = n_blocks();
        double f = 0.0;
        for (int j = 0; j < J; ++j) {
            for (int k = j; k < J; ++k) {
                if (C(j, k)) {
                    const double cov_jk = cov_value_(ws, j, k, eta[j], eta[k]);
                    const double mult = j == k ? 1.0 : 2.0;
                    f += mult * opt_.scheme.g(cov_jk);
                }
            }
        }
        return f;
    }
    double rho_tot_from_correlation_(const Matrix& Corr, const BoolMatrix& C) const {
        const int J = static_cast<int>(Corr.rows());

        double num = 0.0;
        double den = 0.0;

        for (int j = 0; j < J; ++j) {
            for (int k = j + 1; k < J; ++k) {
                if (!C(j, k)) continue;

                const double corr_jk = Corr(j, k);
                if (!std::isfinite(corr_jk)) continue;

                num += opt_.scheme.sign_invariant ? std::abs(corr_jk) : corr_jk;
                den += 1.0;
            }
        }

        return den > 0.0 ? num / den : 0.0;
    }
    double rho_tot_(const BlockRefList& blocks, const BoolMatrix& C) const {
        Matrix Corr;
        correlation_matrix_(blocks, Corr);
        return rho_tot_from_correlation_(Corr, C);
    }
    double rho_tot_with_weights_(const BlockRefList& blocks, const std::vector<Vector>& weights, const BoolMatrix& C) const {
        Matrix Corr;
        correlation_matrix_(eta_with_weights_(blocks, weights), Corr);
        return rho_tot_from_correlation_(Corr, C);
    }
    double criterion_score_with_weights_(const BlockRefList& blocks, const std::vector<Vector>& weights, const BoolMatrix& C) {
        const int J = n_blocks();

        double num = 0.0;
        double den = 0.0;

        const std::vector<Vector> eta = eta_with_weights_for_evaluation_(blocks, weights);
        for (int j = 0; j < J; ++j) {
            for (int k = j + 1; k < J; ++k) {
                if (!C(j, k)) continue;

                const double cov_jk = cov_(eta[j], eta[k]);
                num += opt_.scheme.g(cov_jk);
                den += 1.0;
            }
        }

        return den > 0.0 ? num / den : 0.0;
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

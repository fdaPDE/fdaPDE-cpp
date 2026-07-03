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

namespace fdapde {

// RGCCA owns block wiring, component fitting, bootstrap selection, and result storage.
template <typename SamplingStrategy>
class RGCCA {
public:
    using Block = rgcca::internals::BaseBlock<SamplingStrategy>;
    using BlockPtr = std::unique_ptr<Block>;
    using BlockRefList = std::vector<Block*>;
    using BlockOwnerList = std::vector<std::unique_ptr<Block>>;
    using Matrix = typename Block::Matrix;
    using BoolMatrix = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMatrix = typename Block::SparseMatrix;
    using Vector = typename Block::Vector;
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
    template <typename S = SamplingStrategy>
    requires std::same_as<S, IndependentSampling>
    explicit RGCCA(const int n, const Options& opt = Options(), const int n_comp = 1) : n_(n), opt_(opt), n_comp_(n_comp) {}

    template <typename S = SamplingStrategy>
    requires std::same_as<S, TimeDependentSampling>
    explicit RGCCA(const int n, const Triangulation<1, 1>& T, const Options& opt = Options(), const int n_comp = 1) : n_(n), T_(T), opt_(opt), n_comp_(n_comp) {}

    // blocks management
    int add_block(BlockPtr b) {
        if (!b) throw std::invalid_argument("RGCCA/add_block: null block");
        if constexpr (std::same_as<SamplingStrategy, IndependentSampling>) {
            if (b->n() != n()) throw std::invalid_argument("RGCCA/add_block: n mismatch");
        } else {
            add_times_(b->times());
        }
        b->set_bias(opt_.bias);
        b->set_raw_data_mutable(true);
        b->set_mode(opt_.mode);
        b->set_weight_sign_constraint(opt_.weight_sign_constraint);
        b->set_n_comp(n_comp());
        blocks_.emplace_back(std::move(b));
        initialized_ = false;
        return ++J_;
    }
    template<typename S = SamplingStrategy>
    requires std::same_as<S, IndependentSampling>
    int add_multivariate_block(std::string block_name, Matrix&& X) {
        data_blocks_.push_back(std::make_unique<Matrix>(std::move(X)));
        return add_block(rgcca::internals::make_multivariate_block<SamplingStrategy>(block_name, data_blocks_.back().get()));
    }
    template<typename S = SamplingStrategy>
    requires std::same_as<S, TimeDependentSampling>
    int add_multivariate_block(std::string block_name, const Vector& times, Matrix&& X) {
        data_blocks_.push_back(std::make_unique<Matrix>(std::move(X)));
        return add_block(rgcca::internals::make_multivariate_block<SamplingStrategy>(block_name, T_, times, data_blocks_.back().get()));
    }
    template <typename GeoFrame, typename WeightsPenaltyType>
    requires std::same_as<SamplingStrategy, IndependentSampling>
    int add_functional_block(std::string block_name, const GeoFrame& gf, Matrix&& X, WeightsPenaltyType&& weights_penalty) {
        data_blocks_.push_back(std::make_unique<Matrix>(std::move(X)));
        return add_block(rgcca::internals::make_functional_block<SamplingStrategy>(block_name, gf, data_blocks_.back().get(), std::forward<WeightsPenaltyType>(weights_penalty)));
    }
    template <typename GeoFrame, typename WeightsPenaltyType>
    requires std::same_as<SamplingStrategy, TimeDependentSampling>
    int add_functional_block(std::string block_name, const Vector& times, const GeoFrame& gf, Matrix&& X, WeightsPenaltyType&& weights_penalty) {
        data_blocks_.push_back(std::make_unique<Matrix>(std::move(X)));
        return add_block(rgcca::internals::make_functional_block<SamplingStrategy>(block_name, T_, times, gf, data_blocks_.back().get(), std::forward<WeightsPenaltyType>(weights_penalty)));
    }
    void connect(int j, int k, bool on = true) {
        ensure_design_initialized_();

        check_index_(j);
        check_index_(k);

        if (j == k) {
            fdapde::cout << "RGCCA::connect(): ignoring self-connection for block " << j << '\n';
            return;
        }

        C_(j, k) = on;
        C_(k, j) = on;

        design_mode_ = DesignMode::Custom;
    }

    // initialization
    void init() {
        if (n_blocks() < 2) throw std::runtime_error("RGCCA: need ≥ 2 blocks");

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
        if (lambda_grid.empty())
            throw std::invalid_argument("lambda grid cannot be empty");
        lambda_grid_weights_.assign(n_comp(), lambda_grid);
    }
    void set_lambda_grid_weights(const std::vector<std::vector<double>>& lambda_grid) {
        if (lambda_grid.empty())
            throw std::invalid_argument("lambda grid cannot be empty");

        if (static_cast<int>(lambda_grid.size()) == 1) {
            set_lambda_grid_weights(lambda_grid.front());
            return;
        }

        if (static_cast<int>(lambda_grid.size()) != n_comp())
            throw std::invalid_argument("lambda grid must have size 1 or n_comp");

        for (const auto& grid : lambda_grid) {
            if (grid.empty())
                throw std::invalid_argument("lambda grid contains an empty component grid");
        }

        lambda_grid_weights_ = lambda_grid;
    }

    // setters
    void set_n_comp(const int n_comp) {
        auto blocks = main_blocks_();
        set_n_comp_(blocks, n_comp);
    }
    void set_bootstrap_config(const BootstrapConfig bootstrap_config) {
        bootstrap_config_ = bootstrap_config;
    }

    // fit
    std::vector<Result> fit(ComponentCallback on_component = {}) {
        if (!initialized_) init();

        const int J = n_blocks();
        if (J < 2) throw std::runtime_error("RGCCA: need ≥ 2 blocks");
        const bool run_model_selection = bootstrap_model_selection_requested_();
        fdapde::cout << opt_ << '\n';
        if (run_model_selection || opt_.component_significance) {
            fdapde::cout << bootstrap_config_ << '\n';
        }
        if (opt_.component_significance) {
            validate_bootstrap_support_();
            validate_component_significance_config_();
        }
        if (run_model_selection) {
            validate_bootstrap_support_();
            validate_bootstrap_config_();
        }
        if (opt_.lambda_selection_weights == LambdaSelection::Automatic) {
            validate_lambda_grid_weights_();
        }

        // room for results
        std::vector<Result> results;
        results.reserve(n_comp());
        bootstrap_selection_results_.clear();
        bootstrap_selection_results_.reserve(n_comp());

        // components loop
        for (int hh = 0; hh < n_comp(); ++hh) {
            set_h_(hh);

            // bootstrap model selection
            BoolMatrix C_active = C_;
            if (run_model_selection) {
                auto selection = bootstrap_model_selection_();
                C_active = std::move(selection.C_active);
                if (selection.lambda_selected)
                    set_lambda_weights_all(selection.lambda);
            }

            // final fit
            auto step_start = print_step_start_("Final component fit");
            auto blocks = main_blocks_();
            const auto active_blocks = active_blocks_from_C_(C_active);
            init_comp_(blocks, InitStrategy::None, true, &active_blocks);
            Result component_result = fit_component_(blocks, C_active);
            print_step_end_(step_start);

            if (opt_.component_significance) {
                const auto significance = bootstrap_test_component_significance_(C_active);
                annotate_component_significance_(component_result, significance);

                if (!significance.significant) {
                    results.push_back(std::move(component_result));
                    run_component_callback_(on_component, results.back());
                    append_inactive_components_(results, hh + 1, J);
                    break;
                }
            }

            results.push_back(std::move(component_result));

            step_start = print_step_start_("Deflate blocks");
            deflate_all_();
            print_step_end_(step_start);

            run_component_callback_(on_component, results.back());
        }

        // post-processing weights
        auto step_start = print_step_start_("Compute weights_star");
        compute_weights_star_();
        print_step_end_(step_start);

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

    [[nodiscard]] std::pair<Vector, Vector> bootstrap_weights_ci(
        const int h,
        const int lambda_i,
        const int block_j,
        const SparseMatrix& Psi
    ) const {
        const auto& boot_results = bootstrap_selection_result_(h);
        return bootstrap_weights_ci_(boot_results, lambda_i, block_j, Psi);
    }
    [[nodiscard]] std::pair<Vector, Vector> bootstrap_weights_ci(
        const int h,
        const int block_j,
        const SparseMatrix& Psi
    ) const {
        const auto& boot_results = bootstrap_selection_result_(h);
        return bootstrap_weights_ci_(boot_results, bootstrap_lambda_opt_index_(boot_results), block_j, Psi);
    }
    template <typename DataLocs>
    requires(!std::same_as<std::decay_t<DataLocs>, SparseMatrix>)
    [[nodiscard]] std::pair<Vector, Vector> bootstrap_weights_ci(
        const int h,
        const int lambda_i,
        const int block_j,
        const DataLocs& locs
    ) const {
        check_index_(block_j);
        const SparseMatrix Psi = blocks_[block_j]->Psi_at(locs);
        return bootstrap_weights_ci(h, lambda_i, block_j, Psi);
    }
    template <typename DataLocs>
    requires(!std::same_as<std::decay_t<DataLocs>, SparseMatrix>)
    [[nodiscard]] std::pair<Vector, Vector> bootstrap_weights_ci(
        const int h,
        const int block_j,
        const DataLocs& locs
    ) const {
        const auto& boot_results = bootstrap_selection_result_(h);
        check_index_(block_j);
        const SparseMatrix Psi = blocks_[block_j]->Psi_at(locs);
        return bootstrap_weights_ci_(boot_results, bootstrap_lambda_opt_index_(boot_results), block_j, Psi);
    }

private:
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

    // initialization utils
    void check_index_(int j) const {
        if (j < 0 || j >= static_cast<int>(blocks_.size())) throw std::out_of_range("block index");
    }
    const BootstrapResult& bootstrap_selection_result_(const int h) const {
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
    int bootstrap_lambda_opt_index_(const BootstrapResult& boot_results) const {
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
    std::pair<Vector, Vector> bootstrap_weights_ci_(
        const BootstrapResult& boot_results,
        const int lambda_i,
        const int block_j,
        const SparseMatrix& Psi
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

        const Matrix& w_boot = boot_results.w_boot_by_lambda[lambda_i][block_j];
        if (Psi.rows() <= 0 || Psi.cols() != w_boot.rows()) {
            throw std::invalid_argument(
                "RGCCA: Psi must have one column per bootstrap weight coefficient"
            );
        }

        const double alpha_low = (1.0 - boot_results.ci_level) / 2.0;
        const double alpha_high = 1.0 - alpha_low;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const int B_eff = std::min(
            boot_results.B_used_by_lambda[lambda_i],
            static_cast<int>(w_boot.cols())
        );

        Vector ci_low(Psi.rows());
        Vector ci_high(Psi.rows());

        if (B_eff <= 0) {
            ci_low.setConstant(nan);
            ci_high.setConstant(nan);
            return {ci_low, ci_high};
        }

        const Matrix w_eval = Psi * w_boot.leftCols(B_eff);

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

    // blocks utils
    BlockRefList main_blocks_() const {
        BlockRefList out;
        out.reserve(blocks_.size());
        for (const auto& b : blocks_)
            out.push_back(b.get());
        return out;
    }
    struct BootstrapBlocks {
        BlockOwnerList owners;
        BlockRefList refs;
    };
    BootstrapBlocks clone_blocks_() const {
        BootstrapBlocks out;
        out.owners.reserve(blocks_.size());
        out.refs.reserve(blocks_.size());

        for (const auto& b : blocks_) {
            auto copy = b->clone();
            copy->set_raw_data_mutable(false);

            out.refs.push_back(copy.get());
            out.owners.push_back(std::move(copy));
        }

        return out;
    }
    void copy_weights_snapshot_(const BlockRefList& blocks, const std::vector<Vector>& weights) const {
        if (blocks.size() != weights.size())
            throw std::logic_error("copy_weights_snapshot_: size mismatch");

        for (std::size_t j = 0; j < blocks.size(); ++j) {
            if (blocks[j]->n_dofs_weights() != weights[j].size())
                throw std::logic_error("copy_weights_snapshot_: incompatible weight size");

            blocks[j]->weights().col(h_) = weights[j];
        }
    }
    std::vector<std::string> block_names_(const BlockRefList& blocks) const {
        std::vector<std::string> out;
        out.reserve(blocks.size());

        for (auto* b : blocks)
            out.push_back(b->name());

        return out;
    }
    std::vector<int> block_dims_(const BlockRefList& blocks) const {
        std::vector<int> out;
        out.reserve(blocks.size());

        for (auto* b : blocks)
            out.push_back(b->n_dofs_weights());

        return out;
    }

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

            const double rel_obj_drop = (obj_prev - f_obj) / (1.0 + std::abs(obj_prev));
            if (rel_obj_drop > opt_.tol)
                res.monotone = false;

            // stopping criteria
            const double delta_obj = std::abs(f_obj - obj_prev);
            const double delta_w = weights_variation_(blocks, w_prev);
            if (delta_obj < opt_.tol || delta_w < opt_.tol)
                break;

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

    // fit helpers
    void deflate_all_() const {
        for (auto& b : blocks_) b->deflate(opt_.deflation_mode);
    }
    void compute_weights_star_() {
        for (auto& b : blocks_) b->compute_weights_star();
    }
    std::vector<Vector> snapshot_weights_(const BlockRefList& blocks) const {
        std::vector<Vector> out;
        out.reserve(blocks.size());
        for (auto* b : blocks) out.push_back(b->weights().col(h_));
        return out;
    }
    double weights_variation_(const BlockRefList& blocks, const std::vector<Vector>& w_prev) const {
        double acc = 0.0;

        for (int j = 0; j < n_blocks(); ++j) {
            const auto wj = blocks[j]->weights().col(h_);
            acc += (wj - w_prev[j]).squaredNorm();
        }

        return acc;
    }

    // bootstrap
    struct AdaptiveBootstrapState {

        explicit AdaptiveBootstrapState(
            const BootstrapConfig bootstrap_config,
            const int n_threads_,
            const int h,
            const int J_
        ) : n_threads(n_threads_), J(J_) {
            seed = bootstrap_config.seed + static_cast<unsigned>(h);
            B_min = bootstrap_config.B_min;
            B_max = bootstrap_config.B_max;
            check_every = bootstrap_config.check_every;
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

    std::chrono::high_resolution_clock::time_point print_step_start_(std::string_view label) const {
        fdapde::cout << label << " --> ";
        return std::chrono::high_resolution_clock::now();
    }
    void print_step_end_(const std::chrono::high_resolution_clock::time_point start) const {
        const auto end = std::chrono::high_resolution_clock::now();
        const double elapsed_sec = std::chrono::duration<double>(end - start).count();
        fdapde::cout << "<-- " << std::fixed << std::setprecision(3) << elapsed_sec << std::defaultfloat << "s\n";
    }
    void run_component_callback_(const ComponentCallback& on_component, const Result& result) {
        if (!on_component) return;
        compute_weights_star_();
        const auto step_start = print_step_start_("Component callback");
        on_component(*this, result);
        print_step_end_(step_start);
    }

    ModelSelectionResult bootstrap_model_selection_() {
        fdapde::cout << "\n=========================================\n";
        fdapde::cout << "Bootstrap model selection for component " << h_ + 1 << '\n';
        fdapde::cout << "=========================================\n\n";

        const bool select_lambda = weight_lambda_selection_requested_();
        const std::vector<double> lambda_grid = model_selection_lambda_grid_();

        // set the number of threads
        const int n_threads = bootstrap_n_threads_();

        // original blocks
        auto blocks = main_blocks_();
        const int J = static_cast<int>(blocks.size());
        BoolMatrix C_active = C_;
        BoolMatrix C_best = C_;

        // init bootstrap
        auto step_start = print_step_start_("Init bootstrap");
        AdaptiveBootstrapState bootstrap_state(bootstrap_config_, n_threads, h_, J);
        const auto block_dims = block_dims_(blocks);
        BootstrapResult boot_results(
            h_, bootstrap_state.B_max, lambda_grid,
            block_names_(blocks), block_dims, bootstrap_config_.ci_level
        );
        print_step_end_(step_start);

        // preliminary fit
        step_start = print_step_start_("Preliminary fit");
        if (select_lambda)
            set_lambda_weights_all(lambda_grid.back());
        init_comp_(blocks);
        fit_component_(blocks, C_active);
        auto preliminary_w_fit = snapshot_weights_(blocks);
        print_step_end_(step_start);

        step_start = print_step_start_("Clone worker blocks");
        std::vector<BootstrapBlocks> thread_boot_worker(n_threads);
        for (int t = 0; t < n_threads; ++t) {
            thread_boot_worker[t] = clone_blocks_();
        }
        print_step_end_(step_start);
        auto fixed_weight_lambda = [&]() {
            for (const auto& block : blocks) {
                const double lambda = block->lambda_weights();
                if (std::isfinite(lambda)) return lambda;
            }
            return std::numeric_limits<double>::quiet_NaN();
        };
        auto print_weight_lambda = [](const double lambda) {
            fdapde::cout << "lambda=";
            if (std::isfinite(lambda))
                fdapde::cout << lambda;
            else
                fdapde::cout << "none";
        };

        int n_lambda = static_cast<int>(lambda_grid.size());
        for (int lambda_i = n_lambda - 1; lambda_i >= 0; --lambda_i) {
            ensure_bootstrap_lambda_storage_(boot_results, lambda_i, block_dims, J);
            const bool reuse_preliminary_fit = lambda_i == n_lambda - 1;

            if (select_lambda) {
                // current lambda
                const double lambda = lambda_grid[lambda_i];
                fdapde::cout << "- lambda = " << lambda << '\n';
                if (!reuse_preliminary_fit) {
                    set_lambda_weights_all(lambda);
                    for (auto& worker : thread_boot_worker)
                        set_lambda_weights_all_(worker.refs, lambda);
                }
            } else {
                fdapde::cout << "- fixed weight regularization ";
                print_weight_lambda(fixed_weight_lambda());
                fdapde::cout << '\n';
            }

            // init warm start at lambda
            std::vector<Vector> w_fit;
            if (reuse_preliminary_fit) {
                step_start = print_step_start_("  Warm-start fit (reuse preliminary)");
                w_fit = preliminary_w_fit;
            } else {
                step_start = print_step_start_("  Warm-start fit");
                const auto active_blocks = active_blocks_from_C_(C_active);
                init_comp_(blocks, InitStrategy::WarmStart, true, &active_blocks);
                fit_component_(blocks, C_active);
                w_fit = snapshot_weights_(blocks);
            }
            auto w_min = w_fit;
            if (opt_.block_deactivation)
                threshold_inactive_blocks_(w_min, C_active);
            print_step_end_(step_start);

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

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            const double elapsed_sec = duration.count() / 1000.0;

            boot_results.w_fit_by_lambda[lambda_i] = w_fit;
            boot_results.w_min_by_lambda[lambda_i] = w_min;
            step_start = print_step_start_("  Final lambda correlation");
            correlation_matrix_(blocks, w_min, boot_results.corr_min_by_lambda[lambda_i]);
            boot_results.corr_min_by_lambda[lambda_i].array() *= (C_lambda.cast<double>() + Matrix::Identity(J, J)).array();
            print_step_end_(step_start);
            if (!std::isfinite(bootstrap_state.crit))
                bootstrap_state.crit = criterion_score_with_weights_(blocks, w_min, C_);
            boot_results.criterion[lambda_i] = bootstrap_state.crit;
            boot_results.B_used_by_lambda[lambda_i] = bootstrap_state.B_done;

            if (bootstrap_state.crit > bootstrap_state.best_criterion) {
                C_best = C_lambda;
            }

            const int B_discarded =
                bootstrap_state.B_total -
                bootstrap_state.B_design -
                bootstrap_state.B_stale -
                bootstrap_state.B_cancelled -
                bootstrap_state.B_done;
            fdapde::cout << "  Bootstrap used: total=" << bootstrap_state.B_total
                         << ", design=" << bootstrap_state.B_design
                         << ", stale=" << bootstrap_state.B_stale
                         << ", cancelled=" << bootstrap_state.B_cancelled
                         << ", good=" << bootstrap_state.B_done
                         << ", discarded=" << B_discarded << '\n';
            fdapde::cout << "  Time: total=" << std::fixed << std::setprecision(3) << elapsed_sec
                         << "s, avg_fit=" << bootstrap_timing_summary.avg_fit_time()
                         << " ± " << bootstrap_timing_summary.sd_fit_time()
                         << "s, efficiency=" << std::setprecision(1)
                         << 100.0 * bootstrap_timing_summary.efficiency() << "%\n";
            fdapde::cout << "  Iters: avg=" << std::setprecision(1) << bootstrap_timing_summary.avg_iters()
                         << ", max=" << bootstrap_timing_summary.max_iters
                         << ", capped=" << bootstrap_timing_summary.capped_fits
                         << "/" << bootstrap_timing_summary.n_fits << '\n';
            fdapde::cout << "  Design: active_blocks=" << n_active_blocks
                         << ", active_connections=" << n_active_connections
                         << ", crit=" << std::setprecision(3) << bootstrap_state.crit
                         << std::defaultfloat << ", ";
            print_weight_lambda(select_lambda ? lambda_grid[lambda_i] : fixed_weight_lambda());
            fdapde::cout << std::defaultfloat << "\n\n";

            // early stop
            if (early_stop_lambda_(bootstrap_state, lambda_i, bootstrap_config_)) {
                break;
            }

            // reset connections, but keep fully deactivated blocks off
            C_active = reset_connections_keep_inactive_blocks_(C_, C_active);

        }

        step_start = print_step_start_("Resize bootstrap results");
        resize_bootstrap_results_(boot_results, static_cast<int>(lambda_grid.size()), J, block_dims);
        print_step_end_(step_start);
        step_start = print_step_start_("Compute bootstrap correlation CIs");
        compute_bootstrap_corr_cis_(boot_results, J);
        print_step_end_(step_start);

        if (bootstrap_state.best_i < 0)
            throw std::runtime_error("No model candidate was evaluated during bootstrap selection");

        boot_results.lambda_opt_index = bootstrap_state.best_i;
        if (select_lambda) {
            boot_results.lambda_opt = lambda_grid[bootstrap_state.best_i];
            fdapde::cout << "\nOptimal lambda: " << boot_results.lambda_opt << '\n';
        } else {
            fdapde::cout << "\nNo weight lambda selection requested\n";
        }


        boot_results.active_blocks = active_blocks_from_C_(C_best);
        fdapde::cout << "\nBlock deactivation:\n";
        for (int j = 0; j < J; ++j) {
            const double nrm = boot_results.w_min_by_lambda[bootstrap_state.best_i][j].norm();
            fdapde::cout << "  - block[" << std::setw(2) << j << "]: "
                         << "||w_min|| = " << std::fixed << std::setprecision(4) << nrm
                         << std::defaultfloat
                         << ", active = " << (boot_results.active_blocks[j] ? "yes" : "no")
                         << '\n';
        }
        if (J <= 20) {
            fdapde::cout << "\nUpdated design matrix:\n";
            fdapde::cout << C_best << "\n\n";
        } else {
            fdapde::cout << "\nUpdated design matrix: skipped (" << J << " blocks)\n\n";
        }

        bootstrap_selection_results_.push_back(std::move(boot_results));

        ModelSelectionResult out;
        out.lambda_selected = select_lambda;
        out.lambda = bootstrap_selection_results_.back().lambda_opt;
        out.C_active = C_best;
        return out;
    }

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
    void annotate_component_significance_(
        Result& result,
        const ComponentSignificanceResult& significance
    ) const {
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
            const auto step_start = print_step_start_("Inactive component fit");
            Result inactive_result = fit_component_(C_inactive);
            print_step_end_(step_start);
            annotate_component_significance_(inactive_result, significance);
            results.push_back(std::move(inactive_result));
        }
    }
    ComponentSignificanceResult bootstrap_test_component_significance_(const BoolMatrix& C_active) {
        ComponentSignificanceResult out;

        auto blocks = main_blocks_();
        const int J = static_cast<int>(blocks.size());
        out.rho_tot = rho_tot_(blocks, C_active);

        if (count_active_connections_(C_active) == 0 || !std::isfinite(out.rho_tot)) {
            out.B = 0;
            out.p_value = 1.0;
            out.significant = false;
            return out;
        }

        const int B = bootstrap_config_.component_significance_resamples;
        const int n_threads = bootstrap_n_threads_();
        out.B = B;

        auto step_start = print_step_start_("  Clone significance worker blocks");
        std::vector<BootstrapBlocks> thread_boot_worker(n_threads);
        for (int t = 0; t < n_threads; ++t)
            thread_boot_worker[t] = clone_blocks_();
        print_step_end_(step_start);

        std::vector<int> thread_ge_count(n_threads, 0);
        const unsigned seed = bootstrap_config_.seed + static_cast<unsigned>(1000003 * (h_ + 1));
        const auto active_blocks = active_blocks_from_C_(C_active);

        step_start = print_step_start_("  Significance bootstrap");
        parallel_for(0, B, 1, [&](int b) {
            const int tid = this_thread_id();
            auto& boot_blocks = thread_boot_worker[tid];

            set_permuted_row_index_all_(
                boot_blocks.refs,
                seed + static_cast<unsigned>(7919 * (b + 1))
            );

            init_comp_(boot_blocks.refs, InitStrategy::WarmStart, false, &active_blocks);
            fit_component_(boot_blocks.refs, C_active, false, bootstrap_fit_max_iter_());

            const double rho_star = rho_tot_(boot_blocks.refs, C_active);
            if (std::isfinite(rho_star) && rho_star >= out.rho_tot)
                ++thread_ge_count[tid];

            clear_row_index_all_(boot_blocks.refs);
        });
        print_step_end_(step_start);

        const int ge_count = std::accumulate(thread_ge_count.begin(), thread_ge_count.end(), 0);
        out.p_value = static_cast<double>(ge_count) / static_cast<double>(B);
        out.significant = out.p_value <= bootstrap_config_.component_significance_alpha;

        fdapde::cout << "\nSignificance:\n"
                     << "  - rho_tot = " << out.rho_tot << '\n'
                     << "  - p-value = " << out.p_value << '\n'
                     << "  - significant = " << out.significant << '\n';

        return out;
    }
    BootstrapSampleResult fit_bootstrap_sample_(
        BootstrapBlocks& boot_blocks,
        const int b,
        const unsigned seed,
        const BoolMatrix& C_active,
        const std::vector<Vector>& w_fit,
        const std::function<bool()>& cancelled = {}
    ) {
        BootstrapSampleResult out;
        const auto active_blocks = active_blocks_from_C_(C_active);
        const int fit_max_iter = bootstrap_fit_max_iter_();

        if (cancelled && cancelled()) {
            out.cancelled = true;
            return out;
        }

        copy_weights_snapshot_(boot_blocks.refs, w_fit);
        set_row_index_all_(boot_blocks.refs, bootstrap_index_(n_, seed + static_cast<unsigned>(b)));

        init_comp_(boot_blocks.refs, InitStrategy::WarmStart, true, &active_blocks);
        if (cancelled && cancelled()) {
            out.cancelled = true;
            clear_row_index_all_(boot_blocks.refs);
            return out;
        }

        const auto fit_start = std::chrono::high_resolution_clock::now();
        out.fit_started = true;
        const Result fit_result = fit_component_(boot_blocks.refs, C_active, true, fit_max_iter, cancelled);
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
    bool same_design_(const BoolMatrix& lhs, const BoolMatrix& rhs) const {
        if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) return false;
        for (int i = 0; i < lhs.rows(); ++i)
            for (int j = 0; j < lhs.cols(); ++j)
                if (lhs(i, j) != rhs(i, j)) return false;
        return true;
    }
    void reset_good_bootstrap_(
        AdaptiveBootstrapState& state,
        std::vector<Vector>& w_min,
        const std::vector<Vector>& w_fit,
        const BoolMatrix& C_active
    ) const {
        state.B_design += state.B_done;
        state.reset_good();
        ++state.design_epoch;
        state.design_epoch_signal.store(state.design_epoch, std::memory_order_release);
        w_min = w_fit;

        const auto active_blocks = active_blocks_from_C_(C_active);
        for (int j = 0; j < static_cast<int>(w_min.size()); ++j) {
            if (!active_blocks[j])
                w_min[j].setZero();
        }
    }
    bool connection_deactivation_ready_(const AdaptiveBootstrapState& state) const {
        if (!opt_.connection_deactivation) return false;
        if (state.B_done < bootstrap_config_.min_boots_before_connection_deactivation) return false;
        return bootstrap_config_.aggressive_connection_deactivation || state.B_done >= state.B_min;
    }
    void print_bootstrap_progress_(
        const AdaptiveBootstrapState& state,
        const int n_active_blocks,
        const int n_active_connections,
        const AdaptiveStopInfo& stop_info,
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
    void print_bootstrap_design_reset_(
        const AdaptiveBootstrapState& state,
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
                     << " | epoch=" << state.design_epoch
                     << " | ab=" << std::setw(4) << n_active_blocks
                     << ", ac=" << std::setw(6) << n_active_connections
                     << std::defaultfloat
                     << '\n';
    }
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

                BoolMatrix C_snapshot;
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

                const int tid = this_thread_id();
                auto cancelled = [&]() {
                    return stop.load(std::memory_order_acquire) ||
                        epoch != bootstrap_state.design_epoch_signal.load(std::memory_order_acquire);
                };
                auto sample = fit_bootstrap_sample_(
                    thread_boot_worker[tid], b, seed, C_snapshot, w_fit, cancelled
                );

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

                const int b_good = bootstrap_state.B_done;
                for (int j = 0; j < J; ++j) {
                    boot_results.w_boot_by_lambda[lambda_i][j].col(b_good) = sample.w[j];
                    update_w_min_(w_min[j], w_fit[j], sample.w[j]);
                }
                boot_results.corr_boot_by_lambda[lambda_i].col(b_good) =
                    Eigen::Map<const Vector>(sample.corr.data(), J * J);

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

                if (opt_.block_deactivation) {
                    const BoolMatrix C_before = C_active;
                    threshold_inactive_blocks_(w_min, C_active);
                    if (!same_design_(C_before, C_active)) {
                        reset_good_bootstrap_(bootstrap_state, w_min, w_fit, C_active);
                        const int n_active_blocks = count_active_blocks_(C_active);
                        const int n_active_connections = count_active_connections_(C_active);
                        print_bootstrap_design_reset_(
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

                const bool force_check = bootstrap_state.B_done >= bootstrap_state.B_max;
                const bool due_check =
                    force_check ||
                    bootstrap_state.B_done - bootstrap_state.last_check_B_done >= bootstrap_state.check_every;
                if (!due_check) continue;
                bootstrap_state.last_check_B_done = bootstrap_state.B_done;

                if (connection_deactivation_ready_(bootstrap_state)) {
                    const BoolMatrix C_before = C_active;
                    threshold_inactive_connections_(lambda_i, bootstrap_state, boot_results, C_active);
                    if (!same_design_(C_before, C_active)) {
                        reset_good_bootstrap_(bootstrap_state, w_min, w_fit, C_active);
                        const int n_active_blocks = count_active_blocks_(C_active);
                        const int n_active_connections = count_active_connections_(C_active);
                        print_bootstrap_design_reset_(
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

                bootstrap_state.crit = criterion_score_with_weights_(blocks, w_min, C_);
                const AdaptiveStopInfo stop_info = bootstrap_config_.adaptive ?
                    adaptive_stop_(bootstrap_state, bootstrap_config_) : AdaptiveStopInfo{};
                print_bootstrap_progress_(
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
    void ensure_bootstrap_lambda_storage_(
        BootstrapResult& boot_results,
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

    AdaptiveStopInfo adaptive_stop_(AdaptiveBootstrapState& state, const BootstrapConfig& config) const {
        AdaptiveStopInfo out;
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

    bool early_stop_lambda_(
        AdaptiveBootstrapState& state,
        int lambda_i,
        const BootstrapConfig& config
    ) const {
        if (state.crit > state.best_criterion) {
            state.best_criterion = state.crit;
            state.best_i = lambda_i;
            state.no_improve = 0;
        } else {
            ++state.no_improve;
        }

        if (state.no_improve >= config.patience) {
            fdapde::cout << "  early stop: no improvement for "
                         << config.patience
                         << " consecutive lambdas\n";
            return true;
        }

        return false;
    }

    // bootstrap utils
    int threshold_inactive_blocks_(std::vector<Vector>& w_min, BoolMatrix& C_active) const {
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
    /*
    int threshold_inactive_connections_(BoolMatrix& C_active) const {
        return count_active_connections_(C_active);
    }
    */
    int threshold_inactive_connections_(
        int lambda_i,
        AdaptiveBootstrapState& state,
        const BootstrapResult& boot_results,
        BoolMatrix& C_active,
        const int B_eff_override = -1
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
    BoolMatrix reset_connections_keep_inactive_blocks_(const BoolMatrix& C_full, const BoolMatrix& C_current) const {
        BoolMatrix C_reset = C_full;

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
    std::pair<double, double> fisher_z_corr_ci_(
        const Matrix& corr_boot,
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
    void deactivate_isolated_blocks_(BoolMatrix& C_active) const {
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
    int count_active_connections_(const BoolMatrix& C_active) const {
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
    int count_active_blocks_(const BoolMatrix& C_active) const {
        const auto active_blocks = active_blocks_from_C_(C_active);
        return static_cast<int>(std::count(active_blocks.begin(), active_blocks.end(), true));
    }
    std::vector<bool> active_blocks_from_C_(const BoolMatrix& C_active) const {
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
    void update_w_min_(Vector& w_min_j, const Vector& w_fit_j, const Vector& w_bj) const {
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
    void resize_bootstrap_results_(
        BootstrapResult& boot_results,
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
    void compute_bootstrap_corr_cis_(BootstrapResult& boot_results, int J) const {
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


    void set_row_index_all_(const BlockRefList& blocks, const typename Block::IndexVector& idx) {
        for (auto* b : blocks) b->set_row_index(idx);
    }
    void set_permuted_row_index_all_(const BlockRefList& blocks, const unsigned seed) {
        for (int j = 0; j < static_cast<int>(blocks.size()); ++j) {
            std::mt19937_64 rng(seed + static_cast<unsigned>(104729 * (j + 1)));
            blocks[j]->set_row_index(permutation_indices_(blocks[j]->n_raw(), rng));
        }
    }
    void clear_row_index_all_(const BlockRefList& blocks) {
        for (auto* b : blocks) b->clear_row_index();
    }

    typename Block::IndexVector bootstrap_index_(int n, unsigned seed) const {
        std::mt19937_64 rng(seed);
        return bootstrap_indices_(n, rng);
    }
    typename Block::IndexVector bootstrap_indices_(int n, std::mt19937_64& rng) const {

        switch (bootstrap_config_.resampling_strategy) {
            case ResamplingStrategy::Ordinary:
                return ordinary_bootstrap_indices_(n, rng);
            case ResamplingStrategy::Stationary:
                return stationary_bootstrap_indices_(n, bootstrap_config_.stationary_block_length, rng);
        }

        throw std::logic_error("unsupported resampling strategy");
    }
    typename Block::IndexVector ordinary_bootstrap_indices_(const int n, std::mt19937_64& rng) const {

        if (n <= 0)
            throw std::invalid_argument("n must be positive");


        std::uniform_int_distribution<int> U(0, n - 1);

        typename Block::IndexVector idx(n);
        for (int i = 0; i < n; ++i)
            idx(i) = U(rng);

        return idx;
    }
    typename Block::IndexVector permutation_indices_(const int n, std::mt19937_64& rng) const {
        if (n <= 0)
            throw std::invalid_argument("n must be positive");

        typename Block::IndexVector idx(n);
        std::iota(idx.data(), idx.data() + idx.size(), 0);
        std::shuffle(idx.data(), idx.data() + idx.size(), rng);

        return idx;
    }
    typename Block::IndexVector stationary_bootstrap_indices_(const int n, const double mean_block_length, std::mt19937_64& rng) const {
        if (n <= 0)
            throw std::invalid_argument("n must be positive");

        if (!(mean_block_length > 0.0) || !std::isfinite(mean_block_length))
            throw std::invalid_argument("stationary block length must be positive");

        const double p = std::clamp(1.0 / mean_block_length, 0.0, 1.0);

        std::uniform_int_distribution<int> U_index(0, n - 1);
        std::bernoulli_distribution start_new_block(p);

        typename Block::IndexVector idx(n);

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

    void validate_bootstrap_support_() const {
        if constexpr (std::same_as<SamplingStrategy, TimeDependentSampling>) {
            throw std::runtime_error(
                "RGCCA: bootstrap is not supported "
                "for rgcca::TimeDependentSampling"
            );
        }
    }

    void validate_bootstrap_config_() const {
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
        if (bootstrap_config_.fit_max_iter == 0 || bootstrap_config_.fit_max_iter < -1)
            throw std::invalid_argument("RGCCA: bootstrap fit_max_iter must be positive or -1");
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
        if (bootstrap_config_.min_boots_before_connection_deactivation < 0) {
            throw std::invalid_argument(
                "RGCCA: bootstrap min_boots_before_connection_deactivation must be nonnegative"
            );
        }
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
            bootstrap_config_.resampling_strategy == ResamplingStrategy::Stationary &&
            (!(bootstrap_config_.stationary_block_length > 0.0) ||
             !std::isfinite(bootstrap_config_.stationary_block_length))
        ) {
            throw std::invalid_argument("RGCCA: bootstrap stationary_block_length must be finite and positive");
        }
    }
    void validate_component_significance_config_() const {
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
    }

    int bootstrap_fit_max_iter_() const {
        return bootstrap_config_.fit_max_iter > 0 ? bootstrap_config_.fit_max_iter : opt_.max_iter;
    }

    int bootstrap_n_threads_() const {
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

    void validate_lambda_grid_weights_() const {
        if (static_cast<int>(lambda_grid_weights_.size()) != n_comp_) {
            throw std::runtime_error(
                "RGCCA: automatic weight lambda selection requires set_lambda_grid_weights(...) "
                "with one grid or one grid per component"
            );
        }

        for (int h = 0; h < n_comp_; ++h) {
            const auto& grid = lambda_grid_weights_[h];

            if (grid.empty()) {
                throw std::runtime_error("RGCCA: weight lambda grid contains an empty component grid");
            }

            for (std::size_t i = 0; i < grid.size(); ++i) {
                const double lambda = grid[i];

                if (!(lambda > 0.0) || !std::isfinite(lambda)) {
                    throw std::runtime_error("RGCCA: weight lambda grid values must be finite and positive");
                }

                if (i > 0 && grid[i] < grid[i - 1]) {
                    throw std::runtime_error("RGCCA: weight lambda grid must be sorted in nondecreasing order");
                }
            }
        }
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

            // Evaluation normalizes with M, not Omega, so the penalty does not affect the reported component.
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

                if (std::string_view(opt_.scheme.name) == "Horst") num += corr_jk;
                else num += std::abs(corr_jk);
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

#endif // __FDAPDE_RGCCA_MODEL_H__

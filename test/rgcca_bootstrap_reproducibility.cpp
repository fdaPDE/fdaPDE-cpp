#include <fdaPDE/models.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace fdapde;
using namespace fdapde::rgcca;

namespace {

using DenseMatrix = fdapde::rgcca::Matrix;

DenseMatrix make_block(const int n, const double phase, const double noise_scale) {
    DenseMatrix out(n, 3);
    for (int i = 0; i < n; ++i) {
        const double x = static_cast<double>(i);
        const double signal = std::sin(0.17 * x) + 0.5 * std::cos(0.07 * x);
        out(i, 0) = signal + noise_scale * std::sin(0.53 * x + phase);
        out(i, 1) = 0.7 * signal + noise_scale * std::cos(0.31 * x + phase);
        out(i, 2) = 0.3 * signal + noise_scale * std::sin(0.89 * x + phase);
    }
    return out;
}

template <typename Derived>
void print_matrix(const Eigen::MatrixBase<Derived>& matrix) {
    std::cout << matrix.rows() << ' ' << matrix.cols();
    for (int j = 0; j < matrix.cols(); ++j)
        for (int i = 0; i < matrix.rows(); ++i)
            std::cout << ' ' << matrix(i, j);
    std::cout << '\n';
}

template <typename T>
void print_vector(const std::vector<T>& values) {
    std::cout << values.size();
    for (const auto& value : values)
        std::cout << ' ' << value;
    std::cout << '\n';
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2)
        throw std::invalid_argument("usage: rgcca_bootstrap_reproducibility THREADS");

    const int threads = std::stoi(argv[1]);
    constexpr int n = 64;

    RGCCA<IndependentSampling>::Options options;
    options.mode = Mode::CovMax;
    options.init_strategy = InitStrategy::Uniform;
    options.weight_sign_constraint = WeightSignConstraint::NonNegative;
    options.block_deactivation = true;
    options.connection_deactivation = true;
    options.component_significance = true;
    options.block_importance = true;
    options.max_iter = 40;
    options.tol = 1e-10;

    RGCCA<IndependentSampling> model(n, options, 1);
    model.add_multivariate_block("signal_1", make_block(n, 0.1, 0.08));
    model.add_multivariate_block("signal_2", make_block(n, 0.5, 0.10));
    model.add_multivariate_block("signal_3", make_block(n, 1.1, 0.18));
    model.add_multivariate_block("weak", make_block(n, 2.0, 0.75));

    RGCCA<IndependentSampling>::BootstrapConfig bootstrap;
    bootstrap.seed = 90210;
    bootstrap.max_threads = threads;
    bootstrap.B_min = 12;
    bootstrap.B_max = 24;
    bootstrap.check_every = 6;
    bootstrap.check_every_block_deactivation = 4;
    bootstrap.check_every_connection_deactivation = 6;
    bootstrap.stable_checks_required = 2;
    bootstrap.active_block_tol = 0.95;
    bootstrap.active_connection_sign_stability = 0.8;
    bootstrap.active_connection_min_abs_corr = 0.02;
    bootstrap.fit_max_iter = 20;
    bootstrap.component_significance_resamples = 24;
    bootstrap.component_significance_alpha = 0.5;
    bootstrap.block_importance_resamples = 24;
    bootstrap.block_importance_alpha = 0.05;
    model.set_bootstrap_config(bootstrap);

    const auto results = model.fit();
    const auto& selections = model.bootstrap_selection_results();
    if (results.size() != 1 || selections.size() != 1)
        throw std::runtime_error("unexpected RGCCA result count");

    const auto& result = results.front();
    const auto& selection = selections.front();
    const int lambda = selection.lambda_opt_index;
    if (lambda < 0)
        throw std::runtime_error("bootstrap did not select a candidate");
    if (selection.design_epochs_by_lambda[lambda] < 2)
        throw std::runtime_error("test data did not exercise deterministic design epochs");
    if (selection.B_stale_by_lambda[lambda] != 0 ||
        selection.B_cancelled_by_lambda[lambda] != 0)
        throw std::runtime_error("control-point barriers produced stale or cancelled fits");
    if (selection.B_total_by_lambda[lambda] !=
        selection.B_design_by_lambda[lambda] + selection.B_used_by_lambda[lambda])
        throw std::runtime_error("control-point barrier accounting discarded fitted candidates");
    if (result.rho_tot_bootstrap_count != bootstrap.component_significance_resamples ||
        result.block_importance_bootstrap_count != bootstrap.block_importance_resamples)
        throw std::runtime_error("test data did not exercise significance diagnostics");

    const auto& candidate_ids = selection.candidate_ids_by_lambda[lambda];
    for (std::size_t i = 1; i < candidate_ids.size(); ++i) {
        if (candidate_ids[i] != candidate_ids[i - 1] + 1)
            throw std::runtime_error("accepted bootstrap candidate IDs are not contiguous");
    }

    std::cout << std::setprecision(17);
    print_matrix(result.C.cast<int>());
    std::cout << selection.B_used_by_lambda[lambda] << ' '
              << selection.B_total_by_lambda[lambda] << ' '
              << selection.B_design_by_lambda[lambda] << ' '
              << selection.B_stale_by_lambda[lambda] << ' '
              << selection.B_cancelled_by_lambda[lambda] << ' '
              << selection.design_epochs_by_lambda[lambda] << ' '
              << selection.criterion[lambda] << '\n';
    for (const int id : candidate_ids)
        std::cout << id << ' ';
    std::cout << '\n';
    for (const auto& weights : selection.w_min_by_lambda[lambda])
        print_matrix(weights);
    print_matrix(selection.corr_min_by_lambda[lambda]);
    print_matrix(selection.corr_ci_low_by_lambda[lambda]);
    print_matrix(selection.corr_ci_high_by_lambda[lambda]);
    std::cout << result.rho_tot << ' '
              << result.rho_tot_p_value << ' '
              << result.rho_tot_bootstrap_count << ' '
              << result.component_significant << '\n';
    print_vector(result.block_importance);
    print_vector(result.block_importance_p_values);
    print_vector(result.block_importance_significant);
}

// Standalone FRGCCA CovMax timing driver.
//
// Build:
//   cmake --build test/build --target frgcca_covmax_timing
//
// Run from test/build:
//   ./frgcca_covmax_timing

#include <fdaPDE/models.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "ipopt_options.h"

using namespace fdapde;
using namespace fdapde::rgcca;

namespace {

using DenseMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

constexpr int n_blocks = 4;
constexpr int n_comp = 3;

struct Config {
    int n = 3200;
    int n_locs = 426;
    int max_iter = 1000;
    unsigned seed = 0;
    double lambda_weights = 1e-3;
    double tol = 1e-8;
    double sigma_noise = 2.0;
    InitStrategy init_strategy = InitStrategy::Uniform;
    bool fully_connected = false;
    bool nonnegative = false;
    bool bootstrap_lambda_selection = false;
    bool block_deactivation = false;
    bool connection_deactivation = false;
    bool aggressive_connection_deactivation = false;
    bool verbose = false;
    int bootstrap_threads = 1;
    int bootstrap_B_min = 500;
    int bootstrap_B_max = 1000;
    int bootstrap_check_every = 50;
    int min_boots_before_connection_deactivation = 100;
    int bootstrap_patience = 1;
    int bootstrap_fit_max_iter = -1;
    int inspect_bootstrap = -1;
};

class Timeline {
   public:
    using steady_clock = std::chrono::steady_clock;

    Timeline() : start_(steady_clock::now()), last_(start_) {}

    void mark(std::string_view label) {
        print(label, steady_clock::now(), 0.0);
    }

    class Section {
       public:
        Section(Timeline& timeline, std::string label) :
            timeline_(timeline), label_(std::move(label)), start_(steady_clock::now()) {
            timeline_.print("begin " + label_, start_, 0.0);
        }
        ~Section() {
            const auto end = steady_clock::now();
            timeline_.print("end " + label_, end, std::chrono::duration<double>(end - start_).count());
        }

       private:
        Timeline& timeline_;
        std::string label_;
        steady_clock::time_point start_;
    };

   private:
    void print(std::string_view label, steady_clock::time_point now, double section_seconds) {
        const auto wall_now = std::chrono::system_clock::now();
        const std::time_t wall_time = std::chrono::system_clock::to_time_t(wall_now);
        const double total = std::chrono::duration<double>(now - start_).count();
        const double delta = std::chrono::duration<double>(now - last_).count();

        std::cout << "[" << std::put_time(std::localtime(&wall_time), "%H:%M:%S") << "]"
                  << "[+" << std::fixed << std::setw(9) << std::setprecision(3) << total << "s]"
                  << "[dt " << std::setw(9) << delta << "s]";
        if (section_seconds > 0.0)
            std::cout << "[section " << std::setw(9) << section_seconds << "s]";
        std::cout << " " << label << std::defaultfloat << std::endl;

        last_ = now;
    }

    steady_clock::time_point start_;
    steady_clock::time_point last_;
};

void usage(const char* argv0) {
    std::cout
      << "Usage: " << argv0 << " [options]\n\n"
      << "Options:\n"
      << "  --n N                  observations, default 3200\n"
      << "  --n-locs P             functional locations, default 426\n"
      << "  --lambda VALUE         fixed weight lambda, default 1e-3\n"
      << "  --max-iter N           RGCCA max iterations, default 1000\n"
      << "  --tol VALUE            RGCCA tolerance, default 1e-8\n"
      << "  --init uniform|svd     RGCCA init strategy, default uniform\n"
      << "  --seed N               data seed, default 0\n"
      << "  --sigma-noise VALUE    noise standard deviation, default 2\n"
      << "  --fully-connected      use all pairwise block connections\n"
      << "  --nn                   use non-negative weight constraint\n"
      << "  --bootstrap-lambda     use bootstrap weight-lambda selection\n"
      << "  --block-deactivation   enable bootstrap block deactivation\n"
      << "  --connection-deactivation enable bootstrap connection deactivation\n"
      << "  --aggressive-connection-deactivation apply connection deactivation at streaming checks\n"
      << "  --verbose              print RGCCA diagnostics\n"
      << "  --bootstrap-threads N  bootstrap max threads, default 1\n"
      << "  --bootstrap-B-min N    bootstrap B_min, default 500\n"
      << "  --bootstrap-B-max N    bootstrap B_max, default 1000\n"
      << "  --bootstrap-check-every N bootstrap streaming check interval, default 50\n"
      << "  --min-boots-before-connection-deactivation N minimum good boots before connection deactivation, default 100\n"
      << "  --bootstrap-patience N bootstrap lambda early-stop patience, default 1\n"
      << "  --bootstrap-fit-max-iter N bootstrap refit max iterations, default -1 uses --max-iter\n"
      << "  --inspect-bootstrap B print ordinary-bootstrap multiplicity stats for resample B and exit\n"
      << "  --help                 show this message\n";
}

Config parse_args(int argc, char** argv) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto next = [&]() -> std::string {
            if (++i >= argc) throw std::invalid_argument("missing value for " + arg);
            return argv[i];
        };

        if (arg == "--n") cfg.n = std::stoi(next());
        else if (arg == "--n-locs") cfg.n_locs = std::stoi(next());
        else if (arg == "--lambda") cfg.lambda_weights = std::stod(next());
        else if (arg == "--max-iter") cfg.max_iter = std::stoi(next());
        else if (arg == "--tol") cfg.tol = std::stod(next());
        else if (arg == "--init") {
            const std::string value = next();
            if (value == "uniform") cfg.init_strategy = InitStrategy::Uniform;
            else if (value == "svd") cfg.init_strategy = InitStrategy::SVD;
            else throw std::invalid_argument("--init must be 'uniform' or 'svd'");
        }
        else if (arg == "--seed") cfg.seed = static_cast<unsigned>(std::stoul(next()));
        else if (arg == "--sigma-noise") cfg.sigma_noise = std::stod(next());
        else if (arg == "--fully-connected") cfg.fully_connected = true;
        else if (arg == "--nn") cfg.nonnegative = true;
        else if (arg == "--bootstrap-lambda") cfg.bootstrap_lambda_selection = true;
        else if (arg == "--block-deactivation") cfg.block_deactivation = true;
        else if (arg == "--connection-deactivation") cfg.connection_deactivation = true;
        else if (arg == "--aggressive-connection-deactivation") cfg.aggressive_connection_deactivation = true;
        else if (arg == "--verbose") cfg.verbose = true;
        else if (arg == "--bootstrap-threads") cfg.bootstrap_threads = std::stoi(next());
        else if (arg == "--bootstrap-B-min") cfg.bootstrap_B_min = std::stoi(next());
        else if (arg == "--bootstrap-B-max") cfg.bootstrap_B_max = std::stoi(next());
        else if (arg == "--bootstrap-check-every") cfg.bootstrap_check_every = std::stoi(next());
        else if (arg == "--min-boots-before-connection-deactivation")
            cfg.min_boots_before_connection_deactivation = std::stoi(next());
        else if (arg == "--bootstrap-patience") cfg.bootstrap_patience = std::stoi(next());
        else if (arg == "--bootstrap-fit-max-iter") cfg.bootstrap_fit_max_iter = std::stoi(next());
        else if (arg == "--inspect-bootstrap") cfg.inspect_bootstrap = std::stoi(next());
        else if (arg == "--help") {
            usage(argv[0]);
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown option: " + arg);
        }
    }

    if (cfg.n <= 0) throw std::invalid_argument("--n must be positive");
    if (cfg.n_locs <= 1) throw std::invalid_argument("--n-locs must be greater than 1");
    if (cfg.max_iter <= 0) throw std::invalid_argument("--max-iter must be positive");
    if (!(cfg.lambda_weights > 0.0) || !std::isfinite(cfg.lambda_weights))
        throw std::invalid_argument("--lambda must be finite and positive");
    if (!(cfg.tol > 0.0) || !std::isfinite(cfg.tol))
        throw std::invalid_argument("--tol must be finite and positive");
    if (!(cfg.sigma_noise >= 0.0) || !std::isfinite(cfg.sigma_noise))
        throw std::invalid_argument("--sigma-noise must be finite and non-negative");
    if (cfg.bootstrap_threads <= 0) throw std::invalid_argument("--bootstrap-threads must be positive");
    if (cfg.bootstrap_B_min <= 0) throw std::invalid_argument("--bootstrap-B-min must be positive");
    if (cfg.bootstrap_B_max <= 0) throw std::invalid_argument("--bootstrap-B-max must be positive");
    if (cfg.bootstrap_B_min > cfg.bootstrap_B_max)
        throw std::invalid_argument("--bootstrap-B-min cannot exceed --bootstrap-B-max");
    if (cfg.bootstrap_check_every <= 0) throw std::invalid_argument("--bootstrap-check-every must be positive");
    if (cfg.min_boots_before_connection_deactivation < 0)
        throw std::invalid_argument("--min-boots-before-connection-deactivation must be non-negative");
    if (cfg.bootstrap_patience <= 0) throw std::invalid_argument("--bootstrap-patience must be positive");
    if (cfg.bootstrap_fit_max_iter == 0 || cfg.bootstrap_fit_max_iter < -1)
        throw std::invalid_argument("--bootstrap-fit-max-iter must be positive or -1");
    if (cfg.inspect_bootstrap < -1)
        throw std::invalid_argument("--inspect-bootstrap must be non-negative or -1");
    if (cfg.inspect_bootstrap >= cfg.bootstrap_B_max)
        throw std::invalid_argument("--inspect-bootstrap must be smaller than --bootstrap-B-max");

    return cfg;
}

struct BootstrapIndexStats {
    int unique = 0;
    int omitted = 0;
    int max_count = 0;
    int count_ge_2 = 0;
    int count_ge_5 = 0;
    int count_ge_10 = 0;
};

BootstrapIndexStats ordinary_bootstrap_stats(int n, unsigned seed) {
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int> U(0, n - 1);
    std::vector<int> counts(n, 0);

    for (int i = 0; i < n; ++i)
        ++counts[U(rng)];

    BootstrapIndexStats stats;
    for (const int count : counts) {
        stats.max_count = std::max(stats.max_count, count);
        if (count == 0) ++stats.omitted;
        else ++stats.unique;
        if (count >= 2) ++stats.count_ge_2;
        if (count >= 5) ++stats.count_ge_5;
        if (count >= 10) ++stats.count_ge_10;
    }
    return stats;
}

void inspect_bootstrap_indices(const Config& cfg) {
    constexpr unsigned bootstrap_seed = 12345; // RGCCA::BootstrapConfig default, component 1 adds h=0.
    const auto inspected = ordinary_bootstrap_stats(
        cfg.n,
        bootstrap_seed + static_cast<unsigned>(cfg.inspect_bootstrap)
    );

    double avg_max_count = 0.0;
    int global_max_count = 0;
    int global_max_b = 0;
    for (int b = 0; b < cfg.bootstrap_B_max; ++b) {
        const auto stats = ordinary_bootstrap_stats(cfg.n, bootstrap_seed + static_cast<unsigned>(b));
        avg_max_count += static_cast<double>(stats.max_count);
        if (stats.max_count > global_max_count) {
            global_max_count = stats.max_count;
            global_max_b = b;
        }
    }
    avg_max_count /= static_cast<double>(cfg.bootstrap_B_max);

    std::cout << "Ordinary bootstrap index inspection\n"
              << "  n                  = " << cfg.n << "\n"
              << "  B_max              = " << cfg.bootstrap_B_max << "\n"
              << "  seed               = " << bootstrap_seed << "\n"
              << "  avg max multiplicity = " << std::fixed << std::setprecision(3) << avg_max_count << "\n"
              << "  worst max multiplicity = " << global_max_count << " at b=" << global_max_b << "\n"
              << "  inspected b        = " << cfg.inspect_bootstrap << "\n"
              << "  unique             = " << inspected.unique << "\n"
              << "  omitted            = " << inspected.omitted << "\n"
              << "  max multiplicity   = " << inspected.max_count << "\n"
              << "  rows repeated >=2  = " << inspected.count_ge_2 << "\n"
              << "  rows repeated >=5  = " << inspected.count_ge_5 << "\n"
              << "  rows repeated >=10 = " << inspected.count_ge_10 << std::defaultfloat << std::endl;
}

std::vector<double> bootstrap_lambda_grid() {
    return {1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2};
}

const char* init_strategy_name(InitStrategy init_strategy) {
    switch (init_strategy) {
        case InitStrategy::Uniform: return "Uniform";
        case InitStrategy::SVD: return "SVD";
        default: return "Other";
    }
}

double loading(double x, int id) {
    switch (id) {
        case 0: return 0.0;
        case 1: return std::exp(-80.0 * std::pow(x - 0.25, 2));
        case 2: return std::exp(-80.0 * std::pow(x - 0.75, 2));
        case 3: return std::exp(-80.0 * std::pow(x - 0.50, 2));
        default: throw std::logic_error("unknown loading id");
    }
}

DenseMatrix loading_matrix(int n_locs, std::array<int, n_comp> ids) {
    DenseMatrix A(n_locs, n_comp);
    for (int r = 0; r < n_locs; ++r) {
        const double x = static_cast<double>(r) / static_cast<double>(n_locs - 1);
        for (int h = 0; h < n_comp; ++h)
            A(r, h) = loading(x, ids[h]);
    }
    return A;
}

DenseMatrix centered_noise(int n, int p, double sd, std::mt19937_64& rng) {
    std::normal_distribution<double> normal(0.0, sd);
    DenseMatrix E(n, p);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < p; ++j)
            E(i, j) = normal(rng);

    E.rowwise() -= E.colwise().mean();
    return E;
}

std::vector<DenseMatrix> generate_blocks(const Config& cfg) {
    DenseMatrix scale = DenseMatrix::Identity(12, 12);
    DenseMatrix corr = DenseMatrix::Identity(12, 12);

    for (int i = 0; i < 4; ++i) scale(i, i) = 1.5;
    for (int i = 4; i < 8; ++i) scale(i, i) = 1.2;

    constexpr double rho = 0.9;
    corr(0, 2) = corr(2, 0) = rho;
    corr(1, 3) = corr(3, 1) = -rho;
    corr(4, 5) = corr(5, 4) = -rho * 0.9;
    corr(6, 7) = corr(7, 6) = rho * 0.9;
    corr(8, 11) = corr(11, 8) = rho * 0.8;

    const DenseMatrix sigma = scale * corr * scale;
    Eigen::LLT<DenseMatrix> llt(sigma);
    if (llt.info() != Eigen::Success)
        throw std::runtime_error("latent covariance is not positive definite");

    std::mt19937_64 rng(cfg.seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    DenseMatrix Z(cfg.n, 12);
    for (int i = 0; i < Z.rows(); ++i)
        for (int j = 0; j < Z.cols(); ++j)
            Z(i, j) = normal(rng);

    DenseMatrix HH = Z * DenseMatrix(llt.matrixL()).transpose();
    HH.rowwise() -= HH.colwise().mean();
    HH.col(9).setZero();
    HH.col(10).setZero();

    std::array<DenseMatrix, n_blocks> H {
        DenseMatrix(cfg.n, n_comp), DenseMatrix(cfg.n, n_comp), DenseMatrix(cfg.n, n_comp), DenseMatrix(cfg.n, n_comp)
    };
    for (int g = 0; g < n_blocks; ++g) {
        H[g].col(0) = HH.col(g);
        H[g].col(1) = HH.col(4 + g);
        H[g].col(2) = HH.col(8 + g);
    }

    const std::array<DenseMatrix, n_blocks> A {
        loading_matrix(cfg.n_locs, {1, 2, 3}),
        loading_matrix(cfg.n_locs, {2, 1, 0}),
        loading_matrix(cfg.n_locs, {1, 3, 0}),
        loading_matrix(cfg.n_locs, {1, 3, 2})
    };

    std::mt19937_64 noise_rng(cfg.seed + 1000);
    std::vector<DenseMatrix> X;
    X.reserve(n_blocks);
    for (int g = 0; g < n_blocks; ++g)
        X.push_back(H[g] * A[g].transpose() + centered_noise(cfg.n, cfg.n_locs, cfg.sigma_noise, noise_rng));

    return X;
}

void connect_reference_design(RGCCA<IndependentSampling>& rgcca) {
    rgcca.connect(0, 1);
    rgcca.connect(0, 2);
    rgcca.connect(0, 3);
    rgcca.connect(1, 3);
    rgcca.connect(2, 3);
}

void connect_fully_connected(RGCCA<IndependentSampling>& rgcca) {
    for (int j = 0; j < n_blocks; ++j)
        for (int k = j + 1; k < n_blocks; ++k)
            rgcca.connect(j, k);
}

void add_fem_functional_block(
  RGCCA<IndependentSampling>& rgcca, DenseMatrix&& X, int block_id, int n_locs, Timeline& timeline) {
    const std::string name = "X" + std::to_string(block_id);

    {
        Timeline::Section section(timeline, name + " mesh/geoframe/penalty");
        Triangulation<1, 1> D = Triangulation<1, 1>::UnitInterval(n_locs);
        GeoFrame data(D);
        auto& level = data.insert_scalar_layer<POINT>("data", MESH_NODES);
        level.load_blk(name, X.transpose());

        FeSpace Vh(D, P1<1>);
        TrialFunction f(Vh);
        TestFunction v(Vh);
        auto a = integral(D)(dot(grad(f), grad(v)));
        ZeroField<1> u;
        auto F = integral(D)(u * v);

        Timeline::Section add_section(timeline, name + " add_functional_block");
        rgcca.add_functional_block(name, data, std::move(X), fe_normcovmax_elliptic(a, F));
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Config cfg = parse_args(argc, argv);
        if (cfg.inspect_bootstrap >= 0) {
            inspect_bootstrap_indices(cfg);
            return 0;
        }

        const bool run_bootstrap = cfg.bootstrap_lambda_selection || cfg.block_deactivation || cfg.connection_deactivation;
        if (cfg.nonnegative) write_ipopt_options();
        Timeline timeline;

        std::cout << "FRGCCA CovMax timing run\n"
                  << "  blocks             = " << n_blocks << "\n"
                  << "  components         = " << n_comp << "\n"
                  << "  observations       = " << cfg.n << "\n"
                  << "  functional locs    = " << cfg.n_locs << "\n"
                  << "  lambda_weights     = " << cfg.lambda_weights << "\n"
                  << "  max_iter           = " << cfg.max_iter << "\n"
                  << "  tol                = " << cfg.tol << "\n"
                  << "  init_strategy      = " << init_strategy_name(cfg.init_strategy) << "\n"
                  << "  design             = " << (cfg.fully_connected ? "fully connected" : "test reference") << "\n"
                  << "  weight constraint  = " << (cfg.nonnegative ? "NonNegative" : "None") << "\n"
                  << "  lambda selection   = " << (cfg.bootstrap_lambda_selection ? "bootstrap" : "none") << "\n";
        if (run_bootstrap) {
            std::cout << "  bootstrap threads  = " << cfg.bootstrap_threads << "\n"
                      << "  bootstrap B_min    = " << cfg.bootstrap_B_min << "\n"
                      << "  bootstrap B_max    = " << cfg.bootstrap_B_max << "\n"
                      << "  bootstrap check    = " << cfg.bootstrap_check_every << "\n"
                      << "  conn min boots     = " << cfg.min_boots_before_connection_deactivation << "\n"
                      << "  bootstrap patience = " << cfg.bootstrap_patience << "\n"
                      << "  bootstrap fit iter = " << cfg.bootstrap_fit_max_iter << "\n"
                      << "  verbose            = " << cfg.verbose << "\n"
                      << "  block deactivation = " << cfg.block_deactivation << "\n"
                      << "  conn deactivation  = " << cfg.connection_deactivation << "\n";
            std::cout << "  aggressive conn    = " << cfg.aggressive_connection_deactivation << "\n";
        }
        std::cout << "\n";

        std::vector<DenseMatrix> blocks;
        {
            Timeline::Section section(timeline, "generate synthetic functional blocks");
            blocks = generate_blocks(cfg);
        }

        RGCCA<IndependentSampling>::Options options;
        options.mode = Mode::CovMax;
        options.max_iter = cfg.max_iter;
        options.tol = cfg.tol;
        options.init_strategy = cfg.init_strategy;
        options.weight_sign_constraint = cfg.nonnegative ? WeightSignConstraint::NonNegative : WeightSignConstraint::None;
        options.lambda_selection_weights = cfg.bootstrap_lambda_selection ? LambdaSelection::Automatic : LambdaSelection::Manual;
        options.lambda_selection_components = LambdaSelection::Manual;
        options.block_deactivation = cfg.block_deactivation;
        options.connection_deactivation = cfg.connection_deactivation;
        options.verbose = cfg.verbose;

        RGCCA<IndependentSampling> rgcca(cfg.n, options, n_comp);

        for (int j = 0; j < n_blocks; ++j) {
            Timeline::Section section(timeline, "build/add block X" + std::to_string(j + 1));
            add_fem_functional_block(rgcca, std::move(blocks[j]), j + 1, cfg.n_locs, timeline);
        }

        {
            Timeline::Section section(timeline, "connect design and set fixed lambdas");
            if (cfg.fully_connected)
                connect_fully_connected(rgcca);
            else
                connect_reference_design(rgcca);
            if (cfg.bootstrap_lambda_selection)
                rgcca.set_lambda_grid_weights(bootstrap_lambda_grid());
            else
                rgcca.set_lambda_weights_all(cfg.lambda_weights);
            rgcca.set_lambda_components_all(1e-15);
        }

        if (run_bootstrap) {
            RGCCA<IndependentSampling>::BootstrapConfig bootstrap_config;
            bootstrap_config.max_threads = cfg.bootstrap_threads;
            bootstrap_config.B_min = cfg.bootstrap_B_min;
            bootstrap_config.B_max = cfg.bootstrap_B_max;
            bootstrap_config.check_every = cfg.bootstrap_check_every;
            bootstrap_config.patience = cfg.bootstrap_patience;
            bootstrap_config.fit_max_iter = cfg.bootstrap_fit_max_iter;
            bootstrap_config.active_block_tol = 0.1;
            bootstrap_config.aggressive_connection_deactivation = cfg.aggressive_connection_deactivation;
            bootstrap_config.min_boots_before_connection_deactivation =
                cfg.min_boots_before_connection_deactivation;
            rgcca.set_bootstrap_config(bootstrap_config);
        }

        std::vector<Result> results;
        {
            Timeline::Section section(timeline, "rgcca.fit");
            results = rgcca.fit();
        }

        {
            Timeline::Section section(timeline, "print fit summary");
            for (int h = 0; h < static_cast<int>(results.size()); ++h) {
                double worst_delta = 0.0;
                double worst_rel_delta = 0.0;
                for (std::size_t i = 1; i < results[h].obj_history.size(); ++i) {
                    const double delta = results[h].obj_history[i] - results[h].obj_history[i - 1];
                    if (delta < worst_delta) {
                        worst_delta = delta;
                        worst_rel_delta = delta / (1.0 + std::abs(results[h].obj_history[i - 1]));
                    }
                }
                const double obj0 = results[h].obj_history.front();
                const double obj_final = results[h].obj_history.back();
                std::cout << "component " << h + 1
                          << ": iterations=" << results[h].iters
                          << ", obj0=" << obj0
                          << ", obj_final=" << obj_final
                          << ", rho_tot=" << results[h].rho_tot
                          << ", monotone=" << (results[h].monotone ? "yes" : "no")
                          << ", worst_delta=" << worst_delta
                          << ", worst_rel_delta=" << worst_rel_delta
                          << '\n';
            }
        }

        timeline.mark("done");
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "frgcca_covmax_timing: " << e.what() << std::endl;
        return 1;
    }
}

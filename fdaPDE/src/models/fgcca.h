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

#ifndef __FGCCA_H__
#define __FGCCA_H__

#include "header_check.h"

namespace fdapde {
enum class Init { Random, SVD };
enum class DesignMode {Empty, FullyConnected};
enum class LambdaSelection {Manual, Automatic};
enum class TauSelection {Manual, Automatic};
enum class Deflation { None, Scores, Loadings };

struct IndependentSampling;
struct TimeDependentSampling;

namespace internals {

template <typename SamplingStrategy, typename ComponentsPenaltyType> class BaseBlock;   // forward decl for operator<<
template <typename SamplingStrategy, typename ComponentsPenaltyType>
std::ostream& operator<<(std::ostream& os, const BaseBlock<SamplingStrategy, ComponentsPenaltyType>& b);

struct NullSolver {
    // common aliases used by your code
    using vector_t = Eigen::VectorXd;
    using matrix_t = Eigen::MatrixXd;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;

    template<class... Args>
    explicit NullSolver(Args&&...) noexcept {}

    // --- stubs matching the methods you might call ---
    template<class... Args> void discretize(Args&&...) {}
    template<class... Args> void analyze_data(Args&&...) {}
    template<class... Args> void update_response_and_weights(Args&&...) {}
    template<class... Args> void fit(Args&&...) {}

    // EDF / sizes
    double edf(int = 0, int = 0) const { return 0.0; }
    int n_obs()  const { return 0; }
    int n_covs() const { return 0; }

    // outputs used in GCV paths etc.
    const vector_t& response() const { static vector_t z; return z; }
    vector_t fn() const { return {}; }
    const vector_t& f() const { static vector_t z; return z; }

    // if you ever query Psi() in the time smoother
    const sparse_matrix_t& Psi() const { static sparse_matrix_t Z; return Z; }
};
struct empty_penalty {
    using solver_t = NullSolver;                 // <— key line
    template<class... Args>
    explicit empty_penalty(Args&&...) noexcept {}
};


// GCV utils
template<class Fun>
inline std::pair<double,double> argmin_over_log_grid(Fun&& f, const double log10_min, const double log10_max, int n_grid) {
    if(n_grid<2) n_grid=2;
    double best_log=log10_min;
    double best_val=std::numeric_limits<double>::infinity();
    const double step=(log10_max-log10_min)/(n_grid-1);
    for(int i=0;i<n_grid;++i) {
        const double lg = log10_min+i*step;
        double lam = std::pow(10.0,lg);
        if (const double val = f(lam); val < best_val) {
            best_val=val;
            best_log=lg;
        }
    }
    return{std::pow(10.0,best_log),best_val};
}

struct GCVConfig {
    // log10 λ range (broad defaults; adjust if you know scale)
    double log10_min = -12.0;
    double log10_max = 4.0;
    int grid = 100;

    // edf() stochastic trace settings (if your solver uses Hutch++ etc.)
    int edf_r = 100;
    int edf_seed = 12345;

    // safety
    double eps_dof = 1e-12;  // avoid divide-by-zero in denominator
};

template <class Smoother>
struct GCVEval {
    Smoother* s;     // must expose: fit(λ), edf(r,seed), response(), fn(), n_obs(), n_covs()
    GCVConfig cfg;

    double operator()(double lambda) {
        s->fit(lambda); // update fit for this λ

        const int n = s->n_obs();
        const int q = s->n_covs();
        const double trS = s->edf(cfg.edf_r, cfg.edf_seed);

        const auto& y  = s->response();
        const auto yhat = s->fn();
        const double rss = (yhat - y).squaredNorm();

        const double dor = std::max( cfg.eps_dof, static_cast<double>(n) - (static_cast<double>(q) + trS) ); // residual dof
        return (static_cast<double>(n) / (dor * dor)) * rss;
    }
};

template <typename SamplingStrategy = IndependentSampling, typename ComponentsPenaltyType = empty_penalty>
class BaseBlock {
public:
    using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using SparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
    using SparseSolver = eigen_sparse_solver_movable_wrap<Eigen::SimplicialLDLT<SparseMatrix>>;
    using ComponentsSolverType = typename std::decay_t<ComponentsPenaltyType>::solver_t;

    template<typename S = SamplingStrategy>
    requires std::same_as<S, IndependentSampling>
    BaseBlock(std::string block_name, Matrix data, const int n_nodes_loadings, const double tau = 0.0) :
        block_name_(std::move(block_name)), data_(std::move(data)), n_nodes_loadings_(n_nodes_loadings), tau_(tau) {
        I_.resize(n_obs(), n_obs());
        I_.setIdentity();
        n_nodes_components_ = n_obs();
    }

    template<typename S = SamplingStrategy>
    requires (std::same_as<SamplingStrategy, TimeDependentSampling> && !std::same_as<ComponentsPenaltyType, empty_penalty>)
    BaseBlock(std::string block_name, Matrix data, const int n_nodes_loadings, const Matrix& times, ComponentsPenaltyType penalty, const double tau = 0.0) :
        block_name_(std::move(block_name)), data_(std::move(data)), n_nodes_loadings_(n_nodes_loadings), tau_(tau) {
        I_.resize(n_obs(), n_obs());
        I_.setIdentity();
        components_solver_.discretize(penalty.get());
        components_solver_.analyze_data(times, Vector::Zero(n_obs()), I_);
        n_nodes_components_ = components_solver_.n_dofs();
    }

    virtual ~BaseBlock() = default;

    // ---- Uniform public API ----

    // Initialization
    void init() {
        ensure_sigma_();
        ensure_lc_();
    }

    // Data
    [[nodiscard]] const std::string& name() const { return block_name_; }
    [[nodiscard]] const Matrix& data() const { return data_; }
    Matrix& data() {
        invalidate_sigma_();
        return data_;
    }

    // Dimensions
    [[nodiscard]] int n_obs() const { return static_cast<int>(data_.rows()); }
    [[nodiscard]] int n_covs() const { return static_cast<int>(data_.cols()); }
    [[nodiscard]] int n_nodes_loadings() const { return n_nodes_loadings_; }

    // Components
    [[nodiscard]] int n_comp() const { return n_comp_; }
    void set_n_comp(const int n_comp) {
        if (n_comp <= 0) throw std::invalid_argument("n_comp must be > 0");
        n_comp_ = n_comp;
        loadings_ready_ = false;   // force resize on next access
        components_ready_ = false;
    }

    // Deflation
    [[nodiscard]] int h() const { return h_; }
    void set_h(const int idx) {
        if (idx < 0 || idx >= n_comp_) throw std::out_of_range("h");
        h_ = idx;
    }
    void next_component() { set_h(h_ + 1); }
    void deflate(const Deflation mode) {
        if (h() == n_comp()) throw std::out_of_range("h");
        switch (mode) {
        case Deflation::Scores:   deflate_scores_(); break;
        case Deflation::Loadings: deflate_loadings_(); break;
        case Deflation::None: default: break;
        }
        // data_ changed -> Σ invalid; cached scores/loadings are now stale
        invalidate_sigma_();
    }

    // Shrinkage parameter
    [[nodiscard]] double tau() const { return tau_; }
    void set_tau(const double tau) {
        if (tau > 0) tau_ = tau;
        else select_tau_auto_();
        invalidate_sigma_();
    }

    // Hyperparameters (default no-op). The model can call this on all blocks.
    virtual void set_lambda_loadings(const double) { }  // default: ignored
    [[nodiscard]] virtual double lambda_loadings() const { return std::numeric_limits<double>::quiet_NaN(); }
    void set_lambda_components(const double lambda) { lambda_components_ = lambda; }
    [[nodiscard]] double lambda_components() const {
        if constexpr (std::same_as<SamplingStrategy, IndependentSampling>) return std::numeric_limits<double>::quiet_NaN();
        if (lambda_components_ > 0) return lambda_components_;
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Noise variance
    void set_noise_sigma_sqr(double noise_sigma_sqr) { noise_sigma_sqr_ = std::max(0.0, noise_sigma_sqr); }
    std::optional<double> noise_sigma_sqr() const { return noise_sigma_sqr_; }

    // Sigma
    [[nodiscard]] const SparseMatrix& Sigma() const { ensure_sigma_(); return Sigma_; }
    SparseSolver& invSigma() { ensure_sigma_(); return invSigma_; }

    // Loadings & Components
    Matrix& loadings() { ensure_lc_(); return loadings_; }
    Matrix loadings_m() { ensure_lc_(); return Psi_D() * loadings_; }
    Matrix& components() { ensure_lc_(); return components_; }
    Matrix components_m() { ensure_lc_(); return Psi_T() * components_; }

    // Inner-Component initialization
    struct InitInfo { bool active; Vector nu; double s1; double s1_edge; double frac; };
    InitInfo svd_init(double epsilon = 0.5) const {
        InitInfo out{false, Vector::Zero(n_obs()), 0.0, 0.0, 0.0};
        const Matrix& X = data();
        if (X.size() == 0) return out;

        const Eigen::BDCSVD<Matrix> svd(X, Eigen::ComputeThinU);
        if (svd.singularValues().size() == 0) return out;

        out.s1 = svd.singularValues()(0);
        const double fro2 = X.squaredNorm();
        out.frac = (fro2 > 0.0) ? (out.s1*out.s1)/fro2 : 0.0;

        // If no σ² set, fall back to your energy test
        if (!noise_sigma_sqr_.has_value()) {
            out.active = (out.frac >= 1e-3);
        } else {
            const double sigma = std::sqrt(std::max(0.0, *noise_sigma_sqr_));
            const double n = static_cast<double>(n_obs());
            const double p = static_cast<double>(n_covs());
            out.s1_edge = sigma * (std::sqrt(n) + std::sqrt(p)) * (1.0 - epsilon);
            out.active = (out.s1 > out.s1_edge);
        }

        if (out.active) {
            out.nu = svd.matrixU().col(0);
        } else {
            out.nu.setZero();
        }
        return out;
    }

    [[nodiscard]] const SparseMatrix Psi_T() const {
        if constexpr (std::same_as<SamplingStrategy, IndependentSampling>) return I_;
        return components_solver_.Psi();
    }

    void compute(const Vector& nu) {
        // Compute loadings
        loadings().col(h()) = l_compute_(nu);
        if (a_().mean() < 0) loadings().col(h()) *= -1.0;

        // TODO: add time contribution and mu estimator

        // Compute scores
        components().col(h()) = c_compute_();

        // Normalize
        scale_to_unit_score_variance_();
    }

    // ---- Clean virtual interface ----
    [[nodiscard]] virtual const SparseMatrix& Psi_D() const = 0;

    // Virtual printer
    virtual void print(std::ostream& os) const {
        const Matrix& X = data();   // no copy
        using Index = Eigen::Index;
        const Index max_rows = 3, max_cols = 5;
        const Index rows = std::min<Index>(max_rows, X.rows());
        const Index cols = std::min<Index>(max_cols, X.cols());
        os << name() << " Block preview (" << X.rows() << " x " << X.cols() << "):\n";
        for (Index i = 0; i < rows; ++i) {
            for (Index j = 0; j < cols; ++j) {
                os << X(i, j);
                if (j + 1 < cols) os << '\t';
            }
            if (X.cols() > cols) os << "\t...";
            os << '\n';
        }
        if (X.rows() > rows) os << "...\n";
        os << "tau = " << tau() << std::endl;
    }

    // Optional: expose a config setter (or pass config from RGCCA Options later)
    void set_components_gcv_config(const GCVConfig& cfg) { components_gcv_cfg_ = cfg; }

protected:

    // pick λ by GCV given current response/weights already set in components_solver_
    std::pair<bool, double> select_components_gcv_config_() {
        GCVEval<ComponentsSolverType> gcv{ &components_solver_, components_gcv_cfg_ };
        auto [lambda_opt, gcv_opt] = argmin_over_log_grid(
            [&](double lam){ return gcv(lam); },
            components_gcv_cfg_.log10_min, components_gcv_cfg_.log10_max, components_gcv_cfg_.grid
        );
        return {lambda_opt < std::pow(10.0, components_gcv_cfg_.log10_max) , lambda_opt};
    }

    // Virtual utilities
    virtual Vector l_compute_(const Vector& nu) = 0;

    // Components solver
    Vector c_compute_() {
        Vector z = data()*Psi_D()*a_();
        if constexpr (std::same_as<SamplingStrategy, TimeDependentSampling>) {
            components_solver_.update_response_and_weights(z, I_);
            // lambda selection if required
            if(lambda_components_ < 0.0) {
                auto [success, lambda_opt] = select_components_gcv_config_();
                if (!success) return Vector::Zero(n_nodes_loadings());
                lambda_components_ = lambda_opt; // the optimal lambda is saved for subsequent calls
            }
            components_solver_.fit(lambda_components_);
            return components_solver_.f();
        }
        return z;
    };

    // Current loading and component getters
    Vector a_() { ensure_lc_(); return loadings_.col(h());}
    Vector eta_() { ensure_lc_(); return components_.col(h());}

    // tau estimate using Schäfer–Strimmer analytic shrinkage from correlation
    void select_tau_auto_() {
        const int n = n_obs(), p = n_covs();
        if (n < 2 || p < 1) throw std::runtime_error("tau_auto: need n>=2 and p>=1");

        // xs <- scale(x, center=TRUE, scale=TRUE)  [sample sd with (n-1)]
        Eigen::RowVectorXd mu  = data_.colwise().mean();
        Matrix xs = data_.rowwise() - mu;                                   // center
        Eigen::RowVectorXd var = (xs.array().square().colwise().sum() / static_cast<double>(n - 1)).matrix();
        Eigen::RowVectorXd sd  = var.array().sqrt().matrix();
        for (int j = 0; j < p; ++j) if (!(sd[j] > 0.0) || !std::isfinite(sd[j])) sd[j] = 1.0;
        xs.array().rowwise() /= sd.array();                                 // scale

        // XtX = crossprod(xs) = t(xs) %*% xs
        Matrix XtX = xs.transpose() * xs;                                    // p x p

        // xs2 = xs^2 ; V = (n/(n-1)^3) * (crossprod(xs2) - (1/n)*(crossprod(xs))^2)
        const double c = static_cast<double>(n) / std::pow(static_cast<double>(n - 1), 3.0);
        const Matrix xs2T_xs2 = (xs.array().square().matrix()).transpose()
                           * (xs.array().square().matrix());                 // p x p
        Matrix V = c * (xs2T_xs2 - (1.0 / static_cast<double>(n)) * XtX.array().square().matrix());
        V.diagonal().setZero();
        const double num = V.sum();

        // corm = cor(x) = (1/(n-1)) * crossprod(xs)
        Matrix Corm = XtX / static_cast<double>(n - 1);
        Matrix D = Corm;
        D.diagonal().array() -= 1.0;
        const double den = D.squaredNorm();

        const double tau_hat = (den > 0.0) ? std::clamp(num / den, 0.0, 1.0) : 0.0;
        set_tau(tau_hat);  // invalidates Sigma_; recomputed lazily
    }

    // Sigma
    void compute_sigma_() {
        SparseMatrix I(n_covs(), n_covs());
        I.setIdentity();
        const Matrix dense = ((1.0 - tau_) / static_cast<double>(n_obs())) * (data_.transpose() * data_);
        Sigma_ = dense.sparseView(1e-12);
        if (tau_ != 0.0) Sigma_ += tau_ * I;
        Sigma_.makeCompressed();
        invSigma_.compute(Sigma_);
        sigma_ready_ = true;
    }
    void ensure_sigma_() const {
        if (!sigma_ready_) const_cast<BaseBlock*>(this)->compute_sigma_();
    }
    void invalidate_sigma_() { sigma_ready_ = false; }

    // Loadings and Components
    void ensure_lc_() {
        if (!loadings_ready_) {
            loadings_.setZero(n_nodes_loadings_, n_comp_);
            loadings_ready_ = true;
        }
        if (!components_ready_) {
            components_.setZero(n_nodes_components_, n_comp_);
            components_ready_ = true;
        }
    }
    void scale_to_unit_score_variance_() {
        // Scale factor s so that Var(eta) = 1 where eta has length n_obs()
        const double v = (Psi_T()*eta_()).squaredNorm() / static_cast<double>(n_obs());
        if (v <= 0.0) return;
        const double norm = std::sqrt(v);
        components().col(h()) /= norm;
        loadings().col(h()) /= norm;
    }

    // Deflation
    void deflate_scores_() {

        // --- Scores deflation (uncorrelated scores next)
        // R = I - η η^T / (η^T η)
        // Then X <- R X

        ensure_lc_(); // make sure components() is sized
        const Vector eta_m = Psi_T()*eta_();   // effective components (length n_obs)
        const int n = n_obs();

        // Assemble projection matrix
        Matrix R = Matrix::Identity(n, n);
        const double norm = eta_m.squaredNorm();
        if (norm <= 0.0) return;
        R.noalias() -= (eta_m * eta_m.transpose()) / norm;

        // Apply left projection in scores space
        data_ = R * data_;
    }
    void deflate_loadings_() {

        // --- Loadings deflation (orthogonal loadings next)
        // R = I - a a^T / (a^T a)
        // Then X <- X R

        ensure_lc_(); // make sure loadings() is sized so loadings_m() is OK
        const Vector a_m = Psi_D()*a_();   // effective loading (length n_covs)
        const int p = n_covs();

        // Assemble projection matrix
        Matrix R = Matrix::Identity(p, p);
        const double norm = a_m.squaredNorm();
        if (norm <= 0.0) return;
        R.noalias() -= (a_m * a_m.transpose()) / norm;

        // Apply right projection in variable space
        data_ = data_ * R;
    }

    // Scores solver
    ComponentsSolverType components_solver_;
    double lambda_components_{-1};

    // State
    std::string block_name_;
    Matrix data_; // n_obs x n_covs
    int n_nodes_loadings_ {0};
    int n_nodes_components_ {0};
    double tau_ {0.0};
    int n_comp_ {1};
    int h_ {0};

    std::optional<double> noise_sigma_sqr_;

    Matrix loadings_, components_;
    bool loadings_ready_ {false}, components_ready_ {false};

    SparseMatrix I_; // n_obs x n_obs identity matrix
    SparseMatrix Sigma_;
    SparseSolver invSigma_;
    bool sigma_ready_ {false};
    GCVConfig components_gcv_cfg_;
};

// single non-member operator<< visible to all derived classes
template <typename SamplingStrategy, typename ComponentsPenaltyType>
inline std::ostream& operator<<(std::ostream& os, const BaseBlock<SamplingStrategy, ComponentsPenaltyType>& b) {
    b.print(os);   // virtual dispatch -> works for Multivariate/Functional too
    return os;
}

// ========== MultivariateBlock ==========
template <typename SamplingStrategy = IndependentSampling, typename ComponentsPenaltyType = empty_penalty>
class MultivariateBlock final : public BaseBlock<SamplingStrategy, ComponentsPenaltyType> {
public:
    using Base = BaseBlock<SamplingStrategy, ComponentsPenaltyType>;
    using Matrix = typename Base::Matrix;
    using Vector = typename Base::Vector;
    using SparseMatrix = typename Base::SparseMatrix;

    using Base::init;
    using Base::Sigma;
    using Base::invSigma;
    using Base::n_obs;
    using Base::n_covs;
    using Base::n_nodes_loadings;
    using Base::data;
    using Base::loadings;
    using Base::components;
    using Base::h;
    using Base::scale_to_unit_score_variance_;

    template<typename S = SamplingStrategy>
    requires std::same_as<S, IndependentSampling>
    MultivariateBlock(std::string block_name, const Matrix& X, const double tau = 0.0) :
        Base(std::move(block_name), X, static_cast<int>(X.cols()), tau) {
        init_multivariate();
    }

    template<typename S = SamplingStrategy>
    requires (std::same_as<SamplingStrategy, TimeDependentSampling> && !std::same_as<ComponentsPenaltyType, empty_penalty>)
    MultivariateBlock(std::string block_name, const Matrix& X, const Matrix& times, ComponentsPenaltyType components_penalty, const double tau = 0.0) :
        Base(std::move(block_name), X, static_cast<int>(X.cols()), times, std::move(components_penalty), tau) {
        init_multivariate();
    }

    void init_multivariate() {
        Psi_.resize(n_nodes_loadings(), n_nodes_loadings());
        Psi_.setIdentity();
        init();
    }

    [[nodiscard]] const SparseMatrix& Psi_D() const override { return Psi_; }

    // print override
    void print(std::ostream& os) const override {
        Base::print(os);
        os << "type: MultivariateBlock, n_nodes_loadings = n_covs = " << n_nodes_loadings();
        os << "\n";
    }

protected:
    Vector l_compute_(const Vector& nu) override {
        assert(nu.size() == n_obs() && "nu must have size n_obs (rows of X)");
        init();

        return invSigma().solve(data().transpose() * nu);
    }

private:
    SparseMatrix Psi_;   // identity
};

// ========== FunctionalBlock ==========
template <class LoadingsPenaltyType, typename SamplingStrategy = IndependentSampling, typename ComponentsPenaltyType = empty_penalty>
class FunctionalBlock final : public BaseBlock<SamplingStrategy, ComponentsPenaltyType> {
public:
    using Base = BaseBlock<SamplingStrategy, ComponentsPenaltyType>;
    using Vector = typename Base::Vector;
    using Matrix = typename Base::Matrix;
    using SparseMatrix = typename Base::SparseMatrix;
    using CovariatesSolverType = typename std::decay_t<LoadingsPenaltyType>::solver_t;

    using Base::init;
    using Base::Sigma;
    using Base::invSigma;
    using Base::n_obs;
    using Base::n_covs;
    using Base::n_nodes_loadings;
    using Base::data;
    using Base::components;
    using Base::h;
    using Base::scale_to_unit_score_variance_;

    template <typename GeoFrame>
    requires std::same_as<SamplingStrategy, IndependentSampling>
    FunctionalBlock(std::string block_name, GeoFrame& gf, LoadingsPenaltyType&& loadings_penalty, const double tau = 0.0) :
        Base(block_name, gf[0].template col<double>(block_name).as_matrix().transpose(), gf.template triangulation<0>().n_nodes(), tau) {
        components_solver_.discretize(loadings_penalty.get());
        components_solver_.analyze_data(gf, Sigma());
        init();
    }

    template <typename GeoFrame>
    requires (std::same_as<SamplingStrategy, TimeDependentSampling> && !std::same_as<ComponentsPenaltyType, empty_penalty>)
    FunctionalBlock(std::string block_name, GeoFrame& gf,
                LoadingsPenaltyType&& loadings_penalty,
                const Matrix& times, ComponentsPenaltyType components_penalty, const double tau = 0.0) :
        Base(block_name, gf[0].template col<double>(block_name).as_matrix().transpose(), gf.template triangulation<0>().n_nodes(), times, std::move(components_penalty), tau) {
        components_solver_.discretize(loadings_penalty.get());
        components_solver_.analyze_data(gf, Sigma());
        init();
    }

    [[nodiscard]] const SparseMatrix& Psi_D() const override { return components_solver_.Psi(); }

    // The model must set lambda before calling l_compute
    void set_lambda_loadings(const double lambda) override { lambda_loadings_ = lambda; }
    [[nodiscard]] double lambda_loadings() const override {
        if (lambda_loadings_ > 0) return lambda_loadings_;
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Optional: expose a config setter (or pass config from RGCCA Options later)
    void set_loadings_gcv_config(const GCVConfig& cfg) { loadings_gcv_cfg_ = cfg; }

    // print override
    void print(std::ostream& os) const override {
        Base::print(os);
        os << "type: FunctionalBlock, n_nodes_loadings = " << n_nodes_loadings();
        os << ", lambda = " << lambda_loadings_;
        os << "\n";
    }
protected:
    // pick λ by GCV given current response/weights already set in components_solver_
    std::pair<bool, double> select_loadings_gcv_() {
        GCVEval<CovariatesSolverType> gcv{ &components_solver_, loadings_gcv_cfg_ };
        auto [lambda_opt, gcv_opt] = argmin_over_log_grid(
            [&](double lam){ return gcv(lam); },
            loadings_gcv_cfg_.log10_min, loadings_gcv_cfg_.log10_max, loadings_gcv_cfg_.grid
        );
        return {lambda_opt < std::pow(10.0, loadings_gcv_cfg_.log10_max) , lambda_opt};
    }

    Vector l_compute_(const Vector& nu) override {
        assert(nu.size() == n_obs() && "nu must have size n_obs (rows of X)");
        init();

        const Vector z = invSigma().solve(data().transpose() * nu);
        components_solver_.update_response_and_weights(z, Sigma());

        // lambda selection if required
        if(lambda_loadings_ < 0.0) {
            auto [success, lambda_opt] = select_loadings_gcv_();
            if (!success) return Vector::Zero(n_nodes_loadings());
            lambda_loadings_ = lambda_opt; // the optimal lambda is saved for subsequent calls
        }

        components_solver_.fit(lambda_loadings_);
        return components_solver_.f();
    }
private:
    CovariatesSolverType components_solver_;
    double lambda_loadings_ = -1.0;   // < 0 means "use GCV"
    GCVConfig loadings_gcv_cfg_;
};



template <typename SamplingStrategy = IndependentSampling, typename ComponentsPenaltyType = empty_penalty>
requires std::same_as<SamplingStrategy, IndependentSampling>
inline std::unique_ptr<internals::BaseBlock<SamplingStrategy, ComponentsPenaltyType>>
make_multivariate_block(std::string name, Eigen::Matrix<double, Dynamic, Dynamic>& data, double tau = 0.0) {
    return std::make_unique<internals::MultivariateBlock<SamplingStrategy, ComponentsPenaltyType>>(std::move(name), std::move(data), tau);
}

template <typename SamplingStrategy = IndependentSampling, typename ComponentsPenaltyType = empty_penalty>
requires (std::same_as<SamplingStrategy, TimeDependentSampling> && !std::same_as<ComponentsPenaltyType, empty_penalty>)
inline std::unique_ptr<internals::BaseBlock<SamplingStrategy, ComponentsPenaltyType>>
make_multivariate_block(std::string name, Eigen::Matrix<double, Dynamic, Dynamic>& data, Eigen::Matrix<double, Dynamic, Dynamic>& times,
                        ComponentsPenaltyType components_penalty, double tau = 0.0) {
    return std::make_unique<internals::MultivariateBlock<SamplingStrategy, ComponentsPenaltyType>>(
        std::move(name), std::move(data), times, std::move(components_penalty), tau);
}


template <typename GeoFrame, typename LoadingsPenaltyType, typename SamplingStrategy = IndependentSampling, typename ComponentsPenaltyType = empty_penalty>
requires std::same_as<SamplingStrategy, IndependentSampling>
std::unique_ptr<internals::BaseBlock<SamplingStrategy, ComponentsPenaltyType>>
make_functional_block(std::string name, GeoFrame& gf, LoadingsPenaltyType&& pen, double tau = 0.0) {
    return std::make_unique<internals::FunctionalBlock<LoadingsPenaltyType, SamplingStrategy, ComponentsPenaltyType>>(
      std::move(name), gf, std::forward<LoadingsPenaltyType>(pen), tau);
}

template <typename GeoFrame, typename LoadingsPenaltyType, typename SamplingStrategy = IndependentSampling, typename ComponentsPenaltyType = empty_penalty>
requires (std::same_as<SamplingStrategy, TimeDependentSampling> && !std::same_as<ComponentsPenaltyType, empty_penalty>)
std::unique_ptr<internals::BaseBlock<SamplingStrategy, ComponentsPenaltyType>>
make_functional_block(std::string name, GeoFrame& gf, LoadingsPenaltyType&& loadings_penalty,
                      Eigen::Matrix<double, Dynamic, Dynamic>& times, ComponentsPenaltyType components_penalty, double tau = 0.0) {
    return std::make_unique<internals::FunctionalBlock<LoadingsPenaltyType, SamplingStrategy, ComponentsPenaltyType>>(
        std::move(name), gf, std::forward<LoadingsPenaltyType>(loadings_penalty),
        times, std::move(components_penalty), tau);
}



}   // namespace internals

// ===== Scheme (g, w, phi) =====
struct Scheme {
    std::function<double(double)> g;   // g(t)
    std::function<double(double)> w;   // w(t)
    double phi = 1.0;
    const char* name = "custom";

    static Scheme Horst() {
        return {[](double t) { return t; }, [](double) { return 1.0; }, 1.0, "Horst"};
    }
    static Scheme Centroid() {
        return {
            [](double t) { return std::abs(t); }, [](double t) { return t >= 0 ? 1.0 : -1.0; }, 1.0, "Centroid"};
    }
    static Scheme Factorial() {
        return {[](double t) { return t * t; }, [](double t) { return t; }, 2.0, "Factorial"};
    }
};

// ===== Results =====
struct Result {
    using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

    int h = 0;
    int J = 0;
    std::vector<double> obj_history;
    bool monotone = true;
    int iters = 0;
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> C;
    Matrix covariance_matrix;
    std::vector<double> tau_values;
    std::vector<double> lambda_components_values;
    std::vector<double> lambda_loadings_values;
    std::vector<double> active_blocks;
    std::vector<double> s1_blocks;
    std::vector<double> s1_edge_blocks;

    explicit Result(const int n_blocks) : J(n_blocks), C(J, J), covariance_matrix(J,J),
    tau_values(J), lambda_components_values(J), lambda_loadings_values(J), active_blocks(J), s1_blocks(J), s1_edge_blocks(J)  {}

};

// forward declaration of pretty printers
std::ostream& operator<<(std::ostream& os, const Result& r);
std::ostream& operator<<(std::ostream& os, const std::vector<Result>& results);

// ===== RGCCA =====
template <typename SamplingStrategy = IndependentSampling, typename ComponentsPenaltyType = internals::empty_penalty>
class RGCCA {
public:
    using Block = internals::BaseBlock<SamplingStrategy, ComponentsPenaltyType>;
    using BlockPtr = std::unique_ptr<Block>;
    using Matrix = typename Block::Matrix;
    using Vector = typename Block::Vector;

    struct Options {
        int max_iter;
        double tol;
        unsigned seed;
        bool verbose;
        bool cache_covariances;
        Init init;
        LambdaSelection lambda_selection;
        TauSelection tau_selection;
        Deflation deflation_mode;

        explicit Options(
          const int max_iter_ = 1000, const double tol_ = 1e-8, const unsigned seed_ = 0,
          const Init init_ = Init::SVD, const TauSelection tau_selection_ = TauSelection::Automatic,
          const LambdaSelection lambda_selection_ = LambdaSelection::Automatic,
          const Deflation deflation_mode_ = Deflation::Scores,
          const bool verbose_ = false, const bool cache_ = true) :
            max_iter(max_iter_),
            tol(tol_),
            seed(seed_),
            init(init_),
            tau_selection(tau_selection_),
            lambda_selection(lambda_selection_),
            deflation_mode(deflation_mode_),
            verbose(verbose_),
            cache_covariances(cache_) { }
    };

    explicit RGCCA(const int n_obs, const Scheme& scheme = Scheme::Horst(), const Options& opt = Options(), const int n_comp = 1) :
        n_obs_(n_obs), scheme_(std::move(scheme)), opt_(opt), n_comp_(n_comp) {}

    // ===== Blocks =====
    int add_block(BlockPtr b) {
        if (!b) throw std::invalid_argument("RGCCA/add_block: null block");
        if (b->n_obs() != n_obs()) throw std::invalid_argument("RGCCA/add_block: n_obs mismatch");
        b->set_n_comp(n_comp());
        blocks_.emplace_back(std::move(b));
        initialized_ = false;   // topology/caches need a fresh init later
        return ++J_;
    }

    template<typename S = SamplingStrategy>
    requires std::same_as<S, IndependentSampling>
    int add_multivariate_block(std::string name, Matrix& X, const double tau = 0.0) {
        return add_block(internals::make_multivariate_block(std::move(name), X, tau));
    }
    template <typename GeoFrame, typename LoadingsPenaltyType>
    requires std::same_as<SamplingStrategy, IndependentSampling>
    int add_functional_block(std::string name, GeoFrame& gf, LoadingsPenaltyType&& loadings_penalty, const double tau = 0.0) {
        return add_block(internals::make_functional_block(std::move(name), gf, std::forward<LoadingsPenaltyType>(loadings_penalty), tau));
    }
    template<typename S = SamplingStrategy>
    requires std::same_as<S, TimeDependentSampling>
    int add_multivariate_block(std::string name, Matrix& X, Matrix& times, ComponentsPenaltyType components_penalty, const double tau = 0.0) {
        return add_block(internals::make_multivariate_block<S, ComponentsPenaltyType>(
            std::move(name), X, times, std::move(components_penalty), tau));
    }
    template <typename GeoFrame, typename LoadingsPenaltyType>
    requires std::same_as<SamplingStrategy, TimeDependentSampling>
    int add_functional_block(std::string name, GeoFrame& gf, LoadingsPenaltyType&& loadings_penalty,
                             Matrix& times, ComponentsPenaltyType components_penalty,
                             const double tau = 0.0) {
        return add_block(internals::make_functional_block<GeoFrame, LoadingsPenaltyType, SamplingStrategy, ComponentsPenaltyType>(
            std::move(name), gf, std::forward<LoadingsPenaltyType>(loadings_penalty),
            times, std::move(components_penalty), tau));
    }


    [[nodiscard]] int n_blocks() { return J_; }

    // Connect blocks
    void connect(int j, int k, bool on = true) {
        if (!initialized_) init(DesignMode::Empty);
        check_index_(j);
        check_index_(k);
        if (j == k) return;
        C_(k, j) = C_(j, k) = on;
        user_defined_design_ = true;
    }

    // ===== One-shot init (does all resizes) =====
    void init(const DesignMode mode) {
        const int J = n_blocks();
        if (J < 2) throw std::runtime_error("RGCCA/init: need ≥ 2 blocks");
        // resize design
        C_.resize(J, J);
        C_.setConstant(false);
        if (mode == DesignMode::FullyConnected) {
            for (int j = 0; j < J; ++j)
                for (int k = 0; k < J; ++k)
                    if (k != j) C_(j, k) = true;   // diag remains false
        }
        if (noise_sigma_sqr_.has_value()) set_noise_sigma_sqr_all_();
        initialized_ = true;
        user_defined_design_ = (mode == DesignMode::Empty);   // means user will set edges
    }
    void init_comp() {
        if (opt_.tau_selection == TauSelection::Automatic) { set_tau_auto_all_(); }
        if (opt_.lambda_selection == LambdaSelection::Automatic) { set_lambda_auto_all_(); }
        clear_covariance_cache_();
    }

    // Noise
    void set_noise_sigma_sqr(double noise_sigma_sqr) { noise_sigma_sqr_ = std::max(0.0, noise_sigma_sqr); }
    std::optional<double> noise_sigma_sqr() const { return noise_sigma_sqr_; }

    // Parameters setters
    void set_lambda_loadings_all(const double lambda) const {
        for (auto& b : blocks_) b->set_lambda_loadings(lambda);
    }
    void set_lambda_components_all(const double lambda) const {
        for (auto& b : blocks_) b->set_lambda_components(lambda);
    }

    // Components
    void set_n_comp(const int n_comp) {
        if (n_comp <= 0) throw std::invalid_argument("n_comp must be > 0");
        n_comp_ = n_comp;
        for (auto& b : blocks_) b->set_n_comp(n_comp);
        if (h_ >= n_comp) set_h(n_comp - 1);
    }
    [[nodiscard]] int h() const { return h_; }
    void set_h(const int h) {
        if (h < 0 || h >= n_comp_) throw std::out_of_range("component index");
        h_ = h;
        for (auto& b : blocks_) b->set_h(h_);
    }
    void next_comp() {
        if (!is_last_comp()) set_h(h() + 1);
    }
    bool is_last_comp() {
        if (!is_last_comp_) {
            if (h() + 1 == n_comp()) is_last_comp_ = true;
            return false;
        }
        return true;
    }

    // Deflation
    void deflate_all() const {
        for (auto& b : blocks_) b->deflate(opt_.deflation_mode);
    }

    // Fit
    std::vector<Result> fit() {
        if (!initialized_) {
            // user didn't call init ⇒ assume fully connected (off-diagonal true)
            init(DesignMode::FullyConnected);
        }

        // room for results
        std::vector<Result> results;
        results.reserve(n_comp());

        // components loop
        for (set_h(0); !is_last_comp(); next_comp()) {
            init_comp();
            results.push_back(fit_component());
            deflate_all();
        }

        return results;
    }
    Result fit_component() {

        const int J = n_blocks();
        if (J < 2) throw std::runtime_error("RGCCA: need ≥ 2 blocks");

        // room for results
        Result res(n_blocks());

        // make a component local copy of the connection matrix C
        res.C = C_;

        // random or SVD init for ν_l
        for (int j = 0; j < J; ++j) {
            auto& b = blocks_[j];
            b->set_h(h_);
            if (opt_.init == Init::SVD) {
                const auto info = b->svd_init();
                if (!info.active) {
                    res.active_blocks[j] = false;
                    res.s1_blocks[j] = info.s1;
                    res.s1_edge_blocks[j] = info.s1_edge;
                    // nothing to be done loadings and components are already initialized at 0
                } else {
                    res.active_blocks[j] = true;
                    res.s1_blocks[j] = info.s1;
                    res.s1_edge_blocks[j] = info.s1_edge;
                    b->compute(info.nu);  // block handles normalization
                }
            } else { // Random
                std::mt19937_64 rng(opt_.seed);
                std::uniform_real_distribution<double> U(-1.0, 1.0);
                const Vector nu = Vector::NullaryExpr(n_obs_, [&]{ return U(rng); });
                b->compute(nu);  // block handles normalization
            }
        }

        // disconnect blocks that are not active
        for (int j = 0; j < J; ++j) if (!res.active_blocks[j]) {
            for (int k = 0; k < J; ++k) { res.C(j,k) = false; res.C(k,j) = false; }
        }

        // room for objective function evaluations
        res.obj_history.reserve(opt_.max_iter + 1);
        res.obj_history.push_back(objective_());

        if ( !no_connections_(res.C) ) {
            // require lambda selection also at the first iteration
            if (opt_.lambda_selection == LambdaSelection::Automatic) { set_lambda_auto_all_(); }
            for (int s = 0; s < opt_.max_iter; ++s) {
                for (int l = 0; l < J; ++l) {
                    Vector nu_l = Vector::Zero(n_obs_);
                    const Vector eta_l = eta_(*blocks_[l]);
                    for (int k = 0; k < J; ++k) {
                        if (k == l || !res.C(l,k)) continue;   // <— exclude self
                        const Vector eta_k = eta_(*blocks_[k]);
                        const double cov_lk = cov_value_(l, k, eta_l, eta_k);   // uses/saves cache, marks clean
                        const double w_lk = scheme_.w(cov_lk);
                        nu_l.noalias() += w_lk * eta_k;   // no aliasing with RHS
                    }
                    blocks_[l]->compute(nu_l);   // block handles normalization
                    mark_cov_rowcol_dirty_(l);     // η_l changed → invalidate its row/col
                }

                const double f = objective_();
                const double prev = res.obj_history.back();
                res.obj_history.push_back(f);
                res.iters = s + 1;

                if (res.obj_history.back() + 1e-15 < res.obj_history[res.obj_history.size() - 2]) res.monotone = false;
                const double rel  = std::abs(f - prev) / (std::abs(prev) + 1e-16);
                if (rel < opt_.tol) break;
            }
        }
        compute_covariance_matrix_(res.covariance_matrix);
        get_tau(res.tau_values);
        get_lambdas(res.lambda_components_values, res.lambda_loadings_values);

        return res;
    }

    // ===== Accessors =====
    [[nodiscard]] int n_obs() const { return n_obs_; }
    [[nodiscard]] int n_comp() const { return n_comp_; }
    [[nodiscard]] int n_blocks() const { return static_cast<int>(blocks_.size()); }
    [[nodiscard]] const Scheme& scheme() const { return scheme_; }
    [[nodiscard]] const Options& options() const { return opt_; }
    [[nodiscard]] const std::vector<BlockPtr>& blocks() const { return blocks_; }
    [[nodiscard]] const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& C() const { return C_; }
    [[nodiscard]] bool initialized() const { return initialized_; }
    [[nodiscard]] bool user_defined_design() const { return user_defined_design_; }

private:
    void clear_covariance_cache_() {
        const int J = n_blocks();
        // resize covariance cache + dirty mask
        Cov_.setZero(J, J);
        dirty_.setOnes(J, J);
        for (int i = 0; i < J; ++i) {
            Cov_(i, i) = 1.0;
            dirty_(i, i) = 0;
        }
    }

    void get_tau(std::vector<double> & tau_values) const {
        const int J = n_blocks();
        for (std::size_t j = 0; j < tau_values.size(); ++j) {
            tau_values[j] = blocks_[j]->tau();
        }
    }

    void get_lambdas(std::vector<double> & lambda_components_values, std::vector<double> & lambda_loadings_values) const {
        const int J = n_blocks();
        for (std::size_t j = 0; j < J; ++j) {
            lambda_components_values[j] = blocks_[j]->lambda_components();
            lambda_loadings_values[j] = blocks_[j]->lambda_loadings();
        }
    }

    void set_tau_auto_all_() const { for (auto& b : blocks_) b->set_tau(-1); }
    void set_lambda_auto_all_() const { set_lambda_components_all(-1); set_lambda_loadings_all(-1); }

    void set_noise_sigma_sqr_all_() const { for (auto& b : blocks_) b->set_noise_sigma_sqr(*noise_sigma_sqr_); }

    // ===== Helpers =====
    Vector eta_(Block& b) const {
        // TODO: this will not work in general. Define a Psi matrix to a common set of times across blocks and apply it to components() instead of components_m()
        return b.components_m().col(h()); // η_j = X_j a_mj
    }
    [[nodiscard]] double cov_(const Vector& u, const Vector& v) const {
        return (1.0 / static_cast<double>(n_obs_)) * u.dot(v);
    }

    // objective f = Σ_{j,k} C_jk * g( cov(η_j, η_k) )
    double objective_() {
        const int J = n_blocks();
        double f = 0.0;
        for (int j = 0; j < J; ++j) {
            const Vector eta_j = eta_(*blocks_[j]);
            for (int k = j+1; k < J; ++k){
                if (C_(j, k)) {
                    const double cjk = cov_value_(j, k, eta_j, eta_(*blocks_[k]));
                    f += 2 * scheme_.g(cjk);
                }
            }
        }
        return f;
    }

    // Covariance matrix
    void compute_covariance_matrix_(Matrix & Cov) const {
        const int J = n_blocks();
        for (int j = 0; j < J; ++j) {
            const Vector eta_j = eta_(*blocks_[j]);
            for (int k = 0; k < J; ++k) {
                const Vector eta_k = eta_(*blocks_[k]);
                Cov(j, k) = cov_(eta_j, eta_k);
            }
        }
    }

    // --- covariance cache management ---
    void mark_cov_rowcol_dirty_(int l) {
        if (!opt_.cache_covariances) return;
        for (int k = 0; k < n_blocks(); ++k) {
            dirty_(l, k) = 1;
            dirty_(k, l) = 1;
        }
        dirty_(l, l) = 0;
        Cov_(l, l) = 1.0;
    }

    // compute or reuse cov(l,k); when computed, store & mark clean (both (l,k) and (k,l))
    double cov_value_(int l, int k, const Vector& eta_l, const Vector& eta_k) {
        if (!dirty_(l, k)) return Cov_(l, k);
        const double c = cov_(eta_l, eta_k);
        Cov_(l, k) = Cov_(k, l) = c;
        dirty_(l, k) = dirty_(k, l) = 0;
        return c;
    }
    void ensure_cov_shapes_() {
        const int J = n_blocks();
        if (Cov_.rows() != J) {
            Cov_.setZero(J, J);
            for (int i = 0; i < J; ++i) Cov_(i, i) = 1.0;
        }
        if (dirty_.rows() != J || dirty_.cols() != J) {
            dirty_.setOnes(J, J);
            for (int i = 0; i < J; ++i) dirty_(i, i) = 0;
        }
    }

    // indexes
    void check_index_(int j) const {
        if (j < 0 || j >= static_cast<int>(blocks_.size())) throw std::out_of_range("block index");
    }

    bool no_connections_(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> C) {
        bool flag = true;
        const int J = n_blocks();
        for (int j = 0; j < J; ++j) {
            for (int k = 0; k < J; ++k) {
                if (C(j, k)) return false;
            }
        }
        return true;
    }

private:
    int J_ {0};
    int n_obs_;   // global #observations
    int h_ {0};   // current component index
    Scheme scheme_;
    Options opt_;
    int n_comp_{0};
    bool is_last_comp_{false};

    std::optional<double> noise_sigma_sqr_;

    std::vector<BlockPtr> blocks_;

    // topology & caches (sized in init())
    bool initialized_ {false};
    bool user_defined_design_ {false};
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> C_;

    Matrix Cov_ {0, 0};          // cached covariances between η's (Cov_(j,j)=1)
    Eigen::ArrayXXi dirty_ {0, 0};   // 1=dirty, 0=clean
};


// Pretty printer for a single Result
inline std::ostream& operator<<(std::ostream& os, const Result& r) {
    os << "shrinkage parameters used : " << std::endl;
    for (size_t i = 0; i < r.tau_values.size(); ++i) {
        os << "- Block " << i+1  << ": tau = "<< r.tau_values[i] << "\n";
    }
    os << std::endl;
    os << "active blocks :\n";
    for (size_t i = 0; i < r.active_blocks.size(); ++i) {
        os << "- Block " << i+1  << ": " << (r.active_blocks[i] ? "active    " : "non-active" )
           << " ("<< r.s1_blocks[i]<< (r.active_blocks[i] ? " > " : " < ") << r.s1_edge_blocks[i] << ")" << "\n";
    }
    os << std::endl;
    os << "(updated) connections matrix :\n";
    os << r.C << std::endl;
    os << std::endl;
    os << "regularization parameters used : " << std::endl;
    for (size_t i = 0; i < r.tau_values.size(); ++i) {
        os << "- Block " << i+1  << ": lambda_c = "<< r.lambda_components_values[i]
           << ", lambda_l = "<< r.lambda_loadings_values[i] << "\n";
    }
    os << std::endl;
    os << "n_iters   : " << r.iters << "\n";
    os << "monotone  : " << (r.monotone ? "yes" : "no") << "\n";
    os << std::endl;
    os << "objective :\n";
    double prev = 0.0;
    for (size_t i = 0; i < r.obj_history.size(); ++i) {
        const double val = r.obj_history[i];
        os << "- iter " << std::setw(3) << (i + 1)
           << " | fit = " << std::setw(12) << std::setprecision(8) << val
           << " | diff = " << std::setw(12) << (val - prev) << "\n";
        prev = val;
    }
    os << std::endl;
    os << "covariance matrix :\n";
    os << r.covariance_matrix << std::endl;

    return os;
}

// Pretty printer for a vector of Result (components)
inline std::ostream& operator<<(std::ostream& os, const std::vector<Result>& results) {
    for (size_t h = 0; h < results.size(); ++h) {
        os << "========================================\n";
        os << "Component " << (h + 1) << "\n";
        os << "----------------------------------------\n";
        os << results[h]; // delegate to the single-result printer
        os << "\n";
    }
    return os;
}

}

#endif   // __FGCCA_H__
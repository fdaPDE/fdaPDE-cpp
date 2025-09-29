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
enum class TauSelection {Manual, Automatic};
enum class Deflation { None, Scores, Loadings };

namespace internals {

class BaseBlock;   // forward decl for operator<<
std::ostream& operator<<(std::ostream& os, const BaseBlock& b);

class BaseBlock {
public:
    using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using SparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
    using SparseSolver = eigen_sparse_solver_movable_wrap<Eigen::SimplicialLDLT<SparseMatrix>>;

    BaseBlock(std::string block_name, Matrix data, const int n_nodes, const double tau = 0.0, const int n_comp = 1) :
        block_name_(std::move(block_name)), data_(std::move(data)), n_nodes_(n_nodes), tau_(tau), n_comp_(n_comp) { }

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
    [[nodiscard]] int n_nodes() const { return n_nodes_; }

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
    void set_h(int idx) {
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

    // Sigma
    [[nodiscard]] const SparseMatrix& Sigma() const { ensure_sigma_(); return Sigma_; }
    SparseSolver& invSigma() { ensure_sigma_(); return invSigma_; }

    // Loadings
    Matrix& loadings() { ensure_lc_(); return loadings_; }
    Matrix loadings_m() { ensure_lc_(); return Psi() * loadings_; }
    Matrix& components() { ensure_lc_(); return components_; }

    // ---- Clean virtual interface ----
    [[nodiscard]] virtual const SparseMatrix& Psi() const = 0;
    virtual void l_compute(const Vector& nu) = 0;

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

    // Hyperparameters (default no-op). The model can call this on all blocks.
    virtual void set_lambda(double) { }   // default: ignored

protected:

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
        const double c = double(n) / std::pow(static_cast<double>(n - 1), 3.0);
        Matrix xs2T_xs2 = (xs.array().square().matrix()).transpose()
                           * (xs.array().square().matrix());                 // p x p
        Matrix V = c * (xs2T_xs2 - (1.0 / double(n)) * XtX.array().square().matrix());
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
            loadings_.setZero(n_nodes_, n_comp_);
            loadings_ready_ = true;
        }
        if (!components_ready_) {
            components_.setZero(n_obs(), n_comp_);
            components_ready_ = true;
        }
    }
    [[nodiscard]] double scale_to_unit_score_variance_(Vector& eta) const {
        // Scale factor s so that Var(eta) = 1 where eta has length n_obs()
        const double v = eta.squaredNorm() / static_cast<double>(n_obs());
        if (v <= 0.0) return 1.0;
        const double s = 1.0 / std::sqrt(v);
        eta *= s;
        return s;
    }

    // Deflation
    void deflate_scores_() {

        // --- Scores deflation (uncorrelated scores next)
        // R = I - η η^T / (η^T η)
        // Then X <- R X

        ensure_lc_(); // make sure components() is sized
        const Vector eta = components().col(h_);   // effective components  (length n_obs)
        const int n = n_obs();

        // Assemble projection matrix
        Matrix R = Matrix::Identity(n, n);
        const double denom = eta.squaredNorm();
        if (denom <= 0.0) return;
        R.noalias() -= (eta * eta.transpose()) / denom;

        // Apply left projection in scores space
        data_ = R * data_;
    }
    void deflate_loadings_() {

        // --- Loadings deflation (orthogonal loadings next)
        // R = I - a a^T / (a^T a)
        // Then X <- X R

        ensure_lc_(); // make sure loadings() is sized so loadings_m() is OK
        const Vector a_m = loadings_m().col(h_);   // effective loading a_m = Ψ f  (length n_covs)
        const int p = n_covs();

        // Assemble projection matrix
        Matrix R = Matrix::Identity(p, p);
        const double denom = a_m.squaredNorm();
        if (denom <= 0.0) return;
        R.noalias() -= (a_m * a_m.transpose()) / denom;

        // Apply right projection in variable space
        data_ = data_ * R;
    }

    // State
    std::string block_name_;
    Matrix data_; // n_obs x n_covs
    int n_nodes_ {0};
    double tau_ {0.0};
    int n_comp_ {1};
    int h_ {0};

    Matrix loadings_, components_;
    bool loadings_ready_ {false}, components_ready_ {false};

    SparseMatrix Sigma_;
    SparseSolver invSigma_;
    bool sigma_ready_ {false};
};

// single non-member operator<< visible to all derived classes
inline std::ostream& operator<<(std::ostream& os, const BaseBlock& b) {
    b.print(os);   // virtual dispatch -> works for Multivariate/Functional too
    return os;
}

// ========== MultivariateBlock ==========
class MultivariateBlock final : public BaseBlock {
public:
    using Base = BaseBlock;
    using Matrix = Base::Matrix;
    using Vector = Base::Vector;
    using SparseMatrix = Base::SparseMatrix;

    MultivariateBlock(std::string block_name, const Matrix& X, const double tau = 0.0, const int n_comp = 1) :
        Base(std::move(block_name), X, static_cast<int>(X.cols()), tau, n_comp) {
        Psi_.resize(n_nodes(), n_nodes());
        Psi_.setIdentity();
        init();
    }

    [[nodiscard]] const SparseMatrix& Psi() const override { return Psi_; }

    void l_compute(const Vector& nu) override {
        assert(nu.size() == n_obs() && "nu must have size n_obs (rows of X)");
        init();

        // 1) raw loading (no Σ-normalization)
        Vector a_raw = invSigma().solve(data().transpose() * nu);

        // (optional) orientation for interpretability
        if (a_raw.mean() < 0) a_raw = -a_raw;

        // 2) raw scores (η̃) and unit-variance rescale
        Vector eta = data() * a_raw;
        const double s = scale_to_unit_score_variance_(eta); // scales eta to Var(η)=1

        // 3) store scaled loadings and scores
        loadings().col(h())   = s * a_raw;      // keep η = X (Psi a) consistent
        components().col(h()) = eta;            // already scaled to unit variance
    }

    // print override
    void print(std::ostream& os) const override {
        BaseBlock::print(os);
        os << "type: MultivariateBlock, n_nodes = n_covs = " << n_nodes();
        os << "\n";
    }
private:
    SparseMatrix Psi_;   // identity
};

// ========== FunctionalBlock ==========
template <class PenaltyType>
class FunctionalBlock final : public BaseBlock {
public:
    using Base = BaseBlock;
    using Vector = Base::Vector;
    using Matrix = Base::Matrix;
    using SparseMatrix = Base::SparseMatrix;
    using lSolverType = typename std::decay_t<PenaltyType>::solver_t;

    template <typename GeoFrame>
    FunctionalBlock(std::string block_name, GeoFrame& gf, PenaltyType&& penalty, const double tau = 0.0, const int n_comp = 1) :
        Base(block_name, gf[0].template col<double>(block_name).as_matrix().transpose(), gf.template triangulation<0>().n_nodes(), tau, n_comp) {
        solver_.discretize(penalty.get());
        solver_.analyze_data(gf, Sigma());
        init();
    }

    [[nodiscard]] const SparseMatrix& Psi() const override { return solver_.Psi(); }

    // The model must set lambda before calling l_compute
    void set_lambda(double lambda) override { lambda_ = lambda; }

    void l_compute(const Vector& nu) override {
        assert(nu.size() == n_obs() && "nu must have size n_obs (rows of X)");
        init();

        const Vector z = invSigma().solve(data().transpose() * nu);
        solver_.update_response_and_weights(z, Sigma());
        solver_.fit(lambda_);

        Vector f = solver_.f();
        if (f.mean() < 0) f = -f;

        // a_m_raw = Psi f_raw; scores η̃ = X a_m_raw
        Vector a_m = Psi() * f;
        Vector eta = data() * a_m;

        // unit-variance rescale for scores; apply the same factor to both f and a_m
        const double s = scale_to_unit_score_variance_(eta);

        // store
        loadings().col(h())   = s * f;      // pre-Psi loadings kept coherent
        components().col(h()) = eta;            // unit-variance scores
    }

    // print override
    void print(std::ostream& os) const override {
        BaseBlock::print(os);
        os << "type: FunctionalBlock, n_nodes = " << n_nodes();
        os << ", lambda = " << lambda_;
        os << "\n";
    }
private:
    lSolverType solver_;
    double lambda_ = 1e-12;   // owned by the block (set by the model)
};

inline std::unique_ptr<internals::BaseBlock> make_multivariate_block(
  std::string name, Eigen::Matrix<double, Dynamic, Dynamic>& data, double tau = 0.0, int n_comp = 1) {
    return std::make_unique<internals::MultivariateBlock>(std::move(name), data, tau, n_comp);
}

template <typename GeoFrame, typename PenaltyType>
std::unique_ptr<internals::BaseBlock> make_functional_block(
  std::string name, GeoFrame& gf, PenaltyType&& pen, double tau = 0.0, int n_comp = 1) {
    return std::make_unique<internals::FunctionalBlock<PenaltyType>>(
      std::move(name), gf, std::forward<PenaltyType>(pen), tau, n_comp);
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


struct Result {
    using Block = internals::BaseBlock;
    using Matrix = Block::Matrix;

    int h = 0;
    int J = 0;
    std::vector<double> obj_history;
    bool monotone = true;
    int iters = 0;
    Matrix covariance_matrix;
    std::vector<double> tau_values;

    explicit Result(const int n_blocks) : J(n_blocks), covariance_matrix(J,J), tau_values(J) {}

};

// forward declaration of pretty printers
std::ostream& operator<<(std::ostream& os, const Result& r);
std::ostream& operator<<(std::ostream& os, const std::vector<Result>& results);

class RGCCA {
public:
    using Block = internals::BaseBlock;
    using BlockPtr = std::unique_ptr<Block>;
    using Matrix = Block::Matrix;
    using Vector = Block::Vector;

    struct Options {
        int max_iter;
        double tol;
        unsigned seed;
        bool verbose;
        bool cache_covariances;
        Init init;
        TauSelection tau_selection;
        Deflation deflation_mode;

        explicit Options(
          const int max_iter_ = 1000, const double tol_ = 1e-8, const unsigned seed_ = 0,
          const Init init_ = Init::SVD, const TauSelection tau_selection_ = TauSelection::Automatic,
          const Deflation deflation_mode_ = Deflation::Scores,
          const bool verbose_ = false, const bool cache_ = true) :
            max_iter(max_iter_),
            tol(tol_),
            seed(seed_),
            init(init_),
            tau_selection(tau_selection_),
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
    int add_multivariate_block(std::string name, Matrix& X, const double tau = 0.0) {
        return add_block(internals::make_multivariate_block(std::move(name), X, tau));
    }
    template <class Tri, class Pen>
    int add_functional_block(std::string name, GeoFrame<Tri>& gf, Pen&& pen, const double tau = 0.0) {
        return add_block(internals::make_functional_block(std::move(name), gf, std::forward<Pen>(pen), tau));
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
        initialized_ = true;
        user_defined_design_ = (mode == DesignMode::Empty);   // means user will set edges
    }
    void init_comp() {
        if (opt_.tau_selection == TauSelection::Automatic) { set_tau_auto_all_(); }
        clear_covariance_cache_();
    }

    // Parameters setters
    void set_lambda_all(const double lambda) const { for (auto& b : blocks_) b->set_lambda(lambda); }

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

        // random or SVD init for ν_l
        for (auto& b : blocks_) {
            b->set_h(h_);
            Vector nu;

            if (opt_.init == Init::SVD) {
                // Compute U(:,1) from thin SVD of X (n_obs x n_covs)
                const Matrix& X = b->data();
                // Use BDCSVD for larger problems, JacobiSVD is fine for smaller ones
                Eigen::BDCSVD<Matrix> svd(X, Eigen::ComputeThinU);
                nu = svd.matrixU().col(0);
            } else { // Random
                std::mt19937_64 rng(opt_.seed);
                std::uniform_real_distribution<double> U(-1.0, 1.0);
                nu = Vector::NullaryExpr(n_obs_, [&]{ return U(rng); });
            }
            b->l_compute(nu);  // block handles normalization
        }

        // all cov pairs are dirty; we'll fill on demand
        Result res(n_blocks());
        res.obj_history.reserve(opt_.max_iter + 1);
        res.obj_history.push_back(objective_());

        for (int s = 0; s < opt_.max_iter; ++s) {
            for (int l = 0; l < J; ++l) {
                Vector nu_l = Vector::Zero(n_obs_);
                const Vector eta_l = eta_(*blocks_[l]);
                for (int k = 0; k < J; ++k) {
                    if (k == l || !C_(l,k)) continue;   // <— exclude self
                    const Vector eta_k = eta_(*blocks_[k]);
                    const double cov_lk = cov_value_(l, k, eta_l, eta_k);   // uses/saves cache, marks clean
                    const double w_lk = scheme_.w(cov_lk);
                    nu_l.noalias() += w_lk * eta_k;   // no aliasing with RHS
                }
                blocks_[l]->l_compute(nu_l);   // block handles normalization
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
        compute_covariance_matrix_(res.covariance_matrix);
        get_tau(res.tau_values);

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
        for (std::size_t i = 0; i < tau_values.size(); ++i) {
            tau_values[i] = blocks_[i]->tau();
        }
    }

    void set_tau_auto_all_() const { for (auto& b : blocks_) b->set_tau(-1); }

    // ===== Helpers =====
    Vector eta_(Block& b) const {
        return b.components().col(h()); // η_j = X_j a_mj
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

private:
    int J_ {0};
    int n_obs_;   // global #observations
    int h_ {0};   // current component index
    Scheme scheme_;
    Options opt_;
    int n_comp_{0};
    bool is_last_comp_{false};

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
    os << "n_iters   : " << r.iters << "\n";
    os << "monotone  : " << (r.monotone ? "yes" : "no") << "\n";
    os << std::endl;
    os << "shrinkage parameters used : " << std::endl;
    for (size_t i = 0; i < r.tau_values.size(); ++i) {
        os << "- Block " << i+1  << ": tau = "<< r.tau_values[i] << "\n";
    }
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
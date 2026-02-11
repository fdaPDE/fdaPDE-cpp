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

namespace internals {

// Generic secant root finder for a scalar function f(mu)
double find_root_secant(std::function<double(double)> f,
                        double xL = 0.0,          // left starting point
                        double xR_init = 10.0,    // initial right point
                        double step = 10.0,        // how much to increase xR until sign changes
                        int max_expand = 100,      // max expansions to avoid infinite loop
                        double tol = 1e-8,
                        int max_iter = 100) {

    double fL = f(xL);
    if (!std::isfinite(fL))
        throw std::runtime_error("f(xL) not finite.");
    if (fL < 0)
        throw std::runtime_error("f(xL) < 0.");

    // Expand xR until we get fR < 0 (opposite sign)
    double xR = xR_init;
    double fR = f(xR);
    int expand_count = 0;
    while ((fL * fR > 0.0 || !std::isfinite(fR)) && expand_count < max_expand) {
        xR += step;
        fR = f(xR);
        expand_count++;
    }

    if (expand_count == max_expand)
        throw std::runtime_error("Unable to find a sign change: f(x) stays same sign.");

    // --- Secant iterations ---
    double x_prev = xL;
    double f_prev = fL;
    double x_curr = xR;
    double f_curr = fR;

    for (int iter = 0; iter < max_iter; ++iter) {
        if (std::fabs(f_curr - f_prev) < 1e-15)
            break; // avoid division by zero

        // Secant update
        const double x_next = x_curr - f_curr * (x_curr - x_prev) / (f_curr - f_prev);
        const double f_next = f(x_next);

        if (!std::isfinite(f_next))
            break; // probably diverged

        // Convergence check
        if (std::fabs(f_next) < tol || std::fabs(x_next - x_curr) < tol)
            return x_next;

        // Shift
        x_prev = x_curr;
        f_prev = f_curr;
        x_curr = x_next;
        f_curr = f_next;
    }

    // Return best current estimate
    return x_curr;
}

struct empty_t {
    template<class... Args>
    explicit empty_t(Args&&...) noexcept {}
};

struct identity_ls {
    using vector_t = Eigen::VectorXd;
    using matrix_t = Eigen::MatrixXd;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;

    explicit identity_ls(const int n_dofs) : n_dofs_(n_dofs) {}

    // Shapes / accessors used by GCV path (harmless no-ops here)
    [[nodiscard]] int n_dofs()  const { return n_dofs_; }
    [[nodiscard]] int n_obs()   const { return static_cast<int>(y_.size()); }
    [[nodiscard]] int n_covs()  const { return 0; } // no param covariates in identity model
    [[nodiscard]] double edf(int = 0, int = 0) const { return 0; } // hat-trace is 0 for identity

    // Data flow API
    void analyze_data() {
        Psi_.resize(n_dofs_, n_dofs_);
        Psi_.setIdentity();
    }

    void update_response_and_weights(const vector_t& y, const sparse_matrix_t& /*W*/) { y_ = y; }

    // Identity: fitted values equal response
    void fit(double /*lambda*/) { f_ = y_; }

    // Outputs
    [[nodiscard]] const vector_t& response() const { return y_; }
    [[nodiscard]] vector_t fn() const { return f_; }     // fitted values
    [[nodiscard]] const vector_t& f() const { return f_; }
    [[nodiscard]] const sparse_matrix_t& Psi() const { return Psi_; }
    [[nodiscard]] double ftPf(double) { return 0.; }

private:
    int n_dofs_{0};
    vector_t y_;
    vector_t f_;
    sparse_matrix_t Psi_;
};

}

struct IndependentSampling {
    using solver_t = internals::identity_ls;
};
struct TimeDependentSampling {
    using solver_t = internals::fe_ls_elliptic;
    using Matrix = Eigen::MatrixXd;
    using SparseMatrix = Eigen::SparseMatrix<double>;
    using PointEvalType = std::function<SparseMatrix(const Matrix&)>;

    static void discretize(const Triangulation<1, 1>& T, solver_t& solver_) {
        // define physic in space (same for all the blocks)
        FeSpace Vh(T, P1<1>);
        TrialFunction f_T(Vh);
        TestFunction  v_T(Vh);
        auto a_T = integral(T)(dx(f_T) * dx(v_T));
        ZeroField<1> u;
        auto F_T = integral(T)(u * v_T);
        auto penalty = fdapde::fe_ls_elliptic(a_T, F_T);
        solver_.discretize(penalty.get());
    }

    static void compute_Psi(const Triangulation<1, 1>& T, const Matrix& times, SparseMatrix& Psi) {
        // define physic in space (same for all the blocks)
        FeSpace Vh(T, P1<1>);
        TrialFunction f_T(Vh);
        TestFunction  v_T(Vh);
        auto a_T = integral(T)(dx(f_T) * dx(v_T));
        ZeroField<1> u;
        auto F_T = integral(T)(u * v_T);
        auto penalty = fdapde::fe_ls_elliptic(a_T, F_T);
        // compute point eval functor
        using BilinearForm = typename std::decay_t<decltype(penalty.get())>::BilinearForm;
        const BilinearForm& bilinear_form = penalty.get().bilinear_form();
        Psi = internals::point_basis_eval(bilinear_form.trial_space(), times);
    }

};

namespace internals {

template <typename SamplingStrategy> class BaseBlock;   // forward decl for operator<<
template <typename SamplingStrategy>
std::ostream& operator<<(std::ostream& os, const BaseBlock<SamplingStrategy>& b);

// GCV utils
template<class Fun> inline std::pair<double,double> argmin_over_log_grid(Fun&& f, const double log10_min, const double log10_max, int n_grid) {
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
    double log10_max = 0.0;
    int grid = 100;

    // edf() stochastic trace settings (if your solver uses Hutch++ etc.)
    int edf_r = 100;
    int edf_seed = 12345;

    // safety
    double eps_dof = 1e-12;  // avoid divide-by-zero in denominator
};
template <class Smoother> struct GCVEval {
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
template <typename SolverType> std::pair<bool, double> select_lambda_with_gcv(SolverType& solver, const GCVConfig& gcv_cfg) {
    GCVEval<SolverType> gcv{ &solver, gcv_cfg };
    auto [lambda_opt, gcv_opt] = argmin_over_log_grid(
        [&](double lam){ return gcv(lam); },
        gcv_cfg.log10_min, gcv_cfg.log10_max, gcv_cfg.grid
    );
    return {lambda_opt < std::pow(10.0, gcv_cfg.log10_max) , lambda_opt};
}

template <typename SamplingStrategy>
class BaseBlock {
public:
    using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using SparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
    using SparseSolver = eigen_sparse_solver_movable_wrap<Eigen::SimplicialLDLT<SparseMatrix>>;
    using ComponentsSolverType = typename std::decay_t<SamplingStrategy>::solver_t;

    template<typename S = SamplingStrategy>
    requires std::same_as<S, IndependentSampling>
    BaseBlock(const std::string& block_name, const Matrix& data, const int n_dofs_loadings, const double tau = 0.0) :
        block_name_(block_name), data_(data), components_solver_(data.rows()), n_dofs_loadings_(n_dofs_loadings), tau_(tau) {
        // Init components solver
        components_solver_.analyze_data();
    }

    template<typename S = SamplingStrategy>
    requires std::same_as<SamplingStrategy, TimeDependentSampling>
    BaseBlock(const std::string& block_name, const Triangulation<1, 1>& T, const Vector& times, const Matrix& data, const int n_dofs_loadings, const double tau = 0.0) :
        block_name_(block_name), times_(times), data_(data), n_dofs_loadings_(n_dofs_loadings), tau_(tau) {
        // Init sparse identity
        I_.resize(n_obs(), n_obs());
        I_.setIdentity();
        // Init components solver
        SamplingStrategy::discretize(T, components_solver_);
        components_solver_.analyze_data(Matrix{times}, Vector::Zero(n_obs()), I_);
    }

    virtual ~BaseBlock() = default;

    // ---- Uniform public API ----

    // Initialization
    void init() {
        ensure_M_();
        ensure_lc_();
    }

    // Data
    [[nodiscard]] const std::string& name() const { return block_name_; }
    [[nodiscard]] const Matrix& data() const { return data_; }
    Matrix& data() { invalidate_M_(); return data_; }

    // Dimensions
    [[nodiscard]] int n_obs() const { return static_cast<int>(data_.rows()); }
    [[nodiscard]] int n_covs() const { return static_cast<int>(data_.cols()); }
    [[nodiscard]] int n_dofs_loadings() const { return n_dofs_loadings_; }

    // Components
    [[nodiscard]] int n_comp() const { return n_comp_; }
    void set_n_comp(const int n_comp) {
        if (n_comp <= 0) throw std::invalid_argument("n_comp must be > 0");
        n_comp_ = n_comp;
        // force resize on next access
        loadings_ready_ = false;
        components_ready_ = false;
    }

    // Shrinkage parameter
    [[nodiscard]] double tau() const { return tau_; }
    void set_tau(const double tau) {
        if (tau > 0) tau_ = tau;
        else select_tau_auto_();
        invalidate_M_();
    }

    // Normalization parameter
    [[nodiscard]] double rho() const { return rho_; }

    // Bias flag
    void set_bias(const bool bias) { bias_ = bias; }

    // Components regularization utilities
    void set_lambda_components(const double lambda) {
        *lambda_components_ = lambda;
        if (lambda < 0.0) lambda_components_selection_ = true;
    }
    [[nodiscard]] double lambda_components() const {
        if constexpr (std::same_as<SamplingStrategy, IndependentSampling>) return std::numeric_limits<double>::quiet_NaN();
        if (!lambda_components_.has_value() && *lambda_components_ > 0.0) return *lambda_components_;
        return std::numeric_limits<double>::quiet_NaN();
    }
    void set_components_gcv_config(const GCVConfig& cfg) { components_gcv_cfg_ = cfg; }

    // Loadings regularization utilities
    virtual void set_lambda_loadings(const double) {};
    [[nodiscard]] virtual double lambda_loadings() const { return std::numeric_limits<double>::quiet_NaN(); }

    // Noise variance
    void set_noise_variance(const double noise_variance) { noise_variance_ = std::max(0.0, noise_variance); }
    [[nodiscard]] double noise_variance() const {
        if (!noise_variance_.has_value()) return std::numeric_limits<double>::quiet_NaN();
        return *noise_variance_;
    }

    // Inner-Component initialization
    struct InitInfo { bool active{false}; Vector nu; double s1{0}; double s1_edge{0}; double frac{0}; };
    InitInfo svd_init(const bool allow_block_deactivation, const double relaxation = 0.) const {
        InitInfo out{true, Vector::Zero(n_obs()), 0.0, 0.0, 0.0};
        const Matrix& X = data();
        if (X.size() == 0) return out;

        const Eigen::BDCSVD<Matrix> svd(X, Eigen::ComputeThinU);
        if (svd.singularValues().size() == 0) return out;

        if (allow_block_deactivation){
            out.s1 = svd.singularValues()(0);
            const double fro2 = X.squaredNorm();
            out.frac = (fro2 > 0.0) ? (out.s1*out.s1)/fro2 : 0.0;

            // If no σ² set, fall back to your energy test
            if (!noise_variance_.has_value()) {
                out.active = (out.frac >= 1e-3);
            } else {
                const double sigma = std::sqrt(std::max(0.0, *noise_variance_));
                const auto n = static_cast<double>(n_obs());
                const auto m = static_cast<double>(n_covs());
                out.s1_edge = sigma * (std::sqrt(n) + std::sqrt(m)) * (1.0 + relaxation);
                out.active = (out.s1 > out.s1_edge);
            }
        }

        if (out.active) {
            out.nu = svd.matrixU().col(0);
        } else {
            out.nu.setZero();
        }
        return out;
    }

    // M: normalization matrix
    [[nodiscard]] const SparseMatrix& M() const { ensure_M_(); return M_; }
    SparseSolver& invM() { ensure_M_(); return invM_; }

    // Current component index & Deflation
    [[nodiscard]] int h() const { return h_; }
    void set_h(const int idx) { if (idx < 0 || idx >= n_comp_) throw std::out_of_range("h"); h_ = idx; }
    void next_component() { set_h(h_ + 1); }
    void deflate(const Deflation mode) {
        if (h() == n_comp()) throw std::out_of_range("h");
        switch (mode) {
        case Deflation::Scores:   deflate_scores_(); break;
        case Deflation::Loadings: deflate_loadings_(); break;
        case Deflation::None: default: break;
        }
        // data_ changed -> M invalid; cached scores/loadings are now stale
        invalidate_M_();
    }

    // Main compute method
    void compute(const Vector& nu_D) {

        // Compute the loading
        const Vector a_tilde = l_fit_(nu_D);
        rho_ = compute_multipliers_(a_tilde);

        // Normalize the loading and the non-regularized component
        Vector a = a_tilde / rho_;
        Vector s = data() * Psi_D() * a;

        // Save the loading
        loadings().col(h()) = a;

        // Fit regularized scores and save it
        components().col(h()) = c_fit_(s);
    }

    // Scaling method
    void flip_and_scale_to_unit_score_variance(const SparseMatrix& Psi) {
        // Scale factor so that Var(eta_t) = 1 where eta_t has length that depends on Psi
        const Vector eta_t = Psi*eta_();
        const double den = bias_ ? eta_t.size() : std::max(1, static_cast<int>(eta_t.size()) - 1);
        const double v = (eta_t).squaredNorm() / den;
        if (v <= 0.0) return;
        const double norm = std::sqrt(v);

        // Sign flip
        double sign = 1;
        if (a_().mean() < 0) sign = -1.0;

        // Apply to loadings and scores
        components().col(h()) /= sign * norm;
        loadings().col(h()) /= sign * norm;
    }

    std::pair<double, double> reconstruction_constraint_info() {
        const Vector a_m = Psi_D() * a_();
        const Vector eta_t = Psi_T() * eta_();
        const Vector r = data() * a_m - eta_t;
        const double den = n_obs();
        const double mse = r.squaredNorm() / den;

        if (noise_variance_.has_value()) {
            const double edge = noise_variance() * a_m.squaredNorm();
            return {mse, edge};
        }
        return {mse, std::numeric_limits<double>::quiet_NaN()};
    }

    // Psi matrices
    [[nodiscard]] virtual const SparseMatrix& Psi_D() const = 0; // It depends on the loadings solver
    [[nodiscard]] const SparseMatrix& Psi_T() const { return components_solver_.Psi(); }

    // Penalty evaluation
    [[nodiscard]] double evaluate(const Vector& nu) {
        double a = 1./n_obs() * nu.transpose() * components_m().col(h());
        double b = atPa();
        return a; // - b;
    }
    [[nodiscard]] double ntPn() { return components_solver_.ftPf( lambda_components() ); }
    [[nodiscard]] virtual double atPa() { return 0.; }

    // Loadings & Components
    Matrix& loadings() { ensure_lc_(); return loadings_; }
    Matrix loadings_m() { ensure_lc_(); return Psi_D() * loadings_; }
    Matrix& components() { ensure_lc_(); return components_; }
    Matrix components_m() { ensure_lc_(); return Psi_T() * components_; }

    // Times getter (only active with TimeDependentSampling strategy)
    template <typename S = SamplingStrategy>
    requires std::same_as<S, TimeDependentSampling>
    const Vector& times() { return times_; }

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

protected:

    // Loadings solver
    virtual Vector l_fit_(const Vector& nu) = 0;

    // Component solver
    Vector c_fit_(const Vector& s) {

        // always go through the same solver API
        components_solver_.update_response_and_weights(s, I_);

        double lambda = 1e-15;
        if (!lambda_components_.has_value()){
            if (lambda_components_selection_) {
                auto [success, l] = select_lambda_with_gcv(components_solver_, components_gcv_cfg_);
                // if (!success) return Vector::Zero(n_dofs_loadings());
                *lambda_components_ = l;
                lambda_components_selection_ = false;
            }
            lambda = *lambda_components_;
        }
        components_solver_.fit(lambda);

        return components_solver_.f();
    }

    // Lagrange Multipliers utilities
    [[nodiscard]] double compute_multipliers_(const Vector& a_D) const {
        const Vector a_m_D = Psi_D()*a_D;
        // rho_ = 1.;
        const double norm_sqr = a_m_D.dot(M() * a_m_D) + atPa();
        if (norm_sqr <= 0.0) return 1;
        const double norm = std::sqrt(norm_sqr);
        return norm;
    }

    // Current loading and component getters
    [[nodiscard]] Vector a_() { ensure_lc_(); return loadings_.col(h());}
    [[nodiscard]] Vector eta_() { ensure_lc_(); return components_.col(h());}

    // tau estimate using Schäfer–Strimmer analytic shrinkage from correlation
    void select_tau_auto_() {
        const int n = n_obs(), m = n_covs();
        if (n < 2 || m < 1) throw std::runtime_error("tau_auto: need n>=2 and m>=1");

        // xs <- scale(x, center=TRUE, scale=TRUE)  [sample sd with (n-1)]
        Eigen::RowVectorXd mu  = data_.colwise().mean();
        Matrix xs = data_.rowwise() - mu;                                   // center
        Eigen::RowVectorXd var = (xs.array().square().colwise().sum() / static_cast<double>(n - 1)).matrix();
        Eigen::RowVectorXd sd  = var.array().sqrt().matrix();
        for (int j = 0; j < m; ++j) if (!(sd[j] > 0.0) || !std::isfinite(sd[j])) sd[j] = 1.0;
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
        set_tau(tau_hat);  // invalidates M_; recomputed lazily
    }

    // M
    void compute_M_() {
        SparseMatrix I(n_covs(), n_covs());
        I.setIdentity();
        M_ = tau_ * I;
        if (tau_ < 0.999) {
            const double den = bias_ ? n_obs() : std::max(1, n_obs() - 1);
            const Matrix dense = ((1.0 - tau_) / den) * (data_.transpose() * data_);
            M_ += dense.sparseView(1e-12);
        }
        M_.makeCompressed();
        invM_.compute(M_);
        M_ready_ = true;
    }
    void ensure_M_() const {
        if (!M_ready_) const_cast<BaseBlock*>(this)->compute_M_();
    }
    void invalidate_M_() { M_ready_ = false; }

    // Loadings and Components
    void ensure_lc_() {
        if (!loadings_ready_) {
            loadings_.setZero(n_dofs_loadings_, n_comp_);
            loadings_ready_ = true;
        }
        if (!components_ready_) {
            components_.setZero(components_solver_.n_dofs(), n_comp_);
            components_ready_ = true;
        }
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
        const int m = n_covs();

        // Assemble projection matrix
        Matrix R = Matrix::Identity(m, m);
        const double norm = a_m.squaredNorm();
        if (norm <= 0.0) return;
        R.noalias() -= (a_m * a_m.transpose()) / norm;

        // Apply right projection in variable space
        data_ = data_ * R;
    }

    // Components solver
    ComponentsSolverType components_solver_;

    // State
    Vector times_{0};
    const std::string block_name_;
    Matrix data_; // n_obs x n_covs
    int n_dofs_loadings_ {0};
    double tau_ {0.0};
    int n_comp_ {1};
    int h_ {0};
    bool bias_ = true;
    double rho_ {1.0};

    // Parameters
    GCVConfig components_gcv_cfg_;
    std::optional<double> lambda_components_;
    bool lambda_components_selection_ {false};
    std::optional<double> noise_variance_;

    // Results
    Matrix loadings_, components_;

    // Utilities
    SparseMatrix I_; // n_obs x n_obs identity matrix
    SparseMatrix M_;
    SparseSolver invM_;

    // Flags
    bool M_ready_ {false};
    bool loadings_ready_ {false}, components_ready_ {false};
};

// single non-member operator<< visible to all derived classes
template <typename SamplingStrategy>
inline std::ostream& operator<<(std::ostream& os, const BaseBlock<SamplingStrategy>& b) {
    b.print(os);   // virtual dispatch -> works for Multivariate/Functional too
    return os;
}

// ========== MultivariateBlock ==========
template <typename SamplingStrategy>
class MultivariateBlock final : public BaseBlock<SamplingStrategy> {
public:
    using Base = BaseBlock<SamplingStrategy>;
    using Matrix = typename Base::Matrix;
    using Vector = typename Base::Vector;
    using SparseMatrix = typename Base::SparseMatrix;

    using Base::init;
    using Base::M;
    using Base::invM;
    using Base::n_obs;
    using Base::n_covs;
    using Base::n_dofs_loadings;
    using Base::data;
    using Base::loadings;
    using Base::components;
    using Base::h;

    template<typename S = SamplingStrategy>
    requires std::same_as<S, IndependentSampling>
    MultivariateBlock(const std::string& block_name, const Matrix& X, const double tau = 0.0) :
        Base(block_name, X, static_cast<int>(X.cols()), tau) {
        init_multivariate();
    }

    template<typename S = SamplingStrategy>
    requires std::same_as<SamplingStrategy, TimeDependentSampling>
    MultivariateBlock(const std::string& block_name, const Triangulation<1, 1>& T, const Vector& times, const Matrix& X,  const double tau = 0.0) :
        Base(block_name, T, times, X, static_cast<int>(X.cols()), tau) {
        init_multivariate();
    }

    void init_multivariate() {
        Psi_D_.resize(n_covs(), n_dofs_loadings()); // n_covs == n_dofs_loadings in this case
        Psi_D_.setIdentity();
        init();
    }

    // Psi_D
    [[nodiscard]] const SparseMatrix& Psi_D() const override { return Psi_D_; }

    // Print
    void print(std::ostream& os) const override {
        Base::print(os);
        os << "type: MultivariateBlock, n_dofs_loadings = n_covs = " << n_dofs_loadings();
        os << "\n";
    }

protected:
    Vector l_fit_(const Vector& nu) override {
        assert(nu.size() == n_obs() && "nu must have size n_obs (rows of X)");
        init();
        return invM().solve(data().transpose() * nu);
    }

private:
    SparseMatrix Psi_D_; // n_covs x n_covs sparse identity matrix
};

// ========== FunctionalBlock ==========
template <class LoadingsPenaltyType, typename SamplingStrategy>
class FunctionalBlock final : public BaseBlock<SamplingStrategy> {
public:
    using Base = BaseBlock<SamplingStrategy>;
    using Vector = typename Base::Vector;
    using Matrix = typename Base::Matrix;
    using SparseMatrix = typename Base::SparseMatrix;
    using CovariatesSolverType = typename std::decay_t<LoadingsPenaltyType>::solver_t;

    using Base::init;
    using Base::M;
    using Base::invM;
    using Base::n_obs;
    using Base::n_covs;
    using Base::n_dofs_loadings;
    using Base::data;
    using Base::components;
    using Base::h;
    using Base::rho;

    template <typename GeoFrame>
    requires std::same_as<SamplingStrategy, IndependentSampling>
    FunctionalBlock(const std::string& block_name, GeoFrame& gf, LoadingsPenaltyType&& loadings_penalty, const double tau = 0.0) :
        Base(block_name, gf[0].template col<double>(block_name).as_matrix().transpose(), loadings_penalty.get().bilinear_form().n_dofs(), tau) { // TODO get the correct number of dofs
        init_functional(gf, std::forward<LoadingsPenaltyType>(loadings_penalty));
    }

    template <typename GeoFrame>
    requires std::same_as<SamplingStrategy, TimeDependentSampling>
    FunctionalBlock(const std::string& block_name, const Triangulation<1, 1>& T, const Vector& times, GeoFrame& gf, LoadingsPenaltyType&& loadings_penalty, const double tau = 0.0) :
        Base(block_name, T, times, gf[0].template col<double>(block_name).as_matrix().transpose(), loadings_penalty.get().bilinear_form().n_dofs(), tau) { // TODO get the correct number of dofs
        init_functional(gf, std::forward<LoadingsPenaltyType>(loadings_penalty));
    }

    template <typename GeoFrame>
    void init_functional(GeoFrame& gf, LoadingsPenaltyType&& loadings_penalty) {
        loadings_solver_.discretize(loadings_penalty.get());
        loadings_solver_.analyze_data(gf, M());
        init();
    }

    // Psi_D
    [[nodiscard]] const SparseMatrix& Psi_D() const override { return loadings_solver_.Psi(); }

    // Penalty evaluation
    [[nodiscard]] double atPa() override {
        if (success_) return loadings_solver_.ftPf(lambda_loadings()); //  * 0.5 / (rho()*rho());
        return 0.;
    }

    // Loadings regularization utilities
    void set_lambda_loadings(const double lambda) override { lambda_loadings_ = lambda; success_ = true; }
    [[nodiscard]] double lambda_loadings() const override {
        if (lambda_loadings_ > 0) return lambda_loadings_;
        return std::numeric_limits<double>::quiet_NaN();
    }
    void set_loadings_gcv_config(const GCVConfig& cfg) { loadings_gcv_cfg_ = cfg; }

    // Print
    void print(std::ostream& os) const override {
        Base::print(os);
        os << "type: FunctionalBlock, n_dofs_loadings = " << n_dofs_loadings();
        os << ", lambda = " << lambda_loadings_;
        os << "\n";
    }
protected:
    Vector l_fit_(const Vector& nu) override {
        assert(nu.size() == n_obs() && "nu must have size n_obs (rows of X)");
        init();

        const Vector z = invM().solve(data().transpose() * nu);
        loadings_solver_.update_response_and_weights(z, M()*z.size()); // M()*z.size() because the solver normalizes inside

        // lambda selection if required
        if(lambda_loadings_ < 0.0) {
            if (success_){
                auto [success, lambda_opt] = select_lambda_with_gcv(loadings_solver_, loadings_gcv_cfg_);
                if (!success) {
                    success_ = false;
                    return Vector::Zero(n_dofs_loadings());
                }
                lambda_loadings_ = lambda_opt; // the optimal lambda is saved for subsequent calls
            } else {
                return Vector::Zero(n_dofs_loadings());
            }
        }
        loadings_solver_.fit(lambda_loadings_);
        return loadings_solver_.f();
    }
private:
    bool success_ = true;
    CovariatesSolverType loadings_solver_;
    double lambda_loadings_ = -1.0; // < 0 means "use GCV"
    GCVConfig loadings_gcv_cfg_;
};



template <typename SamplingStrategy, typename Matrix = Eigen::Matrix<double, Dynamic, Dynamic>>
requires std::same_as<SamplingStrategy, IndependentSampling>
inline std::unique_ptr<internals::BaseBlock<SamplingStrategy>>
make_multivariate_block(std::string block_name, const Matrix& data, double tau = 0.0) {
    return std::make_unique<internals::MultivariateBlock<SamplingStrategy>>(block_name, data, tau);
}

template <typename SamplingStrategy, typename Matrix = Eigen::Matrix<double, Dynamic, Dynamic>, typename Vector = Eigen::Matrix<double, Dynamic, 1>>
requires std::same_as<SamplingStrategy, TimeDependentSampling>
inline std::unique_ptr<internals::BaseBlock<SamplingStrategy>>
make_multivariate_block(std::string block_name, const Triangulation<1, 1>& T, const Vector& times, const Matrix& data, double tau = 0.0) {
    return std::make_unique<internals::MultivariateBlock<SamplingStrategy>>(block_name, T, times, data, tau);
}

template <typename SamplingStrategy, typename GeoFrame, typename LoadingsPenaltyType>
requires std::same_as<SamplingStrategy, IndependentSampling>
std::unique_ptr<internals::BaseBlock<SamplingStrategy>>
make_functional_block(std::string block_name, GeoFrame& gf, LoadingsPenaltyType&& loadings_penalty, double tau = 0.0) {
    return std::make_unique<internals::FunctionalBlock<LoadingsPenaltyType, SamplingStrategy>>(
      block_name, gf, std::forward<LoadingsPenaltyType>(loadings_penalty), tau);
}

template <typename SamplingStrategy, typename GeoFrame, typename LoadingsPenaltyType, typename Vector = Eigen::Matrix<double, Dynamic, 1>>
requires std::same_as<SamplingStrategy, TimeDependentSampling>
std::unique_ptr<internals::BaseBlock<SamplingStrategy>>
make_functional_block(std::string block_name, const Triangulation<1, 1>& T, const Vector& times, GeoFrame& gf, LoadingsPenaltyType&& loadings_penalty, double tau = 0.0) {
    return std::make_unique<internals::FunctionalBlock<LoadingsPenaltyType, SamplingStrategy>>(
        block_name, T, times, gf, std::forward<LoadingsPenaltyType>(loadings_penalty), tau);
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
    using BoolMatrix = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;

    int h = 0;
    int J = 0;
    std::vector<double> obj_history;
    std::vector<double> loss_history;
    std::vector<double> space_reg_history;
    std::vector<double> time_reg_history;
    bool monotone = true;
    int iters = 0;
    BoolMatrix C;
    Matrix covariance_matrix;
    double noise_variance = 0.0;
    std::vector<double> tau_values;
    std::vector<double> lambda_components_values;
    std::vector<double> lambda_loadings_values;
    std::vector<bool> active_blocks;
    std::vector<double> s1_blocks;
    std::vector<double> s1_edge_blocks;
    std::vector<double> reconstruction_error;
    std::vector<double> reconstruction_edge;

    explicit Result(const int n_blocks) : J(n_blocks), C(J, J), covariance_matrix(J,J),
    tau_values(J), lambda_components_values(J), lambda_loadings_values(J), active_blocks(J), s1_blocks(J), s1_edge_blocks(J),
    reconstruction_error(J), reconstruction_edge(J) {}

};

// forward declaration of pretty printers
std::ostream& operator<<(std::ostream& os, const Result& r);
std::ostream& operator<<(std::ostream& os, const std::vector<Result>& results);

// ===== RGCCA =====
template <typename SamplingStrategy>
class RGCCA {
public:
    using Block = internals::BaseBlock<SamplingStrategy>;
    using BlockPtr = std::unique_ptr<Block>;
    using Matrix = typename Block::Matrix;
    using BoolMatrix = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMatrix = typename Block::SparseMatrix;
    using Vector = typename Block::Vector;
    using SamplingDomain = std::conditional_t<std::same_as<SamplingStrategy, TimeDependentSampling>, Triangulation<1, 1>, internals::empty_t>;

    struct Options {
        int max_iter;
        double tol;
        unsigned seed;
        bool verbose;
        bool cache_covariances;
        bool flip_and_scale;
        bool allow_blocks_deactivation;
        bool bias;
        bool allow_reconstruction_constraint_compensation;
        Init init;
        LambdaSelection lambda_selection;
        TauSelection tau_selection;
        Deflation deflation_mode;
        Scheme scheme;

        explicit Options(
          const int max_iter_ = 1000, const double tol_ = 1e-8, const unsigned seed_ = 0,
          const bool flip_and_scale_ = true,
          const bool allow_blocks_deactivation_ = true,
          const bool bias_ = true,
          const bool allow_reconstruction_constraint_compensation_ = false,
          const Init init_ = Init::SVD, const TauSelection tau_selection_ = TauSelection::Automatic,
          const LambdaSelection lambda_selection_ = LambdaSelection::Automatic,
          const Deflation deflation_mode_ = Deflation::Scores, const Scheme& scheme_ = Scheme::Factorial(),
          const bool verbose_ = false, const bool cache_ = true) :
            max_iter(max_iter_),
            tol(tol_),
            seed(seed_),
            flip_and_scale(flip_and_scale_),
            allow_blocks_deactivation(allow_blocks_deactivation_),
            bias(bias_),
            allow_reconstruction_constraint_compensation(allow_reconstruction_constraint_compensation_),
            init(init_),
            tau_selection(tau_selection_),
            lambda_selection(lambda_selection_),
            deflation_mode(deflation_mode_),
            scheme(scheme_),
            verbose(verbose_),
            cache_covariances(cache_) { }
    };

    template <typename S = SamplingStrategy>
    requires std::same_as<S, IndependentSampling>
    explicit RGCCA(const int n_obs, const Options& opt = Options(), const int n_comp = 1) :
        n_obs_(n_obs), opt_(opt), n_comp_(n_comp) {}

    template <typename S = SamplingStrategy>
    requires std::same_as<S, TimeDependentSampling>
    explicit RGCCA(const int n_obs, const Triangulation<1, 1>& T, const Options& opt = Options(), const int n_comp = 1) :
        n_obs_(n_obs), T_(T), opt_(opt), n_comp_(n_comp) {}

    // ===== Blocks =====
    int add_block(BlockPtr b) {
        if (!b) throw std::invalid_argument("RGCCA/add_block: null block");
        if constexpr (std::same_as<SamplingStrategy, IndependentSampling>){
            if (b->n_obs() != n_obs()) throw std::invalid_argument("RGCCA/add_block: n_obs mismatch");
        } else { add_times_(b->times()); }
        b->set_n_comp(n_comp());
        b->set_bias(opt_.bias);
        blocks_.emplace_back(std::move(b));
        initialized_ = false;   // topology/caches need a fresh init later
        return ++J_;
    }

    template<typename S = SamplingStrategy>
    requires std::same_as<S, IndependentSampling>
    int add_multivariate_block(std::string block_name, Matrix& X, const double tau = 0.0) {
        return add_block(internals::make_multivariate_block<SamplingStrategy>(block_name, X, tau));
    }
    template<typename S = SamplingStrategy>
    requires std::same_as<S, TimeDependentSampling>
    int add_multivariate_block(std::string block_name, const Vector& times, Matrix& X, const double tau = 0.0) {
        return add_block(internals::make_multivariate_block<SamplingStrategy>(block_name, T_, times, X, tau));
    }
    template <typename GeoFrame, typename LoadingsPenaltyType>
    requires std::same_as<SamplingStrategy, IndependentSampling>
    int add_functional_block(std::string block_name, const GeoFrame& gf, LoadingsPenaltyType&& loadings_penalty, const double tau = 0.0) {
        return add_block(internals::make_functional_block<SamplingStrategy>(block_name, gf, std::forward<LoadingsPenaltyType>(loadings_penalty), tau));
    }
    template <typename GeoFrame, typename LoadingsPenaltyType>
    requires std::same_as<SamplingStrategy, TimeDependentSampling>
    int add_functional_block(std::string block_name, const Vector& times, const GeoFrame& gf, LoadingsPenaltyType&& loadings_penalty, const double tau = 0.0) {
        return add_block(internals::make_functional_block<SamplingStrategy>(block_name, T_, times, gf, std::forward<LoadingsPenaltyType>(loadings_penalty), tau));
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
        if (noise_variance_.has_value()) set_noise_variance_all_();
        compute_Psi_();
        initialized_ = true;
        user_defined_design_ = (mode == DesignMode::Empty);   // means user will set edges
    }
    void init_comp() {
        if (opt_.tau_selection == TauSelection::Automatic) { set_tau_auto_all_(); }
        if (opt_.lambda_selection == LambdaSelection::Automatic) { set_lambda_auto_all_(); }
        clear_covariance_cache_();
    }

    // Noise
    void set_noise_variance(double noise_variance) { noise_variance_ = std::max(0.0, noise_variance); }
    [[nodiscard]] double noise_variance() const {
        if (!noise_variance_.has_value()) return std::numeric_limits<double>::quiet_NaN();
        return *noise_variance_;
    }

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
                const auto info = b->svd_init(opt_.allow_blocks_deactivation);
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
        for (int j = 0; j < J; ++j) if (res.active_blocks[j]) {
            bool alone = true;
            for (int k = 0; k < J; ++k) alone &= !res.C(j,k);
            if (alone) res.C(j,j) = true;
        }

        // room for objective function evaluations
        res.obj_history.reserve(opt_.max_iter);
        res.loss_history.reserve(opt_.max_iter);
        res.space_reg_history.reserve(opt_.max_iter);
        res.time_reg_history.reserve(opt_.max_iter);
        {
            const auto [f_obj, f_loss, f_space_reg, f_time_reg]  = objective_(res.C, res.active_blocks);
            res.obj_history.push_back(f_obj);
            res.loss_history.push_back(f_loss);
            res.space_reg_history.push_back(f_space_reg);
            res.time_reg_history.push_back(f_time_reg);
        }
        if ( !no_connections_(res.C) ) {
            // require lambda selection also at the first iteration
            if (opt_.lambda_selection == LambdaSelection::Automatic) { set_lambda_auto_all_(); }
            auto a_prev = snapshot_loadings_();
            for (int s = 0; s < opt_.max_iter; ++s) {
                for (int l = 0; l < J; ++l) {
                    Vector nu_l = Vector::Zero(blocks_[l]->n_obs());
                    const Vector eta_l = eta_(*blocks_[l]);
                    for (int k = 0; k < J; ++k) {
                        if (!res.C(l,k)) continue;
                        const Vector eta_k = eta_(*blocks_[k]);
                        const double cov_lk = cov_value_(l, k, eta_l, eta_k);   // uses/saves cache, marks clean
                        const double w_lk = opt_.scheme.w(cov_lk);
                        nu_l.noalias() += w_lk * eta_(*blocks_[k], *blocks_[l]);   // no aliasing with RHS
                    }
                    /*double ev_prev = 0;
                    if (s>0) {
                        std::cout << "iter = " << s << ", j = " << l+1 << " : f(a_j^s  ) = ";
                        ev_prev = blocks_[l]->evaluate(nu_l);
                        std::cout << ev_prev << std::endl;
                    }*/
                    blocks_[l]->compute(nu_l);   // block handles normalization
                    /*if (s>0) {
                        std::cout << "                  f(a_j^s+1) = ";
                        double ev_post = blocks_[l]->evaluate(nu_l);
                        std::cout << ev_post << " improv = " << ev_post -  ev_prev << std::endl;
                    }*/
                    mark_cov_rowcol_dirty_(l);     // η_l changed → invalidate its row/col

                    // this is only to emulate the loadings of the R implementation, it could be dropped eventually
                    // bool even_scheme = (opt_.scheme.name == std::string("Centroid") || opt_.scheme.name == std::string("Factorial"));
                    // if (even_scheme && blocks_[l]->loadings().col(h())(0) < 0) {
                    //     blocks_[l]->loadings().col(h()) *= -1.0;
                    //     blocks_[l]->components().col(h()) *= -1.0;
                    // }
                }

                const auto [f_obj, f_loss, f_space_reg, f_time_reg] = objective_(res.C, res.active_blocks);
                const double obj_prev = res.obj_history.back();
                res.obj_history.push_back(f_obj);
                res.loss_history.push_back(f_loss);
                res.space_reg_history.push_back(f_space_reg);
                res.time_reg_history.push_back(f_time_reg);
                res.iters = s + 1;

                // check monotonicity
                if (f_obj + 1e-15 < obj_prev) res.monotone = false;

                // check convergence
                const double delta_obj = std::abs(f_obj - obj_prev);
                const double delta_a2  = loadings_delta2_(a_prev);
                if (delta_obj < opt_.tol || delta_a2 < opt_.tol) break;

                // update snapshot for next iter
                a_prev = snapshot_loadings_();
            }
        }
        if (opt_.flip_and_scale) flip_and_scale_all_to_unit_score_variance_();

        // save information about the iteration in the result struct
        res.noise_variance = noise_variance();
        compute_covariance_matrix_(res.covariance_matrix);
        get_tau(res.tau_values);
        get_lambdas(res.lambda_components_values, res.lambda_loadings_values);
        get_reconstruction_constraint_info(res.reconstruction_error, res.reconstruction_edge);

        return res;
    }

    // ===== Accessors =====
    [[nodiscard]] int n_obs() const { return n_obs_; }
    [[nodiscard]] int n_comp() const { return n_comp_; }
    [[nodiscard]] int n_blocks() const { return static_cast<int>(blocks_.size()); }
    [[nodiscard]] const Options& options() const { return opt_; }
    [[nodiscard]] const Scheme& scheme() const { return opt_.scheme; }
    [[nodiscard]] const std::vector<BlockPtr>& blocks() const { return blocks_; }
    [[nodiscard]] const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& C() const { return C_; }
    [[nodiscard]] bool initialized() const { return initialized_; }
    [[nodiscard]] bool user_defined_design() const { return user_defined_design_; }
    [[nodiscard]] const SparseMatrix& Psi_T() const { return Psi_T_; };

private:

    double loadings_delta2_(const std::vector<typename Block::Vector>& a_prev) const {
        const int J = n_blocks();
        double acc = 0.0;
        for (int j = 0; j < J; ++j) {
            // skip inactive blocks if you want exact R behavior after deactivation
            const auto aj = blocks_[j]->loadings().col(h_);
            const auto dj = aj - a_prev[j];
            acc += dj.squaredNorm();
        }
        return acc;
    }

    std::vector<typename Block::Vector> snapshot_loadings_() const {
        const int J = n_blocks();
        std::vector<typename Block::Vector> out;
        out.reserve(J);
        for (int j = 0; j < J; ++j) out.push_back(blocks_[j]->loadings().col(h_));
        return out;
    }

    void flip_and_scale_all_to_unit_score_variance_() {
        for (auto& b : blocks_) {
            b->flip_and_scale_to_unit_score_variance(Psi_T());
        }
    }

    template <typename S = SamplingStrategy>
    requires std::same_as<S, TimeDependentSampling>
    void add_times_(const Vector& t) {
        times_.reserve(times_.size() + static_cast<size_t>(t.size()));
        times_.insert(times_.end(), t.data(), t.data() + t.size());
    }

    void compute_Psi_() {
        if constexpr (std::same_as<SamplingStrategy, TimeDependentSampling>) {
            std::ranges::sort(times_);
            times_.erase(std::ranges::unique(times_).begin(), times_.end());
            Eigen::VectorXd times_eig = Eigen::Map<Eigen::VectorXd>(times_.data(), times_.size());
            SamplingStrategy::compute_Psi(T_, Matrix{times_eig}, Psi_T_);
        } else {
            Psi_T_.resize(n_obs(), n_obs());
            Psi_T_.setIdentity();
        }
    }

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

    void get_reconstruction_constraint_info(std::vector<double> & reconstruction_error, std::vector<double> & reconstruction_edge) {
        const int J = n_blocks();
        for (std::size_t j = 0; j < J; ++j) {
            auto [error, edge] = blocks_[j]->reconstruction_constraint_info();
            reconstruction_error[j] = error;
            reconstruction_edge[j] = edge;
        }
    }

    void set_tau_auto_all_() const { for (auto& b : blocks_) b->set_tau(-1); }
    void set_lambda_auto_all_() const { set_lambda_components_all(-1); set_lambda_loadings_all(-1); }

    void set_noise_variance_all_() const { for (auto& b : blocks_) b->set_noise_variance(*noise_variance_); }

    // ===== Helpers =====

    // eta using the RGCCA own Psi_T (or components_m for independent)
    Vector eta_(Block& b) const {
        if constexpr (std::same_as<SamplingStrategy, TimeDependentSampling>) {
            return Psi_T() * b.components().col(h());
        } else {
            return b.components_m().col(h());
        }
    }
    // eta using reference block's Psi_T
    Vector eta_(Block& b, const Block& ref) const {
        if constexpr (std::same_as<SamplingStrategy, TimeDependentSampling>) {
            return ref.Psi_T() * b.components().col(h());
        } else {
            return b.components_m().col(h());
        }
    }

    // covariance of two vectors
    double cov_(const Vector& u, const Vector& v) const {
        double den = opt_.bias ? u.size() : u.size()-1;
        return u.dot(v) / den;
    }

    // objective f = Σ_{j,k} C_jk * g( cov(η_j, η_k) )
    std::tuple<double, double, double, double> objective_(const BoolMatrix& C, const std::vector<bool>& active_blocks) {
        // std::cout << "Computing objective --->"<<  std::endl;
        const int J = n_blocks();
        double f = 0.0;
        for (int j = 0; j < J; ++j) {
            const Vector eta_j = eta_(*blocks_[j]);
            for (int k = j; k < J; ++k){
                if (C(j, k)) { // c_jk *
                    const double cov_jk = cov_value_(j, k, eta_j, eta_(*blocks_[k]));
                    // std::cout << "j: " << j << ", k: " << k << " -> C_jk:" << cjk << std::endl;
                    double mult = j==k ? 1.0 : 2.0;
                    f += mult * opt_.scheme.g(cov_jk);
                }
            }

        }
        const double f_loss = f;
        double f_space_reg = 0.;
        double f_time_reg = 0.;
        /* for (int j = 0; j < J; ++j) {
            if (active_blocks[j]) {
                const double atPa = blocks_[j]->atPa();
                // std::cout << "j: "<< j <<" atPa = " << atPa << std::endl;
                f_space_reg += atPa;
                const double ntPn = blocks_[j]->ntPn();
                // f_time_reg += ntPn;
                // std::cout << "j: "<< j <<" ntPn = " << ntPn << std::endl;
                f -= atPa + ntPn; //
            }
        } */
        return {f, f_loss, f_space_reg, f_time_reg};
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
    int n_obs_ {0}; // global number of observations
    std::vector<double> times_;
    SamplingDomain T_; // only used by TimeDependentSampling
    SparseMatrix Psi_T_;
    int h_ {0};   // current component index
    Options opt_;
    int n_comp_{0};
    bool is_last_comp_{false};

    std::optional<double> noise_variance_;

    std::vector<BlockPtr> blocks_;

    // topology & caches (sized in init())
    bool initialized_ {false};
    bool user_defined_design_ {false};
    BoolMatrix C_;

    Matrix Cov_ {0, 0};          // cached covariances between η's (Cov_(j,j)=1)
    Eigen::ArrayXXi dirty_ {0, 0};   // 1=dirty, 0=clean
};


// Pretty printer for a single Result
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
            os << "- Block " << i+1  << ": " << (r.active_blocks[i] ? "active    " : "non-active" )
               << " ("<< r.s1_blocks[i]<< (r.active_blocks[i] ? " > " : " < ") << r.s1_edge_blocks[i] << ")" << "\n";
        }
        os << std::endl;
        os << "(updated) connections matrix :\n";
        os << r.C << std::endl;
        os << std::endl;

        os << "regularization parameters used : " << std::endl;
        os << std::scientific;
        for (size_t i = 0; i < r.tau_values.size(); ++i) {
            os << "- Block " << i+1  << ": lambda_c = "<< r.lambda_components_values[i]
               << ", lambda_l = "<< r.lambda_loadings_values[i] << "\n";
        }
        os << std::fixed;
        os << std::endl;
    }
    os << "n_iters   : " << r.iters << "\n";
    os << "monotone  : " << (r.monotone ? "yes" : "no") << "\n";
    os << std::endl;
    os << "objective :\n";
    double prev_obj = r.obj_history[0];
    double prev_loss = r.loss_history[0];
    double prev_space_reg = r.space_reg_history[0];
    double prev_time_reg = r.time_reg_history[0];
    /*
    os << "- iter " << std::setw(3) << (0)
   << "   |   fit = " << std::setw(12) << std::setprecision(8) << std::fixed << prev_obj << " = "
   << std::setw(6) << std::setprecision(4) << std::fixed << prev_loss << " (" <<  ")" << " - "
   << std::setw(6) << std::setprecision(4) << std::fixed << prev_space_reg << " (" << ")" << " - "
   << std::setw(6) << std::setprecision(4) << std::fixed << prev_time_reg << " (" << ")"
   << "   |   overall diff = " << std::setw(7) << "\n";
    os << std::fixed << std::setprecision(8);
    */
    for (size_t i = 1; i < r.obj_history.size(); ++i) {
        const double obj = r.obj_history[i];
        const double loss = r.loss_history[i];
        const double space_reg = r.space_reg_history[i];
        const double time_reg = r.time_reg_history[i];
        os << "- iter " << std::setw(3) << (i)
           << "   |   fit = " << std::setw(12) << std::setprecision(8) << std::fixed << obj << " = "
           << std::setw(6) << std::setprecision(4) << std::fixed << loss << " (" << ((loss - prev_loss) >= 0 ? "+" : "-") << std::setprecision(1) << std::scientific << std::abs(loss - prev_loss) << ")" << " - "
           << std::setw(6) << std::setprecision(4) << std::fixed << space_reg << " (" << ((space_reg - prev_space_reg) > 0 ? "+" : "-") << std::setprecision(1) << std::scientific << std::abs(space_reg - prev_space_reg) << ")" << " - "
           << std::setw(6) << std::setprecision(4) << std::fixed << time_reg << " (" << ((time_reg - prev_time_reg) > 0 ? "+" : "-") << std::setprecision(1) << std::scientific << std::abs(time_reg - prev_time_reg) << ")"
           << "   |   overall diff = " << std::setw(7) << (obj - prev_obj) << "\n";
        os << std::fixed << std::setprecision(8);
        prev_obj = obj;
        prev_loss = loss;
        prev_space_reg = space_reg;
        prev_time_reg = time_reg;
    }
    os << std::endl;
    if (!minimal) {
        os << "reconstruction constraint :\n";
        for (size_t i = 0; i < r.reconstruction_error.size(); ++i) {
            if (!r.active_blocks[i] || r.reconstruction_edge[i] == 0) {
                os << "- Block " << i+1  << ": " << "non-active" << "\n";
            } else {
                const bool check = r.reconstruction_error[i] <= r.reconstruction_edge[i];
                os << "- Block " << i+1  << ": " << (check ? "satisfied    " : "not-satisfied" )
                   << " ("<< std::setw(10) << r.reconstruction_error[i] << (check ? " ≤ " : " > ") << std::setw(10) << r.reconstruction_edge[i] << ")";
                os << ", equality for σ_noise = "
                   << std::sqrt(r.noise_variance) << " -> "
                   << std::sqrt(r.reconstruction_error[i]/r.reconstruction_edge[i] * r.noise_variance);
                std::cout << "\n";
            }
        }
        os << std::endl;
        os << "covariance matrix :\n";
        os << std::fixed << std::setprecision(2);
        os << r.covariance_matrix << std::endl;
        os << std::fixed << std::setprecision(8);
    }

    return os;
}

// Pretty printer for a vector of Result (components)
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

}

#endif   // __FGCCA_H__
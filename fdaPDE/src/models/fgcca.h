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

#include "fdaPDE/src/solvers/nonnegative_ipopt.h"
#include "header_check.h"

namespace fdapde {

enum class Init { Random, SVD };
enum class DesignMode {Empty, FullyConnected};
enum class LambdaSelection {Manual, Automatic};
enum class Mode { CorMax, Regularized, CovMax };
enum class Deflation { None, Scores };
enum class WeightSignConstraint { None, NonNegative };

namespace internals {

void ginv(const Eigen::MatrixXd& X, Eigen::MatrixXd& ginvX, double tol = std::sqrt(std::numeric_limits<double>::epsilon())){
    // SVD
    Eigen::BDCSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const auto& d = svd.singularValues();
    const double d1 = d(0);
    const double thresh = std::max(tol * d1, 0.0);

    // Identify Positive singular values
    std::vector<int> idx;
    idx.reserve(d.size());
    for (int i = 0; i < d.size(); ++i) {
        if (d(i) > thresh) idx.push_back(i);
    }

    // Compute V_pos * diag(1/d_pos) * U_pos^T * z without forming full matrices
    Eigen::MatrixXd Upos(X.rows(), (int)idx.size());
    Eigen::MatrixXd Vpos(X.cols(), (int)idx.size());
    Eigen::VectorXd invd((int)idx.size());

    for (int k = 0; k < (int)idx.size(); ++k) {
        Upos.col(k) = svd.matrixU().col(idx[k]);
        Vpos.col(k) = svd.matrixV().col(idx[k]);
        invd(k) = 1.0 / d(idx[k]);
    }

    ginvX = Vpos * (invd.asDiagonal() * (Upos.transpose()));
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

    // Identity: fitted values equal z
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
    using solver_t = internals::bs_ls_elliptic;
    using Matrix = Eigen::MatrixXd;
    using SparseMatrix = Eigen::SparseMatrix<double>;
    using PointEvalType = std::function<SparseMatrix(const Matrix&)>;

    static void discretize(const Triangulation<1, 1>& T, solver_t& solver_) {
        // define physic in space (same for all the blocks)
        BsSpace Bh(T, 3);
        TrialFunction f_T(Bh);
        TestFunction  v_T(Bh);
        auto a_T = integral(T)(dxx(f_T) * dxx(v_T));
        ZeroField<1> u_T;
        auto F_T = integral(T)(u_T * v_T);
        auto penalty = fdapde::bs_ls_elliptic(a_T, F_T);
        solver_.discretize(penalty.get());
    }

    static void compute_Psi(const Triangulation<1, 1>& T, const Matrix& times, SparseMatrix& Psi) {
        // define physic in space (same for all the blocks)
        BsSpace Bh(T, 3);
        TrialFunction f_T(Bh);
        TestFunction  v_T(Bh);
        auto a_T = integral(T)(dx(f_T) * dx(v_T));
        ZeroField<1> u_T;
        auto F_T = integral(T)(u_T * v_T);
        auto penalty = fdapde::bs_ls_elliptic(a_T, F_T);
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
        double gcv_index = (static_cast<double>(n) / (dor * dor)) * rss;
        return gcv_index;
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
    BaseBlock(const std::string& block_name, const Matrix& data, const int n_dofs_weights) :
        block_name_(block_name), data_(data), components_solver_(data.rows()), n_dofs_weights_(n_dofs_weights) {
        // Init components solver
        components_solver_.analyze_data();
    }

    template<typename S = SamplingStrategy>
    requires std::same_as<SamplingStrategy, TimeDependentSampling>
    BaseBlock(const std::string& block_name, const Triangulation<1, 1>& T, const Vector& times, const Matrix& data, const int n_dofs_weights) :
        block_name_(block_name), times_(times), data_(data), n_dofs_weights_(n_dofs_weights) {
        // Init sparse identity
        I_.resize(n(), n());
        I_.setIdentity();
        // Init components solver
        SamplingStrategy::discretize(T, components_solver_);
        components_solver_.analyze_data(Matrix{times}, Vector::Zero(n()), I_);
    }

    virtual ~BaseBlock() = default;

    // Initialization
    void init() {
        ensure_M_();
        ensure_lc_();
    }

    // Data
    [[nodiscard]] const std::string& name() const { return block_name_; }
    [[nodiscard]] const Matrix& data() const { return data_; }

    // Dimensions
    [[nodiscard]] int n() const { return static_cast<int>(data_.rows()); }
    [[nodiscard]] int m() const { return static_cast<int>(data_.cols()); }
    [[nodiscard]] int n_dofs_weights() const { return n_dofs_weights_; }

    // Components
    [[nodiscard]] int n_comp() const { return n_comp_; }
    void set_n_comp(const int n_comp) {
        if (n_comp <= 0) throw std::invalid_argument("n_comp must be > 0");
        n_comp_ = n_comp;
        // force resize on next access
        weights_ready_ = false;
        components_ready_ = false;
    }

    // Mode & Shrinkage parameter
    [[nodiscard]] double tau() const { return tau_; }
    void select_tau_auto() { select_tau_auto_(); }
    void set_mode(const Mode mode) {
        mode_ = mode;
        if (mode_ == Mode::CorMax) tau_ = 0.0;
        if (mode_ == Mode::CovMax) tau_ = 1.0;
        if (mode_ == Mode::Regularized) select_tau_auto_();
        invalidate_M_();
    }
    [[nodiscard]] Mode mode() const { return mode_; }

    // Bias flag
    void set_bias(const bool bias) { bias_ = bias; }

    // Weights sign constraint
    void set_weight_sign_constraint(const WeightSignConstraint weight_sign_constraint = WeightSignConstraint::None) {
        weight_sign_constraint_ = weight_sign_constraint;
    }
    [[nodiscard]] WeightSignConstraint weight_sign_constraint() const { return weight_sign_constraint_; }

    // Components regularization utilities
    void set_lambda_components(const double lambda) {
        *lambda_components_ = lambda;
        if (lambda < 0.0) lambda_components_selection_ = true;
    }
    [[nodiscard]] double lambda_components() const {
        if constexpr (std::same_as<SamplingStrategy, IndependentSampling>) return std::numeric_limits<double>::quiet_NaN();
        if (lambda_components_.has_value() && *lambda_components_ > 0.0) return *lambda_components_;
        return std::numeric_limits<double>::quiet_NaN();
    }
    void set_components_gcv_config(const GCVConfig& cfg) { components_gcv_cfg_ = cfg; }

    // Weights regularization utilities
    virtual void set_lambda_weights(const double) {};
    [[nodiscard]] virtual double lambda_weights() const { return std::numeric_limits<double>::quiet_NaN(); }

    // Noise variance
    void set_noise_variance(const double noise_variance) { noise_variance_ = std::max(0.0, noise_variance); }
    [[nodiscard]] double noise_variance() const {
        if (!noise_variance_.has_value()) return std::numeric_limits<double>::quiet_NaN();
        return *noise_variance_;
    }

    // Inner-Component initialization
    struct InitInfo { bool active{false}; Vector nu; };
    InitInfo svd_init() const {
        InitInfo out{true, Vector::Zero(n()) };
        const Matrix& X = data();
        if (X.size() == 0) return out;

        const Eigen::BDCSVD<Matrix> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);

        if (svd.singularValues().size() == 0) return out;

        if (X.rows() >= X.cols()) out.nu = X * svd.matrixV().col(0);
        else out.nu = svd.matrixU().col(0);

        return out;
    }

    // M: normalization matrix
    [[nodiscard]] const SparseMatrix& M() const { ensure_M_(); return M_; }
    [[nodiscard]] const Matrix& ginvM() const { ensure_ginvM_(); return ginvM_; }
    [[nodiscard]] SparseSolver& invM() { ensure_invM_(); return invM_; }

    // Current component index
    [[nodiscard]] int h() const { return h_; }
    void set_h(const int idx) { if (idx < 0 || idx >= n_comp_) throw std::out_of_range("h"); h_ = idx; }
    void next_component() { set_h(h_ + 1); }

    // Main compute method
    void compute(const Vector& nu_D) {

        // Compute the weight
        const Vector a = w_fit_(nu_D);
        weights().col(h()) = a;

        // Compute the component and regularize it (the regularization acts only in the TimeDependent sampling scenrio)
        const Vector s = data() * Psi_D() * a;
        components().col(h()) = c_fit_(s);
    }

    // Deflation
    void deflate(const Deflation mode) {
        if (h() == n_comp()) throw std::out_of_range("h");
        switch (mode) {
            case Deflation::Scores: deflate_scores_(); break;
            case Deflation::None: default: break;
        }
        invalidate_M_();
    }

    // Post-processing weights
    void compute_weights_star() {
        ensure_lc_();
        weights_star_.setZero(weights_.rows(), weights_.cols());
        for (int h = 0; h < n_comp_; ++h) {
            Vector a_star = weights_.col(h);
            if (h > 0) {
                const Matrix A_prev = weights_star_.leftCols(h);
                const Matrix P_prev = deflation_projections_.leftCols(h);
                const Vector coeff = P_prev.transpose() * (Psi_D() * weights_.col(h));
                a_star.noalias() -= A_prev * coeff;
            }
            weights_star_.col(h) = a_star;
        }
    }

    std::pair<double, double> reconstruction_constraint_info() {
        const Vector a_m = Psi_D() * a_();
        const Vector eta_t = Psi_T() * eta_();
        const Vector r = data() * a_m - eta_t;
        const double den = n();
        const double mse = r.squaredNorm() / den;

        if (noise_variance_.has_value()) {
            const double edge = noise_variance() * a_m.squaredNorm();
            return {mse, edge};
        }
        return {mse, std::numeric_limits<double>::quiet_NaN()};
    }

    // Psi matrices
    [[nodiscard]] virtual const SparseMatrix& Psi_D() const = 0; // It depends on the weights' solver
    [[nodiscard]] const SparseMatrix& Psi_T() const { return components_solver_.Psi(); }

    // Weights & Components
    Matrix& weights() { ensure_lc_(); return weights_; }
    Matrix weights_m() { ensure_lc_(); return Psi_D() * weights_; }
    Matrix& weights_star() { ensure_lc_(); return weights_star_; }
    Matrix weights_star_m() { ensure_lc_(); return Psi_D() * weights_star_; }
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

    // Weights solver
    virtual Vector w_fit_(const Vector& nu) = 0;

    // Component solver
    Vector c_fit_(const Vector& s) {

        components_solver_.update_response_and_weights(s, I_);

        double lambda = 1e-15;
        if (lambda_components_selection_) {
            auto [success, l] = select_lambda_with_gcv(components_solver_, components_gcv_cfg_);
            lambda_components_ = l;
            lambda_components_selection_ = false;
        }
        if (lambda_components_.has_value()) {
            lambda = *lambda_components_;
        }
        components_solver_.fit(lambda);

        return components_solver_.f();
    }

    // Current weight and component getters
    [[nodiscard]] Vector a_() { ensure_lc_(); return weights_.col(h());}
    [[nodiscard]] Vector eta_() { ensure_lc_(); return components_.col(h());}

    // tau estimate using Schäfer–Strimmer analytic shrinkage from correlation
    void select_tau_auto_() {
        const int n_obs  = n();
        const int n_vars = m();

        if (n_obs < 2 || n_vars < 1) throw std::runtime_error("tau_auto: need n >= 2 and m >= 1");

        Eigen::RowVectorXd mu = data_.colwise().mean();
        Matrix xs = data_.rowwise() - mu;
        Eigen::RowVectorXd var = (xs.array().square().colwise().sum() / static_cast<double>(n_obs - 1)).matrix();
        Eigen::RowVectorXd sd = var.array().sqrt().matrix();

        for (int j = 0; j < n_vars; ++j) {
            if (!(sd[j] > 0.0) || !std::isfinite(sd[j])) sd[j] = 1.0;
        }

        xs.array().rowwise() /= sd.array();
        const Matrix XtX = xs.transpose() * xs;
        const Matrix xs2 = xs.array().square().matrix();
        const Matrix xs2T_xs2 = xs2.transpose() * xs2;

        const double n_d = static_cast<double>(n_obs);
        const double c = n_d / std::pow(n_d - 1.0, 3.0);

        Matrix V = c * (xs2T_xs2 - (1.0 / n_d) * XtX.array().square().matrix());
        V.diagonal().setZero();

        const double num = V.sum();

        Matrix Corm = XtX / (n_d - 1.0);
        Matrix D = Corm;
        D.diagonal().array() -= 1.0;

        const double den = D.squaredNorm();

        tau_ = (den > 0.0) ? std::clamp(num / den, 0.0, 1.0) : 0.0;
        invalidate_M_();
    }

    // M
    void compute_M_() {
        M_.resize(m(), m());

        if (mode_ == Mode::CovMax) {
            M_.setIdentity();
        } else if (mode_ == Mode::CorMax) {
            const double den = bias_ ? n() : std::max(1, n() - 1);
            M_ = (data_.transpose() * data_ / den).sparseView();
        } else {
            SparseMatrix I(m(), m());
            I.setIdentity();

            const double den = bias_ ? n() : std::max(1, n() - 1);
            const Matrix Sigma = ((1.0 - tau_) / den) * (data_.transpose() * data_);

            M_ = tau_ * I;
            M_ += Sigma.sparseView();
        }

        M_.makeCompressed();
        M_ready_ = true;
        ginvM_ready_ = false;
        invM_ready_ = false;
    }
    void compute_invM_() {
        invM_.compute(M_);
        invM_ready_ = true;
    }
    void compute_ginvM_() {
        if (mode_ == Mode::CovMax) {
            ginvM_.setIdentity(m(), m());
        } else {
            ginvM_.resize(m(), m());
            Eigen::MatrixXd M_dense = Eigen::MatrixXd(M());
            ginv(M_dense, ginvM_);
        }
        ginvM_ready_ = true;
    }
    void ensure_M_() const {
        if (!M_ready_) const_cast<BaseBlock*>(this)->compute_M_();
    }
    void ensure_invM_() const {
        ensure_M_(); // sets ginvM_ready = false
        if (!invM_ready_) const_cast<BaseBlock*>(this)->compute_invM_();
    }
    void ensure_ginvM_() const {
        ensure_M_(); // sets ginvM_ready = false
        if (!ginvM_ready_) const_cast<BaseBlock*>(this)->compute_ginvM_();
    }
    void invalidate_M_() {
        M_ready_ = false;
        invalidate_derived_caches_();
    } // this is enough to invalidate also ginvM and invM
    virtual void invalidate_derived_caches_() {}

    // Weights and Components
    void ensure_lc_() {
        if (!weights_ready_) {
            weights_.setZero(n_dofs_weights_, n_comp_);
            weights_star_.setZero(n_dofs_weights_, n_comp_);
            deflation_projections_.setZero(m(), n_comp_);
            weights_ready_ = true;
        }
        if (!components_ready_) {
            components_.setZero(components_solver_.n_dofs(), n_comp_);
            components_ready_ = true;
        }
    }

    // Deflation
    void deflate_scores_() {
        ensure_lc_();

        const Vector y = Psi_T() * eta_();
        const double yy = y.squaredNorm();

        if (yy > 0) {
            Vector p = data_.transpose() * y / yy;

            // Store p_h for post-processing
            deflation_projections_.col(h()) = p;

            data_.noalias() -= y * p.transpose();
        }
    }

    Vector solve_nonnegative_weight_ipopt_(
        const SparseMatrix& Psi,
        const Matrix& Omega,
        const Vector& z,
        const std::vector<int>& dirichlet_dofs = {}
    ) const {

        const int n_weights = static_cast<int>(Psi.cols());

        auto* raw_problem = new NonNegativeWeightProblem(Psi, Omega, z,  dirichlet_dofs);
        Ipopt::SmartPtr<Ipopt::TNLP> problem = raw_problem;
        Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
        Ipopt::ApplicationReturnStatus status = app->Initialize();

        if (status != Ipopt::Solve_Succeeded) {
            std::cerr << "Ipopt initialization failed.\n";
        }

        status = app->OptimizeTNLP(problem);

        return raw_problem->solution();
    }

    // Components' solver
    ComponentsSolverType components_solver_;

    // Dimensions
    int n_dofs_weights_ {0};
    int n_comp_ {1};

    // State
    Vector times_{0};
    const std::string block_name_;
    Matrix data_;
    int h_ {0};

    // Options
    double tau_ {0.0};
    Mode mode_ = Mode::CorMax;
    WeightSignConstraint weight_sign_constraint_ = WeightSignConstraint::None;
    bool bias_ = true;

    // Parameters
    GCVConfig components_gcv_cfg_;
    std::optional<double> lambda_components_;
    std::optional<double> noise_variance_;

    // Results
    Matrix weights_, weights_star_, components_;
    Matrix deflation_projections_;

    // Utilities
    SparseMatrix I_; // n x n identity matrix
    SparseMatrix M_;
    Matrix ginvM_;
    SparseSolver invM_;

    // Flags
    bool M_ready_ {false}, invM_ready_ {false}, ginvM_ready_ {false};
    bool lambda_components_selection_ {false};
    bool weights_ready_ {false}, components_ready_ {false};
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
    using Base::ginvM;
    using Base::n;
    using Base::m;
    using Base::n_dofs_weights;
    using Base::data;
    using Base::weights;
    using Base::components;
    using Base::h;
    using Base::weight_sign_constraint;
    using Base::solve_nonnegative_weight_ipopt_;

    template<typename S = SamplingStrategy>
    requires std::same_as<S, IndependentSampling>
    MultivariateBlock(const std::string& block_name, const Matrix& X) :
        Base(block_name, X, static_cast<int>(X.cols())) {
        init_multivariate();
    }

    template<typename S = SamplingStrategy>
    requires std::same_as<SamplingStrategy, TimeDependentSampling>
    MultivariateBlock(const std::string& block_name, const Triangulation<1, 1>& T, const Vector& times, const Matrix& X) :
        Base(block_name, T, times, X, static_cast<int>(X.cols())) {
        init_multivariate();
    }

    void init_multivariate() {
        Psi_D_.resize(m(), n_dofs_weights()); // m == n_dofs_weights in this case
        Psi_D_.setIdentity();
        init();
    }

    // Psi_D
    [[nodiscard]] const SparseMatrix& Psi_D() const override { return Psi_D_; }

    // Print
    void print(std::ostream& os) const override {
        Base::print(os);
        os << "type: MultivariateBlock, n_dofs_weights = m = " << n_dofs_weights();
        os << "\n";
    }

protected:
    Vector w_fit_(const Vector& nu) override {
        assert(nu.size() == n() && "nu must have size n (rows of X)");
        init();

        Vector z = data().transpose();

        if (weight_sign_constraint() == WeightSignConstraint::NonNegative) {
            const double s = z.transpose() * weights().col(h());
            if (s*s > 0) z /= s;
            return solve_nonnegative_weight_ipopt_(Psi_D(), M(), z);
        }

        const Vector a_tilde = ginvM() * z; // If mode == Mode::CovMax, ginvM = I
        // Normalization
        double rho = a_tilde.dot(M() * a_tilde);
        if (rho <= 0.0) rho = 1.;
        rho = std::sqrt(rho);
        return a_tilde / rho;
    }

private:
    SparseMatrix Psi_D_; // m x m sparse identity matrix
};

// ========== FunctionalBlock ==========
template <class WeightsPenaltyType, typename SamplingStrategy>
class FunctionalBlock final : public BaseBlock<SamplingStrategy> {
public:
    using Base = BaseBlock<SamplingStrategy>;
    using Vector = typename Base::Vector;
    using Matrix = typename Base::Matrix;
    using SparseMatrix = typename Base::SparseMatrix;
    using WeightsSolverType = typename std::decay_t<WeightsPenaltyType>::solver_t;

    using Base::init;
    using Base::M;
    using Base::tau;
    using Base::n;
    using Base::m;
    using Base::n_dofs_weights;
    using Base::data;
    using Base::components;
    using Base::weights;
    using Base::h;
    using Base::weight_sign_constraint;
    using Base::solve_nonnegative_weight_ipopt_;

    template <typename GeoFrame>
    requires std::same_as<SamplingStrategy, IndependentSampling>
    FunctionalBlock(const std::string& block_name, GeoFrame& gf, WeightsPenaltyType&& weights_penalty) :
        Base(block_name, gf[0].template col<double>(block_name).as_matrix().transpose(), weights_penalty.get().bilinear_form().n_dofs()) {
        init_functional(gf, std::forward<WeightsPenaltyType>(weights_penalty));
    }

    template <typename GeoFrame>
    requires std::same_as<SamplingStrategy, TimeDependentSampling>
    FunctionalBlock(const std::string& block_name, const Triangulation<1, 1>& T, const Vector& times, GeoFrame& gf, WeightsPenaltyType&& weights_penalty) :
        Base(block_name, T, times, gf[0].template col<double>(block_name).as_matrix().transpose(), weights_penalty.get().bilinear_form().n_dofs()) {
        init_functional(gf, std::forward<WeightsPenaltyType>(weights_penalty));
    }

    template <typename GeoFrame>
    void init_functional(GeoFrame& gf, WeightsPenaltyType&& weights_penalty) {
        weights_solver_.discretize(weights_penalty.get());
        weights_solver_.analyze_data(gf, M());
        init();
    }

    // Psi_D
    [[nodiscard]] const SparseMatrix& Psi_D() const override { return weights_solver_.Psi(); }

    // Weights regularization utilities
    void set_lambda_weights(const double lambda) override {
        lambda_weights_ = lambda;
        Omega_ready_ = false;
    }
    [[nodiscard]] double lambda_weights() const override {
        if (lambda_weights_ > 0) return lambda_weights_;
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Omega matrix
    const Matrix& Omega() {
        if (!Omega_ready_) {
            Omega_ = Psi_D().transpose() * M() * Psi_D();
            Omega_ += lambda_weights_ * weights_solver_.P();
            Omega_ready_ = true;
        }

        return Omega_;
    }

    // Print
    void print(std::ostream& os) const override {
        Base::print(os);
        os << "type: FunctionalBlock, n_dofs_weights = " << n_dofs_weights();
        os << ", lambda = " << lambda_weights_;
        os << "\n";
    }
protected:
    Vector w_fit_(const Vector& nu) override {
        assert(nu.size() == n() && "nu must have size n (rows of X)");
        init();

        Vector z = data().transpose() * nu;

        if (weight_sign_constraint() == WeightSignConstraint::NonNegative) {
            const double s = z.transpose() * Psi_D() * weights().col(h());
            if (s*s > 0) z /= s;
            return solve_nonnegative_weight_ipopt_(Psi_D(), Omega(), z);
        }

        weights_solver_.update_z_and_weights(z, M());

        // fit
        weights_solver_.fit(lambda_weights_);
        return weights_solver_.f();
    }
    void invalidate_derived_caches_() override {
        Omega_ready_ = false;
    }
private:
    Matrix Omega_;
    bool Omega_ready_ {false};
    WeightsSolverType weights_solver_;
    double lambda_weights_ = 1e-15;
};



template <typename SamplingStrategy, typename Matrix = Eigen::Matrix<double, Dynamic, Dynamic>>
requires std::same_as<SamplingStrategy, IndependentSampling>
inline std::unique_ptr<internals::BaseBlock<SamplingStrategy>>
make_multivariate_block(std::string block_name, const Matrix& data) {
    return std::make_unique<internals::MultivariateBlock<SamplingStrategy>>(block_name, data);
}

template <typename SamplingStrategy, typename Matrix = Eigen::Matrix<double, Dynamic, Dynamic>, typename Vector = Eigen::Matrix<double, Dynamic, 1>>
requires std::same_as<SamplingStrategy, TimeDependentSampling>
inline std::unique_ptr<internals::BaseBlock<SamplingStrategy>>
make_multivariate_block(std::string block_name, const Triangulation<1, 1>& T, const Vector& times, const Matrix& data) {
    return std::make_unique<internals::MultivariateBlock<SamplingStrategy>>(block_name, T, times, data);
}

template <typename SamplingStrategy, typename GeoFrame, typename WeightsPenaltyType>
requires std::same_as<SamplingStrategy, IndependentSampling>
std::unique_ptr<internals::BaseBlock<SamplingStrategy>>
make_functional_block(std::string block_name, GeoFrame& gf, WeightsPenaltyType&& weights_penalty) {
    return std::make_unique<internals::FunctionalBlock<WeightsPenaltyType, SamplingStrategy>>(
      block_name, gf, std::forward<WeightsPenaltyType>(weights_penalty));
}

template <typename SamplingStrategy, typename GeoFrame, typename WeightsPenaltyType, typename Vector = Eigen::Matrix<double, Dynamic, 1>>
requires std::same_as<SamplingStrategy, TimeDependentSampling>
std::unique_ptr<internals::BaseBlock<SamplingStrategy>>
make_functional_block(std::string block_name, const Triangulation<1, 1>& T, const Vector& times, GeoFrame& gf, WeightsPenaltyType&& weights_penalty) {
    return std::make_unique<internals::FunctionalBlock<WeightsPenaltyType, SamplingStrategy>>(
        block_name, T, times, gf, std::forward<WeightsPenaltyType>(weights_penalty));
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
    bool monotone = true;
    int iters = 0;
    BoolMatrix C;
    Matrix covariance_matrix;
    double noise_variance = 0.0;
    std::vector<double> tau_values;
    std::vector<double> lambda_components_values;
    std::vector<double> lambda_weights_values;
    std::vector<bool> active_blocks;
    std::vector<double> s1_blocks;
    std::vector<double> s1_edge_blocks;
    std::vector<double> reconstruction_error;
    std::vector<double> reconstruction_edge;

    explicit Result(const int n_blocks) : J(n_blocks), C(J, J), covariance_matrix(J,J),
    tau_values(J), lambda_components_values(J), lambda_weights_values(J), active_blocks(J), s1_blocks(J), s1_edge_blocks(J),
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
        bool bias;
        Init init;
        LambdaSelection lambda_selection;
        Mode mode;
        WeightSignConstraint weight_sign_constraint;
        Deflation deflation_mode;
        Scheme scheme;

        explicit Options(
          const int max_iter_ = 1000, const double tol_ = 1e-8, const unsigned seed_ = 0,
          const bool bias_ = true,
          const Init init_ = Init::SVD, const Mode mode_ = Mode::CovMax,
          const WeightSignConstraint weight_sign_constraint_ = WeightSignConstraint::None,
          const LambdaSelection lambda_selection_ = LambdaSelection::Manual,
          const Deflation deflation_mode_ = Deflation::Scores, const Scheme& scheme_ = Scheme::Factorial(),
          const bool verbose_ = false, const bool cache_ = true) :
            max_iter(max_iter_),
            tol(tol_),
            seed(seed_),
            bias(bias_),
            init(init_),
            mode(mode_),
            weight_sign_constraint(weight_sign_constraint_),
            lambda_selection(lambda_selection_),
            deflation_mode(deflation_mode_),
            scheme(scheme_),
            verbose(verbose_),
            cache_covariances(cache_) { }
    };

    template <typename S = SamplingStrategy>
    requires std::same_as<S, IndependentSampling>
    explicit RGCCA(const int n, const Options& opt = Options(), const int n_comp = 1) :
        n_(n), opt_(opt), n_comp_(n_comp) {}

    template <typename S = SamplingStrategy>
    requires std::same_as<S, TimeDependentSampling>
    explicit RGCCA(const int n, const Triangulation<1, 1>& T, const Options& opt = Options(), const int n_comp = 1) :
        n_(n), T_(T), opt_(opt), n_comp_(n_comp) {}

    // ===== Blocks =====
    int add_block(BlockPtr b) {
        if (!b) throw std::invalid_argument("RGCCA/add_block: null block");
        if constexpr (std::same_as<SamplingStrategy, IndependentSampling>){
            if (b->n() != n()) throw std::invalid_argument("RGCCA/add_block: n mismatch");
        } else { add_times_(b->times()); }
        b->set_bias(opt_.bias);
        b->set_mode(opt_.mode);
        b->set_weight_sign_constraint(opt_.weight_sign_constraint);
        b->set_n_comp(n_comp());
        blocks_.emplace_back(std::move(b));
        initialized_ = false;   // topology/caches need a fresh init later
        return ++J_;
    }

    template<typename S = SamplingStrategy>
    requires std::same_as<S, IndependentSampling>
    int add_multivariate_block(std::string block_name, Matrix& X) {
        return add_block(internals::make_multivariate_block<SamplingStrategy>(block_name, X));
    }
    template<typename S = SamplingStrategy>
    requires std::same_as<S, TimeDependentSampling>
    int add_multivariate_block(std::string block_name, const Vector& times, Matrix& X) {
        return add_block(internals::make_multivariate_block<SamplingStrategy>(block_name, T_, times, X));
    }
    template <typename GeoFrame, typename WeightsPenaltyType>
    requires std::same_as<SamplingStrategy, IndependentSampling>
    int add_functional_block(std::string block_name, const GeoFrame& gf, WeightsPenaltyType&& weights_penalty) {
        return add_block(internals::make_functional_block<SamplingStrategy>(block_name, gf, std::forward<WeightsPenaltyType>(weights_penalty)));
    }
    template <typename GeoFrame, typename WeightsPenaltyType>
    requires std::same_as<SamplingStrategy, TimeDependentSampling>
    int add_functional_block(std::string block_name, const Vector& times, const GeoFrame& gf, WeightsPenaltyType&& weights_penalty) {
        return add_block(internals::make_functional_block<SamplingStrategy>(block_name, T_, times, gf, std::forward<WeightsPenaltyType>(weights_penalty)));
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
        if (opt_.mode == Mode::Regularized) { set_tau_auto_all_(); }
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
    void set_lambda_weights_all(const double lambda) const {
        for (auto& b : blocks_) b->set_lambda_weights(lambda);
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

    // Deflation
    void deflate_all() const {
        for (auto& b : blocks_) b->deflate(opt_.deflation_mode);
    }

    // Fit
    std::vector<Result> fit() {
        if (!initialized_) {
            // user didn't call init -> assume fully connected (off-diagonal true)
            init(DesignMode::FullyConnected);
        }

        // room for results
        std::vector<Result> results;
        results.reserve(n_comp());

        // components loop
        for (int hh = 0; hh < n_comp(); ++hh) {
            set_h(hh);
            init_comp();
            results.push_back(fit_component());

            if (hh + 1 < n_comp()) {
                deflate_all();
            }
        }

        // post-processing weights
        for (auto& b : blocks_) {
            b->compute_weights_star();
        }

        return results;
    }
    Result fit_component() {

        const int J = n_blocks();
        if (J < 2) throw std::runtime_error("RGCCA: need ≥ 2 blocks");

        // room for results
        Result res(n_blocks());
        res.obj_history.reserve(opt_.max_iter);

        // make a component local copy of the connection matrix C
        res.C = C_;

        // weights initialization
        for (int j = 0; j < J; ++j) {

            auto& b = blocks_[j];
            b->set_h(h_);

            if (opt_.init == Init::SVD) {
                const auto info = b->svd_init();
                res.active_blocks[j] = true;
                b->compute(info.nu);
            }

            if (opt_.init == Init::Random) {
                std::mt19937_64 rng(opt_.seed);
                std::uniform_real_distribution<double> U(-1.0, 1.0);
                const Vector nu = Vector::NullaryExpr(n_, [&]{ return U(rng); });
                b->compute(nu);
            }
        }

        // room for objective function evaluations
        res.obj_history.push_back(objective_(res.C));
        auto a_prev = snapshot_weights_();

        // require lambda selection also at the first iteration
        if (opt_.lambda_selection == LambdaSelection::Automatic) { set_lambda_auto_all_(); }

        for (int s = 0; s < opt_.max_iter; ++s) {
            for (int l = 0; l < J; ++l) {
                Vector nu_l = Vector::Zero(blocks_[l]->n());
                const Vector eta_l = eta_(*blocks_[l]);
                for (int k = 0; k < J; ++k) {
                    if (!res.C(l,k)) continue;
                    const Vector eta_k = eta_(*blocks_[k]);
                    const double cov_lk = cov_value_(l, k, eta_l, eta_k);
                    const double w_lk = opt_.scheme.w(cov_lk);
                    nu_l.noalias() += w_lk * eta_(*blocks_[k], *blocks_[l]);
                }
                blocks_[l]->compute(nu_l);   // block handles normalization
                mark_cov_rowcol_dirty_(l);   // η_l changed → invalidate its row/col
            }

            const double f_obj = objective_(res.C);
            const double obj_prev = res.obj_history.back();
            res.obj_history.push_back(f_obj);
            res.iters = s + 1;

            // check monotonicity
            if (f_obj + 1e-15 < obj_prev) res.monotone = false;

            // check convergence
            const double delta_obj = std::abs(f_obj - obj_prev);
            const double delta_a  = weights_variation_(a_prev);
            if (delta_obj < opt_.tol || delta_a < opt_.tol) break;

            // update snapshot for next iter
            a_prev = snapshot_weights_();
        }

        // save information about the iteration in the result struct
        res.noise_variance = noise_variance();
        compute_covariance_matrix_(res.covariance_matrix);
        get_tau(res.tau_values);
        get_lambdas(res.lambda_components_values, res.lambda_weights_values);
        get_reconstruction_constraint_info(res.reconstruction_error, res.reconstruction_edge);

        return res;
    }

    // ===== Accessors =====
    [[nodiscard]] int n() const { return n_; }
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

    std::vector<typename Block::Vector> snapshot_weights_() const {
        const int J = n_blocks();
        std::vector<typename Block::Vector> out;
        out.reserve(J);
        for (int j = 0; j < J; ++j) out.push_back(blocks_[j]->weights().col(h_));
        return out;
    }

    double weights_variation_(const std::vector<typename Block::Vector>& a_prev) const {
        const int J = n_blocks();
        double acc = 0.0;
        for (int j = 0; j < J; ++j) {
            // skip inactive blocks if you want exact R behavior after deactivation
            const auto aj = blocks_[j]->weights().col(h_);
            const auto dj = aj - a_prev[j];
            acc += dj.squaredNorm();
        }
        return acc;
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
            Psi_T_.resize(n(), n());
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

    void get_lambdas(std::vector<double> & lambda_components_values, std::vector<double> & lambda_weights_values) const {
        const int J = n_blocks();
        for (std::size_t j = 0; j < J; ++j) {
            lambda_components_values[j] = blocks_[j]->lambda_components();
            lambda_weights_values[j] = blocks_[j]->lambda_weights();
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

    void set_tau_auto_all_() const { for (auto& b : blocks_) b->select_tau_auto(); }
    void set_lambda_auto_all_() const { set_lambda_components_all(-1); }

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
    double objective_(const BoolMatrix& C) {
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
    int n_ {0}; // global number of observations
    std::vector<double> times_;
    SamplingDomain T_; // only used by TimeDependentSampling
    SparseMatrix Psi_T_;
    int h_ {0};   // current component index
    Options opt_;
    int n_comp_{0};

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
               << ", lambda_l = "<< r.lambda_weights_values[i] << "\n";
        }
        os << std::fixed;
        os << std::endl;
    }
    os << "n_iters   : " << r.iters << "\n";
    os << "monotone  : " << (r.monotone ? "yes" : "no") << "\n";
    os << std::endl;
    os << "objective :\n";
    double prev_obj = r.obj_history[0];
    for (size_t i = 1; i < r.obj_history.size(); ++i) {
        const double obj = r.obj_history[i];
        os << "- iter " << std::setw(3) << (i)
           << "   |   fit = " << std::setw(12) << std::setprecision(8) << std::fixed << obj
           << "   |   overall diff = " << std::setw(7) << (obj - prev_obj) << "\n";
        os << std::fixed << std::setprecision(8);
        prev_obj = obj;
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
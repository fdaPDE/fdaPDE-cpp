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
#include "fdaPDE/execution.h"
#include "header_check.h"

namespace fdapde {

enum class InitStrategy { None, SVD, Uniform, WarmStart };
enum class DesignMode {Empty, FullyConnected};
enum class LambdaSelection {Manual, Automatic};
enum class Mode { CorMax, Regularized, CovMax };
enum class Deflation { None, Scores };
enum class WeightSignConstraint { None, NonNegative };
enum class ResamplingStrategy {Ordinary};

namespace internals {

void ginv(const Eigen::MatrixXd& X, Eigen::MatrixXd& ginvX, double tol = std::sqrt(std::numeric_limits<double>::epsilon())){
    // SVD
    Eigen::BDCSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const auto& d = svd.singularValues();
    const double d1 = d(0);
    const double thresh = std::max(tol * d1, 0.0);

    // identify Positive singular values
    std::vector<int> idx;
    idx.reserve(d.size());
    for (int i = 0; i < d.size(); ++i) {
        if (d(i) > thresh) idx.push_back(i);
    }

    // compute V_pos * diag(1/d_pos) * U_pos^T * z without forming full matrices
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

    // shapes / accessors used by GCV path (harmless no-ops here)
    [[nodiscard]] int n_dofs()  const { return n_dofs_; }
    [[nodiscard]] int n_obs()   const { return static_cast<int>(y_.size()); }
    [[nodiscard]] int n_covs()  const { return 0; } // no param covariates in identity model
    [[nodiscard]] double edf(int = 0, int = 0) const { return 0; } // hat-trace is 0 for identity

    // data flow API
    void analyze_data() {
        Psi_.resize(n_dofs_, n_dofs_);
        Psi_.setIdentity();
    }

    void update_response_and_weights(const vector_t& y, const sparse_matrix_t& /*W*/) { y_ = y; }

    // identity: fitted values equal z
    void fit(double /*lambda*/) { f_ = y_; }

    // outputs
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
    double log10_min = -9.0;
    double log10_max = 0.0;
    int grid = 20;

    // edf() stochastic trace settings (if your solver uses Hutch++ etc.)
    int edf_r = 100;
    int edf_seed = 12345;

    // safety
    double eps_dof = 1e-12;  // avoid divide-by-zero in denominator
};
template <class Smoother> struct GCVEval {
    Smoother* s; // must expose: fit(λ), edf(r,seed), response(), fn(), n_obs(), n_covs()
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
    using IndexVector = Eigen::Vector<int, Eigen::Dynamic>;
    using SparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
    using SparseSolver = eigen_sparse_solver_movable_wrap<Eigen::SimplicialLDLT<SparseMatrix>>;
    using ComponentsSolverType = typename std::decay_t<SamplingStrategy>::solver_t;

    template<typename S = SamplingStrategy>
    requires std::same_as<S, IndependentSampling>
    BaseBlock(const std::string& block_name, Matrix* data_ptr, const int n_dofs_weights) :
        block_name_(block_name), data_ptr_(data_ptr), components_solver_(data_ptr->rows()), n_dofs_weights_(n_dofs_weights) {
        // init components solver
        init_identity_row_index_();
        components_solver_.analyze_data();
    }

    template<typename S = SamplingStrategy>
    requires std::same_as<SamplingStrategy, TimeDependentSampling>
    BaseBlock(const std::string& block_name, const Triangulation<1, 1>& T, const Vector& times, Matrix* data_ptr, const int n_dofs_weights) :
        block_name_(block_name), times_(times), data_ptr_(data_ptr), n_dofs_weights_(n_dofs_weights) {
        // init sparse identity
        init_identity_row_index_();
        I_.resize(n(), n());
        I_.setIdentity();
        // init components solver
        SamplingStrategy::discretize(T, components_solver_);
        components_solver_.analyze_data(Matrix{times}, Vector::Zero(n()), I_);
    }

    BaseBlock(const BaseBlock& other) :
        n_dofs_weights_(other.n_dofs_weights_),
        n_comp_(other.n_comp_),
        times_(other.times_),
        block_name_(other.block_name_),
        h_(other.h_),
        data_ptr_(other.data_ptr_),
        row_index_(other.row_index_),
        components_solver_(other.components_solver_),
        tau_(other.tau_),
        mode_(other.mode_),
        weight_sign_constraint_(other.weight_sign_constraint_),
        bias_(other.bias_),
        components_gcv_cfg_(other.components_gcv_cfg_),
        lambda_components_(other.lambda_components_),
        weights_(other.weights_),
        weights_star_(other.weights_star_),
        components_(other.components_),
        deflation_projections_(other.deflation_projections_),
        I_(other.I_),
        lambda_components_selection_(other.lambda_components_selection_),
        weights_ready_(other.weights_ready_),
        components_ready_(other.components_ready_) {

        M_ready_ = false;
        invM_ready_ = false;
        ginvM_ready_ = false;
        raw_data_mutable_ = false;

        nn_solver_.reset();
    }

    virtual ~BaseBlock() = default;

    // block initialization
    void init() {
        ensure_M_();
        ensure_lc_();
    }

    // data
    [[nodiscard]] const std::string& name() const { return block_name_; }
    [[nodiscard]] const Matrix& raw_data() const {
        if (!data_ptr_) throw std::logic_error("BaseBlock: Block has no data pointer");
        return *data_ptr_;
    }
    [[nodiscard]] auto data() const {
        return raw_data()(row_index_, Eigen::all);
    }
    void set_raw_data_mutable(const bool value) { raw_data_mutable_ = value; }

    // dimensions
    [[nodiscard]] int n() const { return static_cast<int>(row_index_.size()); }
    [[nodiscard]] int n_raw() const { return static_cast<int>(raw_data().rows()); }
    [[nodiscard]] int m() const { return static_cast<int>(raw_data().cols()); }
    [[nodiscard]] int n_dofs_weights() const { return n_dofs_weights_; }

    // components
    void set_n_comp(const int n_comp) {
        if (n_comp <= 0)
            throw std::invalid_argument("n_comp must be > 0");

        n_comp_ = n_comp;

        // force resize on next access
        weights_ready_ = false;
        components_ready_ = false;
    }
    [[nodiscard]] int n_comp() const { return n_comp_; }
    void set_h(const int idx) {
        if (idx < 0 || idx >= n_comp_)
            throw std::out_of_range("h");

        h_ = idx;
    }
    [[nodiscard]] int h() const { return h_; }
    void next_component() { set_h(h_ + 1); }

    // mode and shrinkage parameter
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

    // normalization matrix
    [[nodiscard]] const SparseMatrix& M() const { ensure_M_(); return M_; }
    [[nodiscard]] const Matrix& ginvM() const { ensure_ginvM_(); return ginvM_; }
    [[nodiscard]] SparseSolver& invM() { ensure_invM_(); return invM_; }

    // components regularization utilities
    void set_lambda_components(const double lambda) {
        lambda_components_ = lambda;
        lambda_components_selection_ = lambda < 0.0;
    }
    [[nodiscard]] double lambda_components() const {
        if constexpr (std::same_as<SamplingStrategy, IndependentSampling>) return std::numeric_limits<double>::quiet_NaN();
        if (lambda_components_.has_value() && *lambda_components_ > 0.0) return *lambda_components_;
        return std::numeric_limits<double>::quiet_NaN();
    }
    void set_components_gcv_config(const GCVConfig& cfg) {
        components_gcv_cfg_ = cfg;
    }

    // weights regularization utilities
    virtual void set_lambda_weights(const double) {};
    [[nodiscard]] virtual double lambda_weights() const {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // weight and component initialization API
    void init_weight_uniform() {
        ensure_lc_();
        Vector a = Vector::Ones(n_dofs_weights_);
        double norm2 = a.dot(Omega() * a);
        if (norm2 <= 0) norm2 = 1.0;
        weights_.col(h()) = a / std::sqrt(norm2);
    }
    void refresh_component() {
        ensure_lc_();
        const Vector s = data() * Psi_D() * weights_.col(h());
        components_.col(h()) = c_fit_(s);
    }

    // inner-component initialization API
    struct InitInfo {
        bool active{false}; Vector nu;
    };
    InitInfo uniform_init() const {
        InitInfo out{true, Vector::Ones(n()) };
        const Matrix& X = data();
        out.nu = X * Psi_D() * Vector::Ones(n_dofs_weights());

        return out;
    }
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

    // bootstrap API
    void set_row_index(const IndexVector& idx) {
        if (idx.size() == 0)
            throw std::invalid_argument("row index cannot be empty");

        for (int i = 0; i < idx.size(); ++i) {
            if (idx(i) < 0 || idx(i) >= n_raw())
                throw std::out_of_range("invalid row index");
        }

        row_index_ = idx;

        if constexpr (std::same_as<SamplingStrategy, TimeDependentSampling>) {
            I_.resize(n(), n());
            I_.setIdentity();
        }

        invalidate_M_();
        components_ready_ = false;
    }
    void clear_row_index() {
        IndexVector idx(n_raw());
        std::iota(idx.data(), idx.data() + idx.size(), 0);
        set_row_index(idx);
    }

    // main compute method
    void compute(const Vector& nu_D) {

        // compute the weight
        const Vector a = w_fit_(nu_D);
        weights().col(h()) = a;

        // compute the component and regularize it (the regularization acts only in the TimeDependent sampling scenario)
        const Vector s = data() * Psi_D() * a;
        components().col(h()) = c_fit_(s);
    }

    // deflation
    void deflate(const Deflation mode) {
        if (h() == n_comp()) return;
        switch (mode) {
            case Deflation::Scores: deflate_scores_(); break;
            case Deflation::None: default: break;
        }
        invalidate_M_();
    }

    // weights post-processing
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

    // model evaluation
    [[nodiscard]] Vector normalized_component_for_evaluation(const Vector& a, const int hh) {
        if (hh != h())
            throw std::logic_error("normalized_component_for_evaluation: requested component differs from current h; M may refer to current deflated data.");

        // weight at locations
        const Vector am = Psi_D() * a;

        // normalization
        // M(), not Omega()! in the evaluation, the regularization term should not be taken into account
        double nrm2 = am.dot(M() * am);
        if (nrm2 <= 0.0 || !std::isfinite(nrm2)) nrm2 = 1.0;

        // components
        const Vector s = data() * am / std::sqrt(nrm2);
        return c_fit_(s);
    }

    // virtual printer
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

    // setters
    void set_bias(const bool bias) { bias_ = bias; }
    void set_weight_sign_constraint(const WeightSignConstraint weight_sign_constraint = WeightSignConstraint::None) {
        weight_sign_constraint_ = weight_sign_constraint;
    }

    // observers
    [[nodiscard]] const SparseMatrix& Psi_T() const { return components_solver_.Psi(); }
    Matrix& weights() { ensure_lc_(); return weights_; }
    Matrix weights_m() { ensure_lc_(); return Psi_D() * weights_; }
    Matrix& weights_star() { ensure_lc_(); return weights_star_; }
    Matrix weights_star_m() { ensure_lc_(); return Psi_D() * weights_star_; }
    Matrix& components() { ensure_lc_(); return components_; }
    Matrix components_m() { ensure_lc_(); return Psi_T() * components_; }
    [[nodiscard]] WeightSignConstraint weight_sign_constraint() const { return weight_sign_constraint_; }
    template <typename S = SamplingStrategy> requires std::same_as<S, TimeDependentSampling> const Vector& times() { return times_; }

    // abstract methods
    [[nodiscard]] virtual const SparseMatrix& Psi_D() const = 0;
    [[nodiscard]] virtual const SparseMatrix& Omega() = 0; // M + regularization when present
    virtual std::unique_ptr<BaseBlock> clone() const = 0;

protected:

    // data utils
    [[nodiscard]] Matrix& mutable_raw_data_() {
        if (!data_ptr_) throw std::logic_error("BaseBlock: Block has no data pointer");
        if (!raw_data_mutable_)
            throw std::logic_error("BaseBlock: This block is not allowed to mutate raw data");
        return *data_ptr_;
    }
    void init_identity_row_index_() {
        row_index_.resize(raw_data().rows());
        std::iota(row_index_.data(), row_index_.data() + row_index_.size(), 0);
    }
    [[nodiscard]] bool is_identity_row_index_() const {
        if (row_index_.size() != n_raw()) return false;
        for (int i = 0; i < row_index_.size(); ++i)
            if (row_index_(i) != i) return false;
        return true;
    }

    // weights and components
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
    [[nodiscard]] Vector a_() { ensure_lc_(); return weights_.col(h());}
    [[nodiscard]] Vector eta_() { ensure_lc_(); return components_.col(h());}

    // tau estimate using Schäfer–Strimmer analytic shrinkage from correlation
    void select_tau_auto_() {
        const int n_obs  = n();
        const int n_vars = m();

        if (n_obs < 2 || n_vars < 1) throw std::runtime_error("tau_auto: need n >= 2 and m >= 1");

        Eigen::RowVectorXd mu = data().colwise().mean();
        Matrix Xc = data().rowwise() - mu;
        Eigen::RowVectorXd var = (Xc.array().square().colwise().sum() / static_cast<double>(n_obs - 1)).matrix();
        Eigen::RowVectorXd sd = var.array().sqrt().matrix();

        for (int j = 0; j < n_vars; ++j) {
            if (!(sd[j] > 0.0) || !std::isfinite(sd[j])) sd[j] = 1.0;
        }

        Xc.array().rowwise() /= sd.array();
        const Matrix XtX = Xc.transpose() * Xc;
        const Matrix Xc2 = Xc.array().square().matrix();
        const Matrix Xc2T_Xc2 = Xc2.transpose() * Xc2;

        const double n_d = static_cast<double>(n_obs);
        const double c = n_d / std::pow(n_d - 1.0, 3.0);

        Matrix V = c * (Xc2T_Xc2 - (1.0 / n_d) * XtX.array().square().matrix());
        V.diagonal().setZero();

        const double num = V.sum();

        Matrix Corm = XtX / (n_d - 1.0);
        Matrix D = Corm;
        D.diagonal().array() -= 1.0;

        const double den = D.squaredNorm();

        tau_ = (den > 0.0) ? std::clamp(num / den, 0.0, 1.0) : 0.0;
        invalidate_M_();
    }

    // normalization matrix
    void compute_M_() {
        M_.resize(m(), m());

        if (mode_ == Mode::CovMax) {
            M_.setIdentity();
        } else {
            const double n_d = static_cast<double>(n());
            const double den = bias_ ? n_d : std::max(1.0, n_d - 1.0);

            const Matrix XtX = data().transpose() * data();
            const Vector mu = data().colwise().mean();

            const Matrix Sigma = (XtX - n_d * (mu * mu.transpose())) / den;

            if (mode_ == Mode::CorMax) {
                M_ = Sigma.sparseView();
            } else {
                SparseMatrix I(m(), m());
                I.setIdentity();

                M_ = tau_ * I;
                M_ += ((1.0 - tau_) * Sigma).sparseView();
            }
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

    // deflation
    void deflate_scores_() {
        ensure_lc_();

        if (!is_identity_row_index_()) throw std::logic_error("deflate_scores_: raw-data deflation requires identity row_index");

        const Vector y = Psi_T() * eta_();
        const double yy = y.squaredNorm();

        if (yy > 0.0) {
            Vector p = data().transpose() * y / yy;

            // Store p_h for post-processing
            deflation_projections_.col(h()) = p;

            mutable_raw_data_().noalias() -= y * p.transpose();
        }

        invalidate_M_();
    }

    // non-negative solver utils
    Vector solve_nonnegative_weight_ipopt_(const Vector& z) {
        if (!nn_solver_)
            nn_solver_ = std::make_unique<NonNegativeWeightSolver>(
                Psi_D(),
                Omega()
            );
        return nn_solver_->solve(z);
    }
    void reset_nonnegative_weight_solver_() {
        nn_solver_.reset();
    }

    // weights solver
    virtual Vector w_fit_(const Vector& nu) = 0;
    [[nodiscard]] Vector normalize_weight_(const Vector& a) {
        const double rho2 = a.dot(Omega() * a);
        if (rho2 <= 0.0 || !std::isfinite(rho2)) return a;
        return a / std::sqrt(rho2);
    }

    // component solver
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

    // dimensions
    int n_dofs_weights_ {0};
    int n_comp_ {1};

    // state
    Vector times_{0};
    const std::string block_name_;
    int h_ {0};

    // data
    Matrix* data_ptr_ = nullptr;
    IndexVector row_index_;

    // solvers
    ComponentsSolverType components_solver_;
    std::unique_ptr<NonNegativeWeightSolver> nn_solver_;

    // options
    double tau_ {0.0};
    Mode mode_ = Mode::CorMax;
    WeightSignConstraint weight_sign_constraint_ = WeightSignConstraint::None;
    bool bias_ = true;

    // parameters
    GCVConfig components_gcv_cfg_;
    std::optional<double> lambda_components_;

    // results
    Matrix weights_, weights_star_, components_;
    Matrix deflation_projections_;

    // utilities
    SparseMatrix I_;
    SparseMatrix M_;
    Matrix ginvM_;
    SparseSolver invM_;

    // flags
    bool raw_data_mutable_ = false;
    bool weights_ready_ {false}, components_ready_ {false};
    bool M_ready_ {false}, invM_ready_ {false}, ginvM_ready_ {false};
    bool lambda_components_selection_ {false};
};

// single non-member operator<< visible to all derived classes
template <typename SamplingStrategy>
inline std::ostream& operator<<(std::ostream& os, const BaseBlock<SamplingStrategy>& b) {
    b.print(os);   // virtual dispatch -> works for Multivariate/Functional too
    return os;
}

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
    using Base::h;
    using Base::weights;
    using Base::components;
    using Base::weight_sign_constraint;

    template<typename S = SamplingStrategy>
    requires std::same_as<S, IndependentSampling>
    MultivariateBlock(const std::string& block_name, Matrix* data_ptr) :
        Base(block_name, data_ptr, static_cast<int>(data_ptr->cols())) {
        init_multivariate_();
    }

    template<typename S = SamplingStrategy>
    requires std::same_as<SamplingStrategy, TimeDependentSampling>
    MultivariateBlock(const std::string& block_name, const Triangulation<1, 1>& T, const Vector& times, Matrix* data_ptr) :
        Base(block_name, T, times, data_ptr, static_cast<int>(data_ptr->cols())) {
        init_multivariate_();
    }

    MultivariateBlock(const MultivariateBlock& other) : Base(other), Psi_D_(other.Psi_D_) { }

    std::unique_ptr<Base> clone() const override {
        return std::make_unique<MultivariateBlock>(*this);
    }

    // Omega
    [[nodiscard]] const SparseMatrix& Omega() override { return M(); };

    // Psi_D
    [[nodiscard]] const SparseMatrix& Psi_D() const override { return Psi_D_; }

    // print
    void print(std::ostream& os) const override {
        Base::print(os);
        os << "type: MultivariateBlock, n_dofs_weights = m = " << n_dofs_weights();
        os << "\n";
    }

protected:
    using Base::solve_nonnegative_weight_ipopt_;
    using Base::reset_nonnegative_weight_solver_;
    using Base::normalize_weight_;

    void init_multivariate_() {
        Psi_D_.resize(m(), n_dofs_weights()); // n_dofs_weights == m in this case
        Psi_D_.setIdentity();
        init();
    }

    Vector w_fit_(const Vector& nu) override {
        assert(nu.size() == n() && "nu must have size n (rows of X)");
        init();

        Vector z = data().transpose() * nu;

        if (weight_sign_constraint() == WeightSignConstraint::NonNegative) {
            return solve_nonnegative_weight_ipopt_(z); // already normalized
        }

        const Vector a_tilde = ginvM() * z;
        return normalize_weight_(a_tilde);
    }

    void invalidate_derived_caches_() override {
        reset_nonnegative_weight_solver_();
    }

private:
    SparseMatrix Psi_D_; // m x m sparse identity matrix
};

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
    using Base::n;
    using Base::m;
    using Base::n_dofs_weights;
    using Base::data;
    using Base::h;
    using Base::weights;
    using Base::components;
    using Base::weight_sign_constraint;

    template <typename GeoFrame>
    requires std::same_as<SamplingStrategy, IndependentSampling>
    FunctionalBlock(const std::string& block_name, GeoFrame& gf, Matrix* data_ptr, WeightsPenaltyType&& weights_penalty) :
        Base(block_name, data_ptr, weights_penalty.get().bilinear_form().n_dofs()) {
        init_functional_(gf, std::forward<WeightsPenaltyType>(weights_penalty));
    }

    template <typename GeoFrame>
    requires std::same_as<SamplingStrategy, TimeDependentSampling>
    FunctionalBlock(const std::string& block_name, const Triangulation<1, 1>& T, const Vector& times, GeoFrame& gf, Matrix* data_ptr, WeightsPenaltyType&& weights_penalty) :
        Base(block_name, T, times, data_ptr, weights_penalty.get().bilinear_form().n_dofs()) {
        init_functional_(gf, std::forward<WeightsPenaltyType>(weights_penalty));
    }

    FunctionalBlock(const FunctionalBlock& other) :
        Base(other), weights_solver_(other.weights_solver_), lambda_weights_(other.lambda_weights_) {
        Omega_ready_ = false;
    }

    std::unique_ptr<Base> clone() const override {
        return std::make_unique<FunctionalBlock>(*this);
    }

    // Psi_D
    [[nodiscard]] const SparseMatrix& Psi_D() const override { return weights_solver_.Psi(); }

    // weights regularization utils
    void set_lambda_weights(const double lambda) override {
        lambda_weights_ = lambda;
        Omega_ready_ = false;
        reset_nonnegative_weight_solver_();
    }
    [[nodiscard]] double lambda_weights() const override {
        return lambda_weights_;
    }

    // Omega matrix
    [[nodiscard]] const SparseMatrix& Omega() override {
        if (!Omega_ready_) {
            Omega_ = Psi_D().transpose() * M() * Psi_D();
            Omega_ += lambda_weights_ * weights_solver_.P();
            Omega_ready_ = true;
        }
        return Omega_;
    }

    // print
    void print(std::ostream& os) const override {
        Base::print(os);
        os << "type: FunctionalBlock, n_dofs_weights = " << n_dofs_weights();
        os << ", lambda = " << lambda_weights_;
        os << "\n";
    }

protected:
    using Base::solve_nonnegative_weight_ipopt_;
    using Base::reset_nonnegative_weight_solver_;

    template <typename GeoFrame>
    void init_functional_(GeoFrame& gf, WeightsPenaltyType&& weights_penalty) {
        weights_solver_.discretize(weights_penalty.get());
        weights_solver_.analyze_data(gf, M());
        init();
    }

    Vector w_fit_(const Vector& nu) override {
        assert(nu.size() == n() && "nu must have size n (rows of X)");
        init();

        Vector z = data().transpose() * nu;

        if (weight_sign_constraint() == WeightSignConstraint::NonNegative) {
            return solve_nonnegative_weight_ipopt_(z);  // already normalized
        }

        weights_solver_.update_z_and_weights(z, M()); // can be optimized by passing M only when it has changed
        weights_solver_.fit(lambda_weights_);
        return weights_solver_.f(); // already normalized
    }
    void invalidate_derived_caches_() override {
        Omega_ready_ = false;
        reset_nonnegative_weight_solver_();
    }

private:
    SparseMatrix Omega_;
    bool Omega_ready_ {false};
    WeightsSolverType weights_solver_;
    double lambda_weights_ = 1e-15;
};

template <typename SamplingStrategy, typename Matrix = Eigen::Matrix<double, Dynamic, Dynamic>>
requires std::same_as<SamplingStrategy, IndependentSampling>
inline std::unique_ptr<internals::BaseBlock<SamplingStrategy>>
make_multivariate_block(std::string block_name, Matrix* data_ptr) {
    return std::make_unique<internals::MultivariateBlock<SamplingStrategy>>(block_name, data_ptr);
}

template <typename SamplingStrategy, typename Matrix = Eigen::Matrix<double, Dynamic, Dynamic>, typename Vector = Eigen::Matrix<double, Dynamic, 1>>
requires std::same_as<SamplingStrategy, TimeDependentSampling>
inline std::unique_ptr<internals::BaseBlock<SamplingStrategy>>
make_multivariate_block(std::string block_name, const Triangulation<1, 1>& T, const Vector& times, Matrix* data_ptr) {
    return std::make_unique<internals::MultivariateBlock<SamplingStrategy>>(block_name, T, times, data_ptr);
}

template <typename SamplingStrategy, typename GeoFrame, typename WeightsPenaltyType, typename Matrix = Eigen::Matrix<double, Dynamic, Dynamic>>
requires std::same_as<SamplingStrategy, IndependentSampling>
std::unique_ptr<internals::BaseBlock<SamplingStrategy>>
make_functional_block(std::string block_name, GeoFrame& gf, Matrix* data_ptr, WeightsPenaltyType&& weights_penalty) {
    return std::make_unique<internals::FunctionalBlock<WeightsPenaltyType, SamplingStrategy>>(block_name, gf, data_ptr, std::forward<WeightsPenaltyType>(weights_penalty));
}

template <typename SamplingStrategy, typename GeoFrame, typename WeightsPenaltyType, typename Vector = Eigen::Matrix<double, Dynamic, 1>, typename Matrix = Eigen::Matrix<double, Dynamic, Dynamic>>
requires std::same_as<SamplingStrategy, TimeDependentSampling>
std::unique_ptr<internals::BaseBlock<SamplingStrategy>>
make_functional_block(std::string block_name, const Triangulation<1, 1>& T, const Vector& times, GeoFrame& gf, Matrix* data_ptr, WeightsPenaltyType&& weights_penalty) {
    return std::make_unique<internals::FunctionalBlock<WeightsPenaltyType, SamplingStrategy>>(block_name, T, times, gf, data_ptr,  std::forward<WeightsPenaltyType>(weights_penalty));
}

} // namespace internals

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
    using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    using BoolMatrix = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;

    int h = 0;
    int J = 0;
    std::vector<double> obj_history;
    bool monotone = true;
    int iters = 0;
    BoolMatrix C;
    Matrix covariance_matrix;
    std::vector<double> tau_values;
    std::vector<double> lambda_components_values;
    std::vector<double> lambda_weights_values;
    std::vector<bool> active_blocks;

    explicit Result(const int n_blocks) : J(n_blocks), C(J, J), covariance_matrix(J,J),
    tau_values(J), lambda_components_values(J), lambda_weights_values(J), active_blocks(J) {}
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
    using BlockList = std::vector<BlockPtr>;
    using BlockRefList = std::vector<Block*>;
    using BlockOwnerList = std::vector<std::unique_ptr<Block>>;
    using Matrix = typename Block::Matrix;
    using BoolMatrix = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMatrix = typename Block::SparseMatrix;
    using Vector = typename Block::Vector;
    using SamplingDomain = std::conditional_t<std::same_as<SamplingStrategy, TimeDependentSampling>, Triangulation<1, 1>, internals::empty_t>;

    struct Options {
        int max_iter;
        double tol;
        double active_block_tol;
        bool verbose;
        bool cache_covariances;
        bool bias;
        InitStrategy init_strategy;
        LambdaSelection lambda_selection_weights;
        LambdaSelection lambda_selection_components;
        Mode mode;
        WeightSignConstraint weight_sign_constraint;
        Deflation deflation_mode;
        Scheme scheme;

        explicit Options(
          const int max_iter_ = 1000, const double tol_ = 1e-8, const double active_block_tol_ = 1e-8, const bool bias_ = true,
          const InitStrategy init_strategy_ = InitStrategy::SVD, const Mode mode_ = Mode::CovMax,
          const WeightSignConstraint weight_sign_constraint_ = WeightSignConstraint::None,
          const LambdaSelection lambda_selection_weights_ = LambdaSelection::Manual,
          const LambdaSelection lambda_selection_components_ = LambdaSelection::Automatic,
          const Deflation deflation_mode_ = Deflation::Scores, const Scheme& scheme_ = Scheme::Factorial(),
          const bool verbose_ = false, const bool cache_ = true) :
            max_iter(max_iter_),
            tol(tol_),
            active_block_tol(active_block_tol_),
            bias(bias_),
            init_strategy(init_strategy_),
            mode(mode_),
            weight_sign_constraint(weight_sign_constraint_),
            lambda_selection_weights(lambda_selection_weights_),
            lambda_selection_components(lambda_selection_components_),
            deflation_mode(deflation_mode_),
            scheme(scheme_),
            verbose(verbose_),
            cache_covariances(cache_) { }
    };
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
    struct BootstrapConfig {
        int B = 1000;
        unsigned seed = 12345;
        ResamplingStrategy resampling_strategy = ResamplingStrategy::Ordinary;
    };
    struct BootstrapSelectionResult {
        using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

        int h = 0;
        int B = 0;

        std::vector<double> lambda_grid;
        std::vector<double> criterion;

        double lambda_opt = std::numeric_limits<double>::quiet_NaN();

        std::vector<std::string> block_names;
        std::vector<bool> active_blocks;

        // [lambda][block] -> vector/matrix
        std::vector<std::vector<Vector>> w_fit_by_lambda;
        std::vector<std::vector<Matrix>> w_boot_by_lambda;
        std::vector<std::vector<Vector>> w_min_by_lambda;

        BootstrapSelectionResult() = default;

        BootstrapSelectionResult(
            const int h_,
            const int B_,
            const std::vector<double>& lambda_grid_,
            const std::vector<std::string>& block_names_,
            const std::vector<int>& block_dims_
        ) :
            h(h_),
            B(B_),
            lambda_grid(lambda_grid_),
            criterion(lambda_grid_.size(), -std::numeric_limits<double>::infinity()),
            block_names(block_names_)
        {
            const std::size_t n_lambda = lambda_grid.size();
            const std::size_t J = block_dims_.size();

            w_fit_by_lambda.resize(n_lambda);
            w_boot_by_lambda.resize(n_lambda);
            w_min_by_lambda.resize(n_lambda);

            active_blocks.resize(J);

            for (std::size_t i = 0; i < n_lambda; ++i) {
                w_boot_by_lambda[i].resize(J);

                for (std::size_t j = 0; j < J; ++j) {
                    w_boot_by_lambda[i][j].setZero(block_dims_[j], B);
                }
            }
        }
    };

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
        return add_block(internals::make_multivariate_block<SamplingStrategy>(block_name, data_blocks_.back().get()));
    }
    template<typename S = SamplingStrategy>
    requires std::same_as<S, TimeDependentSampling>
    int add_multivariate_block(std::string block_name, const Vector& times, Matrix&& X) {
        data_blocks_.push_back(std::make_unique<Matrix>(std::move(X)));
        return add_block(internals::make_multivariate_block<SamplingStrategy>(block_name, T_, times, data_blocks_.back().get()));
    }
    template <typename GeoFrame, typename WeightsPenaltyType>
    requires std::same_as<SamplingStrategy, IndependentSampling>
    int add_functional_block(std::string block_name, const GeoFrame& gf, Matrix&& X, WeightsPenaltyType&& weights_penalty) {
        data_blocks_.push_back(std::make_unique<Matrix>(std::move(X)));
        return add_block(internals::make_functional_block<SamplingStrategy>(block_name, gf, data_blocks_.back().get(), std::forward<WeightsPenaltyType>(weights_penalty)));
    }
    template <typename GeoFrame, typename WeightsPenaltyType>
    requires std::same_as<SamplingStrategy, TimeDependentSampling>
    int add_functional_block(std::string block_name, const Vector& times, const GeoFrame& gf, Matrix&& X, WeightsPenaltyType&& weights_penalty) {
        data_blocks_.push_back(std::make_unique<Matrix>(std::move(X)));
        return add_block(internals::make_functional_block<SamplingStrategy>(block_name, T_, times, gf, data_blocks_.back().get(), std::forward<WeightsPenaltyType>(weights_penalty)));
    }
    void connect(int j, int k, bool on = true) {
        ensure_design_initialized_();

        check_index_(j);
        check_index_(k);

        if (j == k) {
            std::cerr << "RGCCA::connect(): ignoring self-connection for block " << j << '\n';
            return;
        }

        C_(j, k) = on;
        C_(k, j) = on;
    }

    // initialization
    void init(const DesignMode mode = DesignMode::FullyConnected) {
        if (n_blocks() < 2) throw std::runtime_error("RGCCA: need ≥ 2 blocks");

        // design matrix initialization
        switch (mode) {
            case DesignMode::Empty: clear_design_(); break;
            case DesignMode::FullyConnected: set_fully_connected_design_(); break;
        }

        compute_Psi_();

        initialized_ = true;
    }

    // weights and components regularization utilities
    void set_lambda_weights_all(const double lambda) const {
        for (auto& b : blocks_) b->set_lambda_weights(lambda);
    }
    void set_lambda_components_all(const double lambda) const {
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
    void set_n_bootstrap_samples(const int n_bootstrap_samples) {
        bootstrap_config_.B = n_bootstrap_samples;
    }

    // fit
    std::vector<Result> fit() {
        if (!initialized_) init();

        const int J = n_blocks();
        if (J < 2) throw std::runtime_error("RGCCA: need ≥ 2 blocks");

        // room for results
        std::vector<Result> results;
        results.reserve(n_comp());
        bootstrap_selection_results_.clear();
        bootstrap_selection_results_.reserve(n_comp());

        // components loop
        for (int hh = 0; hh < n_comp(); ++hh) {
            set_h_(hh);

            // weights lambda selection
            std::vector<bool> active_blocks(n_blocks(), true);
            if (opt_.lambda_selection_weights == LambdaSelection::Automatic) {
                auto selection = select_lambda_weights_bootstrap_parallel_();
                const double lambda = selection.first;
                active_blocks = std::move(selection.second);
                set_lambda_weights_all(lambda);
            }

            // final fit
            init_comp_();
            results.push_back(fit_component_(active_blocks));

            deflate_all_();
        }

        // post-processing weights
        compute_weights_star_();

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
    [[nodiscard]] const SparseMatrix& Psi_T() const { return Psi_T_; };
    [[nodiscard]] std::vector<BootstrapSelectionResult> bootstrap_selection_results() const { return bootstrap_selection_results_; }

private:

    // initialization utils
    void check_index_(int j) const {
        if (j < 0 || j >= static_cast<int>(blocks_.size())) throw std::out_of_range("block index");
    }
    void initialize_design_(const DesignMode mode) {
        const int J = n_blocks();

        C_.resize(J, J);
        C_.setConstant(false);

        if (mode == DesignMode::FullyConnected) {
            for (int j = 0; j < J; ++j)
                for (int k = 0; k < J; ++k)
                    C_(j, k) = (j != k);
        }
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
    BootstrapBlocks clone_blocks_from_(const BootstrapBlocks& src) const {
        BootstrapBlocks out;
        out.owners.reserve(src.refs.size());
        out.refs.reserve(src.refs.size());

        for (auto* b : src.refs) {
            auto copy = b->clone();
            copy->set_raw_data_mutable(false);

            out.refs.push_back(copy.get());
            out.owners.push_back(std::move(copy));
        }

        return out;
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
    void init_comp_(const BlockRefList& blocks, InitStrategy init_strategy = InitStrategy::None) {
        if (opt_.mode == Mode::Regularized) set_tau_auto_all_(blocks);
        if (opt_.lambda_selection_components == LambdaSelection::Automatic) set_lambda_components_auto_all_(blocks);

        if (init_strategy == InitStrategy::None) init_strategy = opt_.init_strategy;

        for (int j = 0; j < n_blocks(); ++j) {
            auto* b = blocks[j];
            b->set_h(h_);

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
    Result fit_component_(const BlockRefList& blocks, const std::vector<bool>& active_blocks) {
        const int J = n_blocks();
        FitWorkspace ws(J);

        // room for results
        Result res(J);
        res.obj_history.reserve(opt_.max_iter);

        // design update according to current active blocks
        res.C = C_;
        res.active_blocks = active_blocks;
        for (int j = 0; j < J; ++j) {
            if (!active_blocks[j]) {
                res.C.row(j).setConstant(false);
                res.C.col(j).setConstant(false);
                blocks[j]->weights().col(h_).setZero();
                blocks[j]->components().col(h_).setZero();
            }
        }

        // initialization
        res.obj_history.push_back(objective_(blocks, ws, res.C));
        auto a_prev = snapshot_weights_(blocks);
        if (opt_.lambda_selection_components == LambdaSelection::Automatic)
            set_lambda_components_auto_all_(blocks);

        // main loop
        for (int s = 0; s < opt_.max_iter; ++s) {
            for (int l = 0; l < J; ++l) {

                // skip deactivated blocks
                if (!active_blocks[l]) continue;

                // inner-component assembler
                Vector nu_l = Vector::Zero(blocks[l]->n());
                const Vector eta_l = eta_(*blocks[l]);
                for (int k = 0; k < J; ++k) {
                    if (!res.C(l, k)) continue;
                    const Vector eta_k = eta_(*blocks[k]);
                    const double cov_lk = cov_value_(ws, l, k, eta_l, eta_k);
                    const double w_lk = opt_.scheme.w(cov_lk);
                    nu_l.noalias() += w_lk * eta_(*blocks[k], *blocks[l]);
                }

                // block update
                blocks[l]->compute(nu_l);
                mark_cov_rowcol_dirty_(ws, l);
            }

            // update metrics
            const double f_obj = objective_(blocks, ws, res.C);
            const double obj_prev = res.obj_history.back();
            res.obj_history.push_back(f_obj);
            res.iters = s + 1;

            // chek monotonicity
            if (f_obj + 1e-15 < obj_prev)
                res.monotone = false;

            // stopping criteria
            const double delta_obj = std::abs(f_obj - obj_prev);
            const double delta_a = weights_variation_(blocks, a_prev);
            if (delta_obj < opt_.tol || delta_a < opt_.tol)
                break;

            a_prev = snapshot_weights_(blocks);
        }

        // save results
        covariance_matrix_(blocks, res.covariance_matrix);
        get_tau(blocks, res.tau_values);
        get_lambdas(blocks, res.lambda_components_values, res.lambda_weights_values);

        return res;
    }
    Result fit_component_(const BlockRefList& blocks) {
        std::vector<bool> active_blocks(blocks.size(), true);
        return fit_component_(blocks, active_blocks);
    }
    Result fit_component_(const std::vector<bool>& active_blocks) {
        auto blocks = main_blocks_();
        return fit_component_(blocks, active_blocks);
    }
    Result fit_component_() {
        auto blocks = main_blocks_();
        std::vector<bool> active_blocks(blocks.size(), true);
        return fit_component_(blocks, active_blocks);
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
    double weights_variation_(const BlockRefList& blocks, const std::vector<Vector>& a_prev) const {
        double acc = 0.0;

        for (int j = 0; j < n_blocks(); ++j) {
            const auto aj = blocks[j]->weights().col(h_);
            acc += (aj - a_prev[j]).squaredNorm();
        }

        return acc;
    }

    // bootstrap selectors
    std::pair<double, std::vector<bool>> select_lambda_weights_bootstrap_parallel_() {
        if (lambda_grid_weights_[h_].empty())
            throw std::runtime_error("lambda_grid_weights_ is empty");

        const int n_threads = 12;
        parallel_set_num_threads(n_threads);

        auto blocks = main_blocks_();
        const int J = static_cast<int>(blocks_.size());
        const int B = bootstrap_config_.B;

        const int patience = 2;
        int no_improve = 0;
        double best_criterion = -std::numeric_limits<double>::infinity();
        int best_i = -1;

        BootstrapSelectionResult boot_results(h_, B, lambda_grid_weights_[h_], block_names_(blocks), block_dims_(blocks));

        // same bootstrap resamples for all lambda values
        std::vector<typename Block::IndexVector> bootstrap_idx(B);
        const unsigned seed = bootstrap_config_.seed + static_cast<unsigned>(h_);
        std::mt19937_64 rng(seed);
        for (int b = 0; b < B; ++b) {
            bootstrap_idx[b] = bootstrap_indices_(n_, rng);
        }

        // preliminary fit at largest lambda
        set_lambda_weights_all(lambda_grid_weights_[h_].back());
        init_comp_(blocks);
        fit_component_(blocks);

        for (int i = static_cast<int>(lambda_grid_weights_[h_].size()) - 1; i >= 0; --i) {
            const double lambda = lambda_grid_weights_[h_][i];
            std::cout << "- lambda = " << lambda << std::endl;

            set_lambda_weights_all(lambda);

            // fit using the fit on the previous lambda as warm start
            init_comp_(blocks, InitStrategy::WarmStart);
            fit_component_(blocks);

            auto w_fit = weights_(blocks);
            auto w_min = w_fit;

            // thread-local storage
            std::vector<int> thread_count(n_threads, 0);

            std::vector<std::vector<Matrix>> thread_w_boot(n_threads);
            std::vector<std::vector<int>> thread_sample_id(n_threads);
            for (int tid = 0; tid < n_threads; ++tid) {
                thread_w_boot[tid].resize(J);
                thread_sample_id[tid].resize(B);
                for (int j = 0; j < J; ++j) {
                    thread_w_boot[tid][j].resize(w_fit[j].size(), B);
                }
            }

            std::vector<BootstrapBlocks> thread_boot_template(n_threads);
            for (int t = 0; t < n_threads; ++t) {
                thread_boot_template[t] = clone_blocks_();
            }

            std::cout << "Parallelizzazione su " <<  n_threads <<  " threads --> ";
            parallel_for(0, B,[&](int b) {
                const int tid = this_thread_id();

                auto boot_blocks = clone_blocks_from_(thread_boot_template[tid]);

                set_row_index_all_(boot_blocks.refs, bootstrap_idx[b]);

                init_comp_(boot_blocks.refs, InitStrategy::WarmStart);
                fit_component_(boot_blocks.refs);

                auto w_b = weights_(boot_blocks.refs);

                const int local_col = thread_count[tid]++;
                thread_sample_id[tid][local_col] = b;
                for (int j = 0; j < J; ++j) {
                    thread_w_boot[tid][j].col(local_col) = w_b[j];
                }

            });
            std::cout << "<--";

            // sequential merge
            for (int tid = 0; tid < n_threads; ++tid) {
                for (int local_col = 0; local_col < thread_count[tid]; ++local_col) {
                    const int b = thread_sample_id[tid][local_col];

                    for (int j = 0; j < J; ++j) {
                        Vector w_bj = thread_w_boot[tid][j].col(local_col);

                        // Align bootstrap weight sign with original fit
                        if (w_bj.dot(w_fit[j]) < 0.0) {
                            w_bj *= -1.0;
                        }

                        // Save aligned bootstrap weight
                        boot_results.w_boot_by_lambda[i][j].col(b) = w_bj;

                        // Update w_min toward zero, preserving the sign pattern of w_fit[j]
                        for (int r = 0; r < w_min[j].size(); ++r) {
                            if (w_fit[j][r] > 0.0) {
                                // take minimum, but do not go below zero
                                w_min[j][r] = std::max(0.0, std::min(w_min[j][r], w_bj[r]));
                            } else if (w_fit[j][r] < 0.0) {
                                // take maximum, but do not go above zero
                                w_min[j][r] = std::min(0.0, std::max(w_min[j][r], w_bj[r]));
                            } else {
                                w_min[j][r] = 0.0;
                            }
                        }
                    }
                }
            }

            boot_results.w_fit_by_lambda[i] = w_fit;
            boot_results.w_min_by_lambda[i] = w_min;

            const double crit = criterion_score_with_weights_(blocks, w_min, C_);
            boot_results.criterion[i] = crit;
            std::cout << " -> " << crit << std::endl;

            if (crit > best_criterion) {
                best_criterion = crit;
                best_i = i;
                no_improve = 0;
            } else {
                ++no_improve;
            }

            if (no_improve >= patience) {
                std::cout << "early stop: no improvement for " << patience << " consecutive lambdas" << std::endl;
                break;
            }
        }

        // lambda selection
        boot_results.lambda_opt = lambda_grid_weights_[h_][best_i];

        // blocks deactivation
        boot_results.active_blocks.assign(J, true);
        for (int j = 0; j < J; ++j) {
            const double nrm = boot_results.w_min_by_lambda[best_i][j].norm();
            std::cout << "block " << j << " ||w_min|| = " << nrm << " active = " << (nrm >= opt_.active_block_tol) << std::endl;
            boot_results.active_blocks[j] = nrm >= opt_.active_block_tol;
        }

        bootstrap_selection_results_.push_back(std::move(boot_results));
        return {bootstrap_selection_results_.back().lambda_opt, bootstrap_selection_results_.back().active_blocks};
    }

    // bootstrap utils
    void set_row_index_all_(const BlockRefList& blocks, const typename Block::IndexVector& idx) {
        for (auto* b : blocks) b->set_row_index(idx);
    }
    void clear_row_index_all_(const BlockRefList& blocks) {
        for (auto* b : blocks) b->clear_row_index();
    }
    typename Block::IndexVector bootstrap_indices_(int n, std::mt19937_64& rng) const {
        if (bootstrap_config_.resampling_strategy != ResamplingStrategy::Ordinary)
            throw std::runtime_error("Only ordinary bootstrap is implemented");

        std::uniform_int_distribution<int> U(0, n - 1);

        typename Block::IndexVector idx(n);
        for (int i = 0; i < n; ++i)
            idx(i) = U(rng);

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

    // private setters
    void set_tau_auto_all_(const BlockRefList& blocks) const {
        for (auto* b : blocks)
            b->select_tau_auto();
    }
    void set_lambda_components_auto_all_(const BlockRefList& blocks) const {
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
            return Psi_T() * b.components().col(h_);
        } else {
            return b.components_m().col(h_);
        }
    }
    Vector eta_(Block& b, const Block& ref) const {
        // using reference block's Psi_T
        if constexpr (std::same_as<SamplingStrategy, TimeDependentSampling>) {
            return ref.Psi_T() * b.components().col(h_);
        } else {
            return b.components_m().col(h_);
        }
    }
    Vector eta_eval_(Block& b, const Vector& a) {
        // computed wrt a and normalized wrt M (no regularization)
        const Vector eta = b.normalized_component_for_evaluation(a, h_);
        if constexpr (std::same_as<SamplingStrategy, TimeDependentSampling>) {
            return Psi_T() * eta;
        } else {
            return eta;
        }
    }

    // weights
    std::vector<Vector> weights_(const BlockRefList& blocks) const {
        std::vector<Vector> out;
        out.reserve(blocks.size());

        for (auto* b : blocks) {
            out.push_back(b->weights().col(h_));
        }

        return out;
    }

    // covariance
    double cov_(const Vector& u, const Vector& v) const {
        const double den = opt_.bias ? u.size() : std::max<int>(1, u.size() - 1);
        return (u.dot(v) - static_cast<double>(u.size()) * u.mean() * v.mean()) / den;
    }
    double cov_value_(FitWorkspace& ws, int l, int k, const Vector& eta_l, const Vector& eta_k) const {
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
    void covariance_matrix_(const BlockRefList& blocks, Matrix& Cov) const {
        const int J = n_blocks();
        for (int j = 0; j < J; ++j) {
            const Vector eta_j = eta_(*blocks[j]);
            for (int k = 0; k < J; ++k) {
                const Vector eta_k = eta_(*blocks[k]);
                Cov(j, k) = cov_(eta_j, eta_k);
            }
        }
    }

    // optimization criteria
    double objective_(const BlockRefList& blocks, FitWorkspace& ws, const BoolMatrix& C) const {
        const int J = n_blocks();
        double f = 0.0;
        for (int j = 0; j < J; ++j) {
            const Vector eta_j = eta_(*blocks[j]);
            for (int k = j; k < J; ++k) {
                if (C(j, k)) {
                    const double cov_jk = cov_value_(ws, j, k, eta_j, eta_(*blocks[k]));
                    const double mult = j == k ? 1.0 : 2.0;
                    f += mult * opt_.scheme.g(cov_jk);
                }
            }
        }
        return f;
    }
    double rho_tot_with_weights_(const BlockRefList& blocks, const std::vector<Vector>& weights, const BoolMatrix& C) const {
        const int J = n_blocks();

        if (static_cast<int>(weights.size()) != J)
            throw std::logic_error("rho_tot_with_weights_: size mismatch");

        double num = 0.0;
        double den = 0.0;

        std::vector<Vector> eta(J);

        for (int j = 0; j < J; ++j) {
            if (weights[j].size() != blocks[j]->n_dofs_weights())
                throw std::logic_error("rho_tot_with_weights_: incompatible weight size");

            eta[j] = blocks[j]->data() * blocks[j]->Psi_D() *  weights[j];
        }

        for (int j = 0; j < J; ++j) {
            const double nj = std::sqrt(eta[j].squaredNorm());

            for (int k = j + 1; k < J; ++k) {
                if (!C(j, k)) continue;

                const double nk = std::sqrt(eta[k].squaredNorm());

                if (nj > 0.0 && nk > 0.0) {
                    const double corr_jk = eta[j].dot(eta[k]) / (nj * nk);
                    num += opt_.scheme.g(corr_jk);
                    den += 1.0;
                }
            }
        }

        return den > 0.0 ? num / den : 0.0;
    }
    double criterion_score_with_weights_(const BlockRefList& blocks, const std::vector<Vector>& weights, const BoolMatrix& C) {
        const int J = n_blocks();

        if (static_cast<int>(weights.size()) != J)
            throw std::logic_error("criterion_score_with_weights_: size mismatch");

        double num = 0.0;
        double den = 0.0;

        std::vector<Vector> eta(J);

        for (int j = 0; j < J; ++j) {

            if (weights[j].size() != blocks[j]->n_dofs_weights())
                throw std::logic_error("criterion_score_with_weights_: incompatible weight size");

            eta[j] = eta_eval_(*blocks[j], weights[j]);

        }
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

    int J_ {0};
    int n_ {0}; // global number of observations

    std::vector<double> times_;
    SamplingDomain T_; // only used by TimeDependentSampling
    SparseMatrix Psi_T_;

    int h_ {0};   // current component index
    int n_comp_{1};

    std::vector<std::unique_ptr<Matrix>> data_blocks_;
    std::vector<BlockPtr> blocks_;

    BootstrapConfig bootstrap_config_;
    std::vector<BootstrapSelectionResult> bootstrap_selection_results_;
    std::vector<std::vector<double>> lambda_grid_weights_;

    bool initialized_ {false};
    BoolMatrix C_;
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
            os << "- Block " << i+1  << ": " << (r.active_blocks[i] ? "active    " : "non-active" ) << "\n";
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
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
#include "fdaPDE/src/logging.h"
#include "fdaPDE/execution.h"
#include "header_check.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <random>
#include <string_view>

namespace fdapde {

enum class InitStrategy { None, SVD, Uniform, WarmStart };
enum class DesignMode {Empty, Custom, FullyConnected};
enum class LambdaSelection {Manual, Automatic};
enum class Mode { CorMax, Regularized, CovMax };
enum class Deflation { None, Scores };
enum class WeightSignConstraint { None, NonNegative };
enum class ResamplingStrategy { Ordinary, Stationary };



inline const char* to_string(InitStrategy x) {
    switch (x) {
    case InitStrategy::None:      return "None";
    case InitStrategy::SVD:       return "SVD";
    case InitStrategy::Uniform:   return "Uniform";
    case InitStrategy::WarmStart: return "WarmStart";
    }
    return "Unknown";
}

inline const char* to_string(LambdaSelection x) {
    switch (x) {
    case LambdaSelection::Manual:    return "Manual";
    case LambdaSelection::Automatic: return "Automatic";
    }
    return "Unknown";
}

inline const char* to_string(Mode x) {
    switch (x) {
    case Mode::CorMax:      return "CorMax";
    case Mode::Regularized: return "Regularized";
    case Mode::CovMax:      return "CovMax";
    }
    return "Unknown";
}

inline const char* to_string(WeightSignConstraint x) {
    switch (x) {
    case WeightSignConstraint::None:        return "None";
    case WeightSignConstraint::NonNegative: return "NonNegative";
    }
    return "Unknown";
}

inline const char* to_string(Deflation x) {
    switch (x) {
    case Deflation::None:   return "None";
    case Deflation::Scores: return "Scores";
    }
    return "Unknown";
}

inline const char* to_string(ResamplingStrategy x) {
    switch (x) {
    case ResamplingStrategy::Ordinary:   return "Ordinary";
    case ResamplingStrategy::Stationary: return "Stationary";
    }
    return "Unknown";
}

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

inline void validate_rgcca_positive_regularization_lambda_(const double lambda, const char* name) {
    if (!(lambda > 0.0) || !std::isfinite(lambda)) {
        throw std::invalid_argument(std::string("RGCCA: ") + name + " must be finite and positive");
    }
}

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
    using BinaryMatrixT = BinaryMatrix<Dynamic, Dynamic>;
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
    bool raw_score_cache_for_weight(const Vector& a, Vector& out) const {
        if (!raw_score_cache_ready_ || !weights_ready_ || a.size() != weights_.rows())
            return false;

        const auto current_weight = weights_.col(h_);
        if (a.isApprox(current_weight)) {
            out = raw_score_cache_;
            return true;
        }
        if (a.isApprox(-current_weight)) {
            out = -raw_score_cache_;
            return true;
        }
        return false;
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
        if constexpr (std::same_as<SamplingStrategy, IndependentSampling>) return;
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
        const Vector s = data_times_(Psi_D() * weights_.col(h()));
        components_.col(h()) = c_fit_(s);
    }

    // inner-component initialization API
    struct InitInfo {
        bool active{false}; Vector nu;
    };
    InitInfo uniform_init() const {
        InitInfo out{true, Vector::Ones(n()) };
        out.nu = data_times_(Psi_D() * Vector::Ones(n_dofs_weights()));

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

        invalidate_M_after_data_change_();
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
        const Vector s = data_times_(Psi_D() * a);
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
        const Vector s = data_times_(am) / std::sqrt(nrm2);
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
    [[nodiscard]] SparseMatrix Psi_at(const Matrix& locs) const { return Psi_at_(locs); }
    [[nodiscard]] SparseMatrix Psi_at(const BinaryMatrixT& locs) const { return Psi_at_(locs); }

    // abstract methods
    [[nodiscard]] virtual const SparseMatrix& Psi_D() const = 0;
    [[nodiscard]] virtual const SparseMatrix& Omega() = 0; // M + regularization when present
    virtual std::unique_ptr<BaseBlock> clone() const = 0;
    [[nodiscard]] virtual SparseMatrix Psi_at_(const Matrix& locs) const = 0;
    [[nodiscard]] virtual SparseMatrix Psi_at_(const BinaryMatrixT& locs) const = 0;

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
    Vector data_times_(const Vector& x) const {
        raw_score_cache_.resize(n_raw());
        raw_score_cache_.noalias() = raw_data() * x;
        raw_score_cache_ready_ = true;

        if (is_identity_row_index_()) {
            return raw_score_cache_;
        }

        Vector out(n());
        for (int i = 0; i < row_index_.size(); ++i)
            out[i] = raw_score_cache_[row_index_(i)];
        return out;
    }
    Vector data_transpose_times_(const Vector& x) const {
        Vector out(m());
        if (is_identity_row_index_()) {
            out.noalias() = raw_data().transpose() * x;
        } else {
            Vector raw_x = Vector::Zero(n_raw());
            for (int i = 0; i < row_index_.size(); ++i)
                raw_x[row_index_(i)] += x[i];
            out.noalias() = raw_data().transpose() * raw_x;
        }
        return out;
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

            // const Matrix XtX = data().transpose() * data();
            // const Matrix Sigma = XtX / den;

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
    void invalidate_M_after_data_change_() {
        raw_score_cache_ready_ = false;
        if (mode_ != Mode::CovMax) invalidate_M_();
    }
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

        invalidate_M_after_data_change_();
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
        if constexpr (std::same_as<SamplingStrategy, IndependentSampling>) return s;

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
    mutable Vector raw_score_cache_;

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
    mutable bool raw_score_cache_ready_ {false};
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
    using BinaryMatrixT = typename Base::BinaryMatrixT;

    using Base::init;
    using Base::M;
    using Base::ginvM;
    using Base::n;
    using Base::m;
    using Base::n_dofs_weights;
    using Base::data;
    using Base::data_transpose_times_;
    using Base::h;
    using Base::mode;
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
    [[nodiscard]] SparseMatrix Psi_at_(const Matrix&) const override {
        throw std::logic_error("MultivariateBlock: locations cannot define a weight basis; pass a Psi matrix explicitly");
    }
    [[nodiscard]] SparseMatrix Psi_at_(const BinaryMatrixT&) const override {
        throw std::logic_error("MultivariateBlock: locations cannot define a weight basis; pass a Psi matrix explicitly");
    }

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

        Vector z = data_transpose_times_(nu);

        if (weight_sign_constraint() == WeightSignConstraint::NonNegative) {
            return solve_nonnegative_weight_ipopt_(z); // already normalized
        }

        if (mode() == Mode::CovMax) return normalize_weight_(z);
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
    using BinaryMatrixT = typename Base::BinaryMatrixT;
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
        weights_solver_weights_ready_ = false;
    }

    std::unique_ptr<Base> clone() const override {
        return std::make_unique<FunctionalBlock>(*this);
    }

    // Psi_D
    [[nodiscard]] const SparseMatrix& Psi_D() const override { return weights_solver_.Psi(); }
    [[nodiscard]] SparseMatrix Psi_at_(const Matrix& locs) const override { return weights_solver_.eval_basis_at(locs); }
    [[nodiscard]] SparseMatrix Psi_at_(const BinaryMatrixT& locs) const override { return weights_solver_.eval_basis_at(locs); }

    // weights regularization utils
    void set_lambda_weights(const double lambda) override {
        validate_rgcca_positive_regularization_lambda_(lambda, "weight lambda");
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
    using Base::data_transpose_times_;

    template <typename GeoFrame>
    void init_functional_(GeoFrame& gf, WeightsPenaltyType&& weights_penalty) {
        weights_solver_.discretize(weights_penalty.get());
        weights_solver_.analyze_data(gf, M());
        weights_solver_weights_ready_ = true;
        init();
    }

    Vector w_fit_(const Vector& nu) override {
        assert(nu.size() == n() && "nu must have size n (rows of X)");
        init();

        Vector z = data_transpose_times_(nu);

        if (weight_sign_constraint() == WeightSignConstraint::NonNegative) {
            return solve_nonnegative_weight_ipopt_(z);  // already normalized
        }

        if (!weights_solver_weights_ready_) {
            weights_solver_.update_weights(M());
            weights_solver_weights_ready_ = true;
        }
        weights_solver_.update_z(z);
        weights_solver_.fit(lambda_weights_);
        return weights_solver_.f(); // already normalized
    }
    void invalidate_derived_caches_() override {
        Omega_ready_ = false;
        weights_solver_weights_ready_ = false;
        reset_nonnegative_weight_solver_();
    }

private:
    SparseMatrix Omega_;
    bool Omega_ready_ {false};
    bool weights_solver_weights_ready_ {false};
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
    Matrix correlation_matrix;
    std::vector<double> tau_values;
    std::vector<double> lambda_components_values;
    std::vector<double> lambda_weights_values;
    std::vector<bool> active_blocks;
    double rho_tot = std::numeric_limits<double>::quiet_NaN();
    double rho_tot_p_value = std::numeric_limits<double>::quiet_NaN();
    int rho_tot_bootstrap_count = 0;
    bool component_significant = true;

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
        bool verbose;
        bool cache_covariances;
        bool bias;
        InitStrategy init_strategy;
        LambdaSelection lambda_selection_weights;
        LambdaSelection lambda_selection_components;
        bool component_significance;
        bool block_deactivation;
        bool connection_deactivation;
        Mode mode;
        WeightSignConstraint weight_sign_constraint;
        Deflation deflation_mode;
        Scheme scheme;

        explicit Options(
          const int max_iter_ = 1000, const double tol_ = 1e-8, const bool bias_ = true,
          const InitStrategy init_strategy_ = InitStrategy::SVD, const Mode mode_ = Mode::CovMax,
          const WeightSignConstraint weight_sign_constraint_ = WeightSignConstraint::None,
          const LambdaSelection lambda_selection_weights_ = LambdaSelection::Manual,
          const LambdaSelection lambda_selection_components_ = LambdaSelection::Automatic,
          const bool component_significance_ = false,
          const Deflation deflation_mode_ = Deflation::Scores, const Scheme& scheme_ = Scheme::Factorial(),
          const bool verbose_ = false, const bool cache_ = true,
          const bool block_deactivation_ = false, const bool connection_deactivation_ = false) :
            max_iter(max_iter_),
            tol(tol_),
            bias(bias_),
            init_strategy(init_strategy_),
            mode(mode_),
            weight_sign_constraint(weight_sign_constraint_),
            lambda_selection_weights(lambda_selection_weights_),
            lambda_selection_components(lambda_selection_components_),
            component_significance(component_significance_),
            block_deactivation(block_deactivation_),
            connection_deactivation(connection_deactivation_),
            deflation_mode(deflation_mode_),
            scheme(scheme_),
            verbose(verbose_),
            cache_covariances(cache_) { }

        friend std::ostream& operator<<(std::ostream& os, const Options& opt) {
            os << "RGCCA::Options {\n"
               << "  max_iter                    = " << opt.max_iter << '\n'
               << "  tol                         = " << opt.tol << '\n'
               << "  verbose                     = " << opt.verbose << '\n'
               << "  cache_covariances           = " << opt.cache_covariances << '\n'
               << "  bias                        = " << opt.bias << '\n'
               << "  init_strategy               = " << to_string(opt.init_strategy) << '\n'
               << "  lambda_selection_weights    = " << to_string(opt.lambda_selection_weights) << '\n'
               << "  lambda_selection_components = " << to_string(opt.lambda_selection_components) << '\n'
               << "  component_significance      = " << opt.component_significance << '\n'
               << "  block_deactivation          = " << opt.block_deactivation << '\n'
               << "  connection_deactivation     = " << opt.connection_deactivation << '\n'
               << "  mode                        = " << to_string(opt.mode) << '\n'
               << "  weight_sign_constraint      = " << to_string(opt.weight_sign_constraint) << '\n'
               << "  deflation_mode              = " << to_string(opt.deflation_mode) << '\n'
               << "  scheme                      = " << opt.scheme.name << '\n'
               << "}";
            return os;
        }
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

        unsigned seed = 12345;
        int max_threads = 12;

        int B_min = 500;
        int B_max = 1000;
        int B_per_thread_per_batch = 5;

        // adaptive batch
        bool adaptive = true;
        double adaptive_tol = 1e-3;
        int stable_batches_required = 3;

        // block deactivation
        double active_block_tol = 1e-8;

        // connection deactivation
        double active_connection_sign_stability = 0.95;
        double active_connection_min_abs_corr = 0.05;

        // confidence intervals
        double ci_level = 0.95;

        // early stop
        int patience = 1;

        ResamplingStrategy resampling_strategy = ResamplingStrategy::Ordinary;

        // stationary bootstrap: expected block length = 1 / p
        double stationary_block_length = 10.0;

        // component significance test for H0: rho_tot = 0
        int component_significance_resamples = 100;
        double component_significance_alpha = 0.05;

        friend std::ostream& operator<<(std::ostream& os, const BootstrapConfig& config) {
            os << "RGCCA::BootstrapConfig {\n"
               << "  seed                               = " << config.seed << '\n'
               << "  max_threads                        = " << config.max_threads << '\n'
               << "  B_min                              = " << config.B_min << '\n'
               << "  B_max                              = " << config.B_max << '\n'
               << "  B_per_thread_per_batch             = " << config.B_per_thread_per_batch << '\n'
               << "  adaptive                           = " << config.adaptive << '\n'
               << "  adaptive_tol                       = " << config.adaptive_tol << '\n'
               << "  stable_batches_required            = " << config.stable_batches_required << '\n'
               << "  active_block_tol                   = " << config.active_block_tol << '\n'
               << "  active_connection_sign_stability   = " << config.active_connection_sign_stability << '\n'
               << "  active_connection_min_abs_corr     = " << config.active_connection_min_abs_corr << '\n'
               << "  ci_level                           = " << config.ci_level << '\n'
               << "  patience                           = " << config.patience << '\n'
               << "  resampling_strategy                = " << to_string(config.resampling_strategy) << '\n'
               << "  stationary_block_length            = " << config.stationary_block_length << '\n'
               << "  component_significance_resamples   = " << config.component_significance_resamples << '\n'
               << "  component_significance_alpha       = " << config.component_significance_alpha << '\n'
               << "}";
            return os;
        }
    };
    struct BootstrapSelectionResult {
        using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

        int h = 0;
        int B = 0;

        std::vector<double> lambda_grid;
        std::vector<double> criterion;

        double lambda_opt = std::numeric_limits<double>::quiet_NaN();
        int lambda_opt_index = -1;
        double ci_level = std::numeric_limits<double>::quiet_NaN();

        std::vector<std::string> block_names;
        std::vector<bool> active_blocks;

        // [lambda][block] -> vector/matrix
        std::vector<std::vector<Vector>> w_fit_by_lambda;
        std::vector<std::vector<Matrix>> w_boot_by_lambda;
        std::vector<std::vector<Vector>> w_min_by_lambda;
        std::vector<int> B_used_by_lambda;

        // [lambda] -> vector/matrix
        std::vector<Matrix> corr_boot_by_lambda;
        std::vector<Matrix> corr_min_by_lambda;

        // CI
        std::vector<Matrix> corr_ci_low_by_lambda;    // [lambda] -> J x J
        std::vector<Matrix> corr_ci_high_by_lambda;   // [lambda] -> J x J
        // weight CIs are location-dependent; compute them with RGCCA::bootstrap_weights_ci(...).

        BootstrapSelectionResult() = default;

        BootstrapSelectionResult(
            const int h_,
            const int B_,
            const std::vector<double>& lambda_grid_,
            const std::vector<std::string>& block_names_,
            const std::vector<int>& block_dims_,
            const double ci_level_
        ) :
            h(h_),
            B(B_),
            lambda_grid(lambda_grid_),
            criterion(lambda_grid_.size(), -std::numeric_limits<double>::infinity()),
            ci_level(ci_level_),
            block_names(block_names_)
        {
            const std::size_t n_lambda = lambda_grid.size();
            const std::size_t J = block_dims_.size();

            w_fit_by_lambda.resize(n_lambda);
            w_boot_by_lambda.resize(n_lambda);
            w_min_by_lambda.resize(n_lambda);
            B_used_by_lambda.resize(n_lambda);

            active_blocks.resize(J);

            for (std::size_t i = 0; i < n_lambda; ++i) {
                w_boot_by_lambda[i].resize(J);

                for (std::size_t j = 0; j < J; ++j) {
                    w_boot_by_lambda[i][j].setZero(block_dims_[j], B);
                }
            }

            corr_min_by_lambda.resize(n_lambda);
            corr_boot_by_lambda.resize(n_lambda);
            corr_ci_low_by_lambda.resize(n_lambda);
            corr_ci_high_by_lambda.resize(n_lambda);

            for (std::size_t i = 0; i < n_lambda; ++i) {
                corr_boot_by_lambda[i].setZero(J * J, B);
                corr_min_by_lambda[i].setZero(J , J);
                corr_ci_low_by_lambda[i].setZero(J, J);
                corr_ci_high_by_lambda[i].setZero(J, J);
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

        compute_Psi_();

        initialized_ = true;
    }

    // weights and components regularization utilities
    void set_lambda_weights_all(const double lambda) const {
        internals::validate_rgcca_positive_regularization_lambda_(lambda, "weight lambda");
        for (auto& b : blocks_) b->set_lambda_weights(lambda);
    }
    void set_lambda_components_all(const double lambda) const {
        internals::validate_rgcca_positive_regularization_lambda_(lambda, "component lambda");
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
    std::vector<Result> fit() {
        if (!initialized_) init();

        const int J = n_blocks();
        if (J < 2) throw std::runtime_error("RGCCA: need ≥ 2 blocks");
        const bool run_model_selection = bootstrap_model_selection_requested_();
        fdapde::cout << opt_ << std::endl;
        if (run_model_selection || opt_.component_significance) {
            fdapde::cout << bootstrap_config_ << std::endl;
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
            auto step_start = log_step_start_("Final component fit");
            auto blocks = main_blocks_();
            const auto active_blocks = active_blocks_from_C_(C_active);
            init_comp_(blocks, InitStrategy::None, true, &active_blocks);
            Result component_result = fit_component_(blocks, C_active);
            log_step_end_(step_start);

            if (opt_.component_significance) {
                const auto significance = bootstrap_test_component_significance_(C_active);
                annotate_component_significance_(component_result, significance);

                if (!significance.significant) {
                    results.push_back(std::move(component_result));
                    append_inactive_components_(results, hh + 1, J);
                    break;
                }
            }

            results.push_back(std::move(component_result));

            step_start = log_step_start_("Deflate blocks");
            deflate_all_();
            log_step_end_(step_start);
        }

        // post-processing weights
        auto step_start = log_step_start_("Compute weights_star");
        compute_weights_star_();
        log_step_end_(step_start);

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
    [[nodiscard]] const std::vector<BootstrapSelectionResult>& bootstrap_selection_results() const { return bootstrap_selection_results_; }

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

    // initialization utils
    void check_index_(int j) const {
        if (j < 0 || j >= static_cast<int>(blocks_.size())) throw std::out_of_range("block index");
    }
    const BootstrapSelectionResult& bootstrap_selection_result_(const int h) const {
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
    int bootstrap_lambda_opt_index_(const BootstrapSelectionResult& boot_results) const {
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
        const BootstrapSelectionResult& boot_results,
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

            ci_low[r] = empirical_quantile_(values, alpha_low);
            ci_high[r] = empirical_quantile_(values, alpha_high);
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
        const bool update_component_lambdas = true
    ) {
        const int J = n_blocks();
        FitWorkspace ws(J);

        // room for results
        Result res(J);
        res.obj_history.reserve(opt_.max_iter);

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
        auto a_prev = snapshot_weights_(blocks);
        if (update_component_lambdas && opt_.lambda_selection_components == LambdaSelection::Automatic)
            set_lambda_components_auto_all_(blocks);

        // main loop
        for (int s = 0; s < opt_.max_iter; ++s) {
            for (int l = 0; l < J; ++l) {

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
    double weights_variation_(const BlockRefList& blocks, const std::vector<Vector>& a_prev) const {
        double acc = 0.0;

        for (int j = 0; j < n_blocks(); ++j) {
            const auto aj = blocks[j]->weights().col(h_);
            acc += (aj - a_prev[j]).squaredNorm();
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
            B_batch = bootstrap_config.adaptive ?
                n_threads * bootstrap_config.B_per_thread_per_batch : B_max;
            corr_pos_count.setZero(J, J);
            corr_neg_count.setZero(J, J);
        }

        void reset() {
            B_done = 0;
            stable_batches = 0;
            crit_prev_batch = std::numeric_limits<double>::infinity();
            crit = std::numeric_limits<double>::quiet_NaN();
            corr_pos_count.setZero(J, J);
            corr_neg_count.setZero(J, J);
        }

        // config
        int n_threads;
        int J;
        int seed;
        int B_min;
        int B_max;
        int B_batch;

        // state
        int B_done = 0;
        int B_run = 0;

        Eigen::MatrixXi corr_pos_count;
        Eigen::MatrixXi corr_neg_count;

        // adaptive batch
        int stable_batches = 0;
        double crit_prev_batch = std::numeric_limits<double>::infinity();
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

    std::chrono::high_resolution_clock::time_point log_step_start_(std::string_view label) const {
        fdapde::cout << label << " --> " << std::flush;
        return std::chrono::high_resolution_clock::now();
    }
    void log_step_end_(const std::chrono::high_resolution_clock::time_point start) const {
        const auto end = std::chrono::high_resolution_clock::now();
        const double elapsed_sec = std::chrono::duration<double>(end - start).count();
        fdapde::cout << "<-- " << std::fixed << std::setprecision(3)
                  << elapsed_sec << std::defaultfloat << "s" << std::endl;
    }

    ModelSelectionResult bootstrap_model_selection_() {
        fdapde::cout << "\n=========================================" << std::endl;
        fdapde::cout << "Bootstrap model selection for component " << h_ + 1 << std::endl;
        fdapde::cout << "=========================================\n" << std::endl;

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
        auto step_start = log_step_start_("Init bootstrap");
        AdaptiveBootstrapState bootstrap_state(bootstrap_config_, n_threads, h_, J);
        BootstrapSelectionResult boot_results(
            h_, bootstrap_state.B_max, lambda_grid,
            block_names_(blocks), block_dims_(blocks), bootstrap_config_.ci_level
        );
        log_step_end_(step_start);

        // preliminary fit
        step_start = log_step_start_("Preliminary fit");
        if (select_lambda)
            set_lambda_weights_all(lambda_grid.back());
        init_comp_(blocks);
        fit_component_(blocks, C_active);
        auto preliminary_w_fit = snapshot_weights_(blocks);
        log_step_end_(step_start);

        step_start = log_step_start_("  Clone worker blocks");
        std::vector<BootstrapBlocks> thread_boot_worker(n_threads);
        for (int t = 0; t < n_threads; ++t) {
            thread_boot_worker[t] = clone_blocks_();
        }
        log_step_end_(step_start);

        int n_lambda = static_cast<int>(lambda_grid.size());
        for (int lambda_i = n_lambda - 1; lambda_i >= 0; --lambda_i) {
            const bool reuse_preliminary_fit = lambda_i == n_lambda - 1;

            if (select_lambda) {
                // current lambda
                const double lambda = lambda_grid[lambda_i];
                fdapde::cout << "- lambda = " << lambda << std::endl;
                if (!reuse_preliminary_fit) {
                    set_lambda_weights_all(lambda);
                    for (auto& worker : thread_boot_worker)
                        set_lambda_weights_all_(worker.refs, lambda);
                }
            } else {
                fdapde::cout << "- fixed weight regularization" << std::endl;
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

            auto start = std::chrono::high_resolution_clock::now();

            // bootstrapping by batch
            bootstrap_state.reset();
            while (bootstrap_state.B_done < bootstrap_state.B_max) {
                bootstrap_state.B_run = std::min(bootstrap_state.B_batch, bootstrap_state.B_max - bootstrap_state.B_done);

                const int batch_first = bootstrap_state.B_done;
                const int batch_last = bootstrap_state.B_done + bootstrap_state.B_run - 1;
                fdapde::cout << "  "
                          << (bootstrap_config_.adaptive ? "Adaptive" : "Bootstrap")
                          << " batch [" << std::setw(4) << batch_first
                          << ", " << std::setw(4) << batch_last << "] --> "
                          << std::flush;

                const auto batch_step_start = std::chrono::high_resolution_clock::now();
                auto bootstrap_timing = run_bootstrap_batch_(
                    lambda_i,
                    bootstrap_state,
                    thread_boot_worker,
                    C_active,
                    w_fit,
                    w_min,
                    boot_results
                );
                log_step_end_(batch_step_start);

                step_start = log_step_start_("    Post-batch update");
                int n_active_blocks = count_active_blocks_(C_active);
                if (opt_.block_deactivation)
                    n_active_blocks = threshold_inactive_blocks_(w_min, C_active);
                int n_active_connections = count_active_connections_(C_active);
                bootstrap_state.crit = criterion_score_with_weights_(blocks, w_min, C_);
                log_step_end_(step_start);

                fdapde::cout << "  Batch summary: avg_fit_time = " << std::fixed << std::setprecision(3) << bootstrap_timing.avg_fit_time
                          << " ± " << bootstrap_timing.sd_fit_time << std::defaultfloat << "s";
                fdapde::cout << ", eff = " << std::fixed << std::setprecision(1)
                          << 100.0 * bootstrap_timing.efficiency << "%" << std::defaultfloat;

                fdapde::cout << " | ab = " << n_active_blocks;
                fdapde::cout << ", ac = " << n_active_connections;

                fdapde::cout << " | crit = " << std::fixed << std::setprecision(3) <<  bootstrap_state.crit;
                fdapde::cout << std::defaultfloat;

                bootstrap_state.B_done += bootstrap_state.B_run;

                auto print_bootstrap_timing_debug = [&]() {
                    fdapde::cout << "    timing:" << std::endl;
                    fdapde::cout << "      - wall(s): setup=" << std::fixed << std::setprecision(3)
                              << bootstrap_timing.setup_time
                              << ", parallel=" << bootstrap_timing.parallel_time
                              << ", merge=" << bootstrap_timing.merge_time << std::endl;
                    fdapde::cout << "      - efficiency: fit=" << std::fixed << std::setprecision(1)
                              << 100.0 * bootstrap_timing.efficiency
                              << "%, task=" << 100.0 * bootstrap_timing.task_efficiency
                              << "%" << std::endl;
                    fdapde::cout << "      - avg_task(s): prep=" << std::fixed << std::setprecision(3)
                              << bootstrap_timing.avg_prep_time
                              << ", fit=" << bootstrap_timing.avg_fit_time
                              << ", snapshot=" << bootstrap_timing.avg_snapshot_time
                              << ", corr=" << bootstrap_timing.avg_corr_time
                              << ", total=" << bootstrap_timing.avg_task_time
                              << std::defaultfloat << std::endl;
                };

                if (bootstrap_config_.adaptive) {
                    if (adaptive_stop_(bootstrap_state, bootstrap_config_)) {
                        print_bootstrap_timing_debug();
                        break;
                    }
                    print_bootstrap_timing_debug();
                } else {
                    fdapde::cout << std::endl;
                    print_bootstrap_timing_debug();
                }
            }

            BoolMatrix C_lambda = C_active;
            int n_active_connections = count_active_connections_(C_lambda);
            if (opt_.connection_deactivation) {
                step_start = log_step_start_("  Threshold inactive connections");
                n_active_connections = threshold_inactive_connections_(
                    lambda_i, bootstrap_state, boot_results, C_lambda
                );
                log_step_end_(step_start);
            }

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            const double elapsed_sec = duration.count() / 1000.0;

            boot_results.w_fit_by_lambda[lambda_i] = w_fit;
            boot_results.w_min_by_lambda[lambda_i] = w_min;
            step_start = log_step_start_("  Final lambda correlation");
            correlation_matrix_(blocks, w_min, boot_results.corr_min_by_lambda[lambda_i]);
            boot_results.corr_min_by_lambda[lambda_i].array() *= (C_lambda.cast<double>() + Matrix::Identity(J, J)).array();
            log_step_end_(step_start);
            boot_results.criterion[lambda_i] = bootstrap_state.crit;
            boot_results.B_used_by_lambda[lambda_i] = bootstrap_state.B_done;

            if (bootstrap_state.crit > bootstrap_state.best_criterion) {
                C_best = C_lambda;
            }

            fdapde::cout << "  Bootstrap used: " << bootstrap_state.B_done
                      << ", execution time: " << std::fixed << std::setprecision(3) << elapsed_sec << std::defaultfloat << "s"
                      << ", ac = " << n_active_connections
                      << ", crit = " << bootstrap_state.crit << std::endl;

            // early stop
            if (early_stop_lambda_(bootstrap_state, lambda_i, bootstrap_config_)) {
                break;
            }

            // reset connections, but keep fully deactivated blocks off
            C_active = reset_connections_keep_inactive_blocks_(C_, C_active);

        }

        step_start = log_step_start_("Resize bootstrap results");
        resize_bootstrap_results_(boot_results, static_cast<int>(lambda_grid.size()), J);
        log_step_end_(step_start);
        step_start = log_step_start_("Compute bootstrap correlation CIs");
        compute_bootstrap_corr_cis_(boot_results, J);
        log_step_end_(step_start);

        if (bootstrap_state.best_i < 0)
            throw std::runtime_error("No model candidate was evaluated during bootstrap selection");

        boot_results.lambda_opt_index = bootstrap_state.best_i;
        if (select_lambda) {
            boot_results.lambda_opt = lambda_grid[bootstrap_state.best_i];
            fdapde::cout << "\nOptimal lambda: " << boot_results.lambda_opt << std::endl;
        } else {
            fdapde::cout << "\nNo weight lambda selection requested" << std::endl;
        }


        boot_results.active_blocks = active_blocks_from_C_(C_best);
        fdapde::cout << "\nBlock deactivation:" << std::endl;
        for (int j = 0; j < J; ++j) {
            const double nrm = boot_results.w_min_by_lambda[bootstrap_state.best_i][j].norm();
            fdapde::cout << "  - block[" << std::setw(2) << j << "]: "
                      << "||w_min|| = " << std::fixed << std::setprecision(4) << nrm
                      << std::defaultfloat
                      << ", active = " << (boot_results.active_blocks[j] ? "yes" : "no")
                      << std::endl;
        }
        fdapde::cout << "\nUpdated design matrix:" << std::endl;
        fdapde::cout << C_best << std::endl;

        bootstrap_selection_results_.push_back(std::move(boot_results));

        ModelSelectionResult out;
        out.lambda_selected = select_lambda;
        out.lambda = bootstrap_selection_results_.back().lambda_opt;
        out.C_active = C_best;
        return out;
    }

    struct BootstrapBatchTiming {
        double avg_fit_time = 0.0;
        double sd_fit_time = 0.0;
        double wall_time = 0.0;
        double efficiency = 0.0;
        double task_efficiency = 0.0;
        double setup_time = 0.0;
        double parallel_time = 0.0;
        double merge_time = 0.0;
        double avg_task_time = 0.0;
        double avg_prep_time = 0.0;
        double avg_snapshot_time = 0.0;
        double avg_corr_time = 0.0;
    };
    struct ComponentSignificanceResult {
        double rho_tot = std::numeric_limits<double>::quiet_NaN();
        double p_value = std::numeric_limits<double>::quiet_NaN();
        int B = 0;
        bool significant = true;
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
            const auto step_start = log_step_start_("Inactive component fit");
            Result inactive_result = fit_component_(C_inactive);
            log_step_end_(step_start);
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

        auto step_start = log_step_start_("  Clone significance worker blocks");
        std::vector<BootstrapBlocks> thread_boot_worker(n_threads);
        for (int t = 0; t < n_threads; ++t)
            thread_boot_worker[t] = clone_blocks_();
        log_step_end_(step_start);

        std::vector<int> thread_ge_count(n_threads, 0);
        const unsigned seed = bootstrap_config_.seed + static_cast<unsigned>(1000003 * (h_ + 1));
        const auto active_blocks = active_blocks_from_C_(C_active);

        step_start = log_step_start_("  Significance bootstrap");
        parallel_for(0, B, 1, [&](int b) {
            const int tid = this_thread_id();
            auto& boot_blocks = thread_boot_worker[tid];

            set_permuted_row_index_all_(
                boot_blocks.refs,
                seed + static_cast<unsigned>(7919 * (b + 1))
            );

            init_comp_(boot_blocks.refs, InitStrategy::WarmStart, false, &active_blocks);
            fit_component_(boot_blocks.refs, C_active, false);

            const double rho_star = rho_tot_(boot_blocks.refs, C_active);
            if (std::isfinite(rho_star) && rho_star >= out.rho_tot)
                ++thread_ge_count[tid];

            clear_row_index_all_(boot_blocks.refs);
        });
        log_step_end_(step_start);

        const int ge_count = std::accumulate(thread_ge_count.begin(), thread_ge_count.end(), 0);
        out.p_value = static_cast<double>(ge_count) / static_cast<double>(B);
        out.significant = out.p_value <= bootstrap_config_.component_significance_alpha;

        fdapde::cout << "\nSignificance:" << std::endl
                  << "  - rho_tot = " << out.rho_tot << std::endl
                  << "  - p-value = " << out.p_value << std::endl
                  << "  - significant = " << out.significant << std::endl;

        return out;
    }
    BootstrapBatchTiming run_bootstrap_batch_(
        int lambda_i, AdaptiveBootstrapState& bootstrap_state,
        std::vector<BootstrapBlocks>& thread_boot_worker,
        const BoolMatrix& C_active,
        const std::vector<Vector>& w_fit,
        std::vector<Vector>& w_min,
        BootstrapSelectionResult& boot_results
    ) {

        const int J = static_cast<int>(w_fit.size());
        const int n_threads = bootstrap_state.n_threads;
        const int B_run = bootstrap_state.B_run;
        const int B_offset = bootstrap_state.B_done;
        const int seed = bootstrap_state.seed;
        const auto active_blocks = active_blocks_from_C_(C_active);

        const auto setup_start = std::chrono::high_resolution_clock::now();
        auto bootstrap_idx = make_bootstrap_indices_(B_run, B_offset, seed);

        std::vector<double> fit_times_sec(B_run, 0.0);
        std::vector<double> task_times_sec(B_run, 0.0);
        std::vector<double> prep_times_sec(B_run, 0.0);
        std::vector<double> snapshot_times_sec(B_run, 0.0);
        std::vector<double> corr_times_sec(B_run, 0.0);
        const auto setup_end = std::chrono::high_resolution_clock::now();

        const auto parallel_start = std::chrono::high_resolution_clock::now();

        parallel_for(0, B_run, 1, [&](int b) {

            const auto task_start = std::chrono::high_resolution_clock::now();
            const int tid = this_thread_id();
            auto& boot_blocks = thread_boot_worker[tid];
            copy_weights_snapshot_(boot_blocks.refs, w_fit);

            set_row_index_all_(boot_blocks.refs, bootstrap_idx[b]);
            const auto prep_end = std::chrono::high_resolution_clock::now();

            const auto fit_start = std::chrono::high_resolution_clock::now();
            init_comp_(boot_blocks.refs, InitStrategy::WarmStart, true, &active_blocks);
            fit_component_(boot_blocks.refs, C_active);
            const auto fit_end = std::chrono::high_resolution_clock::now();

            fit_times_sec[b] = std::chrono::duration<double>(fit_end - fit_start).count();

            auto w_b = snapshot_weights_(boot_blocks.refs);
            const int b_global = B_offset + b;

            for (int j = 0; j < J; ++j) {
                if (w_b[j].dot(w_fit[j]) < 0.0) w_b[j] *= -1.0;
                boot_results.w_boot_by_lambda[lambda_i][j].col(b_global) = w_b[j];
            }
            const auto snapshot_end = std::chrono::high_resolution_clock::now();

            Matrix corr_b;
            correlation_matrix_raw_data_(boot_blocks.refs, w_b, corr_b);
            boot_results.corr_boot_by_lambda[lambda_i].col(b_global) = Eigen::Map<const Vector>(corr_b.data(), J * J);
            const auto corr_end = std::chrono::high_resolution_clock::now();

            prep_times_sec[b] = std::chrono::duration<double>(prep_end - task_start).count();
            snapshot_times_sec[b] = std::chrono::duration<double>(snapshot_end - fit_end).count();
            corr_times_sec[b] = std::chrono::duration<double>(corr_end - snapshot_end).count();
            task_times_sec[b] = std::chrono::duration<double>(corr_end - task_start).count();

        });

        const auto parallel_end = std::chrono::high_resolution_clock::now();

        const auto merge_start = std::chrono::high_resolution_clock::now();

        for (auto& boot_worker : thread_boot_worker)
            clear_row_index_all_(boot_worker.refs);

        for (int b = 0; b < B_run; ++b) {
            const int b_global = B_offset + b;
            const auto corr_col = boot_results.corr_boot_by_lambda[lambda_i].col(b_global);

            for (int j = 0; j < J; ++j) {
                for (int k = j + 1; k < J; ++k) {
                    const double c = corr_col(j + k * J);

                    if (c > 0.0) {
                        ++bootstrap_state.corr_pos_count(j, k);
                        ++bootstrap_state.corr_pos_count(k, j);
                    } else if (c < 0.0) {
                        ++bootstrap_state.corr_neg_count(j, k);
                        ++bootstrap_state.corr_neg_count(k, j);
                    }
                }
            }
        }

        parallel_for(0, J, 1, [&](int j) {
            for (int b = 0; b < B_run; ++b) {
                const int b_global = B_offset + b;
                const auto w_bj = boot_results.w_boot_by_lambda[lambda_i][j].col(b_global);
                for (int r = 0; r < w_min[j].size(); ++r) {
                    if (w_fit[j][r] > 0.0) {
                        w_min[j][r] = std::max(0.0, std::min(w_min[j][r], w_bj[r]));
                    } else if (w_fit[j][r] < 0.0) {
                        w_min[j][r] = std::min(0.0, std::max(w_min[j][r], w_bj[r]));
                    } else {
                        w_min[j][r] = 0.0;
                    }
                }
            }
        });

        const auto merge_end = std::chrono::high_resolution_clock::now();

        const double setup_time = std::chrono::duration<double>(setup_end - setup_start).count();
        const double parallel_time = std::chrono::duration<double>(parallel_end - parallel_start).count();
        const double merge_time = std::chrono::duration<double>(merge_end - merge_start).count();
        const double wall_time = setup_time + parallel_time + merge_time;

        double total_fit_time = 0.0;
        for (double t : fit_times_sec)
            total_fit_time += t;
        double avg_fit_time = total_fit_time / static_cast<double>(B_run);
        double var = 0.0;
        for (double t : fit_times_sec) {
            const double d = t - avg_fit_time;
            var += d * d;
        }

        double total_task_time = 0.0;
        double total_prep_time = 0.0;
        double total_snapshot_time = 0.0;
        double total_corr_time = 0.0;
        for (int b = 0; b < B_run; ++b) {
            total_task_time += task_times_sec[b];
            total_prep_time += prep_times_sec[b];
            total_snapshot_time += snapshot_times_sec[b];
            total_corr_time += corr_times_sec[b];
        }

        BootstrapBatchTiming timing;
        timing.avg_fit_time = avg_fit_time;
        timing.sd_fit_time = std::sqrt(var / std::max(1, B_run - 1));
        timing.wall_time = wall_time;
        timing.efficiency = total_fit_time / (parallel_time * static_cast<double>(n_threads));
        timing.task_efficiency = total_task_time / (parallel_time * static_cast<double>(n_threads));
        timing.setup_time = setup_time;
        timing.parallel_time = parallel_time;
        timing.merge_time = merge_time;
        timing.avg_task_time = total_task_time / static_cast<double>(B_run);
        timing.avg_prep_time = total_prep_time / static_cast<double>(B_run);
        timing.avg_snapshot_time = total_snapshot_time / static_cast<double>(B_run);
        timing.avg_corr_time = total_corr_time / static_cast<double>(B_run);

        return timing;
    }

    bool adaptive_stop_(AdaptiveBootstrapState& state, const BootstrapConfig& config) const {
        if (state.crit == 0.0) {
            fdapde::cout << ", adaptive stop (crit = 0)" << std::endl;
            return true;
        }

        if (!std::isfinite(state.crit_prev_batch)) {
            state.crit_prev_batch = state.crit;
            fdapde::cout << std::endl;
            return false;
        }

        const double rel_change =
            std::abs(state.crit_prev_batch - state.crit) /
            (std::abs(state.crit_prev_batch) + 1e-12);

        fdapde::cout << ", rel_change = " << std::fixed << std::setprecision(3) << rel_change;

        if (rel_change < config.adaptive_tol) {
            ++state.stable_batches;
        } else {
            state.stable_batches = 0;
        }

        state.crit_prev_batch = state.crit;

        if (state.stable_batches >= config.stable_batches_required && state.B_done >= state.B_min) {
            fdapde::cout << ", adaptive stop (stable)" << std::endl;
            return true;
        }

        fdapde::cout << std::endl;
        return false;
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
                      << " consecutive lambdas" << std::endl;
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
        const BootstrapSelectionResult& boot_results,
        BoolMatrix& C_active
    ) {
        const int J = static_cast<int>(C_active.rows());
        const int B_eff = state.B_done;

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

                const double med_abs_corr = median_(abs_corr);

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
    double median_(std::vector<double>& x) const {
        if (x.empty())
            return 0.0;

        const std::size_t n = x.size();
        const std::size_t mid = n / 2;

        std::nth_element(x.begin(), x.begin() + mid, x.end());

        if (n % 2 == 1)
            return x[mid];

        const double upper = x[mid];

        std::nth_element(x.begin(), x.begin() + mid - 1, x.end());
        const double lower = x[mid - 1];

        return 0.5 * (lower + upper);
    }
    double empirical_quantile_(std::vector<double>& x, const double p) const {
        x.erase(
            std::remove_if(x.begin(), x.end(), [](const double v) { return !std::isfinite(v); }),
            x.end()
        );

        if (x.empty())
            return std::numeric_limits<double>::quiet_NaN();

        std::sort(x.begin(), x.end());

        if (x.size() == 1)
            return x.front();

        const double pos = std::clamp(p, 0.0, 1.0) * static_cast<double>(x.size() - 1);
        const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
        const std::size_t hi = static_cast<std::size_t>(std::ceil(pos));
        const double frac = pos - static_cast<double>(lo);

        return (1.0 - frac) * x[lo] + frac * x[hi];
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

        const double z_low = empirical_quantile_(z_values, alpha_low);
        const double z_high = empirical_quantile_(z_values, alpha_high);

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
    void resize_bootstrap_results_(BootstrapSelectionResult& boot_results, int n_lambdas, int J) const {
        for (int i = 0; i < n_lambdas; ++i) {
            const int B_eff = boot_results.B_used_by_lambda[i];

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
    void compute_bootstrap_corr_cis_(BootstrapSelectionResult& boot_results, int J) const {
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

    std::vector<typename Block::IndexVector> make_bootstrap_indices_(int B_run, int B_offset, unsigned seed) const {
        std::vector<typename Block::IndexVector> bootstrap_idx(B_run);

        for (int b = 0; b < B_run; ++b) {
            std::mt19937_64 rng(seed + static_cast<unsigned>(B_offset + b));
            bootstrap_idx[b] = bootstrap_indices_(n_, rng);
        }

        return bootstrap_idx;
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
                "for TimeDependentSampling"
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
        if (bootstrap_config_.B_per_thread_per_batch <= 0)
            throw std::invalid_argument("RGCCA: bootstrap B_per_thread_per_batch must be positive");
        if (bootstrap_config_.stable_batches_required <= 0)
            throw std::invalid_argument("RGCCA: bootstrap stable_batches_required must be positive");
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
            return Psi_T() * b.components().col(h_);
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
                eta = Psi_T() * eta;
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


// pretty printer for a single Result
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
    os << "n_iters: " << r.iters << "\n";
    os << "monotone: " << (r.monotone ? "yes" : "no") << "\n";
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
        os << "\ncorrelation matrix :\n";
        os << std::fixed << std::setprecision(2);
        os << r.correlation_matrix << std::endl;
        os << std::fixed << std::setprecision(8);
    }
    os << std::endl;
    if (std::isfinite(r.rho_tot_p_value)) {
        os << "significance:\n";
        os << "- rho_tot: " << std::fixed << r.rho_tot << "\n";
        os << "- p-value: " << std::fixed << r.rho_tot_p_value
           << " (" << r.rho_tot_bootstrap_count << " resamples)\n";
        os << "- signif.: " << (r.component_significant ? "yes" : "no") << "\n";
    }

    return os;
}

// pretty printer for a vector of Result (components)
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

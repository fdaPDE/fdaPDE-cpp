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


#ifndef __FDAPDE_RGCCA_BLOCKS_H__
#define __FDAPDE_RGCCA_BLOCKS_H__

#include "gcv.h"
#include "linear_algebra.h"
#include "sampling.h"
#include "fdaPDE/src/solvers/nonnegative_weight_ipopt.h"

namespace fdapde {
namespace rgcca {
namespace internals {

// blocks adapt raw data matrices to the rgcca weight/component update api
template <typename SamplingStrategy>
class BaseBlock {
public:
    using Matrix = ::fdapde::rgcca::Matrix;
    using Vector = ::fdapde::rgcca::Vector;
    using IndexVector = ::fdapde::rgcca::IndexVector;
    using SparseMatrix = ::fdapde::rgcca::SparseMatrix;
    using BinaryMatrixT = BinaryMatrix<Dynamic, Dynamic>;
    using SparseSolver = ::fdapde::internals::eigen_sparse_solver_movable_wrap<Eigen::SimplicialLDLT<SparseMatrix>>;
    using ComponentsSolverType = typename std::decay_t<SamplingStrategy>::solver_t;

    template<typename S = SamplingStrategy>
    requires std::same_as<S, ::fdapde::rgcca::IndependentSampling>
    BaseBlock(const std::string& block_name, Matrix* data_ptr, const int n_dofs_weights) :
        block_name_(block_name), data_ptr_(data_ptr), components_solver_(data_ptr->rows()), n_dofs_weights_(n_dofs_weights) {
        // init components solver
        init_identity_row_index_();
        components_solver_.analyze_data();
    }

    template<typename S = SamplingStrategy>
    requires std::same_as<SamplingStrategy, ::fdapde::rgcca::TimeDependentSampling>
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
        objective_sign_invariant_(other.objective_sign_invariant_),
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

        nn_weights_solver_.reset();
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
    bool raw_score_cache_for_weight(const Vector& w, Vector& out) const {
        if (!raw_score_cache_ready_ || !weights_ready_ || w.size() != weights_.rows())
            return false;

        const auto current_weight = weights_.col(h_);
        if (w.isApprox(current_weight)) {
            out = raw_score_cache_;
            return true;
        }
        if (w.isApprox(-current_weight)) {
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
    void set_mode(const ::fdapde::rgcca::Mode mode) {
        mode_ = mode;
        if (mode_ == ::fdapde::rgcca::Mode::CorMax) tau_ = 0.0;
        if (mode_ == ::fdapde::rgcca::Mode::CovMax) tau_ = 1.0;
        if (mode_ == ::fdapde::rgcca::Mode::Regularized) select_tau_auto_();
        invalidate_M_();
    }
    [[nodiscard]] ::fdapde::rgcca::Mode mode() const { return mode_; }

    // normalization matrix
    [[nodiscard]] const SparseMatrix& M() const { ensure_M_(); return M_; }
    [[nodiscard]] const Matrix& ginvM() const { ensure_ginvM_(); return ginvM_; }
    [[nodiscard]] SparseSolver& invM() { ensure_invM_(); return invM_; }

    // components regularization utilities
    void set_lambda_components(const double lambda) {
        if constexpr (std::same_as<SamplingStrategy, ::fdapde::rgcca::IndependentSampling>) return;
        lambda_components_ = lambda;
        lambda_components_selection_ = lambda < 0.0;
    }
    [[nodiscard]] double lambda_components() const {
        if constexpr (std::same_as<SamplingStrategy, ::fdapde::rgcca::IndependentSampling>) return std::numeric_limits<double>::quiet_NaN();
        if (lambda_components_.has_value() && *lambda_components_ > 0.0) return *lambda_components_;
        return std::numeric_limits<double>::quiet_NaN();
    }
    void set_components_gcv_config(const ::fdapde::internals::GCVConfig& cfg) {
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
        Vector w = Vector::Ones(n_dofs_weights_);
        double norm2 = w.dot(Omega() * w);
        if (norm2 <= 0) norm2 = 1.0;
        weights_.col(h()) = w / std::sqrt(norm2);
        weights_star_computed_until_ = std::min(weights_star_computed_until_, h() - 1);
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

        if constexpr (std::same_as<SamplingStrategy, ::fdapde::rgcca::TimeDependentSampling>) {
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
        const Vector w = w_fit_(nu_D);
        weights().col(h()) = w;
        weights_star_computed_until_ = std::min(weights_star_computed_until_, h() - 1);

        // compute the component and regularize it (the regularization acts only in the TimeDependent sampling scenario)
        const Vector s = data_times_(Psi_D() * w);
        components().col(h()) = c_fit_(s);
    }

    // deflation
    void deflate(const ::fdapde::rgcca::Deflation mode) {
        if (h() == n_comp()) return;
        switch (mode) {
            case ::fdapde::rgcca::Deflation::Scores: deflate_scores_(); break;
            case ::fdapde::rgcca::Deflation::None: default: break;
        }
        invalidate_M_();
    }

    // weights post-processing
    void compute_weights_star(const int h) {
        ensure_lc_();
        if (h < 0 || h >= n_comp_)
            throw std::out_of_range("weights_star component index");

        if (h > weights_star_computed_until_ + 1) {
            compute_weights_star();
            return;
        }

        Vector w_star = weights_.col(h);
        if (h > 0) {
            const Matrix W_prev = weights_star_.leftCols(h);
            const Matrix P_prev = deflation_projections_.leftCols(h);
            const Vector coeff = P_prev.transpose() * (Psi_D() * weights_.col(h));
            w_star.noalias() -= W_prev * coeff;
        }
        weights_star_.col(h) = w_star;
        weights_star_computed_until_ = std::max(weights_star_computed_until_, h);
    }
    void compute_weights_star() {
        ensure_lc_();
        weights_star_.setZero(weights_.rows(), weights_.cols());
        weights_star_computed_until_ = -1;
        for (int h = 0; h < n_comp_; ++h)
            compute_weights_star(h);
    }

    // model evaluation
    [[nodiscard]] Vector normalized_component_for_evaluation(const Vector& w, const int hh) {
        if (hh != h())
            throw std::logic_error("normalized_component_for_evaluation: requested component differs from current h; M may refer to current deflated data.");

        // weight at locations
        const Vector wm = Psi_D() * w;

        // normalization
        // M(), not Omega()! in the evaluation, the regularization term should not be taken into account
        double nrm2 = wm.dot(M() * wm);
        if (nrm2 <= 0.0 || !std::isfinite(nrm2)) nrm2 = 1.0;

        // components
        const Vector s = data_times_(wm) / std::sqrt(nrm2);
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
    void set_weight_sign_constraint(const ::fdapde::rgcca::WeightSignConstraint weight_sign_constraint = ::fdapde::rgcca::WeightSignConstraint::None) {
        weight_sign_constraint_ = weight_sign_constraint;
        reset_nonnegative_weight_solver_();
    }
    void set_objective_sign_invariant(const bool value) {
        objective_sign_invariant_ = value;
        reset_nonnegative_weight_solver_();
    }

    // observers
    [[nodiscard]] const SparseMatrix& Psi_T() const { return components_solver_.Psi(); }
    Matrix& weights() { ensure_lc_(); return weights_; }
    Matrix weights_m() { ensure_lc_(); return Psi_D() * weights_; }
    Matrix& weights_star() { ensure_lc_(); return weights_star_; }
    Matrix weights_star_m() { ensure_lc_(); return Psi_D() * weights_star_; }
    Matrix& components() { ensure_lc_(); return components_; }
    Matrix components_m() {
        ensure_lc_();
        if constexpr (std::same_as<SamplingStrategy, ::fdapde::rgcca::IndependentSampling>) return components_;
        return Psi_T() * components_;
    }
    [[nodiscard]] ::fdapde::rgcca::WeightSignConstraint weight_sign_constraint() const { return weight_sign_constraint_; }
    template <typename S = SamplingStrategy> requires std::same_as<S, ::fdapde::rgcca::TimeDependentSampling> const Vector& times() { return times_; }
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
            weights_star_computed_until_ = -1;
            weights_ready_ = true;
        }
        if (!components_ready_) {
            components_.setZero(components_solver_.n_dofs(), n_comp_);
            components_ready_ = true;
        }
    }
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

        if (mode_ == ::fdapde::rgcca::Mode::CovMax) {
            M_.setIdentity();
        } else {
            const double n_d = static_cast<double>(n());
            const double den = bias_ ? n_d : std::max(1.0, n_d - 1.0);

            const Matrix XtX = data().transpose() * data();
            const Vector mu = data().colwise().mean();
            const Matrix Sigma = (XtX - n_d * (mu * mu.transpose())) / den;

            // const Matrix XtX = data().transpose() * data();
            // const Matrix Sigma = XtX / den;

            if (mode_ == ::fdapde::rgcca::Mode::CorMax) {
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
        if (mode_ == ::fdapde::rgcca::Mode::CovMax) {
            ginvM_.setIdentity(m(), m());
        } else {
            ginvM_.resize(m(), m());
            Eigen::MatrixXd M_dense = Eigen::MatrixXd(M());
            ::fdapde::internals::ginv(M_dense, ginvM_);
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
        if (mode_ != ::fdapde::rgcca::Mode::CovMax) invalidate_M_();
    }
    virtual void invalidate_derived_caches_() {}

    // deflation
    void deflate_scores_() {
        ensure_lc_();

        if (!is_identity_row_index_()) throw std::logic_error("deflate_scores_: raw-data deflation requires identity row_index");

        const Vector y = [&]() {
            if constexpr (std::same_as<SamplingStrategy, ::fdapde::rgcca::IndependentSampling>) return eta_();
            else return Psi_T() * eta_();
        }();
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
    Vector solve_nonnegative_weight_ipopt_(const Vector& z, const bool use_closed_form_solution = false) {
        if (!nn_weights_solver_)
            nn_weights_solver_ = std::make_unique<::fdapde::internals::NonNegativeWeightSolver>(
                Psi_D(),
                Omega(),
                objective_sign_invariant_,
                use_closed_form_solution
            );
        return nn_weights_solver_->solve(z);
    }
    void reset_nonnegative_weight_solver_() {
        nn_weights_solver_.reset();
    }

    // weights solver
    virtual Vector w_fit_(const Vector& nu) = 0;
    [[nodiscard]] Vector normalize_weight_(const Vector& w) {
        const double rho2 = w.dot(Omega() * w);
        if (rho2 <= 0.0 || !std::isfinite(rho2)) return w;
        return w / std::sqrt(rho2);
    }

    // component solver
    Vector c_fit_(const Vector& s) {
        if constexpr (std::same_as<SamplingStrategy, ::fdapde::rgcca::IndependentSampling>) return s;

        components_solver_.update_response_and_weights(s, I_);

        double lambda = 1e-15;
        if (lambda_components_selection_) {
            auto [success, l] = ::fdapde::internals::select_lambda_with_gcv(components_solver_, components_gcv_cfg_);
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
    std::unique_ptr<::fdapde::internals::NonNegativeWeightSolver> nn_weights_solver_;

    // options
    double tau_ {0.0};
    ::fdapde::rgcca::Mode mode_ = ::fdapde::rgcca::Mode::CorMax;
    ::fdapde::rgcca::WeightSignConstraint weight_sign_constraint_ = ::fdapde::rgcca::WeightSignConstraint::None;
    bool objective_sign_invariant_ = true;
    bool bias_ = true;

    // parameters
    ::fdapde::internals::GCVConfig components_gcv_cfg_;
    std::optional<double> lambda_components_;

    // results
    Matrix weights_, weights_star_, components_;
    Matrix deflation_projections_;
    int weights_star_computed_until_ {-1};

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

    using Base::M;
    using Base::n_dofs_weights;

    template<typename S = SamplingStrategy>
    requires std::same_as<S, ::fdapde::rgcca::IndependentSampling>
    MultivariateBlock(const std::string& block_name, Matrix* data_ptr) :
        Base(block_name, data_ptr, static_cast<int>(data_ptr->cols())) {
        init_multivariate_();
    }

    template<typename S = SamplingStrategy>
    requires std::same_as<SamplingStrategy, ::fdapde::rgcca::TimeDependentSampling>
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
    using Base::init;
    using Base::ginvM;
    using Base::n;
    using Base::m;
    using Base::mode;
    using Base::solve_nonnegative_weight_ipopt_;
    using Base::reset_nonnegative_weight_solver_;
    using Base::normalize_weight_;
    using Base::data_transpose_times_;
    using Base::weight_sign_constraint;

    void init_multivariate_() {
        Psi_D_.resize(m(), n_dofs_weights()); // n_dofs_weights == m in this case
        Psi_D_.setIdentity();
        init();
    }

    Vector w_fit_(const Vector& nu) override {
        assert(nu.size() == n() && "nu must have size n (rows of X)");
        init();

        Vector z = data_transpose_times_(nu);

        if (weight_sign_constraint() == ::fdapde::rgcca::WeightSignConstraint::NonNegative) {
            return solve_nonnegative_weight_ipopt_(z, mode() == ::fdapde::rgcca::Mode::CovMax); // already normalized
        }

        if (mode() == ::fdapde::rgcca::Mode::CovMax) return normalize_weight_(z);
        const Vector w_tilde = ginvM() * z;
        return normalize_weight_(w_tilde);
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

    using Base::M;
    using Base::n_dofs_weights;

    template <typename GeoFrame>
    requires std::same_as<SamplingStrategy, ::fdapde::rgcca::IndependentSampling>
    FunctionalBlock(const std::string& block_name, GeoFrame& gf, Matrix* data_ptr, WeightsPenaltyType&& weights_penalty) :
        Base(block_name, data_ptr, weights_penalty.get().bilinear_form().n_dofs()) {
        init_functional_(gf, std::forward<WeightsPenaltyType>(weights_penalty));
    }

    template <typename GeoFrame>
    requires std::same_as<SamplingStrategy, ::fdapde::rgcca::TimeDependentSampling>
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
        ::fdapde::rgcca::internals::validate_positive_regularization_lambda(lambda, "weight lambda");
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
    using Base::init;
    using Base::n;
    using Base::solve_nonnegative_weight_ipopt_;
    using Base::reset_nonnegative_weight_solver_;
    using Base::data_transpose_times_;
    using Base::weight_sign_constraint;

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

        if (weight_sign_constraint() == ::fdapde::rgcca::WeightSignConstraint::NonNegative) {
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

// factory helpers keep RGCCA::add_* independent from concrete block classes
template <typename SamplingStrategy, typename Matrix = ::fdapde::rgcca::Matrix>
requires std::same_as<SamplingStrategy, ::fdapde::rgcca::IndependentSampling>
inline std::unique_ptr<BaseBlock<SamplingStrategy>>
make_multivariate_block(std::string block_name, Matrix* data_ptr) {
    return std::make_unique<MultivariateBlock<SamplingStrategy>>(block_name, data_ptr);
}

template <typename SamplingStrategy, typename Matrix = ::fdapde::rgcca::Matrix, typename Vector = ::fdapde::rgcca::Vector>
requires std::same_as<SamplingStrategy, ::fdapde::rgcca::TimeDependentSampling>
inline std::unique_ptr<BaseBlock<SamplingStrategy>>
make_multivariate_block(std::string block_name, const Triangulation<1, 1>& T, const Vector& times, Matrix* data_ptr) {
    return std::make_unique<MultivariateBlock<SamplingStrategy>>(block_name, T, times, data_ptr);
}

template <typename SamplingStrategy, typename GeoFrame, typename WeightsPenaltyType, typename Matrix = ::fdapde::rgcca::Matrix>
requires std::same_as<SamplingStrategy, ::fdapde::rgcca::IndependentSampling>
std::unique_ptr<BaseBlock<SamplingStrategy>>
make_functional_block(std::string block_name, GeoFrame& gf, Matrix* data_ptr, WeightsPenaltyType&& weights_penalty) {
    return std::make_unique<FunctionalBlock<WeightsPenaltyType, SamplingStrategy>>(block_name, gf, data_ptr, std::forward<WeightsPenaltyType>(weights_penalty));
}

template <typename SamplingStrategy, typename GeoFrame, typename WeightsPenaltyType, typename Vector = ::fdapde::rgcca::Vector, typename Matrix = ::fdapde::rgcca::Matrix>
requires std::same_as<SamplingStrategy, ::fdapde::rgcca::TimeDependentSampling>
std::unique_ptr<BaseBlock<SamplingStrategy>>
make_functional_block(std::string block_name, const Triangulation<1, 1>& T, const Vector& times, GeoFrame& gf, Matrix* data_ptr, WeightsPenaltyType&& weights_penalty) {
    return std::make_unique<FunctionalBlock<WeightsPenaltyType, SamplingStrategy>>(block_name, T, times, gf, data_ptr,  std::forward<WeightsPenaltyType>(weights_penalty));
}

} // namespace internals
} // namespace rgcca
} // namespace fdapde

#endif // __FDAPDE_RGCCA_BLOCKS_H__

#ifdef __FDAPDE_RGCCA_DEFINE_MODEL_BLOCKS__
#ifndef __FDAPDE_RGCCA_MODEL_BLOCKS_H__
#define __FDAPDE_RGCCA_MODEL_BLOCKS_H__

namespace fdapde {

// registers a block and synchronizes model-level options into it
template <typename SamplingStrategy>
int RGCCA<SamplingStrategy>::add_block(typename RGCCA<SamplingStrategy>::BlockPtr b) {
    if (!b) throw std::invalid_argument("RGCCA/add_block: null block");

    // sampling-specific consistency
    if constexpr (std::same_as<SamplingStrategy, rgcca::IndependentSampling>) {
        if (b->n() != n()) throw std::invalid_argument("RGCCA/add_block: n mismatch");
    } else {
        add_times_(b->times());
    }

    // block state owned by the model
    b->set_bias(opt_.bias);
    b->set_raw_data_mutable(true);
    b->set_mode(opt_.mode);
    b->set_weight_sign_constraint(opt_.weight_sign_constraint);
    b->set_objective_sign_invariant(opt_.scheme.sign_invariant);
    b->set_n_comp(n_comp());
    blocks_.emplace_back(std::move(b));
    initialized_ = false;
    return ++J_;
}

// adds a multivariate block for independent sampling
template <typename SamplingStrategy>
template <typename S>
requires std::same_as<S, rgcca::IndependentSampling>
int RGCCA<SamplingStrategy>::add_multivariate_block(std::string block_name, rgcca::Matrix&& X) {
    data_blocks_.push_back(std::make_unique<rgcca::Matrix>(std::move(X)));
    return add_block(
        rgcca::internals::make_multivariate_block<SamplingStrategy>(
            block_name, data_blocks_.back().get()
        )
    );
}

// adds a multivariate block for time-dependent sampling
template <typename SamplingStrategy>
template <typename S>
requires std::same_as<S, rgcca::TimeDependentSampling>
int RGCCA<SamplingStrategy>::add_multivariate_block(
    std::string block_name,
    const rgcca::Vector& times,
    rgcca::Matrix&& X
) {
    data_blocks_.push_back(std::make_unique<rgcca::Matrix>(std::move(X)));
    return add_block(
        rgcca::internals::make_multivariate_block<SamplingStrategy>(
            block_name, T_, times, data_blocks_.back().get()
        )
    );
}

// adds a functional block for independent sampling
template <typename SamplingStrategy>
template <typename GeoFrame, typename WeightsPenaltyType>
requires std::same_as<SamplingStrategy, rgcca::IndependentSampling>
int RGCCA<SamplingStrategy>::add_functional_block(
    std::string block_name,
    const GeoFrame& gf,
    rgcca::Matrix&& X,
    WeightsPenaltyType&& weights_penalty
) {
    data_blocks_.push_back(std::make_unique<rgcca::Matrix>(std::move(X)));
    return add_block(
        rgcca::internals::make_functional_block<SamplingStrategy>(
            block_name,
            gf,
            data_blocks_.back().get(),
            std::forward<WeightsPenaltyType>(weights_penalty)
        )
    );
}

// adds a functional block for time-dependent sampling
template <typename SamplingStrategy>
template <typename GeoFrame, typename WeightsPenaltyType>
requires std::same_as<SamplingStrategy, rgcca::TimeDependentSampling>
int RGCCA<SamplingStrategy>::add_functional_block(
    std::string block_name,
    const rgcca::Vector& times,
    const GeoFrame& gf,
    rgcca::Matrix&& X,
    WeightsPenaltyType&& weights_penalty
) {
    data_blocks_.push_back(std::make_unique<rgcca::Matrix>(std::move(X)));
    return add_block(
        rgcca::internals::make_functional_block<SamplingStrategy>(
            block_name,
            T_,
            times,
            gf,
            data_blocks_.back().get(),
            std::forward<WeightsPenaltyType>(weights_penalty)
        )
    );
}

// returns non-owning references to the model blocks
template <typename SamplingStrategy>
auto RGCCA<SamplingStrategy>::main_blocks_() const
    -> typename RGCCA<SamplingStrategy>::BlockRefList {
    typename RGCCA<SamplingStrategy>::BlockRefList out;
    out.reserve(blocks_.size());
    for (const auto& b : blocks_)
        out.push_back(b.get());
    return out;
}

// clones blocks for bootstrap workers
template <typename SamplingStrategy>
auto RGCCA<SamplingStrategy>::clone_blocks_() const
    -> typename RGCCA<SamplingStrategy>::BootstrapBlocks {
    typename RGCCA<SamplingStrategy>::BootstrapBlocks out;
    out.owners.reserve(blocks_.size());
    out.refs.reserve(blocks_.size());

    for (const auto& b : blocks_) {
        // worker clones must not mutate shared raw data
        auto copy = b->clone();
        copy->set_raw_data_mutable(false);

        out.refs.push_back(copy.get());
        out.owners.push_back(std::move(copy));
    }

    return out;
}

// returns block names in model order
template <typename SamplingStrategy>
std::vector<std::string> RGCCA<SamplingStrategy>::block_names_(
    const typename RGCCA<SamplingStrategy>::BlockRefList& blocks
) const {
    std::vector<std::string> out;
    out.reserve(blocks.size());

    for (auto* b : blocks)
        out.push_back(b->name());

    return out;
}

// returns block weight dimensions in model order
template <typename SamplingStrategy>
std::vector<int> RGCCA<SamplingStrategy>::block_dims_(
    const typename RGCCA<SamplingStrategy>::BlockRefList& blocks
) const {
    std::vector<int> out;
    out.reserve(blocks.size());

    for (auto* b : blocks)
        out.push_back(b->n_dofs_weights());

    return out;
}

} // namespace fdapde

#endif // __FDAPDE_RGCCA_MODEL_BLOCKS_H__
#endif // __FDAPDE_RGCCA_DEFINE_MODEL_BLOCKS__

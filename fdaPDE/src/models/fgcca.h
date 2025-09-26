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

namespace internals {

class BaseBlock;   // forward decl for operator<<
std::ostream& operator<<(std::ostream& os, const BaseBlock& b);

class BaseBlock {
   public:
    using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using SparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
    using SparseSolver = eigen_sparse_solver_movable_wrap<Eigen::SimplicialLDLT<SparseMatrix>>;

    BaseBlock(std::string block_name, Matrix data, int n_nodes, double tau = 0.0, int n_comp = 1) :
        block_name_(std::move(block_name)), data_(std::move(data)), n_nodes_(n_nodes), tau_(tau), n_comp_(n_comp) { }

    virtual ~BaseBlock() = default;

    // ---- Uniform public API ----
    [[nodiscard]] const std::string& name() const { return block_name_; }
    [[nodiscard]] const Matrix& data() const { return data_; }
    Matrix& data() {
        invalidate_sigma_();
        return data_;
    }

    [[nodiscard]] int n_obs() const { return static_cast<int>(data_.rows()); }
    [[nodiscard]] int n_covs() const { return static_cast<int>(data_.cols()); }
    [[nodiscard]] int n_nodes() const { return n_nodes_; }
    [[nodiscard]] int n_comp() const { return n_comp_; }

    [[nodiscard]] int h() const { return h_; }
    void set_h(int idx) {
        if (idx < 0 || idx >= n_comp_) throw std::out_of_range("h");
        h_ = idx;
    }
    void next_component() { set_h(h_ + 1); }

    [[nodiscard]] double tau() const { return tau_; }
    void set_tau(double tau) {
        tau_ = tau;
        invalidate_sigma_();
    }

    [[nodiscard]] const SparseMatrix& Sigma() const {
        ensure_sigma_();
        return Sigma_;
    }
    SparseSolver& invSigma() {
        ensure_sigma_();
        return invSigma_;
    }

    Matrix& loadings() {
        ensure_lc_();
        return loadings_;
    }
    Matrix loadings_m() {
        ensure_lc_();
        return Psi() * loadings_;
    }
    Matrix& components() {
        ensure_lc_();
        return components_;
    }

    void init() {
        ensure_sigma_();
        ensure_lc_();
    }

    // ---- Clean virtual interface ----
    [[nodiscard]] virtual const SparseMatrix& Psi() const = 0;
    virtual void l_compute(const Vector& nu) = 0;

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
    }

    // Hyperparameters (default no-op). The model can call this on all blocks.
    virtual void set_lambda(double) { }   // default: ignored
   protected:
    void compute_sigma_() {
        const int n = n_covs();
        SparseMatrix I(n, n);
        I.setIdentity();
        const Matrix dense = ((1.0 - tau_) / static_cast<double>(n)) * (data_.transpose() * data_);
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

    // state
    std::string block_name_;
    Matrix data_;   // n_obs x n_covs
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
        const Vector a_m = invSigma().solve(data().transpose() * nu);
        const double norm = std::sqrt(a_m.dot(Sigma() * a_m));
        loadings().col(h()) = a_m / norm;
        components().col(h()) = data() * loadings_m().col(h());
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
template <class TriangulationType, class PenaltyType> class FunctionalBlock final : public BaseBlock {
   public:
    using Base = BaseBlock;
    using Vector = Base::Vector;
    using Matrix = Base::Matrix;
    using SparseMatrix = Base::SparseMatrix;
    using lSolverType = typename std::decay_t<PenaltyType>::solver_t;

    FunctionalBlock(
      std::string block_name, GeoFrame<TriangulationType>& gf, PenaltyType&& penalty, const double tau = 0.0,
      const int n_comp = 1) :
        Base(
          block_name, gf[0].template col<double>(block_name).as_matrix().transpose(),
          gf.template triangulation<0>().n_nodes(), tau, n_comp),
        gf_(gf) {
        gf_[0].add_column("z", Vector::Zero(n_covs()));
        solver_.discretize(penalty.get());
        solver_.analyze_data("z ~ f", gf_, Sigma());
        init();
    }

    [[nodiscard]] const SparseMatrix& Psi() const override { return solver_.Psi(); }

    // The model must set lambda before calling l_compute
    void set_lambda(double lambda) override { lambda_ = lambda; }

    void l_compute(const Vector& nu) override {
        assert(nu.size() == n_obs() && "nu must have size n_obs (rows of X)");
        init();
        if (!lambda_.has_value())
            throw std::logic_error("FunctionalBlock: lambda not set. Call set_lambda() before l_compute().");
        gf_[0].template col<double>("z") = invSigma().solve(data().transpose() * nu);
        solver_.analyze_data("z ~ f", gf_, Sigma());
        solver_.fit(*lambda_);

        const Vector f = solver_.f();
        const Vector a_m = Psi() * f;
        const double norm = std::sqrt(a_m.dot(Sigma() * a_m));

        loadings().col(h()) = f / norm;
        components().col(h()) = data() * (a_m / norm);
    }

    // print override
    void print(std::ostream& os) const override {
        BaseBlock::print(os);
        os << "type: FunctionalBlock, n_nodes = " << n_nodes();
        if (lambda_.has_value()) os << ", lambda = " << *lambda_;
        os << "\n";
    }
   private:
    GeoFrame<TriangulationType>& gf_;
    lSolverType solver_;
    std::optional<double> lambda_;   // owned by the block (set by the model)
};

}   // namespace internals

inline std::unique_ptr<internals::BaseBlock>
make_multivariate_block(std::string name, Eigen::Matrix<double, Dynamic, Dynamic>& data,
                      double tau = 0.0, int n_comp = 1) {
    return std::make_unique<internals::MultivariateBlock>(
        std::move(name), data, tau, n_comp);
}

template<typename TriangulationType, typename PenaltyType>
std::unique_ptr<internals::BaseBlock>
make_functional_block(std::string name, GeoFrame<TriangulationType>& gf, PenaltyType&& pen,
                      double tau = 0.0, int n_comp = 1) {
    return std::make_unique<internals::FunctionalBlock<TriangulationType, PenaltyType>>(
        std::move(name), gf, std::forward<PenaltyType>(pen), tau, n_comp);
}

}


#endif   // __FGCCA_H__
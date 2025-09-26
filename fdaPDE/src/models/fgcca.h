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

template <typename Derived>
class BaseBlock {
public:
    using Matrix       = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector       = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using SparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
    using SparseSolver = eigen_sparse_solver_movable_wrap<Eigen::SimplicialLDLT<SparseMatrix>>;
    // using DenseSolver  = Eigen::PartialPivLU<Matrix>;

    BaseBlock(std::string block_name, const double tau = 0., const int n_comp = 1) : block_name_(std::move(block_name)), tau_(tau), n_comp_(n_comp), loadings_ready_(false) { }

    // unified public API

    // data
    Matrix data() const { return derived().data(); }

    // block name
    std::string name() const { return block_name_; }

    // dimensions
    int n_nodes() const { return derived().n_nodes(); } // (block type dependent)
    int n_obs()   const { return data().rows(); }
    int n_stat()  const { return data().cols(); }
    int n_comp()  const { return n_comp_; }

    // current component index
    int h() const { return h_; }

    // shrinkage coefficient
    double tau() const { return tau_; }

    // Psi matrix (block type dependent)
    const SparseMatrix& Psi() const { return derived().Psi(); }

    // covariance matrix
    const SparseMatrix& Sigma() const { return Sigma_; }
    SparseSolver& invSigma() { return invSigma_; }

    // initialize component
    void init() {
        if (!component_initialized_){
            component_initialized_ = true;
            compute_sigma();
        }
    }

    // compute loadings (block type dependent)
    template <typename... Args>
    void l_compute(Args&&... args) {
        derived().l_compute(std::forward<Args>(args)...);
    }

    // loadings vector/matrix sized to n_nodes() x 1 by default
    Matrix& loadings()   { ensure_l_and_c(); return loadings_; }
    Matrix  loadings_m() { ensure_l_and_c(); return Psi() * loadings_; }
    Matrix& components()   { ensure_l_and_c(); return components_; }

    // debug printing operator
    friend std::ostream& operator<<(std::ostream& os, const BaseBlock& block) {
        const Matrix X = block.data();

        using Index = Eigen::Index;
        const Index max_rows = 3;
        const Index max_cols = 5;

        const Index rows = std::min<Index>(max_rows, X.rows());
        const Index cols = std::min<Index>(max_cols, X.cols());

        os << block.name() << " Block preview (" << X.rows() << " x " << X.cols() << "):\n";
        for (Index i = 0; i < rows; ++i) {
            for (Index j = 0; j < cols; ++j) {
                os << X(i, j);
                if (j + 1 < cols) os << '\t';
            }
            if (X.cols() > cols) os << "\t...";
            os << '\n';
        }
        if (X.rows() > rows) os << "...\n";
        return os;
    }

protected:

    void compute_sigma() {
        const int n = n_obs();
        SparseMatrix I(n, n);
        I.setIdentity();
        Sigma_.resize(n, n);
        Matrix temp((1-tau()) / n * data() * data().transpose());
        Sigma_ = temp.sparseView(1e-12) + tau() * I;
        Sigma_.makeCompressed();
        invSigma_.compute(Sigma_);
    }

private:
    const Derived& derived() const { return static_cast<const Derived&>(*this); }
    Derived& derived() { return static_cast<Derived&>(*this); }

    void ensure_l_and_c() {
        if (!loadings_ready_) {
            loadings_.resize(n_nodes(), n_comp());
            loadings_.setZero();
            loadings_ready_ = true;
        }
        if (!components_ready_) {
            components_.resize(n_stat(), n_comp());
            components_.setZero();
            components_ready_ = true;
        }
    }

    // common state
    std::string block_name_;
    double tau_ = 0.;
    int n_comp_ = 1;
    int h_ = 0; // current component
    Matrix components_;
    Matrix loadings_;
    bool loadings_ready_ = false;
    bool components_ready_ = false;
    bool component_initialized_ = false;
    SparseMatrix Sigma_;
    SparseSolver invSigma_;
};



class MultivariateBlock : public BaseBlock<MultivariateBlock> {
public:
    using Base         = BaseBlock<MultivariateBlock>;
    using Matrix       = typename Base::Matrix;
    using SparseMatrix = typename Base::SparseMatrix;

    explicit MultivariateBlock(std::string block_name, const Matrix &data, const double tau = 0., const int n_comp = 1)
        : Base(std::move(block_name), tau, n_comp), data_(data) {

        // initialize first component
        init();

        // initialize the Psi matrix (at the identity)
        const int n = n_obs();
        Psi_.resize(n, n);
        Psi_.setIdentity();
    }

public:
    using Base::loadings;
    using Base::loadings_m;
    using Base::components;

    // ---- required by Base (CRTP) ----
    const Matrix& data() const { return data_; }
    Matrix& data() { return data_; }
    int n_nodes() const { return Base::n_obs(); }
    const SparseMatrix& Psi() const { return Psi_; }
    void l_compute(const Vector& nu) {
        init();
        const Vector a_m(invSigma().solve(data() * nu));
        loadings().col(h()) = a_m / sqrt(a_m.dot( Sigma() * a_m));
        components().col(h()) = data().transpose() * loadings_m().col(h()) ;
    }

private:
    Matrix data_; // shape: n_obs x n_stat
    SparseMatrix Psi_;
    using Base::init;
    using Base::h;
    using Base::invSigma;
};


template <typename TriangulationType, typename Penalty>
class FunctionalBlock : public BaseBlock<FunctionalBlock<TriangulationType, Penalty>> {
public:
    using Base         = BaseBlock<FunctionalBlock<TriangulationType, Penalty>>;
    using Vector       = typename Base::Vector;
    using Matrix       = typename Base::Matrix;
    using SparseMatrix = typename Base::SparseMatrix;
    using lSolverType  = typename std::decay_t<Penalty>::solver_t;

    FunctionalBlock(std::string block_name, GeoFrame<TriangulationType>& gf, Penalty&& penalty, const double tau = 0., const int n_comp = 1)
        : Base(block_name, tau, n_comp), gf_(gf), l_solver_(){

        // add the inner-component column to gf_
        gf_[0].add_column("z", Vector::Zero(Base::n_obs()) );

        // initialize first component
        init();

        // solver
        l_solver_.discretize(penalty.get());
        l_solver_.analyze_data("z ~ f", gf_, Sigma());
    }

    // Base methods
    using Base::loadings;
    using Base::loadings_m;
    using Base::components;
    using Base::Sigma;

    // ---- required by Base (CRTP) ----
    const Matrix data() const { return gf_[0].template col<double>(name()).as_matrix(); }
    Matrix data() { return gf_[0].template col<double>(name()).as_matrix(); }
    int n_nodes() const { return gf_.template triangulation<0>().n_nodes(); }
    const SparseMatrix& Psi() const { return l_solver_.Psi(); }
    void l_compute(const Vector& nu, const double lambda) {
        init();
        gf_[0].template col<double>("z") = invSigma().solve(data() * nu);
        l_solver_.analyze_data("z ~ f", gf_, Sigma());
        l_solver_.fit(lambda);
        const Vector a_m = Psi() * l_solver_.f();
        const double norm = sqrt(a_m.dot( Sigma() * a_m));
        loadings().col(h()) = l_solver_.f() / norm;
        components().col(Base::h()) = data().transpose() * a_m / norm;
    }

private:
    GeoFrame<TriangulationType>& gf_;
    lSolverType l_solver_;
    using Base::name;
    using Base::h;
    using Base::init;
    using Base::invSigma;
};


}

}


#endif   // __FGCCA_H__
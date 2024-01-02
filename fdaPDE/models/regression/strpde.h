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

#ifndef __STRPDE_H__
#define __STRPDE_H__

#include <fdaPDE/linear_algebra.h>
#include <fdaPDE/pde.h>
#include <fdaPDE/utils.h>

#include <memory>
#include <type_traits>

#include "../model_base.h"
#include "../model_macros.h"
#include "../model_traits.h"
#include "../sampling_design.h"
#include "regression_base.h"
using fdapde::core::BlockVector;
using fdapde::core::Kronecker;
using fdapde::core::PDEBase;
using fdapde::core::SMW;
using fdapde::core::SparseBlockMatrix;
using fdapde::core::SparseKroneckerTensorProduct;
using fdapde::core::SplineBasis;

namespace fdapde {
namespace models {

// base class for STRPDE model
template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver> class STRPDE;

// implementation of STRPDE for separable space-time regularization
template <typename PDE, typename SamplingDesign>
class STRPDE<PDE, SpaceTimeSeparable, SamplingDesign, MonolithicSolver> :
    public RegressionBase<STRPDE<PDE, SpaceTimeSeparable, SamplingDesign, MonolithicSolver>> {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
   private:
    typedef SpaceTimeSeparable RegularizationType;
    typedef RegressionBase<STRPDE<PDE, RegularizationType, SamplingDesign, MonolithicSolver>> Base;
    SparseBlockMatrix<double, 2, 2> A_ {};      // system matrix of non-parametric problem (2N x 2N matrix)
    fdapde::SparseLU<SpMatrix<double>> invA_;   // factorization of matrix A
    DVector<double> b_ {};                      // right hand side of problem's linear system (1 x 2N vector)
    DMatrix<double> T_;                         // T = \Psi^T*Q*\Psi + \lambda*R
    SpMatrix<double> K_;                        // P1 \kron R0
   public:
    // import commonly defined symbols from base
    IMPORT_REGRESSION_SYMBOLS;
    using Base::lambda_D;   // smoothing parameter in space
    using Base::lambda_T;   // smoothing parameter in time
    using Base::P1;         // time penalization matrix: [P1_]_{ij} = \int_{[0,T]} (\phi_i)_tt*(\phi_j)_tt
    using Base::P0;         // time mass matrix: [P0_]_{ij} = \int_{[0,T]} \phi_i*\phi_j
    using Base::P;          // discretized penalty: P = \lambda_D*((R1^T*R0^{-1}*R1) \kron Rt) + \lambda_T*(R0 \kron Pt)
    // constructor
    STRPDE() = default;
    STRPDE(const PDE& pde, const DMatrix<double>& time) : Base(pde, time) {};

    void init_model();          // update model object in case of **structural** changes in its definition
    void update_to_weights();   // update model object in case of changes in the weights matrix
    virtual void solve();       // finds a solution to the smoothing problem

    // GCV support
    const DMatrix<double>& T();  // T = \Psi^T*Q*\Psi + + \lambda_T*(Pt \kron R0) + \lambda_D*(R1^T*R0^{-1}*R1)
    double norm(const DMatrix<double>& op1, const DMatrix<double>& op2) const;   // euclidian norm of op1 - op2

    // getters
    const SparseBlockMatrix<double, 2, 2>& A() const { return A_; }
    const fdapde::SparseLU<SpMatrix<double>>& invA() const { return invA_; }

    virtual ~STRPDE() = default;
};

// model initialization. Called by ModelBase::init() as part of the initialization process.
// NB: a change in the smoothing parameter must trigger a re-initialization of the model
template <typename PDE, typename SamplingDesign>
void STRPDE<PDE, SpaceTimeSeparable, SamplingDesign, MonolithicSolver>::init_model() {
    // assemble system matrix for the nonparameteric part of the model
    if (is_empty(K_)) K_ = Kronecker(P1(), pde().R0());
    A_ = SparseBlockMatrix<double, 2, 2>(
      -PsiTD() * W() * Psi() - lambda_T() * K_, lambda_D() * R1().transpose(),
      lambda_D() * R1(),                        lambda_D() * R0()            );
    // cache system matrix for reuse
    invA_.compute(A_);
    // prepare rhs of linear system
    b_.resize(A_.rows());
    b_.block(A_.rows() / 2, 0, A_.rows() / 2, 1) = lambda_D() * u();
    return;
}

// updates model in case of a change in the weights matrix
template <typename PDE, typename SamplingDesign>
void STRPDE<PDE, SpaceTimeSeparable, SamplingDesign, MonolithicSolver>::update_to_weights() {
    // adjust north-west block of matrix A_ and factorize
    A_.block(0, 0) = -PsiTD() * W() * Psi() - lambda_T() * K_;
    invA_.compute(A_);
    return;
}

// finds a solution to the STR-PDE smoothing problem (separable penalization)
template <typename PDE, typename SamplingDesign>
void STRPDE<PDE, SpaceTimeSeparable, SamplingDesign, MonolithicSolver>::solve() {
    BLOCK_FRAME_SANITY_CHECKS;
    DVector<double> sol;   // room for problem' solution

    if (!Base::has_covariates()) {   // nonparametric case
        // update rhs of STR-PDE linear system
        b_.block(0, 0, A_.rows() / 2, 1) = -PsiTD() * W() * y();
        // solve linear system A_*x = b_
        sol = invA_.solve(b_);
        f_ = sol.head(A_.rows() / 2);
    } else {   // parametric case
        // update rhs of STR-PDE linear system
        b_.block(0, 0, A_.rows() / 2, 1) = -PsiTD() * lmbQ(y());   // -\Psi^T*D*Q*z

        // definition of matrices U and V  for application of woodbury formula
        U_ = DMatrix<double>::Zero(A_.rows(), q());
        U_.block(0, 0, A_.rows() / 2, q()) = PsiTD() * W() * X();
        V_ = DMatrix<double>::Zero(q(), A_.rows());
        V_.block(0, 0, q(), A_.rows() / 2) = X().transpose() * W() * Psi();
        // solve system (A_ + U_*(X^T*W_*X)*V_)x = b using woodbury formula from NLA module
        sol = SMW<>().solve(invA_, U_, XtWX(), V_, b_);
        // store result of smoothing
        f_ = sol.head(A_.rows() / 2);
        beta_ = invXtWX().solve(X().transpose() * W()) * (y() - Psi() * f_);
    }
    // store PDE misfit
    g_ = sol.tail(A_.rows() / 2);
    return;
}
  
template <typename PDE, typename SamplingDesign>
const DMatrix<double>& STRPDE<PDE, SpaceTimeSeparable, SamplingDesign, MonolithicSolver>::T() {
    if (!has_covariates())   // case without covariates, Q is the identity matrix
        T_ = PsiTD() * W() * Psi() + P();
    else   // general case with covariates
        T_ = PsiTD() * lmbQ(Psi()) + P();
    return T_;
};

template <typename PDE, typename SamplingDesign>
double STRPDE<PDE, SpaceTimeSeparable, SamplingDesign, MonolithicSolver>::norm(
  const DMatrix<double>& op1, const DMatrix<double>& op2) const {
    return (op1 - op2).squaredNorm();
}

template <typename PDE_, typename SamplingDesign_>
struct model_traits<STRPDE<PDE_, SpaceTimeSeparable, SamplingDesign_, MonolithicSolver>> {
    typedef PDE_ PDE;
    typedef SpaceTimeSeparable regularization;
    typedef SplineBasis<3> TimeBasis;   // use cubic B-splines
    typedef SamplingDesign_ sampling;
    typedef MonolithicSolver solver;
    enum { N = PDE::N, M = PDE::M, n_lambda = 2 };
};

// implementation of STRPDE for parabolic space-time regularization, monolithic solver
template <typename PDE, typename SamplingDesign>
class STRPDE<PDE, SpaceTimeParabolic, SamplingDesign, MonolithicSolver> :
    public RegressionBase<STRPDE<PDE, SpaceTimeParabolic, SamplingDesign, MonolithicSolver>> {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
   private:
    typedef SpaceTimeParabolic RegularizationType;
    typedef RegressionBase<STRPDE<PDE, SpaceTimeParabolic, SamplingDesign, MonolithicSolver>> Base;
    SparseBlockMatrix<double, 2, 2> A_ {};      // system matrix of non-parametric problem (2N x 2N matrix)
    fdapde::SparseLU<SpMatrix<double>> invA_;   // factorization of matrix A
    DVector<double> b_ {};                      // right hand side of problem's linear system (1 x 2N vector)

    SpMatrix<double> L_;   // L \kron R0
   public:
    // import commonly defined symbols from base
    IMPORT_REGRESSION_SYMBOLS;
    using Base::L;         // [L]_{ii} = 1/DeltaT for i \in {1 ... m} and [L]_{i,i-1} = -1/DeltaT for i \in {1 ... m-1}
    using Base::lambda_D;   // smoothing parameter in space
    using Base::lambda_T;   // smoothing parameter in time
    using Base::n_temporal_locs;   // number of time instants m defined over [0,T]
    using Base::s;                 // initial condition
    // constructor
    STRPDE() = default;
    STRPDE(const PDE& pde, const DMatrix<double>& time) : Base(pde, time) {};

    void init_model();          // update model object in case of **structural** changes in its definition
    void update_to_weights();   // update model object in case of changes in the weights matrix
    virtual void solve();       // finds a solution to the smoothing problem

    // getters
    const SparseBlockMatrix<double, 2, 2>& A() const { return A_; }
    const fdapde::SparseLU<SpMatrix<double>>& invA() const { return invA_; }

    virtual ~STRPDE() = default;
};

// perform proper initialization and update of model. Computes quantites which can be reused
// across many calls to solve() and are **not affected by a change in the data**.
// It is implicitly called by ModelBase::init() as part of the initialization process.
// NB: a change in the smoothing parameter must trigger a re-initialization of the model
template <typename PDE, typename SamplingDesign>
void STRPDE<PDE, SpaceTimeParabolic, SamplingDesign, MonolithicSolver>::init_model() {
    // assemble system matrix for the nonparameteric part of the model
    std::cout << "STRPDE init model pt 1" << std::endl;
    if (is_empty(L_)) L_ = Kronecker(L(), pde().R0());
    A_ = SparseBlockMatrix<double, 2, 2>(
      -PsiTD() * W() * Psi(),                lambda_D() * (R1() + lambda_T() * L_).transpose(),
      lambda_D() * (R1() + lambda_T() * L_), lambda_D() * R0()                                );
    // cache system matrix for reuse
    invA_.compute(A_);
    std::cout << "STRPDE init model pt 2" << std::endl;
    // prepare rhs of linear system
    b_.resize(A_.rows());
    b_.block(A_.rows() / 2, 0, A_.rows() / 2, 1) = lambda_D() * u();
    std::cout << "STRPDE init model pt 3" << std::endl;
    return;
}

// updates model in case of a change in the weights matrix
template <typename PDE, typename SamplingDesign>
void STRPDE<PDE, SpaceTimeParabolic, SamplingDesign, MonolithicSolver>::update_to_weights() {
    // adjust north-west block of matrix A_ and factorize
    std::cout << "STRPDE update weights pt 1" << std::endl;
    A_.block(0, 0) = -PsiTD() * W() * Psi();
    invA_.compute(A_);
    std::cout << "STRPDE update weights pt 2" << std::endl;
    return;
}

// finds a solution to the STR-PDE smoothing problem (parabolic penalization, monolithic solution)
template <typename PDE, typename SamplingDesign>
void STRPDE<PDE, SpaceTimeParabolic, SamplingDesign, MonolithicSolver>::solve() {
    BLOCK_FRAME_SANITY_CHECKS;
    DVector<double> sol;   // room for problem' solution

    std::cout << "STRPDE solve pt 1" << std::endl; 

    if (!Base::has_covariates()) {   // nonparametric case
        // update rhs of STR-PDE linear system
        b_.block(0, 0, A_.rows() / 2, 1) = -PsiTD() * W() * y();
        std::cout << "STRPDE solve pt 2.0" << std::endl;
        // solve linear system A_*x = b_
        sol = invA_.solve(b_);
        std::cout << "STRPDE solve pt 2.1" << std::endl;
        f_ = sol.head(A_.rows() / 2);
    } else {   // parametric case
        // rhs of STR-PDE linear system
        b_.block(0, 0, A_.rows() / 2, 1) = -PsiTD() * lmbQ(y());   // -\Psi^T*D*Q*z

        // definition of matrices U and V  for application of woodbury formula
        U_ = DMatrix<double>::Zero(A_.rows(), q());
        U_.block(0, 0, A_.rows() / 2, q()) = PsiTD() * W() * X();
        V_ = DMatrix<double>::Zero(q(), A_.rows());
        V_.block(0, 0, q(), A_.rows() / 2) = X().transpose() * W() * Psi();
        // solve system (A_ + U_*(X^T*W_*X)*V_)x = b using woodbury formula from NLA module
        sol = SMW<>().solve(invA_, U_, XtWX(), V_, b_);
        // store result of smoothing
        f_ = sol.head(A_.rows() / 2);
        beta_ = invXtWX().solve(X().transpose() * W()) * (y() - Psi() * f_);
    }
    std::cout << "STRPDE solve pt 3" << std::endl;
    // store PDE misfit
    g_ = sol.tail(A_.rows() / 2);
    return;
}

template <typename PDE_, typename SamplingDesign_>
struct model_traits<STRPDE<PDE_, SpaceTimeParabolic, SamplingDesign_, MonolithicSolver>> {
    typedef PDE_ PDE;
    typedef SpaceTimeParabolic regularization;
    typedef SamplingDesign_ sampling;
    typedef MonolithicSolver solver;
    enum { N = PDE::N, M = PDE::M, n_lambda = 2 };
};

// implementation of STRPDE for parabolic space-time regularization, iterative solver
template <typename PDE, typename SamplingDesign>
class STRPDE<PDE, SpaceTimeParabolic, SamplingDesign, IterativeSolver> :
    public RegressionBase<STRPDE<PDE, SpaceTimeParabolic, SamplingDesign, IterativeSolver>> {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
   private:
    typedef SpaceTimeParabolic RegularizationType;
    typedef RegressionBase<STRPDE<PDE, SpaceTimeParabolic, SamplingDesign, IterativeSolver>> Base;
    SparseBlockMatrix<double, 2, 2> A_ {};      // system matrix of non-parametric problem (2N x 2N matrix)
    fdapde::SparseLU<SpMatrix<double>> invA_;   // factorization of matrix A
    DVector<double> b_ {};                      // right hand side of problem's linear system (1 x 2N vector)

    // the functional minimized by the iterative scheme
    // J(f,g) = \sum_{k=1}^m (z^k - \Psi*f^k)^T*(z^k - \Psi*f^k) + \lambda_S*(g^k)^T*(g^k)
    double J(const DMatrix<double>& f, const DMatrix<double>& g) const;
    // internal solve routine used by the iterative method
    void solve(std::size_t t, BlockVector<double>& f_new, BlockVector<double>& g_new) const;
    // vector of input data points at time k
    DMatrix<double> y(std::size_t k) const {
        return Base::y().block(Base::n_spatial_locs() * k, 0, Base::n_spatial_locs(), 1);
    }
   public:
    // import commonly defined symbols from base
    IMPORT_REGRESSION_SYMBOLS;
    using Base::DeltaT;             // distance between two time instants
    using Base::lambda_D;           // smoothing parameter in space
    using Base::lambda_T;           // smoothing parameter in time
    using Base::max_iter_;          // maximum number of allowed iterations before forced stop
    using Base::n_temporal_locs;    // number of time instants m defined over [0,T]
    using Base::tol_;               // tolerance on std::abs((Jnew - Jold)/Jnew)
    // constructor
    STRPDE() = default;
    STRPDE(const PDE& pde, const DMatrix<double>& time) : Base(pde, time) {};

    void init_model() { return; };          // update model object in case of **structural** changes in its definition
    void update_to_weights() { return; };   // update model object in case of changes in the weights matrix
    virtual void solve();                   // finds a solution to the smoothing problem

    virtual ~STRPDE() = default;
};

// J(f,g) = \sum_{k=1}^m (z^k - \Psi*f^k)^T*(z^k - \Psi*f^k) + \lambda_S*(g^k)^T*(g^k)
template <typename PDE, typename SamplingDesign>
double STRPDE<PDE, SpaceTimeParabolic, SamplingDesign, IterativeSolver>::J(
  const DMatrix<double>& f, const DMatrix<double>& g) const {
    double SSE = 0;
    // SSE = \sum_{k=1}^m (z^k - \Psi*f^k)^T*(z^k - \Psi*f^k)
    for (std::size_t t = 0; t < n_temporal_locs(); ++t) {
        SSE += (y(t) - Psi() * f.block(n_spatial_basis() * t, 0, n_spatial_basis(), 1)).squaredNorm();
    }
    return SSE + lambda_D() * g.squaredNorm();
}

// internal solve routine used by the iterative method
template <typename PDE, typename SamplingDesign>
void STRPDE<PDE, SpaceTimeParabolic, SamplingDesign, IterativeSolver>::solve(
  std::size_t t, BlockVector<double>& f_new, BlockVector<double>& g_new) const {
    DVector<double> x = invA_.solve(b_);
    f_new(t) = x.topRows(n_spatial_basis());
    g_new(t) = x.bottomRows(n_spatial_basis());
    return;
}

// finds a solution to the STR-PDE smoothing problem (parabolic penalization, iterative solution)
template <typename PDE, typename SamplingDesign>
void STRPDE<PDE, SpaceTimeParabolic, SamplingDesign, IterativeSolver>::solve() {
    // compute starting point (f^(k,0), g^(k,0)) k = 1 ... m for iterative minimization of functional J(f,g)
    A_ = SparseBlockMatrix<double, 2, 2>(
      PsiTD() * Psi(),    lambda_D() * R1().transpose(),
      lambda_D() * R1(), -lambda_D() * R0()            );
    // cache system matrix and its factorization
    invA_.compute(A_);
    b_.resize(A_.rows());

    // compute f^(k,0), k = 1 ... m as solution of Ax = b_(k)
    BlockVector<double> f_old(n_temporal_locs(), n_spatial_basis());
    // solve n_temporal_locs() space only linear systems
    for (std::size_t t = 0; t < n_temporal_locs(); ++t) {
        // right hand side at time step t
        b_ << PsiTD() * y(t),   // should put W()
          lambda_D() * lambda_T() * u(t);
        // solve linear system Ax = b_(t) and store estimate of spatial field
        f_old(t) = invA_.solve(b_).head(A_.rows() / 2);
    }

    // compute g^(k,0), k = 1 ... m as solution of the system
    //    G0 = [(\lambda_S*\lambda_T)/DeltaT * R_0 + \lambda_S*R_1^T]
    //    G0*g^(k,0) = \Psi^T*y^k + (\lambda_S*\lambda_T/DeltaT*R_0)*g^(k+1,0) - \Psi^T*\Psi*f^(k,0)
    SpMatrix<double> G0 =
      (lambda_D() * lambda_T() / DeltaT()) * R0() + SpMatrix<double>((lambda_D() * R1()).transpose());
    Eigen::SparseLU<SpMatrix<double>, Eigen::COLAMDOrdering<int>> invG0;
    invG0.compute(G0);   // compute factorization of matrix G0

    BlockVector<double> g_old(n_temporal_locs(), n_spatial_basis());
    // solve n_temporal_locs() distinct problems (in backward order)
    // at last step g^(t+1,0) is zero
    b_ = PsiTD() * (y(n_temporal_locs() - 1) - Psi() * f_old(n_temporal_locs() - 1));
    g_old(n_temporal_locs() - 1) = invG0.solve(b_);
    // general step
    for (int t = n_temporal_locs() - 2; t >= 0; --t) {
        // compute rhs at time t: \Psi^T*y^t + (\lambda_S*\lambda_T/DeltaT*R_0)*g^(t+1,0) - \Psi^T*\Psi*f^(t,0)
        b_ = PsiTD() * (y(t) - Psi() * f_old(t)) + (lambda_D() * lambda_T() / DeltaT()) * R0() * g_old(t + 1);
        // solve linear system G0*g^(t,1) = b_t and store estimate of PDE misfit
        g_old(t) = invG0.solve(b_);
    }

    // initialize value of functional J to minimize
    double Jold = std::numeric_limits<double>::max();
    double Jnew = J(f_old.get(), g_old.get());
    std::size_t i = 1;   // iteration number

    // build system matrix for the iterative scheme
    A_.block(0, 1) += lambda_D() * lambda_T() / DeltaT() * R0();
    A_.block(1, 0) += lambda_D() * lambda_T() / DeltaT() * R0();
    invA_.compute(A_);
    b_.resize(A_.rows());

    // internal iteration variables
    BlockVector<double> f_new(n_temporal_locs(), n_spatial_basis()), g_new(n_temporal_locs(), n_spatial_basis());
    // iterative scheme for minimization of functional J
    while (i < max_iter_ && std::abs((Jnew - Jold) / Jnew) > tol_) {
        // at step 0 f^(k-1,i-1) is zero
        b_ << PsiTD() * y(0) + (lambda_D() * lambda_T() / DeltaT()) * R0() * g_old(1), lambda_D() * u(0);
        // solve linear system
        solve(0, f_new, g_new);

        // general step
        for (std::size_t t = 1; t < n_temporal_locs() - 1; ++t) {
            // \Psi^T*y^k   + (\lambda_D*\lambda_T/DeltaT)*R_0*g^(k+1,i-1),
            // \lambda_D*u^k + (\lambda_D*\lambda_T/DeltaT)*R_0*f^(k-1,i-1)
            b_ << PsiTD() * y(t) + (lambda_D() * lambda_T() / DeltaT()) * R0() * g_old(t + 1),
	      lambda_D() * (lambda_T() / DeltaT() * R0() * f_old(t - 1) + u(t));
            // solve linear system
            solve(t, f_new, g_new);
        }

        // at last step g^(k+1,i-1) is zero
        b_ << PsiTD() * y(n_temporal_locs() - 1),
	  lambda_D() * (lambda_T() / DeltaT() * R0() * f_old(n_temporal_locs() - 2) + u(n_temporal_locs() - 1));
        // solve linear system
        solve(n_temporal_locs() - 1, f_new, g_new);

        // prepare for next iteration
        Jold = Jnew;
        f_old = f_new;
        g_old = g_new;
        Jnew = J(f_old.get(), g_old.get());
        i++;
    }

    // store solution
    f_ = f_old.get();
    g_ = g_old.get();
    return;
}
  
template <typename PDE_, typename SamplingDesign_>
struct model_traits<STRPDE<PDE_, SpaceTimeParabolic, SamplingDesign_, IterativeSolver>> {
    typedef PDE_ PDE;
    typedef SpaceTimeParabolic regularization;
    typedef SamplingDesign_ sampling;
    typedef IterativeSolver solver;
    enum { N = PDE::N, M = PDE::M, n_lambda = 2 };
};

// strpde trait
template <typename Model> struct is_strpde {
    static constexpr bool value = is_instance_of<Model, STRPDE>::value;
};

}   // namespace models
}   // namespace fdapde

#endif   // __STRPDE_H__

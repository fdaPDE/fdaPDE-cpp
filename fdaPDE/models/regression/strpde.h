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
using fdapde::core::SMW;
using fdapde::core::SparseBlockMatrix;
using fdapde::core::KroneckerTensorProduct;
using fdapde::core::SplineBasis;

namespace fdapde {
namespace models {

// STRPDE model signature
template <typename RegularizationType, typename SolutionPolicy> class STRPDE;

// implementation of STRPDE for separable space-time regularization
template <>
class STRPDE<SpaceTimeSeparable, monolithic> :
    public RegressionBase<STRPDE<SpaceTimeSeparable, monolithic>, SpaceTimeSeparable> {
   private:
    SparseBlockMatrix<double, 2, 2> A_ {};      // system matrix of non-parametric problem (2N x 2N matrix)
    fdapde::SparseLU<SpMatrix<double>> invA_;   // factorization of matrix A
    DVector<double> b_ {};                      // right hand side of problem's linear system (1 x 2N vector)
    DMatrix<double> T_;                         // T = \Psi^T*Q*\Psi + \lambda*R
    SpMatrix<double> K_;                        // P1 \kron R0
   public:
    using RegularizationType = SpaceTimeSeparable;
    using Base = RegressionBase<STRPDE<RegularizationType, monolithic>, RegularizationType>;
    // import commonly defined symbols from base
    IMPORT_REGRESSION_SYMBOLS
    using Base::lambda_D;   // smoothing parameter in space
    using Base::lambda_T;   // smoothing parameter in time
    using Base::P;          // discretized penalty: P = \lambda_D*((R1^T*R0^{-1}*R1) \kron Rt) + \lambda_T*(R0 \kron Pt)
    using Base::P0;         // time mass matrix: [P0_]_{ij} = \int_{[0,T]} \phi_i*\phi_j
    using Base::P1;         // time penalization matrix: [P1_]_{ij} = \int_{[0,T]} (\phi_i)_tt*(\phi_j)_tt
    // constructor
    STRPDE() = default;
    STRPDE(const Base::PDE& space_penalty, const Base::PDE& time_penalty, Sampling s) :
        Base(space_penalty, time_penalty, s) {};

    void init_model() {
        if (runtime().query(runtime_status::is_lambda_changed)) {
            // assemble system matrix for the nonparameteric part
            if (is_empty(K_)) K_ = Kronecker(P1(), pde().mass());
            A_ = SparseBlockMatrix<double, 2, 2>(
              -PsiTD() * W() * Psi() - lambda_T() * K_, lambda_D() * R1().transpose(),
	      lambda_D() * R1(),                        lambda_D() * R0()            );
            invA_.compute(A_);
            // prepare rhs of linear system
            b_.resize(A_.rows());
            b_.block(A_.rows() / 2, 0, A_.rows() / 2, 1) = lambda_D() * u();
            return;
        }
        if (runtime().query(runtime_status::require_W_update)) {
            // adjust north-west block of matrix A_ and factorize
            A_.block(0, 0) = -PsiTD() * W() * Psi() - lambda_T() * K_;
            invA_.compute(A_);
            return;
        }
    }
    void solve() {
        fdapde_assert(y().rows() != 0);

        // std::cout << "strpde solve: n tot " << y().size() << std::endl;
        // std::cout << "strpde solve: n " << n_obs() << std::endl;
        // std::cout << "strpde solve: n missing " << Base::masked_obs().count() << std::endl;

        DVector<double> sol;             // room for problem' solution
        if (!Base::has_covariates()) {   // nonparametric case
            // update rhs of STR-PDE linear system
            b_.block(0, 0, A_.rows() / 2, 1) = -PsiTD() * W() * y();
            // solve linear system A_*x = b_
            sol = invA_.solve(b_);
            f_ = sol.head(A_.rows() / 2);
        } else {   // parametric case
            // update rhs of STR-PDE linear system
            b_.block(0, 0, A_.rows() / 2, 1) = -PsiTD() * lmbQ(y());   // -\Psi^T*D*Q*z
            // matrices U and V for application of woodbury formula
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
    // GCV support
    double norm(const DMatrix<double>& op1, const DMatrix<double>& op2) const { return (op1 - op2).squaredNorm(); }
    // getters
    const SparseBlockMatrix<double, 2, 2>& A() const { return A_; }
    const fdapde::SparseLU<SpMatrix<double>>& invA() const { return invA_; }
};

// implementation of STRPDE for parabolic space-time regularization, monolithic approach
template <>
class STRPDE<SpaceTimeParabolic, monolithic> :
    public RegressionBase<STRPDE<SpaceTimeParabolic, monolithic>, SpaceTimeParabolic> {
   private:
    SparseBlockMatrix<double, 2, 2> A_ {};      // system matrix of non-parametric problem (2N x 2N matrix)
    fdapde::SparseLU<SpMatrix<double>> invA_;   // factorization of matrix A
    DVector<double> b_ {};                      // right hand side of problem's linear system (1 x 2N vector)
    SpMatrix<double> L_;                        // L \kron R0
   public:
    using RegularizationType = SpaceTimeParabolic;
    using Base = RegressionBase<STRPDE<RegularizationType, monolithic>, RegularizationType>;
    // import commonly defined symbols from base
    IMPORT_REGRESSION_SYMBOLS
    using Base::L;          // [L]_{ii} = 1/DeltaT for i \in {1 ... m} and [L]_{i,i-1} = -1/DeltaT for i \in {1 ... m-1}
    using Base::lambda_D;   // smoothing parameter in space
    using Base::lambda_T;   // smoothing parameter in time
    using Base::n_temporal_locs;   // number of time instants m defined over [0,T]
    using Base::s;                 // initial condition
    // constructor
    STRPDE() = default;
    STRPDE(const Base::PDE& pde, Sampling s) : Base(pde, s) {};

    void init_model() {   // update model object in case of **structural** changes in its definition
        // assemble system matrix for the nonparameteric part of the model
        if (is_empty(L_)) L_ = Kronecker(L(), pde().mass());
        A_ = SparseBlockMatrix<double, 2, 2>(
          -PsiTD() * W() * Psi(), lambda_D() * (R1() + lambda_T() * L_).transpose(),
          lambda_D() * (R1() + lambda_T() * L_), lambda_D() * R0());
        // cache system matrix for reuse
        invA_.compute(A_);
        // prepare rhs of linear system
        b_.resize(A_.rows());
        b_.block(A_.rows() / 2, 0, A_.rows() / 2, 1) = lambda_D() * u();
        return;
    }
    void update_to_weights() {   // update model object in case of changes in the weights matrix
        // adjust north-west block of matrix A_ and factorize
        A_.block(0, 0) = -PsiTD() * W() * Psi();
        invA_.compute(A_);
        return;
    }
    void solve() {
        fdapde_assert(y().rows() != 0);
        DVector<double> sol;             // room for problem' solution
        if (!Base::has_covariates()) {   // nonparametric case
            // update rhs of STR-PDE linear system
            b_.block(0, 0, A_.rows() / 2, 1) = -PsiTD() * W() * y();
            // solve linear system A_*x = b_
            sol = invA_.solve(b_);
            f_ = sol.head(A_.rows() / 2);
        } else {   // parametric case
            // rhs of STR-PDE linear system
            b_.block(0, 0, A_.rows() / 2, 1) = -PsiTD() * lmbQ(y());   // -\Psi^T*D*Q*z
            // matrices U and V for application of woodbury formula
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
    // getters
    const SparseBlockMatrix<double, 2, 2>& A() const { return A_; }
    const fdapde::SparseLU<SpMatrix<double>>& invA() const { return invA_; }
    double norm(const DMatrix<double>& op1, const DMatrix<double>& op2) const {   // euclidian norm of op1 - op2
        return (op1 - op2).squaredNorm(); // NB: to check, defined just for compiler
    }
};

// implementation of STRPDE for parabolic space-time regularization, iterative approach
template <>
class STRPDE<SpaceTimeParabolic, iterative> :
    public RegressionBase<STRPDE<SpaceTimeParabolic, iterative>, SpaceTimeParabolic> {
   private:
    SparseBlockMatrix<double, 2, 2> A_ {};      // system matrix of non-parametric problem (2N x 2N matrix)
    fdapde::SparseLU<SpMatrix<double>> invA_;   // factorization of matrix A
    DVector<double> b_ {};                      // right hand side of problem's linear system (1 x 2N vector)

    // the functional minimized by the iterative scheme
    // J(f,g) = \sum_{k=1}^m (z^k - \Psi*f^k)^T*(z^k - \Psi*f^k) + \lambda_S*(g^k)^T*(g^k)
    double J(const DMatrix<double>& f, const DMatrix<double>& g) const {
        double SSE = 0;
        // SSE = \sum_{k=1}^m (z^k - \Psi*f^k)^T*(z^k - \Psi*f^k)
        for (int t = 0; t < n_temporal_locs(); ++t) {
            SSE += (y(t) - Psi() * f.block(n_spatial_basis() * t, 0, n_spatial_basis(), 1)).squaredNorm();
        }
        return SSE + lambda_D() * g.squaredNorm();
    }
    // internal solve routine used by the iterative method
    void solve(int t, BlockVector<double>& f_new, BlockVector<double>& g_new) const {
        DVector<double> x = invA_.solve(b_);
        f_new(t) = x.topRows(n_spatial_basis());
        g_new(t) = x.bottomRows(n_spatial_basis());
        return;
    }
    // internal utilities
    DMatrix<double> y(int k) const { return y().block(n_spatial_locs() * k, 0, n_spatial_locs(), 1); }
    DMatrix<double> u(int k) const { return u_.block(n_basis() * k, 0, n_basis(), 1); }

    // quantities related to iterative scheme
    double tol_ = 1e-4;           // tolerance used as stopping criterion
    int max_iter_ = 50;   // maximum number of allowed iterations
   public:
    using RegularizationType = SpaceTimeParabolic;
    using Base = RegressionBase<STRPDE<RegularizationType, iterative>, RegularizationType>;
    // import commonly defined symbols from base
    IMPORT_REGRESSION_SYMBOLS
    using Base::DeltaT;            // distance between two time instants
    using Base::lambda_D;          // smoothing parameter in space
    using Base::lambda_T;          // smoothing parameter in time
    using Base::n_temporal_locs;   // number of time instants m defined over [0,T]
    using Base::pde_;              // parabolic differential operator df/dt + Lf - u
    // constructor
    STRPDE() = default;
    STRPDE(const Base::PDE& pde, Sampling s) : Base(pde, s) { pde_.init(); };

    // redefine SpaceTimeParabolicBase properties affected by iterative approach
    void tensorize_psi() { return; } // avoid tensorization of \Psi matrix
    void init_regularization() {
        pde_.init();
	s_ = pde_.initial_condition();
        // compute time step (assuming equidistant points)
        DeltaT_ = time_[1] - time_[0];
        u_ = pde_.force();   // compute forcing term
        // correct first n rows of discretized force as (u_1 + R0*s/DeltaT)
        u_.block(0, 0, n_basis(), 1) += (1.0 / DeltaT_) * (pde_.mass() * s_);
    }
    // getters
    const SpMatrix<double>& R0() const { return pde_.mass(); }    // mass matrix in space
    const SpMatrix<double>& R1() const { return pde_.stiff(); }   // discretization of differential operator L
    int n_basis() const { return pde_.n_dofs(); }         // number of basis functions

    void init_model() { return; };
    void solve() {
        fdapde_assert(y().rows() != 0);
        // compute starting point (f^(k,0), g^(k,0)) k = 1 ... m for iterative minimization of functional J(f,g)
        A_ = SparseBlockMatrix<double, 2, 2>(
          PsiTD() * Psi(),   lambda_D() * R1().transpose(),
	  lambda_D() * R1(), -lambda_D() * R0()           );
        // cache system matrix and its factorization
        invA_.compute(A_);
        b_.resize(A_.rows());

        // compute f^(k,0), k = 1 ... m as solution of Ax = b_(k)
        BlockVector<double> f_old(n_temporal_locs(), n_spatial_basis());
        // solve n_temporal_locs() space only linear systems
        for (int t = 0; t < n_temporal_locs(); ++t) {
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
        int i = 1;   // iteration number
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
            for (int t = 1; t < n_temporal_locs() - 1; ++t) {
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
    // setters
    void set_tolerance(double tol) { tol_ = tol; }
    void set_max_iter(int max_iter) { max_iter_ = max_iter; }
};

}   // namespace models
}   // namespace fdapde

#endif   // __STRPDE_H__

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

#ifndef __SRPDE_H__
#define __SRPDE_H__

#include <fdaPDE/linear_algebra.h>
#include <fdaPDE/pde.h>
#include <fdaPDE/utils.h>

#include <memory>
#include <type_traits>

#include "../model_base.h"
#include "../model_macros.h"
#include "../sampling_design.h"
#include "regression_base.h"
using fdapde::core::PDEBase;
using fdapde::core::SMW;
using fdapde::core::SparseBlockMatrix;

namespace fdapde {
namespace models {

template <typename PDE, typename SamplingDesign> class SRPDE : public RegressionBase<SRPDE<PDE, SamplingDesign>> {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
   private:
    typedef RegressionBase<SRPDE<PDE, SamplingDesign>> Base;
    SparseBlockMatrix<double, 2, 2> A_ {};         // system matrix of non-parametric problem (2N x 2N matrix)
    fdapde::SparseLU<SpMatrix<double>> invA_ {};   // factorization of matrix A
    DVector<double> b_ {};                         // right hand side of problem's linear system (1 x 2N vector)
   public:
    IMPORT_REGRESSION_SYMBOLS;
    using Base::lambda_D;   // smoothing parameter in space
    using Base::n_basis;    // number of spatial basis
    // constructor
    SRPDE() = default;
    SRPDE(const PDE& pde) : Base(pde) {};

    void init_model();          // update model object in case of **structural** changes in its definition
    void update_to_weights();   // update model object in case of changes in the weights matrix
    virtual void solve();       // finds a solution to the smoothing problem  
    double norm(const DMatrix<double>& op1, const DMatrix<double>& op2) const;   // euclidian norm of op1 - op2

    // getters
    const SparseBlockMatrix<double, 2, 2>& A() const { return A_; }
    const fdapde::SparseLU<SpMatrix<double>>& invA() const { return invA_; }

    virtual ~SRPDE() = default;
};

// implementative details

// model initialization. Called by ModelBase::init() as part of the initialization process.
// NB: a change in the smoothing parameter must trigger a re-initialization of the model
template <typename PDE, typename SamplingDesign> void SRPDE<PDE, SamplingDesign>::init_model() {
    // assemble system matrix for nonparameteric part
    A_ = SparseBlockMatrix<double, 2, 2>(
      -PsiTD() * W() * Psi(), lambda_D() * R1().transpose(),
      lambda_D() * R1(),      lambda_D() * R0()            );
    // cache non-parametric matrix factorization for reuse
    invA_.compute(A_);
    // prepare rhs of linear system
    b_.resize(A_.rows());
    b_.block(n_basis(), 0, n_basis(), 1) = lambda_D() * u();
    return;
}

// updates model in case of a change in the weights matrix
template <typename PDE, typename SamplingDesign> void SRPDE<PDE, SamplingDesign>::update_to_weights() {
    // adjust north-west block of matrix A_ and factorize
    A_.block(0, 0) = -PsiTD() * W() * Psi();
    invA_.compute(A_);
    return;
}

// finds a solution to the SR-PDE smoothing problem
template <typename PDE, typename SamplingDesign> void SRPDE<PDE, SamplingDesign>::solve() {
    BLOCK_FRAME_SANITY_CHECKS;
    DVector<double> sol;   // room for problem' solution

    if (!has_covariates()) {   // nonparametric case
        // update rhs of SR-PDE linear system
        b_.block(0, 0, n_basis(), 1) = -PsiTD() * W() * y();
        // solve linear system A_*x = b_
        sol = invA_.solve(b_);
        f_ = sol.head(n_basis());
    } else {   // parametric case
        // update rhs of SR-PDE linear system
        b_.block(0, 0, n_basis(), 1) = -PsiTD() * lmbQ(y());   // -\Psi^T*D*Q*z

        // definition of matrices U and V  for application of woodbury formula
        U_ = DMatrix<double>::Zero(2 * n_basis(), q());
        U_.block(0, 0, n_basis(), q()) = PsiTD() * W() * X();
        V_ = DMatrix<double>::Zero(q(), 2 * n_basis());
        V_.block(0, 0, q(), n_basis()) = X().transpose() * W() * Psi();
        // solve system (A_ + U_*(X^T*W_*X)*V_)x = b using woodbury formula from NLA module
        sol = SMW<>().solve(invA_, U_, XtWX(), V_, b_);
        // store result of smoothing
        f_ = sol.head(n_basis());
        beta_ = invXtWX().solve(X().transpose() * W()) * (y() - Psi() * f_);
    }
    // store PDE misfit
    g_ = sol.tail(n_basis());
    return;
}

// returns the euclidean norm of op1 - op2
template <typename PDE, typename SamplingDesign>
double SRPDE<PDE, SamplingDesign>::norm(const DMatrix<double>& op1, const DMatrix<double>& op2) const {
    return (op1 - op2).squaredNorm();
}

template <typename PDE_, typename SamplingDesign_> struct model_traits<SRPDE<PDE_, SamplingDesign_>> {
    typedef PDE_ PDE;
    typedef SpaceOnly regularization;
    typedef SamplingDesign_ sampling;
    typedef MonolithicSolver solver;
    enum { N = PDE::N, M = PDE::M, n_lambda = 1 };
};

// srpde trait
template <typename Model> struct is_srpde {
    static constexpr bool value = is_instance_of<Model, SRPDE>::value;
};

}   // namespace models
}   // namespace fdapde

#endif   // __SRPDE_H__

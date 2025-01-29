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

#ifndef __FE_PARABOLIC_DRIVER_H__
#define __FE_PARABOLIC_DRIVER_H__

#include "header_check.h"

namespace fdapde {
namespace internals {

struct fe_parabolic_driver_base {
    using VectorType = Eigen::Matrix<double, Dynamic, 1>;
    using MatrixType = Eigen::Matrix<double, Dynamic, Dynamic>;
    using SparseMatrixType = Eigen::SparseMatrix<double>;
    using DiagonalMatrixType = Eigen::DiagonalMatrix<double, Dynamic, Dynamic>;
    using SparseSolverType = Eigen::SparseLU<SparseMatrixType>;
    using DenseSolverType  = Eigen::PartialPivLU<MatrixType>;
  
    fe_parabolic_driver_base() noexcept = default;

  // need setters for initial condition
  // need to handle boundary conditions
  
    template <typename BilinearForm_, typename LinearForm_, typename GeoFrame, typename InitialCondition_>
    fe_parabolic_driver_base(
      const GeoFrame& gf, BilinearForm_&& bilinear_form, LinearForm_&& linear_form, const InitialCondition_& s) :
      R1_(bilinear_form.assemble()), u_(linear_form.assemble()), s_(s) {
        using BilinearForm = std::decay_t<BilinearForm_>;
        using LinearForm   = std::decay_t<LinearForm_>;
        using FeSpace = typename BilinearForm::TrialSpace;

	internals::fe_mass_assembly_loop<FeSpace> mass_assembler(bilinear_form.trial_space());
        R0_ = mass_assembler.assemble();     // mass matrix
        n_dofs_ = bilinear_form.n_dofs();    // number of basis functions over physical domain

	// extract spatial and temporal locations (only SpaceMajor format support)
	std::vector<double> space_locs;
	std::vector<double> time_locs;

	// get n_spatial_locations, n_temporal_locations
	
	// then we can derive the observational mask (for each time instant, which spatial location is observed)

	// then we can assemble the Psi (which, for fully observed, same number of locations, is exactly the tensorized Psi with the identity)

	// then derive delta (assert locations are at same time distance)

	// assemble matrix associated with derivation in time L_
        // [L_]_{ii} = 1/DeltaT for i \in {1 ... m} and [L_]_{i,i-1} = -1/DeltaT for i \in {1 ... m-1}
        std::vector<fdapde::Triplet<double>> triplet_list;
        triplet_list.reserve(2 * m_);
        // start assembly loop
        double invDeltaT = 1.0 / DeltaT_;
        triplet_list.emplace_back(0, 0, invDeltaT);
        for (int i = 1; i < m_; ++i) {
            triplet_list.emplace_back(i, i, invDeltaT);
            triplet_list.emplace_back(i, i - 1, -invDeltaT);
        }
        L_.resize(m_, m_);
        L_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        L_.makeCompressed();

	// correct first n rows of discretized force as (u_1 + R0*s/DeltaT)
        u_.block(0, 0, n_dofs_, 1) += (1.0 / DeltaT_) * (R0_ * s_); // ------------------------- this might be ok only for the monolithic

        // evaluate basis system on physical domain
        switch (gf.layer_category(0).value()) {
        case ltype::point: {
            const auto& layer = gf.get_as(layer_t::point, 0);
            if (layer.locs_at_mesh_nodes()) {
                // locations at mesh nodes
                Psi_.resize(n_dofs_, n_dofs_);
                Psi_.setIdentity();   // \psi_i(p_j) = 1 \iff i == j, otherwise \psi_i(p_j) = 0
            } else {
                Psi_ = internals::point_basis_eval(bilinear_form.trial_space(), layer.coordinates());
            }
            D_ = VectorType::Ones(Psi_.rows()).asDiagonal();
            break;
        }
        case ltype::areal: {
            const auto& layer = gf.get_as(layer_t::areal, 0);
            const auto& [psi, measure_vect] =
              internals::areal_basis_eval(bilinear_form.trial_space(), layer.incidence_matrix());
            Psi_ = psi;
            D_ = measure_vect.asDiagonal();   // regions' measure
            break;
        }
        }
	// we build only one Psi on the full grid, and a vector of observational masks for each time instant
    }
    // observers
    int n_dofs() const { return n_dofs_; }
    const SparseMatrixType& R0() const { return R0_; }
    const SparseMatrixType& R1() const { return R1_; }
    const SparseMatrixType& Psi() const { return Psi_; }
    const VectorType& u() const { return u_; }
    // penalty matrix
    // MatrixType P(double lambda) const {
    //     if (invR0_.info() != Eigen::Success) { invR0_.compute(R0_); }
    //     return lambda * R1_.transpose() * invR0_.solve(R1_);
    // }
   protected:
    int n_dofs_ = 0;         // number of basis on physical domain
    SparseMatrixType R0_;    // mass matrix: [R0]_{ij} = \int_D \psi_i * \psi_j
    SparseMatrixType R1_;    // discretization of bilinear form a: [R1]_{ij} = \int_D a(\psi_i, \psi_j)
    SparseMatrixType Psi_;   // evaluation of basis system \psi_1, ..., \psi_{n_dofs_} on physical domain
    SparseMatrixType K_;
    VectorType u_;           // discretized forcing term u_i = \int_D u * \psi_i
    DiagonalMatrixType D_;   // vector of regions' measures (areal sampling)
    mutable SparseSolverType invR0_;
};

template <typename SolutionStrategy> struct fe_parabolic_driver_impl;

// solves \min_{f, \beta} \| W^{1/2} * (y_i - x_i^\top * \beta - f(p_i, t_j)) \|_2^2 + \int_D \int_T (L_D(f) - u_D)^2 +
// \int_T \int_D (L_T(f) - u_T)^2
template <> struct fe_parabolic_driver_impl<monolithic_tag> : fe_parabolic_driver_base {
   private:
    template <typename GeoFrame, typename WeightMatrix>
    void init_(const std::string& formula, const GeoFrame& gf, const WeightMatrix& W) {
        // parse formula
        Formula formula_(formula);
        std::vector<std::string> covs;
        for (const std::string& token : formula_.rhs()) {
            if (gf.contains(token)) { covs.push_back(token); }
        }
        q_ = covs.size();
        // extract data from geoframe
	n_obs_ = gf[0].rows();
        y_.resize(n_obs_);
        gf[0].template col<double>(formula_.lhs()).data().assign_to(y_);
        if (q_ != 0) {
            // assemble design matrix
            X_.resize(gf[0].rows(), q_);
            for (int i = 0; i < q_; ++i) { gf[0].template col<double>(covs[i]).data().assign_to(X_.col(i)); }
            XtX_ = X_.transpose() * W * X_;
            invXtX_ = XtX_.partialPivLu();
            invXtXXt_ = invXtX_.solve(X_.transpose() * W);   // (X^\top * X)^{-1} * X^\top * W
            // woodbury decomposition matrices
            U_ = MatrixType::Zero(2 * n_dofs_, q_);
            U_.block(0, 0, n_dofs_, q_) = Psi_.transpose() * D_ * W * X_;
            V_ = MatrixType::Zero(q_, 2 * n_dofs_);
            V_.block(0, 0, q_, n_dofs_) = X_.transpose() * W * Psi_;
        }
	// tensorize
	L_ = kronecker(L_, R0_); // ------------------------ this L_ is the one of the basis, need to correct
	R0_ = kronecker(Im_, R0_); // ----------------------- Im_ identity matrix
        R1_ = kronecker(Im_, R1_);
	return;
    }
   public:
    fe_parabolic_driver_impl() noexcept = default;
    template <typename BilinearForm, typename LinearForm, typename GeoFrame, typename WeightMatrix>
    fe_parabolic_driver_impl(
      const std::string& formula, const GeoFrame& gf, BilinearForm&& bilinear_form, LinearForm&& linear_form,
      const WeightMatrix& W) :
        fe_parabolic_driver_base(gf, bilinear_form, linear_form) {
        init_(formula, gf, W);
    }
    template <typename BilinearForm, typename LinearForm, typename GeoFrame>
    fe_parabolic_driver_impl(
      const std::string& formula, const GeoFrame& gf, BilinearForm&& bilinear_form, LinearForm&& linear_form) :
        fe_parabolic_driver_base(gf, bilinear_form, linear_form) {
        init_(formula, gf, Eigen::Matrix<double, Dynamic, 1>::Ones(y_.rows()).asDiagonal());
    }

    void operator()(double lambda_D, double lambda_T) {
        // assemble system matrix for the nonparameteric part
        A_ = SparseBlockMatrix<double, 2, 2>(
          -Psi_.transpose() * D_ * Psi_, lambda_D * (R1_ + lambda_T * L_).transpose(), lambda_D * (R1_ + lambda_T * L_),
          lambda_D * R0_);
        invA_.compute(A_);
        // linear system rhs
        b_.resize(2 * n_dofs_);
        b_.block(n_dofs_, 0, n_dofs_, 1) = lambda_D * u_;
	
        VectorType x;
        if (q_ == 0) {   // nonparametric case
            b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * y_;
            x = invA_.solve(b_);
            f_ = x.head(n_dofs_);
        } else {   // parametric case
            b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * internals::lmbQ(X_, invXtX_, y_);
            x = woodbury_system_solve(invA_, U_, XtX_, V_, b_);
            f_ = x.head(n_dofs_);
            beta_ = invXtXXt_ * (y_ - Psi_ * f_);
        } 
        g_ = x.tail(n_dofs_);   // PDE misfit
        return;
    }
    template <typename WeightMatrix> void operator()(double lambda_D, double lambda_T, WeightMatrix&& W) {
        // assemble system matrix for the nonparameteric part
        A_ = SparseBlockMatrix<double, 2, 2>(
          -Psi_.transpose() * D_ * W * Psi_, lambda_D * (R1_ + lambda_T * L_).transpose(),
          lambda_D * (R1_ + lambda_T * L_), lambda_D * R0_);
        invA_.compute(A_);
        // linear system rhs
        VectorType b_(2 * n_dofs_);
        b_.block(n_dofs_, 0, n_dofs_, 1) = lambda_D * u_;

        VectorType x;
        if (q_ == 0) {   // nonparametric case
            b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * W * y_;
            x = invA_.solve(b_);
            f_ = x.head(n_dofs_);
        } else {   // parametric case
            XtX_ = X_.transpose() * W * X_;
            invXtX_ = XtX_.partialPivLu();
            b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * internals::lmbQ(W, X_, invXtX_, y_);
            // woodbury matrices
            U_.block(0, 0, n_dofs_, q_) = Psi_.transpose() * D_ * W * X_;
            V_.block(0, 0, q_, n_dofs_) = X_.transpose() * W * Psi_;
            // solve A * x = (A_ + U_ * (X^\top*W*X) * V_) * x = b
            x = woodbury_system_solve(invA_, U_, XtX_, V_, b_);
            f_ = x.head(n_dofs_);
            beta_ = invXtX_.solve(X_.transpose() * W) * (y_ - Psi_ * f_);
        }
        g_ = x.tail(n_dofs_);   // PDE misfit
        return;
    }
    // observers
    const VectorType& f() const { return f_; }
    const VectorType& beta() const { return beta_; }
    const VectorType& g() const { return g_; }
   private:
    int q_;                       // number of covariates
    int n_obs_;
    MatrixType X_;                // n_obs x q design matrix
    VectorType y_;                // n_obs x 1 observation vector
    MatrixType U_, V_;            // (2 * n_dofs) x q matrices [\Psi^\top * D * y, 0] and [X^\top * \Psi, 0]
    MatrixType XtX_, invXtXXt_;   // q x q matrix X^\top * X and q x n_obs matrix (X^\top * X)^{-1} * X^\top
    DenseSolverType invXtX_;      // factorization of q x q matrix X^\top * X
    SparseSolverType invA_;
    VectorType f_, beta_, g_;
};

template <> struct fe_parabolic_driver_impl<iterative_tag> : fe_parabolic_driver_base {

};
  
}   // namespace internals
}   // namespace fdapde

#endif // __FE_PARABOLIC_DRIVER_H__

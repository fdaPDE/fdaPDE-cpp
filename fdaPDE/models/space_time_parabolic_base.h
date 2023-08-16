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

#ifndef __SPACE_TIME_PARABOLIC_BASE_H__
#define __SPACE_TIME_PARABOLIC_BASE_H__

#include <fdaPDE/finite_elements.h>
#include <fdaPDE/linear_algebra.h>
#include <fdaPDE/utils.h>
using fdapde::core::is_parabolic;
using fdapde::core::Kronecker;
using fdapde::core::SparseKroneckerProduct;

#include "space_time_base.h"

namespace fdapde {
namespace models {

// base class for parabolic regularization solved using either a moholitic or iterative solution strategy
template <typename Model, typename Solver> class SpaceTimeParabolicBase;

// base class for parabolic regularization, monholitic solver
template <typename Model> class SpaceTimeParabolicBase<Model, MonolithicSolver> : public SpaceTimeBase<Model> {
    static_assert(
      is_parabolic<typename model_traits<Model>::PDE::OperatorType>::value,
      "you have asked for parabolic regularization but using a non-parabolic differential operator");
   private:
    // let m the number of time points
    DMatrix<double> s_;     // N x 1 initial condition vector
    DMatrix<double> u_;     // discretized forcing [1/DeltaT * (u_1 + R_0*s) \ldots u_n]
    SpMatrix<double> Im_;   // m x m sparse identity matrix (assembled once and cached for reuse)
    SpMatrix<double> L_;    // m x m matrix associated with the derivation in time
    double DeltaT_;         // time step (assumes equidistant points in time)
    SpMatrix<double> R0_;   // Im \kron R0 (R0: discretization of the identity operator)
    SpMatrix<double> R1_;   // Im \kron R1 (R1: discretization of the differential operator L in the regularizing PDE)

    SpMatrix<double> penT_;                      // discretization of the time derivative: L \kron R0
    fdapde::SparseLU<SpMatrix<double>> invR0_;   // factorization of Im \kron R0
    // discretized penalty: (Im \kron R1 + L \kron R0)^T*(I_m \kron R0)^{-1}*(Im \kron R1 + L \kron R0)
    SpMatrix<double> pen_;
   public:
    typedef typename model_traits<Model>::PDE PDE;                             // PDE used in space
    typedef typename model_traits<Model>::regularization TimeRegularization;   // regularization in time
    typedef SpaceTimeBase<Model> Base;
    using Base::lambda_D;   // smoothing parameter in space
    using Base::lambda_T;   // smoothing parameter in time
    using Base::model;      // underlying model object
    using Base::pde_;       // regularizing term in space
    using Base::time_;      // time interval [0,T]

    // constructor
    SpaceTimeParabolicBase() = default;
    SpaceTimeParabolicBase(const PDE& pde, const DVector<double>& time) : SpaceTimeBase<Model>(pde, time) { }
    // init data structure related to parabolic regularization
    void init_regularization() {
        std::size_t m_ = time_.rows();   // number of time points
        DeltaT_ = time_[1] - time_[0];   // time step (assuming equidistant points)

        // assemble once the m x m identity matrix and cache for fast access
        Im_.resize(m_, m_);
        Im_.setIdentity();

        // assemble matrix associated with derivation in time L_
        // [L_]_{ii} = 1/DeltaT for i \in {1 ... m} and [L_]_{i,i-1} = -1/DeltaT for i \in {1 ... m-1}
        std::vector<fdapde::Triplet<double>> triplet_list;
        triplet_list.reserve(2 * m_);
        // start assembly loop
        double invDeltaT = 1.0 / DeltaT_;
        triplet_list.emplace_back(0, 0, invDeltaT);
        for (std::size_t i = 1; i < m_; ++i) {
            triplet_list.emplace_back(i, i, invDeltaT);
            triplet_list.emplace_back(i, i - 1, -invDeltaT);
        }
        // finalize construction
        L_.resize(m_, m_);
        L_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        L_.makeCompressed();
        // compute tensorized matrices
        R0_ = Kronecker(Im_, pde_->R0());
        R1_ = Kronecker(Im_, pde_->R1());
    }

    // getters
    const SpMatrix<double>& R0() const { return R0_; }
    const SpMatrix<double>& R1() const { return R1_; }
    std::size_t n_basis() const { return pde_->domain().dof(); }   // number of basis functions
    // matrices proper of separable regularization
    const SpMatrix<double>& L() const { return L_; }
    // return discretized force corrected by initial conditions
    const DMatrix<double>& u() {
        if (is_empty(u_)) {   // compute once and cache
            u_ = pde_->force();
            // correct first n rows of discretized force as (u_1 + R0*s/DeltaT)
            u_.block(0, 0, model().n_basis(), 1) += (1.0 / DeltaT_) * (pde_->R0() * s_);
        }
        return u_;
    }
    const DMatrix<double>& s() { return s_; }   // initial condition
    double DeltaT() const { return DeltaT_; }

    // computes and cache matrices (Im \kron R0)^{-1} and L \kron R0, returns the discretized penalty P =
    // \lambda_D*((Im \kron R1 + \lambda_T*(L \kron R0))^T*(I_m \kron R0)^{-1}*(Im \kron R1 + \lambda_T*(L \kron R0)))
    auto P() {
        if (is_empty(pen_)) {   // compute once and cache result
            invR0_.compute(R0());
            penT_ = Kronecker(L_, pde_->R0());
        }
        return lambda_D() * (R1() + lambda_T() * penT_).transpose() * invR0_.solve(R1() + lambda_T() * penT_);
    }
    // setters
    void set_initial_condition(const DMatrix<double>& s) { s_ = s; }

    // destructor
    virtual ~SpaceTimeParabolicBase() = default;
};

// base class for parabolic regularization, iterative solver
template <typename Model> class SpaceTimeParabolicBase<Model, IterativeSolver> : public SpaceTimeBase<Model> {
    static_assert(
      is_parabolic<typename model_traits<Model>::PDE::OperatorType>::value,
      "you have asked for parabolic regularization but using a non-parabolic differential operator");
   protected:
    typedef typename model_traits<Model>::PDE PDE;                             // PDE used for regularization in space
    typedef typename model_traits<Model>::regularization TimeRegularization;   // regularization in time
    typedef SpaceTimeBase<Model> Base;
    using Base::model;   // underlying model object
    using Base::pde_;    // regularizing term in space
    using Base::time_;   // time interval [0,T]

    DMatrix<double> s_;   // N x 1 initial condition vector
    DMatrix<double> u_;   // discretized forcing [1/DeltaT * (u_1 + R_0*s) \ldots u_n]
    double DeltaT_;       // time step (assumes equidistant points in time)

    // quantities related to iterative scheme
    double tol_ = 1e-4;           // tolerance used as stopping criterion
    std::size_t max_iter_ = 50;   // maximum number of allowed iterations
   public:
    // constructor
    SpaceTimeParabolicBase() = default;
    SpaceTimeParabolicBase(const PDE& pde, const DVector<double>& time) : SpaceTimeBase<Model>(pde, time) { }
    // init data required for iterative solution of parabolic regularization
    void init_regularization() {
        // compute time step (assuming equidistant points)
        DeltaT_ = time_[1] - time_[0];
        u_ = pde_->force();   // compute forcing term
        // correct first n rows of discretized force as (u_1 + R0*s/DeltaT)
        u_.block(0, 0, model().n_basis(), 1) += (1.0 / DeltaT_) * (pde_->R0() * s_);
    }

    // getters
    const SpMatrix<double>& R0() const { return pde_->R0(); }   // mass matrix in space
    const SpMatrix<double>& R1() const { return pde_->R1(); }   // discretization of differential operator L
    DMatrix<double> u(std::size_t k) const {                    // discretization of forcing term u at time k
        return u_.block(model().n_basis() * k, 0, model().n_basis(), 1);
    }
    double DeltaT() const { return DeltaT_; }
    const DMatrix<double>& s() const { return s_; }                // initial condition
    std::size_t n_basis() const { return pde_->domain().dof(); }   // number of basis functions

    // setters
    void set_initial_condition(const DMatrix<double>& s) { s_ = s; }
    void set_tolerance(double tol) { tol_ = tol; }
    void set_max_iter(std::size_t max_iter) { max_iter_ = max_iter; }

    // destructor
    virtual ~SpaceTimeParabolicBase() = default;
};

}   // namespace models
}   // namespace fdapde

#endif   // __SPACE_TIME_PARABOLIC_BASE_H__

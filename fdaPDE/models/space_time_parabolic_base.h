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

#include <fdaPDE/pde.h>
#include <fdaPDE/linear_algebra.h>
#include <fdaPDE/utils.h>
using fdapde::core::Kronecker;

#include "space_time_base.h"

namespace fdapde {
namespace models {

// base class for parabolic regularization
template <typename Model>
class SpaceTimeParabolicBase : public SpaceTimeBase<Model, SpaceTimeParabolic> {
   public:
    using PDE = erase<heap_storage, core::PDE__>;
    using Base = SpaceTimeBase<Model, SpaceTimeParabolic>;
    using Base::lambda_D;   // smoothing parameter in space
    using Base::lambda_T;   // smoothing parameter in time
    using Base::model;      // underlying model object
    using Base::time_;      // time interval [0,T]
    using Base::df_;        // model's data
    // constructor
    SpaceTimeParabolicBase() = default;
    SpaceTimeParabolicBase(const PDE& parabolic_penalty) :
        Base(parabolic_penalty.time_domain()), pde_(parabolic_penalty) { }
    // init data structure related to parabolic regularization
    void init_regularization() {
        pde_.init();
        s_ = pde_.initial_condition();   // derive initial condition from parabolic problem
        int m_ = time_.rows();           // number of time points
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
        for (int i = 1; i < m_; ++i) {
            triplet_list.emplace_back(i, i, invDeltaT);
            triplet_list.emplace_back(i, i - 1, -invDeltaT);
        }
        L_.resize(m_, m_);
        L_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        L_.makeCompressed();
        // compute tensorized matrices
        R0_ = Kronecker(Im_, pde_.mass());
        R1_ = Kronecker(Im_, pde_.stiff());
	// correct first n rows of discretized force as (u_1 + R0*s/DeltaT)
	u_ = pde_.force();
	u_.block(0, 0, model().n_basis(), 1) += (1.0 / DeltaT_) * (pde_.mass() * s_);
    }
    // setters
    void set_penalty(const PDE& pde) {
        pde_ = pde;
        model().runtime().set(runtime_status::require_penalty_init);
    }
    // getters
    const PDE& pde() const { return pde_; }   // regularizing term df/dt + Lf - u
    const SpMatrix<double>& R0() const { return R0_; }
    const SpMatrix<double>& R1() const { return R1_; }
    int n_basis() const { return pde_.n_dofs(); }   // number of basis functions
    int n_spatial_basis() const { return pde_.n_dofs(); }
    const SpMatrix<double>& L() const { return L_; }
    const DMatrix<double>& u() const { return u_; }   // discretized force corrected by initial conditions
    const DMatrix<double>& s() { return s_; }         // initial condition
    double DeltaT() const { return DeltaT_; }
    // computes and cache matrices (Im \kron R0)^{-1} and L \kron R0, returns the discretized penalty P =
    // \lambda_D*((Im \kron R1 + \lambda_T*(L \kron R0))^T*(I_m \kron R0)^{-1}*(Im \kron R1 + \lambda_T*(L \kron R0)))
    auto P() {
        if (is_empty(pen_)) {   // compute once and cache result
            invR0_.compute(R0());
            penT_ = Kronecker(L_, pde_.mass());
        }
        return lambda_D() * (R1() + lambda_T() * penT_).transpose() * invR0_.solve(R1() + lambda_T() * penT_);
    }
    // destructor
    virtual ~SpaceTimeParabolicBase() = default;
   protected:
    PDE pde_ {};   // parabolic differential penalty df/dt + Lf - u
    // let m the number of time points
    DMatrix<double> s_;     // N x 1 initial condition vector
    DMatrix<double> u_;     // discretized forcing [1/DeltaT * (u_1 + R_0*s) \ldots u_n]
    SpMatrix<double> Im_;   // m x m sparse identity matrix (assembled once and cached for reuse)
    SpMatrix<double> L_;    // m x m matrix associated with the derivation in time
    double DeltaT_;         // time step (assumes equidistant points in time)
    SpMatrix<double> R0_;   // Im \kron R0 (R0: spatial mass matrix)
    SpMatrix<double> R1_;   // Im \kron R1 (R1: spatial penalty discretization)
    SpMatrix<double> penT_;                      // L_ \kron pde.R0
    fdapde::SparseLU<SpMatrix<double>> invR0_;   // factorization of Im \kron R0
    // discretized penalty: (Im \kron R1 + L \kron R0)^T*(I_m \kron R0)^{-1}*(Im \kron R1 + L \kron R0)
    SpMatrix<double> pen_;
};
  
}   // namespace models
}   // namespace fdapde

#endif   // __SPACE_TIME_PARABOLIC_BASE_H__

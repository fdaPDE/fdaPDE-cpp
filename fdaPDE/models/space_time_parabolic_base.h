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
using fdapde::core::SparseKroneckerTensorProduct;

#include "space_time_base.h"

namespace fdapde {
namespace models {

// base class for parabolic regularization
template <typename Model>
class SpaceTimeParabolicBase : public SpaceTimeBase<Model, SpaceTimeParabolic> {
   protected:
    // let m the number of time points
    DMatrix<double> s_;     // N x 1 initial condition vector
    DMatrix<double> u_;     // discretized forcing [1/DeltaT * (u_1 + R_0*s) \ldots u_n]
    SpMatrix<double> Im_;   // m x m sparse identity matrix (assembled once and cached for reuse)
    SpMatrix<double> L_;    // m x m matrix associated with the derivation in time
    double DeltaT_;         // time step (assumes equidistant points in time)
    SpMatrix<double> R0_;   // Im \kron R0 (R0: spatial mass matrix)
    SpMatrix<double> R1_;   // Im \kron R1 (R1: spatial penalty discretization)

    SpMatrix<double> penT_;                      // discretization of the time derivative: L \kron R0
    fdapde::SparseLU<SpMatrix<double>> invR0_;   // factorization of Im \kron R0
    // discretized penalty: (Im \kron R1 + L \kron R0)^T*(I_m \kron R0)^{-1}*(Im \kron R1 + L \kron R0)
    SpMatrix<double> pen_;
   public:
    typedef SpaceTimeBase<Model, SpaceTimeParabolic> Base;
    using Base::lambda_D;   // smoothing parameter in space
    using Base::lambda_T;   // smoothing parameter in time
    using Base::model;      // underlying model object
    using Base::pde_;       // regularizing term in space
    using Base::time_;      // time interval [0,T]
    using Base::df_;        // model's data

    // constructor
    SpaceTimeParabolicBase() = default;
    SpaceTimeParabolicBase(const pde_ptr& pde, const DVector<double>& time) : Base(pde, time) { }
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
        R0_ = Kronecker(Im_, pde_.mass());
        R1_ = Kronecker(Im_, pde_.stiff());
	// correct first n rows of discretized force as (u_1 + R0*s/DeltaT)
	u_ = pde_.force();
	u_.block(0, 0, model().n_basis(), 1) += (1.0 / DeltaT_) * (pde_.mass() * s_);
    }

    // getters
    const SpMatrix<double>& R0() const { return R0_; }
    const SpMatrix<double>& R1() const { return R1_; }
    std::size_t n_basis() const { return pde_.n_dofs(); }   // number of basis functions
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
    // setters
    // shift = true, cause the removal of the first time instant of data, in case it has been used to estimate the IC
    void set_initial_condition(const DMatrix<double>& s, bool shift = true) {
        s_ = s;
        if (shift) { // left shrink time domain by one step
            std::size_t m = time_.rows();       // number of time instants
            time_ = time_.tail(m - 1).eval();   // correct time interval [0,T] (eval() to avoid aliasing)
            pde_.set_forcing(pde_.forcing_data().rightCols(m - 1));
            model().runtime().set(runtime_status::require_pde_init);   // force pde (re-)initialization
            // remove from data the first time instant, reindex points
            model().set_data(df_.tail(model().n_spatial_locs()).extract(), true);
        }
    }
  
    // destructor
    virtual ~SpaceTimeParabolicBase() = default;
};
  
}   // namespace models
}   // namespace fdapde

#endif   // __SPACE_TIME_PARABOLIC_BASE_H__

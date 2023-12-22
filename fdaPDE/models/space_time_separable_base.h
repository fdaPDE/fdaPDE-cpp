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

#ifndef __SPACE_TIME_SEPARABLE_BASE_H__
#define __SPACE_TIME_SEPARABLE_BASE_H__

#include <fdaPDE/linear_algebra.h>
#include <fdaPDE/splines.h>
#include <fdaPDE/utils.h>
using fdapde::core::Assembler;
using fdapde::core::bilaplacian;
using fdapde::core::Kronecker;
using fdapde::core::reaction;
using fdapde::core::SparseKroneckerTensorProduct;
using fdapde::core::SPLINE;
using fdapde::core::SplineBasis;

namespace fdapde {
namespace models {

// base class for separable regularization
template <typename Model>
class SpaceTimeSeparableBase : public SpaceTimeBase<Model, SpaceTimeSeparable> {
   protected:
    // let \phi_i the i-th basis function in time
    SpMatrix<double> Phi_;   // [Phi_]_{ij} = \phi_i(t_j)
    DVector<double> u_;      // stacked discretized forcing [u_1 \ldots u_n, \ldots, u_1 \ldots u_n]
    SpMatrix<double> R0_;    // P0 \kron R0 (R0: mass matrix in space)
    SpMatrix<double> R1_;    // P0 \kron R1 (R1: stiff matrix in space)

    pde_ptr time_penalty_ {};     // time regularizing term
    DVector<double> time_locs_;   // time instants t_1, ..., t_m
    SpMatrix<double> P_D_;        // discretization of space regularization: (R1^T*R0^{-1}*R1) \kron P0
    SpMatrix<double> P_T_;        // discretization of time regularization:  (R0 \kron P1)
   public:
    typedef SpaceTimeBase<Model, SpaceTimeSeparable> Base;
    using Base::lambda_D;   // smoothing parameter in space
    using Base::lambda_T;   // smoothing parameter in time
    using Base::model;      // underlying model object
    using Base::pde_;       // regularizing term in space Lf = u
    using Base::time_;      // time interval [0,T]

    // constructor
    SpaceTimeSeparableBase() = default;
    SpaceTimeSeparableBase(const pde_ptr& space_penalty, const pde_ptr& time_penalty) :
        Base(space_penalty, time_penalty.dof_coords()), time_penalty_(time_penalty) { }
    // init data structure related to separable regularization
    void init_regularization() {
        time_penalty_.init();   // initialize time-penalty operator
        // compute \Phi = [\Phi]_{ij} = \phi_i(t_j)
        // assume time instants equal to time nodes if not provided
        DVector<double> time_locs = is_empty(time_locs_) ? time_ : time_locs_;
        Phi_ = time_penalty_.eval_basis(core::eval::pointwise, time_locs)->Psi;
        // compute tensorized matrices
        R0_ = Kronecker(time_penalty_.mass(), pde_.mass());
        R1_ = Kronecker(time_penalty_.mass(), pde_.stiff());
    }
    // setters
    void set_temporal_locations(const DVector<double>& time_locations) { time_locs_ = time_locations; }
    void set_time_penalty(const pde_ptr& time_penalty) {
        time_penalty_ = time_penalty;
        model().runtime().set(runtime_status::require_penalty_init);
    }
    // getters
    const SpMatrix<double>& R0() const { return R0_; }
    const SpMatrix<double>& R1() const { return R1_; }
    std::size_t n_basis() const { return pde_.n_dofs() * time_penalty_.n_dofs(); }   // number of basis functions
    std::size_t n_temporal_basis() const { return time_penalty_.n_dofs(); }          // number of time basis functions
    // matrices proper of separable regularization
    const SpMatrix<double>& P0() const { return time_penalty_.mass(); }
    const SpMatrix<double>& P1() const { return time_penalty_.stiff(); }
    const SpMatrix<double>& Phi() const { return Phi_; }
    inline std::size_t n_temporal_locs() const { return is_empty(time_locs_) ? time_.rows() : time_locs_.rows(); }
    const DVector<double>& time_locs() const { return is_empty(time_locs_) ? time_ : time_locs_; }
    const pde_ptr& time_penalty() const { return time_penalty_; }
  
    // return stacked version of discretized forcing field
    const DVector<double>& u() {
        if (is_empty(u_)) {   // compute once and cache
            std::size_t N = Base::n_spatial_basis();
            u_.resize(n_basis());
            // in separable regularization PDE doesn't depend on time. stack forcing term m times
            for (std::size_t i = 0; i < n_temporal_basis(); ++i) { u_.segment(i * N, N) = pde_.force(); }
        }
        return u_;
    }

    // computes and cache matrices (R1^T*R0^{-1}*R1) \kron Rt and R0 \kron Pt.
    // returns the discretized penalty P = \lambda_D*((R1^T*R0^{-1}*R1) \kron Rt) + \lambda_T*(R0 \kron Pt)
    auto P() {
        if (is_empty(P_D_)) {   // compute once and cache result
            fdapde::SparseLU<SpMatrix<double>> invR0_;
            invR0_.compute(pde_.mass());
            P_D_ =
              Kronecker(pde_.stiff().transpose() * invR0_.solve(pde_.stiff()), P0());   // (R1^T*R0^{-1}*R1) \kron P0
            P_T_ = Kronecker(pde_.mass(), P1());                                        // (R0 \kron P1)
        }
        return lambda_D() * P_D_ + lambda_T() * P_T_;
    }
  
    // destructor
    virtual ~SpaceTimeSeparableBase() = default;
};

}   // namespace models
}   // namespace fdapde

#endif   // __SPACE_TIME_SEPARABLE_BASE_H__

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

#include <fdaPDE/pde.h>
#include <fdaPDE/linear_algebra.h>
#include <fdaPDE/utils.h>
using fdapde::core::pde_ptr;
using fdapde::core::Kronecker;

namespace fdapde {
namespace models {

// base class for separable regularization
template <typename Model>
class SpaceTimeSeparableBase : public SpaceTimeBase<Model, SpaceTimeSeparable> {
   protected:
    pde_ptr space_pde_ {};   // regularizing term in space Lf - u
    pde_ptr time_pde_ {};    // regularizing term in time
    // let \phi_i the i-th basis function in time
    SpMatrix<double> Phi_;          // [Phi_]_{ij} = \phi_i(t_j)
    DVector<double> u_;             // stacked discretized forcing [u_1 \ldots u_n, \ldots, u_1 \ldots u_n]
    SpMatrix<double> R0_;           // P0 \kron R0 (R0: mass matrix in space)
    SpMatrix<double> R1_;           // P0 \kron R1 (R1: stiff matrix in space)
    DVector<double> time_locs_;     // time instants t_1, ..., t_m
    mutable SpMatrix<double> PD_;   // discretization of space regularization: (R1^T*R0^{-1}*R1) \kron P0
    mutable SpMatrix<double> PT_;   // discretization of time regularization:  (R0 \kron P1)
   public:
    using Base = SpaceTimeBase<Model, SpaceTimeSeparable>;
    static constexpr int n_lambda = n_smoothing_parameters<SpaceTimeSeparable>::value;
    using Base::lambda_D;   // smoothing parameter in space
    using Base::lambda_T;   // smoothing parameter in time
    using Base::model;      // underlying model object
    using Base::time_;      // time interval [0,T]

    // constructor
    SpaceTimeSeparableBase() = default;
    SpaceTimeSeparableBase(const pde_ptr& space_pde, const pde_ptr& time_pde) :
        Base(time_pde.dof_coords()), space_pde_(space_pde), time_pde_(time_pde) { }
    // init data structure related to separable regularization
    void init_regularization() {
        space_pde_.init();   // initialize penalty in space
        time_pde_.init();    // initialize penalty in time
        // compute \Phi matrix [\Phi]_{ij} = \phi_i(t_j)
        DVector<double> time_locs = is_empty(time_locs_) ? time_ : time_locs_;
        Phi_ = time_pde_.eval_basis(core::eval::pointwise, time_locs)->Psi;
        // compute tensorized matrices
        R0_ = Kronecker(time_pde_.mass(), space_pde_.mass());    // P0 \kron R0
        R1_ = Kronecker(time_pde_.mass(), space_pde_.stiff());   // P0 \kron R1
    }
    // setters
    void set_temporal_locations(const DVector<double>& time_locations) { time_locs_ = time_locations; }
    void set_time_pde(const pde_ptr& time_pde) {
        time_pde_ = time_pde;
        model().runtime().set(runtime_status::require_penalty_init);
    }
    void set_space_pde(const pde_ptr& space_pde) {
        space_pde_ = space_pde;
        model().runtime().set(runtime_status::require_penalty_init);
    }
    // getters
    const SpMatrix<double>& R0() const { return R0_; }
    const SpMatrix<double>& R1() const { return R1_; }
    std::size_t n_basis() const { return space_pde_.n_dofs() * time_pde_.n_dofs(); }
    std::size_t n_spatial_basis() const { return space_pde_.n_dofs(); }   // number of space basis functions
    std::size_t n_temporal_basis() const { return time_pde_.n_dofs(); }   // number of time basis functions
    // matrices proper of separable regularization
    const SpMatrix<double>& P0() const { return time_pde_.mass(); }
    const SpMatrix<double>& P1() const { return time_pde_.stiff(); }
    const SpMatrix<double>& Phi() const { return Phi_; }
    inline std::size_t n_temporal_locs() const { return is_empty(time_locs_) ? time_.rows() : time_locs_.rows(); }
    const DVector<double>& time_locs() const { return is_empty(time_locs_) ? time_ : time_locs_; }
    const pde_ptr& time_pde() const { return time_pde_; }
    const pde_ptr& pde() const { return space_pde_; }
    // return stacked version of discretized forcing field
    const DVector<double>& u() {
        if (is_empty(u_)) {   // compute once and cache
            std::size_t N = n_spatial_basis();
            u_.resize(n_basis());
            // in separable regularization PDE doesn't depend on time. stack forcing term m times
            for (std::size_t i = 0; i < n_temporal_basis(); ++i) { u_.segment(i * N, N) = space_pde_.force(); }
        }
        return u_;
    }
    const SpMatrix<double>& PT() const {   // time-penalty component (P1 \kron R0)
        if (is_empty(PT_)) { PT_ = Kronecker(P1(), space_pde_.mass()); }
        return PT_;
    }
    const SpMatrix<double>& PD() const {   // space-penalty component (P0 \kron (R1^T*R0^{-1}*R1))
        if (is_empty(PD_)) {
            fdapde::SparseLU<SpMatrix<double>> invR0_;
            invR0_.compute(space_pde_.mass());
            PD_ = Kronecker(P0(), space_pde_.stiff().transpose() * invR0_.solve(space_pde_.stiff()));
        }
        return PD_;
    }
    // discretized penalty P = \lambda_D*(P0 \kron (R1^T*R0^{-1}*R1)) + \lambda_T*(P1 \kron R0)
    auto P(const SVector<n_lambda>& lambda) const { return lambda[0] * PD() + lambda[1] * PT(); }
    auto P() const { return P(Base::lambda()); }

    // destructor
    virtual ~SpaceTimeSeparableBase() = default;
};

}   // namespace models
}   // namespace fdapde

#endif   // __SPACE_TIME_SEPARABLE_BASE_H__

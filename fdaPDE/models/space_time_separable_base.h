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
    SpMatrix<double> P0_;    // mass matrix in time: [P0_]_{ij} = \int_{[0,T]} \phi_i*\phi_j
    SpMatrix<double> P1_;    // penalty matrix in time: [P1_]_{ij} = \int_{[0,T]} (\phi_i)_tt*(\phi_j)_tt
    SpMatrix<double> Phi_;   // [Phi_]_{ij} = \phi_i(t_j)
    DVector<double> u_;      // stacked discretized forcing [u_1 \ldots u_n, \ldots, u_1 \ldots u_n]
    SpMatrix<double> R0_;    // P0_ \kron R0 (R0: discretization of the identity operator)
    SpMatrix<double> R1_;    // P0_ \kron R1 (R1: discretization of the differential operator L in the regularizing PDE)

    using TimeBasis = SplineBasis<3>;   // basis used for discretization in time
    TimeBasis basis_;
    DVector<double> time_locs_;   // time instants t_1, ..., t_m
    SpMatrix<double> P_D_;        // discretization of space regularization: (R1^T*R0^{-1}*R1) \kron Rt
    SpMatrix<double> P_T_;        // discretization of time regularization:  (R0 \kron Pt)
   public:
    typedef SpaceTimeBase<Model, SpaceTimeSeparable> Base;
    using Base::lambda_D;   // smoothing parameter in space
    using Base::lambda_T;   // smoothing parameter in time
    using Base::model;      // underlying model object
    using Base::pde_;       // regularizing term in space Lf = u
    using Base::time_;      // time interval [0,T]

    // constructor
    SpaceTimeSeparableBase() = default;
    SpaceTimeSeparableBase(const pde_ptr& pde, const DVector<double>& time) : Base(pde, time) { }
    // init data structure related to separable regularization
    void init_regularization() {
        basis_ = TimeBasis(time_);

        // compute \Phi = [\Phi]_{ij} = \phi_i(t_j)
        // assume time instants equal to time nodes if not provided
        DVector<double> time_locs = is_empty(time_locs_) ? time_ : time_locs_;
        int m = time_locs.rows();
        int M = basis_.size();
        Phi_.resize(m, M);

        // start matrix assembly
        std::vector<fdapde::Triplet<double>> triplet_list;
        triplet_list.reserve(m * M);
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < m; ++j) {   // evaluate basis function at given m time locations
                triplet_list.emplace_back(j, i, basis_[i](SVector<1>(time_locs[j])));
            }
        }
        // finalize construction
        Phi_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        Phi_.prune(0.0);   // remove zeros
        Phi_.makeCompressed();

        // discretize operators in time using spline-approximation
        core::IntegratorTable<1, 3, core::GaussLegendre> quadrature;
        Assembler<SPLINE, DVector<double>, TimeBasis, decltype(quadrature)> assembler(time_, quadrature);
        P0_ = assembler.discretize_operator(reaction<SPLINE>(1.0));   // mass matrix
        P1_ = assembler.discretize_operator(-bilaplacian<SPLINE>());

        // compute tensorized matrices
        R0_ = Kronecker(P0_, pde_.R0());
        R1_ = Kronecker(P0_, pde_.R1());
    }
    // setters
    void set_temporal_locations(const DVector<double>& time_locations) { time_locs_ = time_locations; }
    // getters
    const SpMatrix<double>& R0() const { return R0_; }
    const SpMatrix<double>& R1() const { return R1_; }
    std::size_t n_basis() const { return pde_.n_dofs() * basis_.size(); }   // number of basis functions
    std::size_t n_temporal_basis() const { return basis_.size(); }           // number of time basis functions
    // matrices proper of separable regularization
    const SpMatrix<double>& P0() const { return P0_; }
    const SpMatrix<double>& P1() const { return P1_; }
    const SpMatrix<double>& Phi() const { return Phi_; }
    inline std::size_t n_temporal_locs() const { return is_empty(time_locs_) ? time_.rows() : time_locs_.rows(); }
    const DVector<double>& time_locs() const { return is_empty(time_locs_) ? time_ : time_locs_; }

    // return stacked version of discretized forcing field
    const DVector<double>& u() {
        if (is_empty(u_)) {   // compute once and cache
            std::size_t N = Base::n_spatial_basis();
            u_.resize(n_basis());
            // in separable regularization PDE doesn't depend on time. stack forcing term m times
            for (std::size_t i = 0; i < basis_.size(); ++i) { u_.segment(i * N, N) = pde_.force(); }
        }
        return u_;
    }

    // computes and cache matrices (R1^T*R0^{-1}*R1) \kron Rt and R0 \kron Pt.
    // returns the discretized penalty P = \lambda_D*((R1^T*R0^{-1}*R1) \kron Rt) + \lambda_T*(R0 \kron Pt)
    auto P() {
        if (is_empty(P_D_)) {   // compute once and cache result
            fdapde::SparseLU<SpMatrix<double>> invR0_;
            invR0_.compute(pde_.R0());
            P_D_ = Kronecker(pde_.R1().transpose() * invR0_.solve(pde_.R1()), P0_);   // (R1^T*R0^{-1}*R1) \kron P0
            P_T_ = Kronecker(pde_.R0(), P1_);                                          // (R0 \kron P1)
        }
        return lambda_D() * P_D_ + lambda_T() * P_T_;
    }
  
    // destructor
    virtual ~SpaceTimeSeparableBase() = default;
};

}   // namespace models
}   // namespace fdapde

#endif   // __SPACE_TIME_SEPARABLE_BASE_H__

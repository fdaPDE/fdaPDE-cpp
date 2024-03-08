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

#ifndef __MODEL_TYPE_ERASURE_H__
#define __MODEL_TYPE_ERASURE_H__

#include <fdaPDE/utils.h>
using fdapde::core::BlockFrame;
#include "model_traits.h"
#include "sampling_design.h"

namespace fdapde {
namespace models {

// minimal type erased wrappers for models  
template <typename RegularizationType> struct StatisticalModel__ { };

// base statistical models interface (do not permute BASE_MODEL_FN_PTRS members!)
#define BASE_MODEL_ERASED_INTERFACE                                                                                    \
    void init() { invoke<void, 0>(*this); }                                                                            \
    void solve() { invoke<void, 1>(*this); }                                                                           \
    void set_data(const BlockFrame<double, int>& data, bool reindex = false) {                                         \
        invoke<void, 2>(*this, data, reindex);                                                                         \
    }                                                                                                                  \
    decltype(auto) n_locs()  const { return invoke<int, 3>(*this); }                                                   \
    decltype(auto) n_basis() const { return invoke<int, 4>(*this); }                                                   \
    decltype(auto) data()          { return invoke<BlockFrame<double, int>&, 5>(*this); }                              \
    decltype(auto) data()    const { return invoke<const BlockFrame<double, int>&, 6>(*this); }                        \
    decltype(auto) R0()      const { return invoke<const SpMatrix<double>&, 7>(*this); }                               \
    decltype(auto) R1()      const { return invoke<const SpMatrix<double>&, 8>(*this); }                               \
    decltype(auto) u()       const { return invoke<const DMatrix<double>&, 9>(*this); }                                \
    decltype(auto) Psi()     const { return invoke<const SpMatrix<double>&, 10>(*this, not_nan()); }                   \
    decltype(auto) PsiTD()   const { return invoke<const SpMatrix<double>&, 11>(*this, not_nan()); }

#define BASE_MODEL_FN_PTRS                                                                                             \
    &M::init,     /* full model-stack initialization */                                                                \
    &M::solve,    /* model fit function */                                                                             \
    &M::set_data, /* given a BlockFrame instance, sets model's data and request data stack initialization */           \
    &M::n_locs,   /* number of locations where data are observed, given as n_spatial_basis() * n_temporal_basis() */   \
    &M::n_basis,  /* number of basis functions used in discretization */                                               \
    /* const and non-const access to model's data */                                                                   \
    static_cast<BlockFrame<double, int>& (ModelBase<M>::*)()>(&M::data),                                               \
    static_cast<const BlockFrame<double, int>& (ModelBase<M>::*)() const>(&M::data),                                   \
    &M::R0, /* mass matrix, possibly tensorized */                                                                     \
    &M::R1, /* discretization matrix of the differential operator L, possibly tensorized */                            \
    &M::u,  /* discretization of forcing term, possibly tensorized */                                                  \
    static_cast<const SpMatrix<double>& (SamplingBase<M>::*)(not_nan) const>(&M::Psi),  /* [\Psi]_ij = \psi_j(p_i) */  \
    static_cast<const SpMatrix<double>& (SamplingBase<M>::*)(not_nan) const>(&M::PsiTD) /* \Psi^\top * D */

// penalty independent model interface
template <> struct StatisticalModel__<void> {
    using RegularizationType = void;
    template <typename M> using Base = ModelBase<M>;  
    static constexpr int n_lambda = fdapde::Dynamic;
    // interface implementation
    template <typename M>
    using fn_ptrs = fdapde::mem_fn_ptrs<
      BASE_MODEL_FN_PTRS,
      static_cast<void (Base<M>::*)(const DVector<double>&)>(&M::set_lambda),   // dynamic-sized \lambda setter
      static_cast<DVector<double> (Base<M>::*)(int) const>(&M::lambda)>;        // dynamic-sized access to \lambda
    // interface implementation
    BASE_MODEL_ERASED_INTERFACE
    void set_lambda(const DVector<double>& lambda) { invoke<void, 12>(*this, lambda); }
    decltype(auto) lambda() const { return invoke<DVector<double>, 13>(*this, fdapde::Dynamic); }
};

// space-only statistical model interface
template <> struct StatisticalModel__<SpaceOnly> {
    using RegularizationType = SpaceOnly;
    template <typename M> using Base = SpaceOnlyBase<M>;
    static constexpr int n_lambda = 1;
    // interface implementation
    template <typename M>
    using fn_ptrs = fdapde::mem_fn_ptrs<
      BASE_MODEL_FN_PTRS,
      &M::set_lambda_D,   // sets model's smoothing parameter in space \lambd_D
      static_cast<void (Base<M>::*)(const SVector<1>&)>(&M::set_lambda),   // static-sized \lambda setter
      static_cast<SVector<1> (Base<M>::*)() const>(&M::lambda),            // static-sized access to \lambda
      &M::set_spatial_locations>;   // sets spatial locations p_1, \ldots, p_n over \mathcal{D}
    BASE_MODEL_ERASED_INTERFACE
    void set_lambda_D(double lambda_D) { invoke<void, 12>(*this, lambda_D); }
    void set_lambda(const SVector<1>& lambda) { invoke<void, 13>(*this, lambda); }
    SVector<1> lambda() const { return invoke<SVector<1>, 14>(*this); }
    void set_spatial_locations(const DMatrix<double>& locs) { invoke<void, 15>(*this, locs); }
};

// space-time separable interface
template <> struct StatisticalModel__<SpaceTimeSeparable> {
    using RegularizationType = SpaceTimeSeparable;
    template <typename M> using Base = SpaceTimeBase<M, SpaceTimeSeparable>;
    static constexpr int n_lambda = 2;
    // interface implementation
    template <typename M>
    using fn_ptrs = fdapde::mem_fn_ptrs<
      BASE_MODEL_FN_PTRS,
      &M::set_lambda_D,   // sets model's smoothing level in space \lambda_D
      &M::set_lambda_T,   // sets model's smoothing level in time \lambda_T
      static_cast<void (Base<M>::*)(const SVector<2>&)>(&M::set_lambda),   // static-sized \lambda vector setter
      static_cast<SVector<2> (Base<M>::*)() const>(&M::lambda),            // static-sized access to \lambda vector
      &M::set_spatial_locations,     // sets spatial locations p_1, \ldots, p_n over \mathcal{D}
      &M::set_temporal_locations>;   // sets temporal locations t_1, \ldots, t_m over T
    BASE_MODEL_ERASED_INTERFACE
    void set_lambda_D(double lambda_D) { invoke<void, 12>(*this, lambda_D); }
    void set_lambda_T(double lambda_T) { invoke<void, 13>(*this, lambda_T); }
    void set_lambda(const SVector<2>& lambda) { invoke<void, 14>(*this, lambda); }
    SVector<2> lambda() const { return invoke<SVector<2>, 15>(*this); }
    void set_spatial_locations(const DMatrix<double>& locs) { invoke<void, 16>(*this, locs); }
    void set_temporal_locations(const DMatrix<double>& locs) { invoke<void, 17>(*this, locs); }
};

// space-time parabolic interface
template <> struct StatisticalModel__<SpaceTimeParabolic> {
    using RegularizationType = SpaceTimeParabolic;
    template <typename M> using Base = SpaceTimeBase<M, SpaceTimeParabolic>;
    static constexpr int n_lambda = 2;
    // interface implementation
    template <typename M>
    using fn_ptrs = fdapde::mem_fn_ptrs<
      BASE_MODEL_FN_PTRS,
      &M::set_lambda_D,   // sets model's smoothing level in space \lambda_D
      &M::set_lambda_T,   // sets model's smoothing level in time \lambda_T
      static_cast<void (Base<M>::*)(const SVector<2>&)>(&M::set_lambda),   // static-sized \lambda vector setter
      static_cast<SVector<2> (Base<M>::*)() const>(&M::lambda),            // static-sized access to \lambda vector
      &M::set_spatial_locations>;   // sets spatial locations p_1, \ldots, p_n over \mathcal{D}
    BASE_MODEL_ERASED_INTERFACE
    void set_lambda_D(double lambda_D) { invoke<void, 12>(*this, lambda_D); }
    void set_lambda_T(double lambda_T) { invoke<void, 13>(*this, lambda_T); }
    void set_lambda(const SVector<2>& lambda) { invoke<void, 14>(*this, lambda); }
    SVector<2> lambda() const { return invoke<SVector<2>, 15>(*this); }
    void set_spatial_locations(const DMatrix<double>& locs) { invoke<void, 16>(*this, locs); }
};

}   // namespace models
}   // namespace fdapde

#endif   // __MODEL_TYPE_ERASURE_H__

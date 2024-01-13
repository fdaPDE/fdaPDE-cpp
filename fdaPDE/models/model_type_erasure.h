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

namespace fdapde {
namespace models {

// minimal type erased wrappers for models  
template <typename RegularizationType> struct IStatModel { };

// base statistical models interface (do not permute fn_ptrs members!)
#define DEFINE_BASE_STAT_MODEL_INTERFACE                                                                               \
    void init() { fdapde::invoke<void, 0>(*this); }                                                                    \
    void solve() { fdapde::invoke<void, 1>(*this); }                                                                   \
    void set_lambda_D(double lambda) { fdapde::invoke<void, 2>(*this, lambda); }                                       \
    void set_data(const BlockFrame<double, int>& data, bool reindex = false) {                                         \
        fdapde::invoke<void, 3>(*this, data, reindex);                                                                 \
    }                                                                                                                  \
    void set_spatial_locations(const DMatrix<double>& locs) { fdapde::invoke<void, 4>(*this, locs); }                  \
    BlockFrame<double, int>& data() { return fdapde::invoke<BlockFrame<double, int>&, 5>(*this); }                     \
    decltype(auto) R0() const { return fdapde::invoke<const SpMatrix<double>&, 6>(*this); }                            \
    std::size_t n_locs() const { return fdapde::invoke<std::size_t, 7>(*this); }                                       \
    std::size_t n_basis() const { return fdapde::invoke<std::size_t, 8>(*this); }

// space-only statistical model interface
template <> struct IStatModel<SpaceOnly> {
    using RegularizationType = SpaceOnly;
    template <typename M>
    using fn_ptrs = fdapde::mem_fn_ptrs<
      &M::init, &M::solve,   // initializtion and solution of the regression problem
      &M::set_lambda_D, &M::set_data, &M::set_spatial_locations,
      static_cast<BlockFrame<double, int>& (ModelBase<M>::*)()>(&M::data), &M::R0, &M::n_locs, &M::n_basis,
      static_cast<void (SpaceOnlyBase<M>::*)(const SVector<1>&)>(&M::set_lambda)>;
    // interface implementation
    DEFINE_BASE_STAT_MODEL_INTERFACE;
    void set_lambda(const SVector<1>& lambda) { fdapde::invoke<void, 9>(*this, lambda); }
};

// space-time separable statistical model interface
template <> struct IStatModel<SpaceTimeSeparable> {
    using RegularizationType = SpaceTimeSeparable;
    template <typename M>
    using fn_ptrs = fdapde::mem_fn_ptrs<
      &M::init, &M::solve,   // initializtion and solution of the regression problem
      &M::set_lambda_D, &M::set_data, &M::set_spatial_locations,
      static_cast<BlockFrame<double, int>& (ModelBase<M>::*)()>(&M::data), &M::R0, &M::n_locs, &M::n_basis,
      static_cast<void (SpaceTimeBase<M, SpaceTimeSeparable>::*)(const SVector<2>&)>(&M::set_lambda), &M::set_lambda_T,
      &M::set_temporal_locations>;
    // interface implementation
    DEFINE_BASE_STAT_MODEL_INTERFACE;
    void set_lambda(const SVector<2>& lambda) { fdapde::invoke<void, 9>(*this, lambda); }
    void set_lambda_T(double lambda_T) { fdapde::invoke<void, 10>(*this, lambda_T); }
    void set_temporal_locations(const DMatrix<double>& locs) { fdapde::invoke<void, 11>(*this, locs); }
};

// space-time parabolic statistical model interface
template <> struct IStatModel<SpaceTimeParabolic> {
    using RegularizationType = SpaceTimeParabolic;
    template <typename M>
    using fn_ptrs = fdapde::mem_fn_ptrs<
      &M::init, &M::solve,   // initializtion and solution of the regression problem
      &M::set_lambda_D, &M::set_data, &M::set_spatial_locations,
      static_cast<BlockFrame<double, int>& (ModelBase<M>::*)()>(&M::data), &M::R0, &M::n_locs, &M::n_basis,
      static_cast<void (SpaceTimeBase<M, SpaceTimeParabolic>::*)(const SVector<2>&)>(&M::set_lambda), &M::set_lambda_T>;
    // interface implementation
    DEFINE_BASE_STAT_MODEL_INTERFACE;
    void set_lambda(const SVector<2>& lambda) { fdapde::invoke<void, 9>(*this, lambda); }
    void set_lambda_T(double lambda_T) { fdapde::invoke<void, 10>(*this, lambda_T); }
};

// regularization-independent model interface
template <> struct IStatModel<void> {
    template <typename M>
    using fn_ptrs = fdapde::mem_fn_ptrs<
      &M::init, &M::solve, static_cast<void (ModelBase<M>::*)(const DVector<double>&)>(&M::set_lambda),
      static_cast<DVector<double> (ModelBase<M>::*)(int) const>(&M::lambda),
      static_cast<const SpMatrix<double>& (SamplingBase<M>::*)(not_nan) const>(&M::Psi),
      static_cast<const SpMatrix<double>& (SamplingBase<M>::*)(not_nan) const>(&M::PsiTD), &M::set_data,
      &M::n_locs, &M::n_basis>;
    // interface implementation
    void init()  { fdapde::invoke<void, 0>(*this); }
    void solve() { fdapde::invoke<void, 1>(*this); }
    void set_lambda(const DVector<double>& lambda) { fdapde::invoke<void, 2>(*this, lambda); }
    decltype(auto) lambda() const { return fdapde::invoke<DVector<double>, 3>(*this, fdapde::Dynamic); }
    decltype(auto) Psi()    const { return fdapde::invoke<const SpMatrix<double>&, 4>(*this, not_nan()); }
    decltype(auto) PsiTD()  const { return fdapde::invoke<const SpMatrix<double>&, 5>(*this, not_nan()); }
    void set_data(const BlockFrame<double, int>& data, bool reindex = false) {
        fdapde::invoke<void, 6>(*this, data, reindex);
    }
    std::size_t n_locs() const { return fdapde::invoke<std::size_t, 7>(*this); }
    std::size_t n_basis() const { return fdapde::invoke<std::size_t, 8>(*this); }
};

}   // namespace models
}   // namespace fdapde

#endif   // __MODEL_TYPE_ERASURE_H__

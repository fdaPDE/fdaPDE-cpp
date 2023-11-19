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

#ifndef __MODEL_WRAPPERS_H__
#define __MODEL_WRAPPERS_H__

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
    BlockFrame<double, int>& data() { return fdapde::invoke<BlockFrame<double, int>&, 5>(*this); }

// space-only statistical model interface
template <> struct IStatModel<SpaceOnly> {
    template <typename M>
    using fn_ptrs = fdapde::mem_fn_ptrs<
      &M::init, &M::solve,   // initializtion and solution of the regression problem
      &M::set_lambda_D, &M::set_data, &M::set_spatial_locations,
      static_cast<BlockFrame<double, int>& (ModelBase<M>::*)()>(&M::data),
      &M::set_lambda>;
    // interface implementation
    DEFINE_BASE_STAT_MODEL_INTERFACE;
    void set_lambda(const SVector<1>& lambda) { fdapde::invoke<void, 6>(*this, lambda); }
};

// space-time separable statistical model interface
template <> struct IStatModel<SpaceTimeSeparable> {
    template <typename M>
    using fn_ptrs = fdapde::mem_fn_ptrs<
      &M::init, &M::solve,   // initializtion and solution of the regression problem
      &M::set_lambda_D, &M::set_data, &M::set_spatial_locations,
      static_cast<BlockFrame<double, int>& (ModelBase<M>::*)()>(&M::data),
      &M::set_lambda, &M::set_lambda_T, &M::set_temporal_locations>;
    // interface implementation
    DEFINE_BASE_STAT_MODEL_INTERFACE;
    void set_lambda(const SVector<2>& lambda) { fdapde::invoke<void, 6>(*this, lambda); }
    void set_lambda_T(double lambda_T) { fdapde::invoke<void, 7>(*this, lambda_T); }
    void set_temporal_locations(const DMatrix<double>& locs) { fdapde::invoke<void, 8>(*this, locs); }
};

// space-time parabolic statistical model interface
template <> struct IStatModel<SpaceTimeParabolic> {
    template <typename M>
    using fn_ptrs = fdapde::mem_fn_ptrs<
      &M::init, &M::solve,   // initializtion and solution of the regression problem
      &M::set_lambda_D, &M::set_data, &M::set_spatial_locations,
      static_cast<BlockFrame<double, int>& (ModelBase<M>::*)()>(&M::data),
      &M::set_lambda, &M::set_lambda_T, &M::set_initial_condition>;
    // interface implementation
    DEFINE_BASE_STAT_MODEL_INTERFACE;
    void set_lambda(const SVector<2>& lambda) { fdapde::invoke<void, 6>(*this, lambda); }
    void set_lambda_T(double lambda_T) { fdapde::invoke<void, 7>(*this, lambda_T); }
    void set_initial_condition(const DVector<double>& s, bool shift = true) {
        fdapde::invoke<void, 8>(*this, s, shift);
    }
};
  
}   // namespace models
}   // namespace fdapde

#endif   // __MODEL_WRAPPERS_H__

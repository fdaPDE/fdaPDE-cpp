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

#ifndef __MODEL_TRAITS_H__
#define __MODEL_TRAITS_H__

#include <fdaPDE/utils.h>
using fdapde::is_base_of_template;

namespace fdapde {
namespace models {
  
// supported regularization strategies, forward declarations
template <typename Model> class ModelBase;                             // base class for any fdaPDE statistical model
template <typename Model> class SpaceOnlyBase;                         // base class for any spatial model
template <typename Model, typename PenaltyType> class SpaceTimeBase;   // base class for any spatio-temporal model
template <typename Model> class SpaceTimeSeparableBase;                // spatio-temporal, separable regularization
template <typename Model> class SpaceTimeParabolicBase;                // spatio-temporal, parabolic regularization

// empty classes for tagging regularization types
struct SpaceOnly { };
struct SpaceTime { };   // generic space-time tag
struct SpaceTimeSeparable { };
struct SpaceTimeParabolic { };

// trait to detect if a type inherits ModelBase (can be regarded as a statistical model).
template <typename T> struct is_stat_model {
    static constexpr bool value = fdapde::is_base_of_template<ModelBase, T>::value;
};

// space-only regularization trait
template <typename Model> struct is_space_only {
    static constexpr bool value = std::is_same<typename std::decay<Model>::type::RegularizationType, SpaceOnly>::value;
};
// traits for space-time regularizations
template <typename Model> struct is_space_time {
    static constexpr bool value = !is_space_only<Model>::value;
};
template <typename Model> struct is_space_time_separable {
    static constexpr bool value =
      std::is_same<typename std::decay<Model>::type::RegularizationType, SpaceTimeSeparable>::value;
};
template <typename Model> struct is_space_time_parabolic {
    static constexpr bool value =
      std::is_same<typename std::decay<Model>::type::RegularizationType, SpaceTimeParabolic>::value;
};
  
template <typename Model_> struct has_single_penalty {
    using Model = typename std::decay<Model_>::type;
    static constexpr bool value = (is_space_only<Model>::value || is_space_time_parabolic<Model>::value);
};
template <typename Model_> struct has_double_penalty {
    static constexpr bool value = is_space_time_separable<typename std::decay<Model_>::type>::value;
};

// selects the regularization base depending on the Model regularization type
template <typename Model, typename RegularizationType> struct select_regularization_base {
    using type = typename std::conditional<
      std::is_same<RegularizationType, SpaceOnly>::value, SpaceOnlyBase<Model>,
      typename std::conditional<
        std::is_same<RegularizationType, SpaceTimeSeparable>::value, SpaceTimeSeparableBase<Model>,
        SpaceTimeParabolicBase<Model>>::type>::type;
};

// selects the number of smoothing parameters
template <typename RegularizationType> struct n_smoothing_parameters {
    static constexpr int value = []() constexpr -> int {
        if constexpr (std::is_same<RegularizationType, SpaceOnly>::value) {
            return 1;
        } else {
            return 2;
        }    
    }();
};

}   // namespace models
}   // namespace fdapde

#endif   // __MODEL_TRAITS_H__

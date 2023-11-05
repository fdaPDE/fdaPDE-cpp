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

// base class for model traits
template <typename B> struct model_traits;

// supported resolution strategies for the smoothing problem
struct MonolithicSolver { };
template <typename Model> struct is_solver_monolithic {
    static constexpr bool value = std::is_same<typename model_traits<Model>::solver, MonolithicSolver>::value;
};
struct IterativeSolver { };
template <typename Model> struct is_solver_iterative {
    static constexpr bool value = std::is_same<typename model_traits<Model>::solver, IterativeSolver>::value;
};

// supported regularization strategies, forward declarations
template <typename Model> class ModelBase;       // base class for any fdaPDE statistical model
template <typename Model> class SpaceOnlyBase;   // base class for any spatial model
template <typename Model> class SpaceTimeBase;   // base class for any spatio-temporal model
template <typename Model, typename Solver> class SpaceTimeSeparableBase;   // spatio-temporal, separable regularization
template <typename Model, typename Solver> class SpaceTimeParabolicBase;   // spatio-temporal, parabolic regularization

// empty classes for tagging regularization types
struct SpaceOnly { };
struct SpaceTime { };   // generic space-time tag
struct SpaceTimeSeparable { };
struct SpaceTimeParabolic { };

// trait to detect if a type inherits ModelBase (can be regarded as a statistical model).
template <typename T> struct is_stat_model {
    static constexpr bool value = fdapde::is_base_of_template<ModelBase, T>::value;
};

template <typename Model> struct is_space_only {
    static constexpr bool value =
      std::is_same<typename model_traits<typename std::decay<Model>::type>::regularization, SpaceOnly>::value;
};
template <typename Model> struct is_space_time {
    static constexpr bool value = !is_space_only<Model>::value;
};
// specific traits for space-time regularizations
template <typename Model> struct is_space_time_separable {   // separable regularization
    static constexpr bool value =
      std::is_same<typename model_traits<typename std::decay<Model>::type>::regularization, SpaceTimeSeparable>::value;
};
template <typename Model> struct is_space_time_parabolic {   // parabolic regularization
    static constexpr bool value =
      std::is_same<typename model_traits<typename std::decay<Model>::type>::regularization, SpaceTimeParabolic>::value;
};

// selects the regularization base depending on the Model regularization type
template <typename Model> struct select_regularization_type {
    using type = typename std::conditional<
      is_space_only<Model>::value, SpaceOnlyBase<Model>,
      typename std::conditional<
        is_space_time_separable<Model>::value, SpaceTimeSeparableBase<Model, typename model_traits<Model>::solver>,
        SpaceTimeParabolicBase<Model, typename model_traits<Model>::solver>>::type>::type;
};

// selects the number of smoothing parameters
template <typename Regularization> class n_smoothing_parameters {
    static constexpr int compute() {
        if constexpr (std::is_same<typename std::decay<Regularization>::type, SpaceOnly>::value) {
            return 1;
        } else {
            return 2;
        }
    }
   public:
    static constexpr int value = n_smoothing_parameters<Regularization>::compute();
};

// allowed sampling strategies
struct GeoStatLocations { };   // general sampling locations p_1, ..., p_n
struct GeoStatMeshNodes { };   // sampling locations p_1, ..., p_n coincide with mesh nodes
struct Areal { };              // areal sampling at subdomains D_1, ..., D_n
// traits for sampling design in space
template <typename Model> struct is_sampling_areal {
    static constexpr bool value = std::is_same<typename model_traits<Model>::sampling, Areal>::value;
};
template <typename Model> struct is_sampling_pointwise_at_mesh {
    static constexpr bool value = std::is_same<typename model_traits<Model>::sampling, GeoStatMeshNodes>::value;
};
template <typename Model> struct is_sampling_pointwise_at_locs {
    static constexpr bool value = std::is_same<typename model_traits<Model>::sampling, GeoStatLocations>::value;
};

}   // namespace models
}   // namespace fdapde

#endif   // __MODEL_TRAITS_H__

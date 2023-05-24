#ifndef __MODEL_TRAITS_H__
#define __MODEL_TRAITS_H__

#include "../core/utils/Traits.h"
using fdaPDE::is_base_of_template;

namespace fdaPDE{
namespace models{

  // base class for model traits
  template <typename B> struct model_traits;
  
  // supported resolution strategies for the smoothing problem
  struct MonolithicSolver {};
  template <typename Model>
  struct is_solver_monolithic { 
    static constexpr bool value = std::is_same<
      typename model_traits<Model>::solver,
      MonolithicSolver>::value;
  };
  struct IterativeSolver{};
  template <typename Model>
  struct is_solver_iterative { 
    static constexpr bool value = std::is_same<
      typename model_traits<Model>::solver,
      IterativeSolver>::value;
  };

  // supported regularization strategies, forward declarations
  template <typename Model> class ModelBase; // base class for any fdaPDE statistical model
  template <typename Model> class SpaceOnlyBase; // base class for any spatial model
  template <typename Model> class SpaceTimeBase; // base class for any spatio-temporal model
  template <typename Model, typename Solver> class SpaceTimeSeparableBase; // spatio-temporal, separable regularization
  template <typename Model, typename Solver> class SpaceTimeParabolicBase; // spatio-temporal, parabolic regularization
  
  // empty classes for tagging regularization types
  struct SpaceOnly {};
  struct SpaceTimeSeparable {};
  struct SpaceTimeParabolic {};
    
  // trait to detect if a type implements ModelBase (can be regarded as a statistical model). This is the least
  // constraining requirement an algorithm can require on a type
  template <typename T>
  struct is_stat_model { static constexpr bool value = fdaPDE::is_base_of_template<ModelBase, T>::value; };  

  // traits to detect if a model is space-only or space-time 
  template <typename Model>
  struct is_space_only {
    static constexpr bool value = std::is_same<
      typename model_traits<typename std::decay<Model>::type>::regularization,
      SpaceOnly>::value;
  };
  template <typename Model>
  struct is_space_time { static constexpr bool value = !is_space_only<Model>::value; };
  // specific traits for space-time regularizations
  template <typename Model>
  struct is_space_time_separable { // separable regularization
    static constexpr bool value = std::is_same<
      typename model_traits<typename std::decay<Model>::type>::regularization,
      SpaceTimeSeparable>::value;
  };
  template <typename Model>
  struct is_space_time_parabolic { // parabolic regularization
    static constexpr bool value = std::is_same<
      typename model_traits<typename std::decay<Model>::type>::regularization,
      SpaceTimeParabolic>::value;
  };
  
  // selects the regularization type for Model
  template <typename Model>
  struct select_regularization_type {
    using type = typename std::conditional<
      is_space_only<Model>::value, SpaceOnlyBase<Model>,
      typename std::conditional<
	is_space_time_separable<Model>::value,				
	SpaceTimeSeparableBase<Model, typename model_traits<Model>::solver>,
	SpaceTimeParabolicBase<Model, typename model_traits<Model>::solver>>::type
      >::type;
  };

  // selects the number of smoothing parameters
  template <typename Regularization>
  class n_smoothing_parameters {
    static constexpr int compute() {
      // space only model have only one smoothness level in space
      if constexpr(std::is_same<typename std::decay<Regularization>::type, SpaceOnly>::value) return 1;
      else return 2; // space-time models regularize both in time and space
    }
  public:
    static constexpr int value = n_smoothing_parameters<Regularization>::compute();
  };

  // allowed sampling strategies
  struct GeoStatLocations {}; // sampling locations p_1, ..., p_n are provided with no particular structure in the domain
  struct GeoStatMeshNodes {}; // sampling locations p_1, ..., p_n coincide with mesh nodes
  struct Areal {};            // subdomains D_1, ..., D_n are provided
  // traits for sampling design in space
  template <typename Model>
  struct is_sampling_areal {
    static constexpr bool value = std::is_same<
      typename model_traits<Model>::sampling, Areal
      >::value;
  };
  template <typename Model>
  struct is_sampling_pointwise_at_mesh {
    static constexpr bool value = std::is_same<
      typename model_traits<Model>::sampling, GeoStatMeshNodes
      >::value;
  };
  template <typename Model>
  struct is_sampling_pointwise_at_locs {
    static constexpr bool value = std::is_same<
      typename model_traits<Model>::sampling, GeoStatLocations
      >::value;
  };
  
}}

#endif // __MODEL_TRAITS_H__

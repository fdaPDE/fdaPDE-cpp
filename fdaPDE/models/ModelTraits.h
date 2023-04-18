#ifndef __MODEL_TRAITS_H__
#define __MODEL_TRAITS_H__

#include "../core/utils/Traits.h"
using fdaPDE::is_base_of_template;

namespace fdaPDE{
namespace models{

  // base class for model traits
  template <typename B> struct model_traits;
  
  // possible solution strategies for a model
  enum SolverType { Monolithic, Iterative };
  // traits for detection of solution strategy
  template <typename Model>
  struct is_solver_monolithic { 
    static constexpr bool value = model_traits<Model>::solver == SolverType::Monolithic; };
  template <typename Model>
  struct is_solver_iterative { 
    static constexpr bool value = model_traits<Model>::solver == SolverType::Iterative; };

  // supported regularization strategies, forward declarations
  template <typename Model> class ModelBase; // base class for any fdaPDE statistical model
  template <typename Model> class SpaceOnlyBase; // base class for any spatial model
  template <typename Model> class SpaceTimeBase; // base class for any spatio-temporal model
  template <typename Model, SolverType solver> class SpaceTimeSeparableBase; // spatio-temporal, separable regularization
  template <typename Model, SolverType Solver> class SpaceTimeParabolicBase; // spatio-temporal, parabolic regularization
  
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
      typename model_traits<typename std::decay<Model>::type>::RegularizationType,
      SpaceOnly>::value;
  };
  template <typename Model>
  struct is_space_time { static constexpr bool value = !is_space_only<Model>::value; };
  // specific traits for space-time regularizations
  template <typename Model>
  struct is_space_time_separable { // separable regularization
    static constexpr bool value = std::is_same<
      typename model_traits<typename std::decay<Model>::type>::RegularizationType,
      SpaceTimeSeparable>::value;
  };
  template <typename Model>
  struct is_space_time_parabolic { // parabolic regularization
    static constexpr bool value = std::is_same<
      typename model_traits<typename std::decay<Model>::type>::RegularizationType,
      SpaceTimeParabolic>::value;
  };
  
  // selects the regularization type for Model
  template <typename Model>
  struct select_regularization_type {
    using type = typename std::conditional<
      is_space_only<Model>::value, SpaceOnlyBase<Model>,
      typename std::conditional<
	is_space_time_separable<Model>::value,				
	SpaceTimeSeparableBase<Model, model_traits<Model>::solver>,
	SpaceTimeParabolicBase<Model, model_traits<Model>::solver>>::type
      >::type;
  };

  // selects the number of smoothing parameters given a regularization
  template <typename Regularization>
  class n_smoothing_parameters {
    static constexpr int compute() {
      if constexpr(std::is_same<typename std::decay<Regularization>::type, SpaceOnly>::value) return 1;
      else return 2;
    }
  public:
    static constexpr int value = n_smoothing_parameters<Regularization>::compute();
  };

  // allowed sampling strategies
  enum Sampling { GeoStatLocations, GeoStatMeshNodes, Areal };
  // traits for sampling design in space  
  template <typename Model>
  struct is_sampling_areal { 
    static constexpr bool value = model_traits<Model>::sampling == Sampling::Areal; };
  template <typename Model>
  struct is_sampling_pointwise_at_mesh { 
    static constexpr bool value = model_traits<Model>::sampling == Sampling::GeoStatMeshNodes; };
  template <typename Model>
  struct is_sampling_pointwise_at_locs { 
    static constexpr bool value = model_traits<Model>::sampling == Sampling::GeoStatLocations; };
      
}}

#endif // __MODEL_TRAITS_H__

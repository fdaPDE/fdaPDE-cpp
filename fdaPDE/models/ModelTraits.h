#ifndef __MODEL_TRAITS_H__
#define __MODEL_TRAITS_H__

#include "../core/utils/Traits.h"
using fdaPDE::is_base_of_template;

namespace fdaPDE{
namespace models{

  // collection of traits used across the models module

  // possible solution strategies for a model
  enum SolverType { Monolithic, Iterative };
  
  template <typename Model> class ModelBase; // base class for any fdaPDE statistical model
  template <typename Model> class SpaceOnlyBase; // base class for any spatial model
  template <typename Model> class SpaceTimeBase; // base class for any spatio-temporal model
  template <typename Model> class SpaceTimeSeparableBase; // spatio-temporal model with separable regularization
  template <typename Model, SolverType Solver> class SpaceTimeParabolicBase; // spatio-temporal model with parabolic regularization
  
  // empty classes for tagging regularization types
  struct SpaceOnlyTag {};
  struct SpaceTimeSeparableTag {};
  struct SpaceTimeParabolicTag {};
  
  // base class for model traits
  template <typename B> struct model_traits;
  
  // trait to detect if a type implements ModelBase (can be regarded as a statistical model)
  template <typename T>
  struct is_stat_model {
    static constexpr bool value = fdaPDE::is_base_of_template<ModelBase, T>::value;
  };  
  
  // trait for the selection of the type of regularization on the base of the property of a model
  template <typename Model>
  struct select_regularization_type {
    using type = typename std::conditional<
      std::is_same<typename model_traits<Model>::RegularizationType, SpaceOnlyTag>::value,
      SpaceOnlyBase<Model>,
      typename std::conditional<
        std::is_same<typename model_traits<Model>::RegularizationType, SpaceTimeSeparableTag>::value,
        SpaceTimeSeparableBase<Model>,
        SpaceTimeParabolicBase<Model, model_traits<Model>::solver>>::type
      >::type;
  };

  // trait to detect if a model is a space-time model
  template <typename Model>
  struct is_space_time {
    static constexpr bool value = !std::is_same<
      typename model_traits<typename std::decay<Model>::type>::RegularizationType, SpaceOnlyTag>::value;
  };

  // trait to detect if a model has a non-gaussian error distribution
  class Gaussian; // tag used for distinguish a generalized model from a non-generalized one
  template <typename Model>
  struct is_generalized {
    static constexpr bool value = !std::is_same<
      typename model_traits<typename std::decay<Model>::type>::DistributionType, Gaussian>::value;
  };
  
  // trait to select the number of smoothing parameters
  template <typename Model>
  class n_smoothing_parameters {
    static constexpr int compute() {
      if constexpr(is_space_time<Model>::value) return 2;
      else return 1;
    }
  public:
    static constexpr int value = n_smoothing_parameters<Model>::compute();
  };
  
  // allowed sampling strategies
  enum Sampling { GeoStatLocations, GeoStatMeshNodes, Areal };
  // traits for sampling design in space  
  template <typename Model>
  struct is_sampling_areal { 
    static constexpr bool value = model_traits<Model>::sampling == Sampling::Areal;
  };
  template <typename Model>
  struct is_sampling_pointwise_at_mesh { 
    static constexpr bool value = model_traits<Model>::sampling == Sampling::GeoStatMeshNodes;
  };
  template <typename Model>
  struct is_sampling_pointwise_at_locs { 
    static constexpr bool value = model_traits<Model>::sampling == Sampling::GeoStatLocations;
  };
  
  // macros for the import of common symbols to avoid long annoying lists of using declarations in model implemetations

  // this macro is intended to import all **common** symbols a model type can expect from its parent classes
#define IMPORT_MODEL_SYMBOLS				                                 \
  using Base::y;       /* vector of observations y = [y_1 ... y_n] */                    \
  using Base::n_obs;   /* number of observations n */                                    \
  using Base::n_basis; /* number of basis function for discretization in space N */      \
  using Base::n_locs;  /* number of locations p_1 ... p_n where data are observed */     \
  using Base::Psi;     /* n x N matrix of spatial basis evaluations at p_1 ... p_n */    \
  using Base::PsiTD;   /* block P^T*D, being D the matrix of subdomains' measure */      \
                       /* returns P^T if sampling is not areal */			 \
  using Base::R1;      /* discretization of differential operator L (tensorized for */   \
                       /* space-time problems) */					 \
  using Base::R0;      /* mass matrix in space (tensorized for space-time problems) */   \
  using Base::u;       /* discretization of forcing term */			         \
  using Base::pde;     /* differential operator L (regularizing term) */                 \
  
  // this macro is intended to import all **common** symbols a model can expect from a Regression base
  // symbols specific for the regularization type used need to be imported via dedicated using declaration
#define IMPORT_REGRESSION_SYMBOLS					                 \
  IMPORT_MODEL_SYMBOLS;			              			                 \
  /* data access */							                 \
  using Base::W;             /* matrix of observation weights W_ = diag[W_1 ... W_n] */  \
  using Base::q;	     /* number of covariates */			                 \
  using Base::X;	     /* n x q design matrix X = [X_1 ... X_q] */                 \
  using Base::XtWX;  	     /* q x q dense matrix X^T*W*X */		                 \
  using Base::invXtWX;	     /* partialPivLU factorization of X^T*W*X */                 \
  /* utilities */							                 \
  using Base::hasCovariates; /* true if the model is semi-parametric */	                 \
  using Base::hasWeights;    /* true if heteroscedastic observations are assumed */      \
  using Base::lmbQ;	     /* efficient left multiplication by Q */	                 \
  /* room for problem solution */					                 \
  using Base::f_;            /* estimate of the nonparametric part of the model */       \
  using Base::g_;            /* PDE misfit */				                 \
  using Base::beta_;         /* estimate of coefficient vector for parametric part */    \

  // standardized definitions for stat model BlockFrame. layers below will make heavy assumptions on
  // the layout of the BlockFrame, use these instead of manually typing the block name when accessing df_
#define OBSERVATIONS_BLK "y" // matrix of observations
#define INDEXES_BLK "i"      // vector of observation indices
  
}}

#endif // __MODEL_TRAITS_H__

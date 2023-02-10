#ifndef __FUNCTIONAL_BASE_H__
#define __FUNCTIONAL_BASE_H__

#include "../../core/utils/Symbols.h"
#include "../ModelBase.h"
#include "../SamplingDesign.h"
#include "../ModelTraits.h"
using fdaPDE::models::select_regularization_type;

namespace fdaPDE {
namespace models {

  // base class for any *functional* fdaPDE model
  template <typename Model>
  class FunctionalBase
    : public select_regularization_type<Model>::type, public SamplingDesign<Model, model_traits<Model>::sampling> {
  public:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    typedef typename select_regularization_type<Model>::type Base;
    using Base::pde_; // differential operator L 
    
    FunctionalBase() = default;
    // space-only constructor
    template <typename... SamplingData,
	      typename U = Model, // fake type to enable substitution in SFINAE
	      typename std::enable_if<
		std::is_same<typename model_traits<U>::RegularizationType, SpaceOnlyTag>::value,
		int>::type = 0> 
    FunctionalBase(const PDE& pde, const SamplingData&... s) 
      : select_regularization_type<Model>::type(pde),
      SamplingDesign<Model, model_traits<Model>::sampling>(s...) {};

    // space-time constructor
    template <typename... SamplingData,
	      typename U = Model, // fake type to enable substitution in SFINAE
	      typename std::enable_if<!
		std::is_same<typename model_traits<U>::RegularizationType, SpaceOnlyTag>::value,
		int>::type = 0> 
    FunctionalBase(const PDE& pde, const DVector<double>& time, const SamplingData&... s)
      : select_regularization_type<Model>::type(pde, time),
      SamplingDesign<Model, model_traits<Model>::sampling>(s...) {};

    // copy constructor, copy only pde object (as a consequence also the problem domain)
    FunctionalBase(const FunctionalBase& rhs) { pde_ = rhs.pde_; }

  };
  
}}

#endif // __FUNCTIONAL_BASE_H__

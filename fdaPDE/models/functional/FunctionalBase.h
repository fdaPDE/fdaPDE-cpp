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
    : public select_regularization_type<Model>::type,
      public SamplingDesign<Model, typename model_traits<Model>::sampling> {
  protected:
    // vector of smoothing parameters
    std::vector<SVector<model_traits<Model>::n_lambda>> lambdas_;
  public:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    typedef typename select_regularization_type<Model>::type Base;
    using Base::pde_;  // differential operator L 
    
    FunctionalBase() = default;
    // space-only constructor
    template <typename U = Model, // fake type to enable substitution in SFINAE
	      typename std::enable_if<
		std::is_same<typename model_traits<U>::regularization, SpaceOnly>::value,
		int>::type = 0> 
    FunctionalBase(const PDE& pde) 
      : select_regularization_type<Model>::type(pde),
      SamplingDesign<Model, typename model_traits<Model>::sampling>() {};

    // space-time constructor
    template <typename U = Model, // fake type to enable substitution in SFINAE
	      typename std::enable_if<!
		std::is_same<typename model_traits<U>::regularization, SpaceOnly>::value,
		int>::type = 0> 
    FunctionalBase(const PDE& pde, const DVector<double>& time)
      : select_regularization_type<Model>::type(pde, time),
      SamplingDesign<Model, typename model_traits<Model>::sampling>() {};

    // copy constructor, copy only pde object (as a consequence also the problem domain)
    FunctionalBase(const FunctionalBase& rhs) { pde_ = rhs.pde_; }

    // accepts a collection of \lambda parameters if a not fixed_lambda method is selected
    void setLambda(const std::vector<SVector<model_traits<Model>::n_lambda>>& lambdas) { lambdas_ = lambdas; }
    const std::vector<SVector<model_traits<Model>::n_lambda>>& lambdas() const { return lambdas_; }
  };
  
}}

#endif // __FUNCTIONAL_BASE_H__

#ifndef __I_REGRESSION_MODEL_H__
#define __I_REGRESSION_MODEL_H__

#include "../../core/utils/Symbols.h"
#include "../../core/utils/Traits.h"
#include "../iStatModel.h"
using fdaPDE::models::select_regularization_type;
#include "../iSpaceOnlyModel.h"
#include "../../core/utils/DataStructures/BlockFrame.h"
#include <memory>

namespace fdaPDE {
namespace models {

#define DESIGN_MATRIX_BLK "X" // design matrix
#define WEIGHTS_BLK "W" // weights for heteroscedastic observations
  
  // base class for any regression model
  template <typename Model>
  struct iRegressionModel : public select_regularization_type<Model>::type {
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    typedef typename select_regularization_type<Model>::type Base;
    IMPORT_STAT_MODEL_SYMBOLS;
    
    iRegressionModel() = default;
    // space-only constructor
    template <typename U = Model, // fake type to enable substitution
	      typename std::enable_if<
		std::is_same<typename model_traits<U>::RegularizationType, SpaceOnly>::value,
		int>::type = 0> 
    iRegressionModel(const PDE& pde) : select_regularization_type<Model>::type(pde) {};
    // space-time constructor
    template <typename U = Model, // fake type to enable substitution
	      typename std::enable_if<!
		std::is_same<typename model_traits<U>::RegularizationType, SpaceOnly>::value,
		int>::type = 0> 
    iRegressionModel(const PDE& pde, const DVector<double>& time) : select_regularization_type<Model>::type(pde, time) {};
    
    // copy constructor, copy only pde object (as a consequence also the problem domain)
    iRegressionModel(const iRegressionModel& rhs) { pde_ = rhs.pde_; }

    // getters
    std::size_t q() const { return df_.hasBlock(DESIGN_MATRIX_BLK) ?
	df_.template get<double>(DESIGN_MATRIX_BLK).cols() : 0; }
    const DMatrix<double>& X() const { return df_.template get<double>(DESIGN_MATRIX_BLK); } // covariates
    auto W() const { return df_.template col<double>(WEIGHTS_BLK, 0); } // observations' weights
    
    // utilities
    bool hasCovariates() const { return q() != 0; } // true if the model has a parametric part
    bool hasWeights() const { return df_.hasBlock(WEIGHTS_BLK); }
    
    // abstract part of the interface, must be implemented by concrete models   
    // getters to problem's solution (estimate of spatial field, PDE misfit and parametric vector)
    virtual const DMatrix<double>& f() const = 0;
    virtual const DMatrix<double>& g() const = 0;
    virtual const DMatrix<double>& beta() const = 0;
    // efficient multiplication by matrix Q
    virtual DMatrix<double> lmbQ(const DMatrix<double>& x) = 0;
    virtual DMatrix<double> fitted() = 0; // computes fitted values at observations' locations
    // compute prediction at new unseen datapoint
    virtual double predict(const DVector<double>& covs, const std::size_t loc) const = 0;    
  };

  // this macro is intended to import all **common** symbols a model can expect from a Regression base
  // symbols specific for the regularization type used need to be imported via explicit using declaration
#define IMPORT_REGRESSION_SYMBOLS		\
  IMPORT_STAT_MODEL_SYMBOLS;			\
  using Base::q;				\
  using Base::X;				\
  using Base::W;				\
  using Base::hasCovariates;			\
  using Base::hasWeights;			\
  
  // trait to detect if a type implements iRegressionModel
  template <typename T>
  struct is_regression_model {
    static constexpr bool value = fdaPDE::is_base_of_template<iRegressionModel, T>::value;
  };

}}

#endif // __I_REGRESSION_MODEL_H__

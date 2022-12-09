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
    IMPORT_STAT_MODEL_SYMBOLS(Model);
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    
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

#define IMPORT_REGRESSION_MODEL_SYMBOLS( ... )				\
  using iRegressionModel<__VA_ARGS__>::q;				\
  using iRegressionModel<__VA_ARGS__>::X;				\
  using iRegressionModel<__VA_ARGS__>::W;				\
  using iRegressionModel<__VA_ARGS__>::hasCovariates;			\
  using iRegressionModel<__VA_ARGS__>::hasWeights;			\
  
#define IMPORT_SPACE_ONLY_REGRESSION_MODEL_SYMBOLS( ... )	\
  IMPORT_SPACE_ONLY_MODEL_SYMBOLS(__VA_ARGS__)			\
  IMPORT_REGRESSION_MODEL_SYMBOLS(__VA_ARGS__)			\

#define IMPORT_SPACE_TIME_REGRESSION_MODEL_SYMBOLS( ... )	\
  IMPORT_SPACE_TIME_MODEL_SYMBOLS(__VA_ARGS__)			\
  IMPORT_REGRESSION_MODEL_SYMBOLS(__VA_ARGS__)			\
  
  // trait to detect if a type implements iRegressionModel
  template <typename T>
  struct is_regression_model {
    static constexpr bool value = fdaPDE::is_base_of_template<iRegressionModel, T>::value;
  };

}}

#endif // __I_REGRESSION_MODEL_H__

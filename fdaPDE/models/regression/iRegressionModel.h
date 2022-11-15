#ifndef __I_REGRESSION_MODEL_H__
#define __I_REGRESSION_MODEL_H__

#include "../../core/utils/Symbols.h"
#include "../../core/utils/Traits.h"
#include "../iStatModel.h"
using fdaPDE::models::iStatModel;
#include "../../core/utils/DataStructures/BlockFrame.h"
#include <memory>

namespace fdaPDE {
namespace models {

#define STAT_MODEL_X_BLK "X" // design matrix
  
  // base class for any regression model
  template <typename PDE>
  struct iRegressionModel : public iStatModel<PDE> {
    IMPORT_STAT_MODEL_SYMBOLS(PDE);

    iRegressionModel() = default;
    iRegressionModel(const PDE& pde, double lambda)
      : iStatModel<PDE>(pde, lambda) {};
    // copy constructor, copy only pde object (as a consequence also the problem domain)
    iRegressionModel(const iRegressionModel& rhs) { pde_ = rhs.pde_; }

    // getters
    std::size_t q() const { return df_.hasBlock(STAT_MODEL_X_BLK) ?
	df_.template get<double>(STAT_MODEL_X_BLK).cols() : 0; }
    const DMatrix<double>& X() const { return df_.template get<double>(STAT_MODEL_X_BLK); } // covariates

    // utilities
    bool hasCovariates() const { return q() != 0; } // true if the model has a parametric part
    
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

#define IMPORT_REGRESSION_MODEL_SYMBOLS(E)	   \
  IMPORT_STAT_MODEL_SYMBOLS(E)			   \
  using iRegressionModel<E>::q;			   \
  using iRegressionModel<E>::X;			   \
  using iRegressionModel<E>::hasCovariates;	   \

  // trait to detect if a type implements iRegressionModel
  template <typename T>
  struct is_regression_model {
    static constexpr bool value = fdaPDE::is_base_of_template<iRegressionModel, T>::value;
  };

}}

#endif // __I_REGRESSION_MODEL_H__

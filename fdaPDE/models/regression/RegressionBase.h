#ifndef __REGRESSION_BASE_H__
#define __REGRESSION_BASE_H__

#include "../../core/utils/Symbols.h"
#include "../ModelBase.h"
#include "../SamplingDesign.h"
#include "../ModelTraits.h"
using fdaPDE::models::select_regularization_type;
#include "../SpaceOnlyBase.h"
#include "../SpaceTimeBase.h"
#include "../space_time/SpaceTimeSeparableBase.h"

namespace fdaPDE {
namespace models {

#define DESIGN_MATRIX_BLK "X" // design matrix
#define WEIGHTS_BLK "W" // weights for heteroscedastic observations
  
  // base class for any *regression* fdaPDE model
  template <typename Model>
  class RegressionBase : public select_regularization_type<Model>::type, public SamplingDesign<Model, model_traits<Model>::sampling> {
  protected:
    DiagMatrix<double> W_; // diagonal matrix of weights (implements possible heteroscedasticity)
    DMatrix<double> XtWX_{}; // q x q dense matrix X^T*W*X
    Eigen::PartialPivLU<DMatrix<double>> invXtWX_{}; // factorization of the dense q x q matrix XtWX_.
  public:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    typedef typename select_regularization_type<Model>::type Base;
    using Base::pde_;
    using Base::df_;
    
    RegressionBase() = default;
    // space-only constructor
    template <typename... SamplingData,
	      typename U = Model, // fake type to enable substitution in SFINAE
	      typename std::enable_if<
		std::is_same<typename model_traits<U>::RegularizationType, SpaceOnlyTag>::value,
		int>::type = 0> 
    RegressionBase(const PDE& pde, const SamplingData&... s) 
      : select_regularization_type<Model>::type(pde),
      SamplingDesign<Model, model_traits<Model>::sampling>(s...) {}; // s can either be DMatrix<double> or DMatrix<int>

    // space-time constructor
    template <typename... SamplingData,
	      typename U = Model, // fake type to enable substitution in SFINAE
	      typename std::enable_if<!
		std::is_same<typename model_traits<U>::RegularizationType, SpaceOnlyTag>::value,
		int>::type = 0> 
    RegressionBase(const PDE& pde, const DVector<double>& time, const SamplingData&... s)
      : select_regularization_type<Model>::type(pde, time),
      SamplingDesign<Model, model_traits<Model>::sampling>(s...) {}; // s can either be DMatrix<double> or DMatrix<int>
    
    // copy constructor, copy only pde object (as a consequence also the problem domain)
    RegressionBase(const RegressionBase& rhs) { pde_ = rhs.pde_; }

    // getters
    std::size_t q() const { return df_.hasBlock(DESIGN_MATRIX_BLK) ?
	df_.template get<double>(DESIGN_MATRIX_BLK).cols() : 0; }
    const DMatrix<double>& X() const { return df_.template get<double>(DESIGN_MATRIX_BLK); } // covariates
    const DiagMatrix<double>& W() const { return W_; } // observations' weights
    const DMatrix<double>& XtWX() const { return XtWX_; } 
    const Eigen::PartialPivLU<DMatrix<double>>& invXtWX() const { return invXtWX_; }
    
    // utilities
    bool hasCovariates() const { return q() != 0; } // true if the model has a parametric part
    bool hasWeights() const { return df_.hasBlock(WEIGHTS_BLK); } // true if heteroscedastic observation are assumed

    // an efficient way to perform a left multiplication by Q implementing the following
    //  given the design matrix X, the weight matrix W and x
    //    compute v = X^T*W*x
    //    solve Yz = v
    //    return Wx - WXz = W(I-H)x = Qx
    // it is required to having assigned a design matrix X to the model before calling this method
    DMatrix<double> lmbQ(const DMatrix<double>& x) {
      DMatrix<double> v = X().transpose()*W_*x; // X^T*W*x
      DMatrix<double> z = invXtWX_.solve(v);  // (X^T*W*X)^{-1}*X^T*W*x
      // compute W*x - W*X*z = W*x - (W*X*(X^T*W*X)^{-1}*X^T*W)*x = W(I - H)*x = Q*x
      return W_*x - W_*X()*z;
    }

    // perform proper preprocessing of input data and initialization of regression base.
    // This is called in ModelBase::setData() and executed after initialization of the block frame
    void preprocess() {
      // if heteroscedastic observations are provided as datum, prepare weights matrix
      if(hasWeights()) W_ = df_.template get<double>(WEIGHTS_BLK).asDiagonal();
      else // othersise set W_ to identity and assume homoscedastic observations
	W_ = DVector<double>::Ones(Base::n_obs()).asDiagonal();
      
      if(hasCovariates()){ // parametric model
	// compute q x q dense matrix X^T*W*X and its factorization
	XtWX_ = X().transpose()*W_*X();
	invXtWX_ = XtWX_.partialPivLu();
      }
    }
    
    // abstract part of the interface, must be implemented by concrete models   
    virtual DMatrix<double> fitted() = 0; // computes fitted values at observations' locations
    // compute prediction at new unseen datapoint
    virtual double predict(const DVector<double>& covs, const std::size_t loc) const = 0;    
  };
  
  // trait to detect if a type is a regression model
  template <typename T>
  struct is_regression_model {
    static constexpr bool value = fdaPDE::is_base_of_template<RegressionBase, T>::value;
  };

}}

#endif // __REGRESSION_BASE_H__

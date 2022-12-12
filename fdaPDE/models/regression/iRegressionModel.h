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
  class iRegressionModel : public select_regularization_type<Model>::type {
  protected:
    // diagonal matrix of weights (implements possible heteroscedasticity)
    DiagMatrix<double> W_;
    // q x q dense matrix X^T*W*X
    DMatrix<double> XTX_{};
    // partial LU (with pivoting) factorization of the dense (square invertible) q x q matrix XTX_.
    Eigen::PartialPivLU<DMatrix<double>> invXTX_{};

    // perform proper preprocessing of input data and initialization of regression base
    void analyzeData() {
      // if heteroscedastic observations are provided as datum, prepare weights matrix
      if(hasWeights()) W_ = W().asDiagonal();
      else W_ = DVector<double>::Ones(obs()).asDiagonal(); // homoscedastic observations
      
      if(hasCovariates()){ // parametric model
	// compute q x q dense matrix X^T*W*X and its factorization
	XTX_ = X().transpose()*W_*X();
	invXTX_ = XTX_.partialPivLu();
      }
    }
    
  public:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    typedef typename select_regularization_type<Model>::type Base;
    IMPORT_STAT_MODEL_SYMBOLS;
    
    iRegressionModel() = default;
    // space-only constructor
    template <typename U = Model, // fake type to enable substitution in SFINAE
	      typename std::enable_if<
		std::is_same<typename model_traits<U>::RegularizationType, SpaceOnly>::value,
		int>::type = 0> 
    iRegressionModel(const PDE& pde) : select_regularization_type<Model>::type(pde) {};
    // space-time constructor
    template <typename U = Model, // fake type to enable substitution in SFINAE
	      typename std::enable_if<!
		std::is_same<typename model_traits<U>::RegularizationType, SpaceOnly>::value,
		int>::type = 0> 
    iRegressionModel(const PDE& pde, const DVector<double>& time)
      : select_regularization_type<Model>::type(pde, time) {};
    
    // copy constructor, copy only pde object (as a consequence also the problem domain)
    iRegressionModel(const iRegressionModel& rhs) { pde_ = rhs.pde_; }

    // getters
    std::size_t q() const { return df_.hasBlock(DESIGN_MATRIX_BLK) ?
	df_.template get<double>(DESIGN_MATRIX_BLK).cols() : 0; }
    const DMatrix<double>& X() const { return df_.template get<double>(DESIGN_MATRIX_BLK); } // covariates
    auto W() const { return df_.template col<double>(WEIGHTS_BLK, 0); } // observations' weights
   
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
      DMatrix<double> z = invXTX_.solve(v);  // (X^T*W*X)^{-1}*X^T*W*x
      // compute W*x - W*X*z = W*x - (W*X*(X^T*W*X)^{-1}*X^T*W)*x = W(I - H)*x = Q*x
      return W_*x - W_*X()*z;
    }
    
    // abstract part of the interface, must be implemented by concrete models   
    virtual DMatrix<double> fitted() = 0; // computes fitted values at observations' locations
    // compute prediction at new unseen datapoint
    virtual double predict(const DVector<double>& covs, const std::size_t loc) const = 0;    
  };

  // this macro is intended to import all **common** symbols a model can expect from a Regression base
  // symbols specific for the regularization type used need to be imported via explicit using declaration
#define IMPORT_REGRESSION_SYMBOLS		\
  IMPORT_STAT_MODEL_SYMBOLS;			\
  /* private members */				\
  using Base::W_;				\
  using Base::XTX_;				\
  using Base::invXTX_;				\
  /* getters */					\
  using Base::q;				\
  using Base::X;				\
  using Base::W;				\
  /* utilities */				\
  using Base::hasCovariates;			\
  using Base::hasWeights;			\
  using Base::lmbQ;				\
  
  // trait to detect if a type is a regression model
  template <typename T>
  struct is_regression_model {
    static constexpr bool value = fdaPDE::is_base_of_template<iRegressionModel, T>::value;
  };

}}

#endif // __I_REGRESSION_MODEL_H__

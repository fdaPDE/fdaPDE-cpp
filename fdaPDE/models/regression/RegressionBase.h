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
#include "../space_time/SpaceTimeParabolicBase.h"

namespace fdaPDE {
namespace models {

#define DESIGN_MATRIX_BLK "X" // design matrix
#define WEIGHTS_BLK "W" // weights for heteroscedastic observations
  
  // base class for any *regression* fdaPDE model
  template <typename Model>
  class RegressionBase :
    public select_regularization_type<Model>::type, public SamplingDesign<Model, model_traits<Model>::sampling> {
  protected:
    DiagMatrix<double> W_{}; // diagonal matrix of weights (implements possible heteroscedasticity)
    DMatrix<double> XtWX_{}; // q x q dense matrix X^T*W*X
    Eigen::PartialPivLU<DMatrix<double>> invXtWX_{}; // factorization of the dense q x q matrix XtWX_.

    // room for problem solution
    DMatrix<double> f_{};    // estimate of the spatial field (1 x N vector)
    DMatrix<double> g_{};    // PDE misfit
    DMatrix<double> beta_{}; // estimate of the coefficient vector (1 x q vector)
  public:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    typedef typename select_regularization_type<Model>::type Base;
    using Base::pde_; // differential operator L 
    using Base::df_;  // BlockFrame for problem's data storage
    using Base::Psi;  // matrix of spatial basis evaluation at locations p_1 ... p_n
    using Base::n_basis; // number of basis function over domain D
    
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
    DMatrix<double> lmbQ(const DMatrix<double>& x) const {
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

    // getters to problem solution
    const DMatrix<double>& f() const    { return f_; };
    const DMatrix<double>& g() const    { return g_; };
    const DMatrix<double>& beta() const { return beta_; };
    
    // computes fitted values \hat y = \Psi*f_ + X*beta_
    DMatrix<double> fitted() const {
      DMatrix<double> hat_y = Psi()*f_;
      if(hasCovariates()) hat_y += X()*beta_;
      return hat_y;
    }
    
    // compute prediction of model at new unseen data location: \hat y(p_{n+1}) = X_{n+1}^T*\beta + f_*\psi(p_{n+1})
    // TODO: divide from space and space-time cases
    double predict(const DVector<double>& p, std::size_t t, const DVector<double>& covs) {
      // vector \psi(p_{n+1}) = (\psi_1(p_{n+1}), \ldots, \psi_N(p_{n+1}))
      DVector<double> psi_new = DVector<double>::Zero(n_basis());

      auto gse = this->gse(); // require pointer to geometric search engine
      auto e = gse.search(p); // search element containing point
      // fill entries in psi_new vector
      for(std::size_t j = 0; j < pde_->basis()[e->ID()].size(); ++j){
	// extract \phi_h from basis
	std::size_t h = e->nodeIDs()[j]; // row index of psi_new vector
	auto psi_h = pde_->basis()[e->ID()][j];
	// store value of basis function \psi_h evaluated at query point p
	psi_new[h] = psi_h(p);
      }
      
      double result = DVector<double>(f_.block(n_basis()*t,0, n_basis(),1)).dot(psi_new);
      if(hasCovariates()) result += covs.dot(DVector<double>(beta_));
      return result;
    }
  };
  
  // trait to detect if a type is a regression model
  template <typename T>
  struct is_regression_model {
    static constexpr bool value = fdaPDE::is_base_of_template<RegressionBase, T>::value;
  };

}}

#endif // __REGRESSION_BASE_H__

#ifndef __REGRESSION_BASE_H__
#define __REGRESSION_BASE_H__

#include "../../core/utils/Symbols.h"
#include "../ModelBase.h"
#include "../SamplingDesign.h"
#include "../ModelTraits.h"
#include "../ModelMacros.h"
#include "../SpaceOnlyBase.h"
#include "../SpaceTimeBase.h"
#include "../space_time/SpaceTimeSeparableBase.h"
#include "../space_time/SpaceTimeParabolicBase.h"

namespace fdaPDE {
namespace models {
  
  // base class for any *regression* fdaPDE model
  template <typename Model>
  class RegressionBase :
    public select_regularization_type<Model>::type,
    public SamplingDesign<Model, typename model_traits<Model>::sampling> {
  protected:
    DiagMatrix<double> W_{}; // diagonal matrix of weights (implements possible heteroscedasticity)
    DMatrix<double> XtWX_{}; // q x q dense matrix X^T*W*X
    Eigen::PartialPivLU<DMatrix<double>> invXtWX_{}; // factorization of the dense q x q matrix XtWX_.

    // matrices required for Woodbury decomposition
    DMatrix<double> U_;  // [\Psi^T*D*W*X, 0] 
    DMatrix<double> V_;  // [X^T*W*\Psi,   0]

    // quantites related to missing data setting
    std::unordered_set<std::size_t> nan_idxs_; // indexes of missing observations
    SpMatrix<double> B_; // matrix \Psi where rows corresponding to NaN observations are zeroed
    
    // room for problem solution
    DVector<double> f_{};    // estimate of the spatial field (1 x N vector)
    DVector<double> g_{};    // PDE misfit
    DVector<double> beta_{}; // estimate of the coefficient vector (1 x q vector)
  public:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    typedef typename select_regularization_type<Model>::type Base;
    typedef SamplingDesign<Model, typename model_traits<Model>::sampling> SamplingBase;
    using Base::pde_;        // differential operator L 
    using Base::df_;         // BlockFrame for problem's data storage
    using Base::n_basis;     // number of basis function over domain D
    using Base::idx;         // indices of observations
    using SamplingBase::Psi; // matrix of spatial basis evaluation at locations p_1 ... p_n
    using SamplingBase::D;   // matrix of subdomains measures (for areal sampling)
    
    RegressionBase() = default;
    // space-only constructor
    template <typename U = Model, // fake type to enable substitution in SFINAE
	      typename std::enable_if<
		std::is_same<typename model_traits<U>::regularization, SpaceOnly>::value,
		int>::type = 0> 
    RegressionBase(const PDE& pde) 
      : select_regularization_type<Model>::type(pde),
      SamplingDesign<Model, typename model_traits<Model>::sampling>() {};
    // space-time constructor
    template <typename U = Model, // fake type to enable substitution in SFINAE
	      typename std::enable_if<!
		std::is_same<typename model_traits<U>::regularization, SpaceOnly>::value,
		int>::type = 0> 
    RegressionBase(const PDE& pde, const DVector<double>& time)
      : select_regularization_type<Model>::type(pde, time),
      SamplingDesign<Model, typename model_traits<Model>::sampling>() {};
    // copy constructor, copy only pde object (as a consequence also the problem domain)
    RegressionBase(const RegressionBase& rhs) { pde_ = rhs.pde_; }

    // getters
    const DMatrix<double>& y() const { return df_.template get<double>(OBSERVATIONS_BLK); } // observation vector y
    std::size_t q() const { return df_.hasBlock(DESIGN_MATRIX_BLK) ? df_.template get<double>(DESIGN_MATRIX_BLK).cols() : 0; }
    const DMatrix<double>& X() const { return df_.template get<double>(DESIGN_MATRIX_BLK); } // covariates
    const DiagMatrix<double>& W() const { return W_; } // observations' weights
    const DMatrix<double>& XtWX() const { return XtWX_; } 
    const Eigen::PartialPivLU<DMatrix<double>>& invXtWX() const { return invXtWX_; }
    const DVector<double>& f() const { return f_; }; // estimate of spatial field
    const DVector<double>& g() const { return g_; }; // PDE misfit
    const DVector<double>& beta() const { return beta_; }; // estimate of regression coefficients
    const std::unordered_set<std::size_t>& nan_idxs() const { return nan_idxs_; } // missing data indexes
    std::size_t n_obs() const { return y().rows(); } // number of observations
    // getters to Woodbury decomposition matrices
    const DMatrix<double>& U() const { return U_; }
    const DMatrix<double>& V() const { return V_; }
    // access to NaN corrected \Psi and \Psi^T*D matrices
    const SpMatrix<double>& Psi() const { return has_nan() ? B_ : Psi(not_nan()); }
    auto PsiTD() const { return has_nan() ? B_.transpose()*D() : Psi(not_nan()).transpose()*D();}

    // utilities
    bool hasCovariates() const { return q() != 0; } // true if the model has a parametric part
    bool hasWeights() const { return df_.hasBlock(WEIGHTS_BLK); } // true if heteroscedastic observation are assumed
    bool has_nan() const { return nan_idxs_.size() != 0; } // true if there are missing data
    DMatrix<double> lmbQ(const DMatrix<double>& x) const; // efficient multiplication by matrix Q
    DMatrix<double> fitted() const; // computes fitted values \hat y = \Psi*f_ + X*beta_

    // initialization methods 
    void init_data(); // update model's status to data (called by ModelBase::setData())
    void init_nan();  // regression models' missing data logic (called by SamplingBase::init_sampling())
  };

  # include "RegressionBase.tpp"
  
  // trait to detect if a type is a regression model
  template <typename T>
  struct is_regression_model {
    static constexpr bool value = fdaPDE::is_base_of_template<RegressionBase, T>::value;
  };

}}

#endif // __REGRESSION_BASE_H__

#ifndef __FUNCTIONAL_BASE_H__
#define __FUNCTIONAL_BASE_H__

#include <algorithm>
#include "../../core/utils/Symbols.h"
#include "../ModelBase.h"
#include "../SamplingDesign.h"
#include "../ModelTraits.h"
using fdaPDE::models::select_regularization_type;

namespace fdaPDE {
namespace models {

  // ** to be moved **
  struct fixed_lambda {};
  struct gcv_lambda_selection {};
  struct kcv_lambda_selection {};
  
  // base class for any *functional* fdaPDE model
  template <typename Model>
  class FunctionalBase
    : public select_regularization_type<Model>::type,
      public SamplingDesign<Model, typename model_traits<Model>::sampling> {
  protected:
    // vector of smoothing parameters
    std::vector<SVector<model_traits<Model>::n_lambda>> lambdas_;

    // quantites related to missing data setting
    std::vector<std::unordered_set<std::size_t>> nan_idxs_; // for each statistical unit, the indexes of missing observations
    std::vector<SpMatrix<double>> B_; // matrix \Psi of the i-th statistical unit where rows corresponding to NaN observations are zeroed
  public:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    typedef typename select_regularization_type<Model>::type Base;
    typedef SamplingDesign<Model, typename model_traits<Model>::sampling> SamplingBase;
    using Base::pde_;          // differential operator L
    using Base::df_;           // BlockFrame for problem's data storage
    using Base::n_locs;        // number of spatial (spatio-temporal) data locations
    using SamplingBase::Psi;   // matrix of basis evaluations at locations p_1 ... p_n
    using SamplingBase::PsiTD; // block \Psi^T*D
    using SamplingBase::D;     // matrix of subdomains measures (for areal sampling)
    
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

    // getters
    const DMatrix<double>& X() const { return df_.template get<double>(OBSERVATIONS_BLK); } // observation matrix y
    std::size_t n_stat_units() const { return X().rows(); }
    const std::vector<std::unordered_set<std::size_t>>& nan_idxs() const { return nan_idxs_; } // missing data indexes
    std::size_t n_nan() const; // total number of missing data points
    std::size_t n_obs() const { return X().size() - n_nan(); }; // number of not missing observations
    // access to NaN corrected \Psi and \Psi^T*D matrices of i-th unit
    const SpMatrix<double>& Psi(std::size_t i) const { return has_nan(i) ? B_[i] : Psi(not_nan()); }
    auto PsiTD(std::size_t i) const { return has_nan(i) ? B_[i].transpose()*D() : Psi(not_nan()).transpose()*D();}

    // setters
    // accepts a collection of \lambda parameters if a not fixed_lambda method is selected ** not stable **
    void setLambda(const std::vector<SVector<model_traits<Model>::n_lambda>>& lambdas) { lambdas_ = lambdas; }
    const std::vector<SVector<model_traits<Model>::n_lambda>>& lambdas() const { return lambdas_; } // ** not stable **
    
    // utilities
    bool has_nan(std::size_t i) const { return !nan_idxs_[i].empty(); } // true if the i-th unit has missing data
    bool has_nan() const; // true if any of the statistical unit has missing data
    
    // initialization methods
    void update_data() { return; }   
    void init_nan(); // functional models' missing data logic (called by SamplingBase::init_sampling())
  };

  template <typename Model>
  std::size_t FunctionalBase<Model>::n_nan() const {
    return std::accumulate(nan_idxs_.begin(), nan_idxs_.end(), 0,
			   [] (std::size_t n, const std::unordered_set<std::size_t>& stat_unit_nan) {
			     return n + stat_unit_nan.size();
			   });
  }
  
  #include "FunctionalBase.tpp"
  
}}

#endif // __FUNCTIONAL_BASE_H__

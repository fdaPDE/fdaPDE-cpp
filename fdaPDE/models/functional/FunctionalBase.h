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
    using SamplingBase::Psi;   // matrix of spatial basis evaluation at locations p_1 ... p_n
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

    // utilities
    bool has_nan(std::size_t i) const { return !nan_idxs_[i].empty(); } // test for i-th statistical unit missing data
    bool has_nan() const { // true if there are missing data in any of the statistical unit
      for(auto s : nan_idxs_) { if(!s.empty()) return true; }
      return false;
    }

    void init_data() { return; }
    
    // functional models' missing data logic
    void init_nan() {
      nan_idxs_.clear(); // clean previous missingness structure
      nan_idxs_.resize(n_stat_units()); B_.resize(n_stat_units());
      // \Psi matrix dimensions
      std::size_t n = Psi(not_nan()).rows();
      std::size_t N = Psi(not_nan()).cols();
      // for i-th statistical unit, analyze missingness structure and set \Psi_i
      for(std::size_t i = 0; i < n_stat_units(); ++i){
	// derive missingness pattern for i-th statistical unit
	for(std::size_t j = 0; j < n_locs(); ++j){
	  if(std::isnan(X()(i,j))) // requires -ffast-math compiler flag to be disabled
	    nan_idxs_[i].insert(j);
	}

	// NaN detected for this unit, start assembly
	if(!nan_idxs_[i].empty()){
	  for(std::size_t i = 0; i < n_stat_units(); ++i){
	    B_[i].resize(n, N); // reserve space
	    std::vector<fdaPDE::Triplet<double>> tripletList;
	    tripletList.reserve(n*N);
	    for(int k = 0; k < Psi(not_nan()).outerSize(); ++k)
	      for(SpMatrix<double>::InnerIterator it(Psi(not_nan()),k); it; ++it){
		if(nan_idxs_[i].find(it.row()) == nan_idxs_[i].end())
		  // no missing data at this location for i-th statistical unit
		  tripletList.emplace_back(it.row(), it.col(), it.value());
	      }
	    // finalize construction
	    B_[i].setFromTriplets(tripletList.begin(), tripletList.end());
	    B_[i].makeCompressed();
	  }	  
	}
	// otherwise no matrix is assembled, full \Psi is returned by Psi(std::size_t) getter
      }
      return;
    }
    // access to NaN corrected \Psi and \Psi^T*D matrices of i-th subject
    const SpMatrix<double>& Psi(std::size_t i) const { return has_nan(i) ? B_[i] : Psi(not_nan()); }
    auto PsiTD(std::size_t i) const { return has_nan(i) ? B_[i].transpose()*D() : Psi(not_nan()).transpose()*D();}
    
    // accepts a collection of \lambda parameters if a not fixed_lambda method is selected
    void setLambda(const std::vector<SVector<model_traits<Model>::n_lambda>>& lambdas) { lambdas_ = lambdas; }
    const std::vector<SVector<model_traits<Model>::n_lambda>>& lambdas() const { return lambdas_; }
  };
  
}}

#endif // __FUNCTIONAL_BASE_H__

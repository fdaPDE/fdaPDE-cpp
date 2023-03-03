#ifndef __SPACE_TIME_BASE_H__
#define __SPACE_TIME_BASE_H__

#include <memory>
#include <type_traits>

#include "../core/utils/Symbols.h"
#include "../core/utils/Traits.h"
#include "ModelBase.h"
using fdaPDE::models::ModelBase;

namespace fdaPDE {
namespace models {
  
  // abstract base interface for any *space-time* fdaPDE statistical model. This class is not directly usable from
  // models since the type of time regularization is not yet defined at this point
  template <typename Model>
  class SpaceTimeBase : public ModelBase<Model> {
    // check Model refers to the space-time case
    static_assert(!std::is_same<typename model_traits<Model>::RegularizationType, SpaceOnlyTag>::value);
  protected:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    typedef ModelBase<Model> Base;
    using Base::pde_; // regularizing PDE
    using Base::lambda_; // vector of smoothing parameters

    DVector<double> time_; // time domain [0, T]
  public:
    // constructor
    SpaceTimeBase() = default;
    SpaceTimeBase(const PDE& pde, const DVector<double>& time) : ModelBase<Model>(pde), time_(time) {};
    // copy constructor
    SpaceTimeBase(const SpaceTimeBase& rhs) { pde_ = rhs.pde_; time_ = rhs.time_; }

    // setters
    void setLambdaS(double lambdaS) { lambda_[0] = lambdaS; }
    void setLambdaT(double lambdaT) { lambda_[1] = lambdaT; }    
    void setTimeDomain(const DVector<double>& time) { time_ = time; }
    // getters
    inline double lambdaS() const { return lambda_[0]; }
    inline double lambdaT() const { return lambda_[1]; }
    inline const DVector<double>& time_domain() const { return time_; }
    inline std::size_t n_time() const { return time_.rows(); } // number of time instants

    // remove first n time instants from the problem
    void shift_time(std::size_t n) {
      std::size_t m = time_.rows(); // number of time instants
      time_ = time_.bottomRows(m - n); // correct time interval [0,T]
      // correct regularization PDE informations
      pde_->setForcing(pde_->forcingData().rightCols(m - n));
    }

    // impose missing data by setting the (j*n_basis + i)-th row of \Psi to zero if no data
    // is observed at space-time point (p_i, t_j)
    // void setNaN() {
    //   // compute locations indexes where data are not observed
    //   std::unordered_set<std::size_t> missing_idx;
    //   for(std::size_t i = 0; i < model().n_locs()*n_time(); ++i) missing_idx.emplace(i);
    //   for(std::size_t i = 0; i < model().idx().rows(); ++i) missing_idx.erase(model().idx()(i,0));
    //   // impose NaN at space-time location (p_i, t_j) by setting the (j*n_basis + i)-th row of \Psi to zero
    //   for(int k = 0; k < model().Psi_.outerSize(); ++k){
    // 	for(SpMatrix<double>::InnerIterator it(model().Psi_,k); it; ++it){
    // 	  if(missing_idx.find(it.row()) != missing_idx.end())
    // 	    it.valueRef() = 0;
    // 	}
    //   }
    //   return;
    // }
    
    // destructor
    virtual ~SpaceTimeBase() = default;  
  };

}}

#endif // __SPACE_TIME_BASE_H__

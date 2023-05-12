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
    static_assert(is_space_time<Model>::value);
  protected:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    typedef ModelBase<Model> Base;
    using Base::pde_;    // regularizing PDE
    using Base::lambda_; // vector of smoothing parameters
    using Base::df_;     // model's data
    using Base::model;   // underlying model object
    
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
    const DVector<double>& time_domain() const { return time_; }
    const DVector<double>& time_locs() const { return time_; } // for space-time separable regularization it might be different
    inline std::size_t n_temporal_locs() const { return time_.rows(); } // number of time instants

    // remove first n time instants from the problem
    void shift_time(std::size_t n) {
      std::size_t m = time_.rows(); // number of time instants
      time_ = time_.bottomRows(m - n); // correct time interval [0,T]
      // correct regularization PDE informations
      pde_->setForcing(pde_->forcingData().rightCols(m - n));
      // remove from data first n time instants, reindex points
      model().setData(df_.tail(n * model().n_spatial_locs()).extract(), true);
    }
    
    // destructor
    virtual ~SpaceTimeBase() = default;  
  };

}}

#endif // __SPACE_TIME_BASE_H__

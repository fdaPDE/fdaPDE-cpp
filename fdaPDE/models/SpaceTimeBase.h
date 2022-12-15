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
    using Base::pde_;

    double lambdaS_; // smoothing parameter in space
    double lambdaT_; // smoothing parameter in time
    DVector<double> time_; // time domain [0, T] 
  public:
    // constructor
    SpaceTimeBase() = default;
    SpaceTimeBase(const PDE& pde, const DVector<double>& time) : ModelBase<Model>(pde), time_(time) {};
    // copy constructor
    SpaceTimeBase(const SpaceTimeBase& rhs) { pde_ = rhs.pde_; time_ = rhs.time_; }

    // setters
    void setLambdaS(double lambdaS) { lambdaS_ = lambdaS; }
    void setLambdaT(double lambdaT) { lambdaT_ = lambdaT; }    
    // getters
    inline double lambdaS() const { return lambdaS_; }
    inline double lambdaT() const { return lambdaT_; }
    inline const DVector<double>& time_domain() const { return time_; }
    
    // destructor
    virtual ~SpaceTimeBase() = default;  
  };

}}

#endif // __SPACE_TIME_BASE_H__

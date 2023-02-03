#ifndef __SPACE_ONLY_BASE_H__
#define __SPACE_ONLY_BASE_H__

#include <memory>
#include <type_traits>

#include "../core/utils/Symbols.h"
#include "ModelBase.h"
using fdaPDE::models::ModelBase;

namespace fdaPDE {
namespace models {

  // abstract base interface for any *space-only* fdaPDE statistical model.
  template <typename Model>
  class SpaceOnlyBase : public ModelBase<Model> {
  protected:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    typedef ModelBase<Model> Base;
    using Base::pde_;  // regularizing PDE
    using Base::model; // underlying model object
    using Base::lambda_; // vector of smoothing parameters
  public:  
    // constructor
    SpaceOnlyBase() = default;
    SpaceOnlyBase(const PDE& pde) : ModelBase<Model>(pde) {};
    void init_regularization() { return; } // do nothing
    
    // setters
    void setLambdaS(double lambda) { lambda_[0] = lambda; } 
    // getters
    double lambdaS() const { return lambda_[0]; } // smoothing parameter
    const SpMatrix<double>& R0()  const { return pde_->R0(); }    // mass matrix in space
    const SpMatrix<double>& R1()  const { return pde_->R1(); }    // discretization of differential operator L
    const DMatrix<double>&  u()   const { return pde_->force(); } // discretization of forcing term u
    const SpMatrix<double>& Psi() const { return model().Psi_; }  // matrix of spatial basis evaluation
    
    // destructor
    virtual ~SpaceOnlyBase() = default;  
  };
  
}}

#endif // __SPACE_ONLY_BASE_H__

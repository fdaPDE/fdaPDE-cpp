#ifndef __I_GCV_H__
#define __I_GCV_H__

#include <memory>
#include <Eigen/SparseLU>
#include "../core/utils/Symbols.h"

namespace fdaPDE{
namespace calibration{

  // abstract base class for models capable to support selection of smoothing parameters via GCV optimization
  class iGCV {
  protected:
    fdaPDE::SparseLU<SpMatrix<double>> invR0_{};
    DMatrix<double> R_{}; // R = R1^T*R0^{-1}*R1
    DMatrix<double> T_{}; // T = \Psi^T*Q*\Psi + \lambda*R
    DMatrix<double> Q_{}; // Q_ = I - H, whatever H is for the model
  public:
    // constructor
    iGCV() {};

    // the following methods should compute matrices once and cache them for reuse using the provided data members.
    // Subsequent calls should immediately return the cached result.
    virtual const DMatrix<double>& Q() = 0; // returns matrix Q = I - H
    virtual const DMatrix<double>& T() = 0; // returns matrix T = \Psi^T*Q*\Psi + \lambda*R, stores also matrix R_
    // computes the norm of y - \hat y, according to the choosen distribution
    virtual double norm(const DMatrix<double>& obs, const DMatrix<double>& fitted) const = 0;

    // utilities
    fdaPDE::SparseLU<SpMatrix<double>>& invR0() { return invR0_; };
    const DMatrix<double>& R() { return R_; } 
    virtual ~iGCV() = default;
  };
}}

#endif // __I_GCV_H__

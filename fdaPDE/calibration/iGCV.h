#ifndef __I_GCV_H__
#define __I_GCV_H__

#include <memory>
#include <Eigen/SparseLU>

#include "../core/utils/Symbols.h"

// abstract base class for models capable to support selection of smoothing parameters via GCV optimization
class iGCV {
protected:
  // SparseLU has a deleted copy construcor, need to wrap it in a movable object to allow copy construction of derived models
  std::shared_ptr<Eigen::SparseLU<SpMatrix<double>>> invR0_{};
  std::shared_ptr<DMatrix<double>> R_{}; // R = R1^T*R0^{-1}*R1
  std::shared_ptr<DMatrix<double>> T_{}; // T = \Psi^T*Q*\Psi + \lambda*K
  
public:
  // constructor
  iGCV() {
    // initialize pointer to SparseLU solver
    invR0_ = std::make_shared<Eigen::SparseLU<SpMatrix<double>>>();
  };
  // performs computation of matrix T
  virtual std::shared_ptr<DMatrix<double>> T() = 0;
  // getters
  std::shared_ptr<DMatrix<double>> R() const { return R_; }
  std::shared_ptr<DMatrix<double>> T() const { return T_; }
  Eigen::SparseLU<SpMatrix<double>>& invR0() { return *invR0_; }
  
  virtual ~iGCV() = default;
};

#endif // __I_GCV_H__

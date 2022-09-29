#ifndef __I_GCV_H__
#define __I_GCV_H__

// abstract base class for models capable to support selection of smoothing parameters via GCV optimization
class iGCV {
protected:
  std::shared_ptr<DMatrix<double>> R_{}; // R = R1^T*R0^{-1}*R1
  std::shared_ptr<DMatrix<double>> T_{}; // T = \Psi^T*Q*\Psi + \lambda*K

public:
  // constructor
  iGCV() = default;
  // performs computation of matrix T
  virtual std::shared_ptr<DMatrix<double>> T() = 0;
  // getters
  std::shared_ptr<DMatrix<double>> R() const { return R_; }
  std::shared_ptr<DMatrix<double>> T() const { return T_; }

  virtual ~iGCV() = default;
};

#endif // __I_GCV_H__

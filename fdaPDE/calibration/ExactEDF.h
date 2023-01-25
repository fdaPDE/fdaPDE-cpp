#ifndef __EXACT_EDF_H__
#define __EXACT_EDF_H__

#include <Eigen/LU>
#include "../core/utils/Symbols.h"
#include "../models/regression/RegressionBase.h"
using fdaPDE::models::is_regression_model;

namespace fdaPDE {
namespace calibration{

  // Evaluates exactly the trace of matrix S = \Psi*T^{-1}*\Psi^T*Q. This is obtained by explicitly computing S and extracting its trace.
  template <typename Model>
  class ExactEDF {
    static_assert(is_regression_model<Model>::value);
  private:
    Model& model_;

    // compute matrix S = \Psi*T^{-1}*\Psi^T*Q
    const DMatrix<double>& S() {
      // compute \Psi^T*D*Q (take into account of areal sampling)
      if(model_.hasCovariates())
	E_ = model_.PsiTD()*model_.Q();
      else E_ = model_.PsiTD();
      
      // factorize matrix T
      invT_ = model_.T().partialPivLu();
      V_ = invT_.solve(E_); // V = invT*E = T^{-1}*\Psi^T*Q
      S_ = model_.Psi()*V_; // S = \Psi*V
      return S_;
    };    

  public:
    // let public access since those might be required for the computation of GCV derivatives
    // in the following R = R1^T*R0^{-1}*R1 should be correctly exposed by model M
    Eigen::PartialPivLU<DMatrix<double>> invT_{}; // T^{-1}
    DMatrix<double> E_{}; // \Psi^T*Q
    DMatrix<double> V_{}; // T^{-1}*\Psi^T*Q
    DMatrix<double> S_{}; // \Psi*T^{-1}*\Psi^T*Q = \Psi*V_

    // constructor
    explicit ExactEDF(Model& model) : model_(model) {};
    double compute() { return S().trace(); } // computes Tr[S]
  };

}}
#endif // __EXACT_EDF_H__

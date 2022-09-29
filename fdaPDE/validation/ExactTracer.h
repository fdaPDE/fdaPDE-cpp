#ifndef __EXACT_TRACER_H__
#define __EXACT_TRACER_H__

#include <memory>
#include <Eigen/LU>

#include "../core/utils/Symbols.h"

// Evaluates exactly the trace of matrix S = \Psi*T^{-1}*\Psi^T*Q. This is obtained by explicitly computing S and extracting its trace.
// Exact evaluation of GCV requires also additional quantites which are evaluated here since strictly releated with partial results obtained
// during computation of Tr[S]
class ExactTracer {
private:
  Eigen::PartialPivLU<DMatrix<double>> invT_{}; // T^{-1}

  // in the following R = R1^T*R0^{-1}*R1 should be correctly exposed by model M
  DMatrix<double> E_; // \Psi^T*Q
  DMatrix<double> V_; // T^{-1}*\Psi^T*Q
  DMatrix<double> L_; // T^{-1}*R

public:
  // constructor
  ExactTracer() = default;
  
  // compute matrix S = \Psi*T^{-1}*\Psi^T*Q
  template <typename M>
  std::unique_ptr<DMatrix<double>> S(M& model) {
    // compute \Psi^T*Q
    // an optimization can be put in place if \Psi is a permutation matrix (locations are a subset of nodes)
    if(model.q() != 0)
      E_ = model.Psi()->transpose()*(*model.Q());
    else
      E_ = model.Psi()->transpose();
    // factorize matrix T
    invT_ = model.T()->partialPivLu();

    // compute V = invT*E = T^{-1}*\Psi^T*Q
    V_ = invT_.solve(E_);
    // S = \Psi*T^{-1}*\Psi^T*Q
    // if locations are equal to mesh nodes then premultiply by \Psi reduces to the permutation of the rows of T^{-1}*\Psi^T*Q = V
    std::unique_ptr<DMatrix<double>> S = std::make_unique<DMatrix<double>>(model.loc(), model.loc());
    for(std::size_t k = 0; k < model.Psi()->outerSize(); ++k){
      for (SpMatrix<double>::InnerIterator it(*model.Psi(),k); it; ++it){
	S->row(it.row()) = V_.row(it.col());
      }
    }
    return S;
  };
  
  // evaluate the trace of S in an exact way
  template <typename M>
  double compute(M& model) {
    return S(model)->trace();
  }

  // the following methods are employed when the GCV field is optimized without numerical approximation of its derivatives
  // evaluate the trace of the derivative of S with respect to \lambda
  // compute first derivative of matrix S (used during optimization of gcv field when analytical derivatives are required).
  // Matrix S should have been already computed before call to this method
  template <typename M>
  std::shared_ptr<DMatrix<double>> dS(const M& model) {
    // \frac{dS}{d\lambda} = -\Psi*L*V = -\Psi*(T^{-1}*R)*(T^{-1}*E)
    L_ = invT_.solve(model.R());
    DMatrix<double> F = -L_*invT_.solve(E_);

    // if locations are equal to mesh nodes then premultiply by \Psi reduces to the permutation of the rows of (T^{-1}*R)*(T^{-1}*E)
    std::shared_ptr<DMatrix<double>> dS = std::make_shared<DMatrix<double>>(model.loc(), model.loc());
    for(std::size_t k = 0; k < model.Psi()->outerSize(); ++k){
      for (SpMatrix<double>::InnerIterator it(*model.Psi(),k); it; ++it){
	dS->row(it.row()) = F.row(it.col());
      }
    }
    return dS;
  }
  
  template <typename M>
  double derive(M& model) {
    return dS(model)->trace();
  }

  DMatrix<double> L() const { return L_; }
  
};

#endif // __EXACT_TRACER_H__

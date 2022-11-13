#ifndef __EXACT_GCV_ENGINE_H__
#define __EXACT_GCV_ENGINE_H__

#include <memory>
#include <Eigen/LU>
#include "../core/utils/Symbols.h"

namespace fdaPDE{
namespace calibration{

  // Evaluates exactly the trace of matrix S = \Psi*T^{-1}*\Psi^T*Q. This is obtained by explicitly computing S and extracting its trace.
  // Exact optimization of GCV requires also additional quantites which are evaluated here for efficiency reasons (some costly computations performed
  // during evaluation of Tr[S] can be cached and reused for the computation of other quantites)
  class ExactGCVEngine {
  private:
    Eigen::PartialPivLU<DMatrix<double>> invT_{}; // T^{-1}

    // in the following R = R1^T*R0^{-1}*R1 should be correctly exposed by model M
    DMatrix<double> E_{}; // \Psi^T*Q
    DMatrix<double> V_{}; // T^{-1}*\Psi^T*Q

    // first and second GCV derivatives related quantites (used only for exact GCV optimization)
    DMatrix<double> L_{}; // T^{-1}*R
    DMatrix<double> F_{}; // (T^{-1}*R)*(T^{-1}*E)
    DVector<double> h_{}; // h = (\lambda*L - I)*T^{-1}*R1^T*R0^{-1}*u
    DVector<double> p_{}; // p = \Psi*h - dS*z

    DMatrix<double> S_  {}; // S = \Psi*T^{-1}*\Psi^T*Q
    DMatrix<double> dS_ {}; // dS = -\Psi*(T^{-1}*R)*(T^{-1}*E)
    DMatrix<double> ddS_{}; // ddS = 2*\Psi*L*F

    // compute matrix S = \Psi*T^{-1}*\Psi^T*Q
    template <typename M>
    void S(M& model) {
      // compute \Psi^T*Q
      if(model.hasCovariates())
	E_ = model.Psi().transpose()*model.Q();
      else
	E_ = model.Psi().transpose();
    
      // factorize matrix T
      invT_ = model.T().partialPivLu();
      V_ = invT_.solve(E_); // V = invT*E = T^{-1}*\Psi^T*Q
      S_ = model.lmbPsi(V_);
      return;
    };

    // the following methods are employed when the GCV field is optimized without numerical approximation of its derivatives
    // Matrix S should have been already computed before call to these methods

    // compute first derivative of matrix S: dS = -\Psi*(T^{-1}*R)*(T^{-1}*E)
    template <typename M>
    void dS(M& model) {
      L_ = invT_.solve(model.R()); // T^{-1}*R
      F_ = L_*invT_.solve(E_);      // (T^{-1}*R)*(T^{-1}*E)
      dS_ = model.lmbPsi(-F_); // this takes into account of sampling strategy
      return;
    }
  
    // compute second derivative of matrix S: ddS = 2*\Psi*L*F
    // dS must have already been computed
    template <typename M>
    void ddS(M& model){
      DMatrix<double> C = 2*L_*F_; // compute temporary 2*L*F
      ddS_ = model.lmbPsi(C); // this takes into account of sampling strategy
      return;
    }

  public:
    // constructor
    ExactGCVEngine() = default;

    // computes Tr[S]
    template <typename M> double compute(M& model) {
      S(model); // compute matrix S
      return S_.trace();
    }
    
    // expression of GCV first derivative is
    // dGCV(\lambda) = \frac{2n}{(n - (q + Tr[S]))^2}[ \sigma^2 * Tr[dS] + a ]

    // computes the a term in the dGCV expression, given by
    // a = p.dot(z - \hat z)
    //   p = \Psi*h - t
    //     h = (\lambda*L - I)*T^{-1}*g
    //       g = R1^T*R0^{-1}*u
    //     t = dS*z
    template <typename M>
    double a(M& model){
      // NB: dS_ must already contain valid data
      DMatrix<double> g = model.R1().transpose()*model.invR0().solve(model.u());
      // cache h and p since needed for computation of second derivative
      h_ = (model.lambda()*L_ - DMatrix<double>::Identity(model.loc(), model.loc()))*invT_.solve(g);
      p_ = model.Psi()*h_ - dS_*model.z();
      // return a = p.dot(z - \hat z)
      return (( model.z() - model.fitted() ).transpose() * p_).coeff(0,0);
    }

    // computes Tr[dS]
    template <typename M> double derive(M& model) {
      dS(model); // compute first derivative of S
      return dS_.trace();
    }
  
    // expression of GCV second derivative is (edf = n - (q + Tr[S]))
    // ddGCV(\lambda) = \frac{2n}{edf^2}[ \frac{1}{edf}(3*\sigma^2*Tr[dS] + 4*a)*Tr[dS] + \sigma^2*Tr[ddS] + b ]

    // computes the b term in the ddGCV expression, given by
    // b = p.dot(Q*p) + (-ddS*z - 2*\Psi*L*h).dot(z - \hat z) 
       //   p = \Psi*h - t
       //     h = (\lambda*L - I)*T^{-1}*g
       //       g = R1^T*R0^{-1}*u
       //     t = dS*z
    template <typename M>
    double b(M& model){
      // NB: ddS_ must already contain valid data
      DMatrix<double> C = 2*L_*h_;
      // perform efficient multiplication by permutation matrix Psi
      DMatrix<double> D(model.loc(), 1); // 2*\Psi*L*h
      for(std::size_t k = 0; k < model.Psi().outerSize(); ++k){
	for(SpMatrix<double>::InnerIterator it(model.Psi(),k); it; ++it){
	  D.row(it.row()) = C.row(it.col());
	}
      }
      DVector<double> Qp_ = model.lmbQ(p_); // efficient computation of Q*p
      // return b = p.dot(Q*p) + (-ddS*z - 2*\Psi*L*h).dot(z - \hat z)
      return (( model.z() - model.fitted() ).transpose() * ( -ddS_*model.z() - D )).coeff(0,0) + p_.dot(Qp_);
    }

    // computes Tr[ddS]
    template <typename M> double deriveTwice(M& model) {
      ddS(model); // compute second derivative of S
      return ddS_.trace();
    }
  };
}}
  
#endif // __EXACT_GCV_ENGINE_H__

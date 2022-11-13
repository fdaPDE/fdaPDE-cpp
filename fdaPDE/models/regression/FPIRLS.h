#ifndef __FPIRLS_H__
#define __FPIRLS_H__

#include "../../core/utils/Symbols.h"
#include "../../core/NLA/SparseBlockMatrix.h"
using fdaPDE::core::NLA::SparseBlockMatrix;
#include "../../core/NLA/SMW.h"
using fdaPDE::core::NLA::SMW;
#include "Distributions.h"
#include "iRegressionModel.h"
using fdaPDE::models::is_regression_model;
#include <cstddef>

#include <chrono>

namespace fdaPDE{
namespace models{

  // a general implementation of the Functional Penalized Iterative Reweighted Least Square (FPIRLS) algorithm
  template <typename Distribution>
  class FPIRLS {
  private:
    // data characterizing the behaviour of the algorithm
    Distribution distribution_{};
    
    double tolerance_;
    std::size_t max_iter_;
    
    // let g() the link function of the considered distribution and y = (y_1, y_2, ..., y_n) the (1 x n) vector of observations
    DVector<double> mu_{};     // \mu^k = [ \mu^k_1, ..., \mu^k_n ] : mean vector at step k
    DVector<double> theta_{};  // \theta^k = [ g(\mu^k_1), ..., g(\mu^k_n) ]
    DVector<double> G_{};      // G^k = diag(g'(\mu^k_1), ..., g'(\mu^k_n))
    DVector<double> py_{};     // \tilde z^k = G^k(y-u^k) + \theta^k : pseudo-observations vector at step k
    DVector<double> V_{};      // V^k = diag(v(\mu^k_1), ..., v(\mu^k_n)) : variance matrix at step k
    DVector<double> W_{};      // W^k = ((G^k)^{-2})*((V^k)^{-1})
    
    DVector<double> f_{};
    DVector<double> g_{};
    DVector<double> beta_{};
    
  public:
    // constructor
    FPIRLS(double tolerance, std::size_t max_iter)
      : tolerance_(tolerance), max_iter_(max_iter) {};
    
    // in general observations are not known at construction time
    template <typename M>
    void compute(M& m_) {
      static_assert(is_regression_model<M>::value);
      // get number of data and preallocate space
      std::size_t n = m_.z().rows();
      theta_.resize(n); G_.resize(n); py_.resize(n); V_.resize(n); W_.resize(n);

      // algorithm initialization
      mu_ = m_.z();
      distribution_.preprocess(mu_);
      std::size_t k = 0; // current iteration

      // assemble system matrix for the nonparameteric part of the model. Develop here once,
      // the scheme will modifiy at each iteration only the nord-west block
      std::size_t N = m_.loc();
      SparseBlockMatrix<double,2,2>
	A(SpMatrix<double>(N,N), m_.lambda() * m_.R1().transpose(),
	  m_.lambda() * m_.R1(), m_.lambda() * m_.R0()            );
      SpMatrix<double> A_ = A.derived();

      DMatrix<double> b_(A_.rows(),1);
      b_ << DMatrix<double>::Zero(N,1),
	m_.lambda()*m_.u();

      // algorithm stops when an enought small difference between two consecutive values of the functional to
      // minimize are recordered
      double J_old = tolerance_+1; double J_new = 0;

      while(k < max_iter_ && std::abs(J_new - J_old) > tolerance_){
        theta_ = distribution_.link(mu_);
	G_ = distribution_.der_link(mu_);
	V_ = distribution_.variance(mu_);
	W_ = ((G_.array().pow(2)*V_.array()).inverse()).matrix();
	// compute pseudo observations
	py_ = G_.asDiagonal()*(m_.z() - mu_) + theta_;

	// solve weighted least square problem
	// \argmin_{\beta, f} [ \norm(W^{1/2}(y - X\beta - f_n))^2 + \lambda \int_D (Lf - u)^2 ]

	// assemble system matrix for the nonparameteric part of the model
	SparseBlockMatrix<double,2,2>
	  B_(-m_.PsiTD()*W_.asDiagonal()*m_.Psi(), SpMatrix<double>(N,N),
	     SpMatrix<double>(N,N),                SpMatrix<double>(N,N));
	SpMatrix<double> C_ = A_ + B_.derived();
	C_.makeCompressed();

        DVector<double> sol; // room for problem' solution
  
	if(!m_.hasCovariates()){ // nonparametric case
          // define rhs for this iteration
          b_.topRows(N) = -m_.PsiTD()*W_.asDiagonal()*py_;
	  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
	  solver.compute(C_); 
	  
	  // solve linear system C_*x = b_
          sol = solver.solve(b_);
	  f_ = sol.head(A_.rows()/2); // spatial-field estimate	
	}else{ // semi-parametric case
	  // compute q x q dense matrix X^T*W_*X and its factorization
	  DMatrix<double> XTX_ = m_.W().transpose()*W_.asDiagonal()*m_.W();
	  Eigen::PartialPivLU<DMatrix<double>> invXTX_ = XTX_.partialPivLu();
 
	  // perform efficient multiplication of pseudo_y_ by Q
	  // compute P*x - P*W*z = P*x - (P*W*(W^T*P*W)^{-1}*W^T*P)*x = P(I - H)*x = Q*x
	  DMatrix<double> Qy = W_.asDiagonal()*(py_ - m_.W()*invXTX_.solve(m_.W().transpose()*W_.asDiagonal()*py_));

	  // rhs of SR-PDE linear system
	  b_.topRows(N) = -m_.PsiTD()*Qy; // -\Psi^T*D*Q*z

	  // definition of matrices U and V for application of woodbury formula
	  DMatrix<double> U = DMatrix<double>::Zero(A_.rows(), m_.q());
	  U.block(0,0, A_.rows()/2, m_.q()) = m_.PsiTD()*W_.asDiagonal()*m_.W();
	  DMatrix<double> V = DMatrix<double>::Zero(m_.q(), A_.rows());
	  V.block(0,0, m_.q(), A_.rows()/2) = m_.W().transpose()*W_.asDiagonal()*m_.Psi();

	  // Define system solver. Use SMW solver from NLA module
	  SMW<> solver{};
	  solver.compute(C_);
	  // solve system Mx = b
	  sol = solver.solve(U, XTX_, V, b_);
	  // extract solution estimates
	  f_ = sol.head(A_.rows()/2); // spatial-field estimate
	  beta_ = invXTX_.solve(m_.W().transpose()*W_.asDiagonal())*(py_ - m_.Psi()*f_);
	}
	// extract estimates for the non-parametric part of the model
	g_ = sol.tail(A_.rows()/2); // PDE misfit
	
	// update value of \mu_
	DVector<double> fitted = m_.hasCovariates() ?
	  DVector<double>(m_.Psi()*f_ + m_.W()*beta_) : DVector<double>(m_.Psi()*f_);
	for(std::size_t i = 0; i < n; ++i)
	  mu_[i] = distribution_.inv_link(fitted[i]);
	
	// compute value of functional J for this pair (\beta, f): \norm{V^{-1/2}(y - \mu)}^2 + \int_D (Lf-u)^2
	DVector<double> V = distribution_.variance(mu_).array().sqrt().inverse().matrix();
	double J = (V.asDiagonal()*(m_.z() - mu_)).squaredNorm() + m_.lambda()*g_.dot(m_.R0()*g_); // \int_D (Lf-u)^2

	// prepare for next iteration
	k++;
	J_old = J_new; J_new = J;
      }
      return;
    }

    DVector<double> weights() const { return W_; }
    DVector<double> beta() const { return beta_; }
    DVector<double> f() const { return f_; }
    
  };
  
}}



#endif // __FPIRLS_H__

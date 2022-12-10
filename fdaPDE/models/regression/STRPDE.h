#ifndef __STRPDE_H__
#define __STRPDE_H__

#include <memory>
#include <type_traits>
#include <unsupported/Eigen/src/SparseExtra/MarketIO.h>
// CORE imports
#include "../../core/utils/Symbols.h"
#include "../../core/FEM/PDE.h"
#include "../space_time/SpaceTime.h"
#include "../space_time/SplineBasis.h"
#include "../iStatModel.h"
using fdaPDE::core::FEM::PDEBase;
#include "../../core/NLA/SparseBlockMatrix.h"
using fdaPDE::core::NLA::SparseBlockMatrix;
#include "../../core/NLA/SMW.h"
using fdaPDE::core::NLA::SMW;
// calibration module imports
#include "../../calibration/iGCV.h"
using fdaPDE::calibration::iGCV;
// regression module imports
#include "iRegressionModel.h"
using fdaPDE::models::iRegressionModel;
#include "../iSpaceTimeModel.h"
#include "../space_time/SpaceTime.h"
using fdaPDE::models::PhiAssembler;
#include "../space_time/TimeAssembler.h"
using fdaPDE::models::TimeAssembler;

#include "../../core/NLA/KroneckerProduct.h"
using fdaPDE::core::NLA::Kronecker;

namespace fdaPDE{
namespace models{

  // base class for STRPDE model
  template <typename PDE, typename TimeRegularization> class STRPDE;

  template <typename PDE>
  class STRPDE<PDE, SpaceTimeSeparable> : public iRegressionModel<STRPDE<PDE, SpaceTimeSeparable>>/*, public iGCV*/ {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
  private:
    typedef SpaceTimeSeparable TimeRegularization;
    typedef iRegressionModel<STRPDE<PDE, TimeRegularization>> Base;
    
    // diagonal matrix of weights (implements possible heteroscedasticity)
    DiagMatrix<double> W_;
    // system matrix of non-parametric problem (2N x 2N matrix)
    //     | -\Psi^T*D*W*\Psi  \lambda*R1^T |
    // A = |                                |
    //     |    \lambda*R1     \lambda*R0   |
    SpMatrix<double> A_{};
    // right hand side of problem's linear system (1 x 2N vector)
    //     | -\Psi^T*D*Q*y |
    // b = |               |, Q = W(I-H), H = X*(X^T*W*X)^{-1}*X^T*W
    //     |   \lambda*u   |
    DVector<double> b_{};
    // q x q dense matrix X^T*W*X
    DMatrix<double> XTX_{};
    // partial LU (with pivoting) factorization of the dense (square invertible) q x q matrix XTX_.
    Eigen::PartialPivLU<DMatrix<double>> invXTX_{};

    // problem solution
    DMatrix<double> f_{};    // estimate of the spatial field (1 x N vector)
    DMatrix<double> g_{};    // PDE misfit
    DMatrix<double> beta_{}; // estimate of the coefficient vector (1 x q vector)

    // perform proper initialization of model
    void init();
  public:
    // import commonly defined symbols from base
    IMPORT_REGRESSION_SYMBOLS;
    using Base::lambdaS;
    using Base::lambdaT;
    // import symbols needed for separable regularization
    using Base::Pt; // time penalization matrix
    using Base::Rt; // time mass matrix
    using Base::time_;
    
    // constructor
    STRPDE() = default;
    STRPDE(const PDE& pde, const DMatrix<double>& time)
      : iRegressionModel<STRPDE<PDE, SpaceTimeSeparable>>(pde, time) {};
    
    // iStatModel interface implementation
    virtual void solve(); // finds a solution to the smoothing problem

    // iRegressionModel interface implementation
    virtual DMatrix<double> lmbQ(const DMatrix<double>& x);
    virtual DMatrix<double> fitted();
    virtual double predict(const DVector<double>& covs, const std::size_t loc) const;
    // getters to problem solution
    virtual const DMatrix<double>& f() const { return f_; };
    virtual const DMatrix<double>& g() const { return g_; };
    virtual const DMatrix<double>& beta() const { return beta_; };
    
    // iGCV interface implementation
    // virtual const DMatrix<double>& T(); // T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
    // virtual const DMatrix<double>& Q(); // Q = W(I - H) = W - W*X*(X^T*W*X)^{-1}X^T*W
    // returns the euclidian norm of y - \hat y
    // virtual double norm(const DMatrix<double>& obs, const DMatrix<double>& fitted) const;
    
    virtual ~STRPDE() = default;
  };

  // compile time informations related to the model
  template <typename PDE_>
  struct model_traits<STRPDE<PDE_, SpaceTimeSeparable>> {
    typedef PDE_ PDE;
    typedef SpaceTimeSeparable RegularizationType;
    typedef SplineBasis<3> TimeBasis;
  };

  template<typename PDE>
  DMatrix<double> STRPDE<PDE, SpaceTimeSeparable>::fitted() {
    return DMatrix<double>::Zero(1,1);
  }

  template<typename PDE>
  DMatrix<double> STRPDE<PDE, SpaceTimeSeparable>::lmbQ(const DMatrix<double>& x) {
    return DMatrix<double>::Zero(1,1);
  }

  template<typename PDE>
  double STRPDE<PDE, SpaceTimeSeparable>::predict(const DVector<double>& covs, const std::size_t loc) const {
    return 0;
  }

  // finds a solution to the STR-PDE smoothing problem (separable penalization)
  template <typename PDE>
  void STRPDE<PDE, SpaceTimeSeparable>::solve() {
    // assemble system matrix for the nonparameteric part of the model
    SparseKroneckerProduct<> P = Kronecker(Pt(), pde_->R0());
    SparseBlockMatrix<double,2,2>
      A(-Psi().transpose()*Psi()-lambdaT()*P, lambdaS()*R1().transpose(),
	lambdaS()*R1(),                       lambdaS()*R0()            );
    // cache system matrix for reuse
    A_ = A.derived();
    b_.resize(A_.rows());
    DVector<double> sol; // room for problem' solution

    // prepare input data
    DVector<double> stacked_y(obs()*time_.rows());
    DVector<double> u_bestiale = DVector<double>::Zero(A_.rows()/2);
    for(std::size_t i = 0; i < time_.rows(); ++i){
      stacked_y.segment(i*obs(), obs()) = y().col(i);
    }

    if(!Base::hasCovariates()){ // nonparametric case
      // rhs of STR-PDE linear system
      b_ << -Psi().transpose()*stacked_y,
	    lambdaS()*u_bestiale;
    
      // define system solver. Use a sparse solver
      Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver{};
      solver.compute(A_);
      // solve linear system A_*x = b_
      sol = solver.solve(b_);

      // store result of smoothing
      f_ = sol.head(A_.rows()/2);
    }
    // }else{ // parametric case
    //   if(!Base::isAlloc(XTX_)){
    // 	// compute q x q dense matrix W^T*W and its factorization
    // 	XTX_ = X().transpose()*W_*X();
    // 	invXTX_ = XTX_.partialPivLu();
    //   }
    //   // rhs of SR-PDE linear system
    //   b_ << -PsiTD()*lmbQ(y()), // -\Psi^T*D*Q*z
    // 	lambda()*u();
    
    //   std::size_t q_ = X().cols(); // number of covariates
    //   // definition of matrices U and V  for application of woodbury formula
    //   DMatrix<double> U = DMatrix<double>::Zero(A_.rows(), q_);
    //   U.block(0,0, A_.rows()/2, q_) = PsiTD()*W_*X();
    //   DMatrix<double> V = DMatrix<double>::Zero(q_, A_.rows());
    //   V.block(0,0, q_, A_.rows()/2) = X().transpose()*W_*Psi();

    //   // Define system solver. Use SMW solver from NLA module
    //   SMW<> solver{};
    //   solver.compute(A_);
    //   // solve system Mx = b
    //   sol = solver.solve(U, XTX_, V, b_);
    //   // store result of smoothing 
    // 	   f_    = sol.head(A_.rows()/2);
    //   beta_ = invXTX_.solve(X().transpose()*W_)*(y() - Psi()*f_);
    // }
    // store PDE misfit
    g_ = sol.tail(A_.rows()/2);
    return;
  }
  
#include "STRPDE.tpp"
}}


#endif // __STRPDE_H__

#ifndef __SRPDE_H__
#define __SRPDE_H__

#include <memory>
#include <type_traits>
// CORE imports
#include "../../core/utils/Symbols.h"
#include "../../core/FEM/PDE.h"
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

namespace fdaPDE{
namespace models{
  
  template <typename PDE>
  class SRPDE : public iRegressionModel<PDE>, public iGCV {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
  private:
    // system matrix of non-parametric problem (2N x 2N matrix)
    //     | -\Psi^T*D*\Psi  \lambda*R1^T |
    // A = |                              |
    //     |   \lambda*R1    \lambda*R0   |
    SpMatrix<double> A_{};
    // right hand side of problem's linear system (1 x 2N vector)
    //     | -\Psi^T*D*Q*y |
    // b = |               |
    //     |   \lambda*u   |
    DVector<double> b_{};
    // q x q dense matrix X^T*X
    DMatrix<double> XTX_{};
    // partial LU (with pivoting) factorization of the dense (square invertible) q x q matrix XTX_.
    Eigen::PartialPivLU<DMatrix<double>> invXTX_{};

    // problem solution
    DMatrix<double> f_{};    // estimate of the spatial field (1 x N vector)
    DMatrix<double> g_{};    // PDE misfit
    DMatrix<double> beta_{}; // estimate of the coefficient vector (1 x q vector)
  public:
    IMPORT_REGRESSION_MODEL_SYMBOLS(PDE);
    
    // constructor
    SRPDE() = default;
    SRPDE(const PDE& pde, double lambda)
      : iRegressionModel<PDE>(pde, lambda) {};

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
    virtual const DMatrix<double>& T(); // T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
    virtual const DMatrix<double>& Q(); // Q = I - H = I - X*(X^T*X)^{-1}X^T
    // returns the euclidian norm of y - \hat y
    virtual double norm(const DMatrix<double>& obs, const DMatrix<double>& fitted) const;
    
    virtual ~SRPDE() = default;
  };

#include "SRPDE.tpp"
}}
    
#endif // __SRPDE_H__

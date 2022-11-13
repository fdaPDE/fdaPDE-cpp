#ifndef __GSRPDE_H__
#define __GSRPDE_H__

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

#include "FPIRLS.h"
using fdaPDE::models::FPIRLS;

#include <chrono>

namespace fdaPDE{
namespace models{
  
  template <typename PDE, typename Distribution>
  class GSRPDE : public iRegressionModel<PDE>/*, public iGCV*/ {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
  private:
    // FPIRLS engine
    FPIRLS<Distribution> fpirls;
    // weight matrix obtained at FPIRLS convergence
    DiagMatrix<double> P_;
    // q x q dense matrix W^T*P*W
    DMatrix<double> WTW_{};
    // n x n projection matrix onto the orthogonal space of Im(W), Q_ = I - H_ = I - W*(W*T*P*W)^{-1}*W^T*P
    DMatrix<double> Q_{};
    // partial LU (with pivoting) factorization of the dense (square invertible) q x q matrix WTW_.
    Eigen::PartialPivLU<DMatrix<double>> invWTW_{};

    // problem solution
    DMatrix<double> f_{};    // estimate of the spatial field (1 x N vector)
    DMatrix<double> beta_{}; // estimate of the coefficient vector (1 x q vector)
  public:
    IMPORT_REGRESSION_MODEL_SYMBOLS(PDE);
    
    // constructor
    GSRPDE() = default;
    GSRPDE(const PDE& pde, double lambda, double tolerance, std::size_t max_iter)
      : iRegressionModel<PDE>(pde, lambda), fpirls(tolerance, max_iter) {};

    // iStatModel interface implementation
    virtual void solve(); // finds a solution to the smoothing problem

    // iRegressionModel interface implementation
    virtual DMatrix<double> lmbQ(const DMatrix<double>& x);
    virtual DMatrix<double> fitted() const;
    virtual double predict(const DVector<double>& covs, const std::size_t loc) const;
    // getters to problem solution
    virtual const DMatrix<double>& f() const { return f_; };
    virtual const DMatrix<double>& beta() const { return beta_; };
    
    // iGCV interface implementation
    virtual const DMatrix<double>& T(); // T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
    virtual const DMatrix<double>& Q();
    // returns the total deviance of the model as \sum dev(y - \hat \mu)
    virtual double norm(const DMatrix<double>& obs, const DMatrix<double>& fitted) const;

    virtual ~GSRPDE() = default;
  };
  
  #include "GSRPDE.tpp"
}}

#endif // __GSRPDE_H__

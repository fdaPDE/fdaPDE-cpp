#ifndef __SRPDE_H__
#define __SRPDE_H__

#include <memory>
#include <type_traits>
// CORE imports
#include "../../core/utils/Symbols.h"
#include "../../core/FEM/PDE.h"
using fdaPDE::core::FEM::PDEBase;
#include "../ModelBase.h"
#include "../../core/NLA/SparseBlockMatrix.h"
using fdaPDE::core::NLA::SparseBlockMatrix;
#include "../../core/NLA/SMW.h"
using fdaPDE::core::NLA::SMW;
// calibration module imports
#include "../../calibration/iGCV.h"
using fdaPDE::calibration::iGCV;
// regression module imports
#include "../SamplingDesign.h"
#include "RegressionBase.h"
using fdaPDE::models::RegressionBase;

namespace fdaPDE{
namespace models{
  
  template <typename PDE, Sampling SamplingDesign>
  class SRPDE : public RegressionBase<SRPDE<PDE, SamplingDesign>>, public iGCV {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
  private:
    typedef RegressionBase<SRPDE<PDE, SamplingDesign>> Base;

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

    // problem solution
    DMatrix<double> f_{};    // estimate of the spatial field (1 x N vector)
    DMatrix<double> g_{};    // PDE misfit
    DMatrix<double> beta_{}; // estimate of the coefficient vector (1 x q vector)
  public:
    IMPORT_REGRESSION_SYMBOLS;
    using Base::lambda; // smoothing parameter in space
    // constructor
    SRPDE() = default;
    template <typename... SamplingData>
    SRPDE(const PDE& pde, const SamplingData&... s) : RegressionBase<SRPDE<PDE, SamplingDesign>>(pde, s...) {};
    
    // iStatModel interface implementation
    virtual void solve(); // finds a solution to the smoothing problem

    // RegressionBase interface implementation
    virtual DMatrix<double> fitted();
    virtual double predict(const DVector<double>& covs, const std::size_t loc) const;
    // getters to problem solution
    virtual const DMatrix<double>& f() const { return f_; };
    virtual const DMatrix<double>& g() const { return g_; };
    virtual const DMatrix<double>& beta() const { return beta_; };
    
    // iGCV interface implementation
    virtual const DMatrix<double>& T(); // T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
    virtual const DMatrix<double>& Q(); // Q = W(I - H) = W - W*X*(X^T*W*X)^{-1}X^T*W
    // returns the euclidian norm of y - \hat y
    virtual double norm(const DMatrix<double>& obs, const DMatrix<double>& fitted) const;
    
    virtual ~SRPDE() = default;
  };

  // compile time informations related to the model
  template <typename PDE_, Sampling SamplingDesign>
  struct model_traits<SRPDE<PDE_, SamplingDesign>> {
    typedef PDE_ PDE;
    typedef SpaceOnlyTag RegularizationType;
    static constexpr Sampling sampling = SamplingDesign;
  };

#include "SRPDE.tpp"
}}
    
#endif // __SRPDE_H__

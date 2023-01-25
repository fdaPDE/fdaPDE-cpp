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
    SpMatrix<double> A_{}; // system matrix of non-parametric problem (2N x 2N matrix)
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> invA_; // factorization of matrix A
    DVector<double> b_{}; // right hand side of problem's linear system (1 x 2N vector)

    // matrices related to woodbury decomposition
    DMatrix<double> U_{};
    DMatrix<double> V_{};
    
  public:
    IMPORT_REGRESSION_SYMBOLS;
    using Base::lambda; // smoothing parameter in space
    // constructor
    SRPDE() = default;
    template <typename... SamplingData>
    SRPDE(const PDE& pde, const SamplingData&... s) : RegressionBase<SRPDE<PDE, SamplingDesign>>(pde, s...) {};
    
    // iStatModel interface implementation
    virtual void solve(); // finds a solution to the smoothing problem
    
    // iGCV interface implementation
    virtual const DMatrix<double>& T(); // T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
    virtual const DMatrix<double>& Q(); // Q = W(I - H) = W - W*X*(X^T*W*X)^{-1}X^T*W
    // returns the euclidian norm of y - \hat y
    virtual double norm(const DMatrix<double>& obs, const DMatrix<double>& fitted) const;

    // getters
    const SpMatrix<double>& A() const { return A_; }
    const Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>& invA() const { return invA_; }
    const DMatrix<double>& U() const { return U_; }
    const DMatrix<double>& V() const { return V_; }
    
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

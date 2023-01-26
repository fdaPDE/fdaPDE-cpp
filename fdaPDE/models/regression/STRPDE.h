#ifndef __STRPDE_H__
#define __STRPDE_H__

#include <memory>
#include <type_traits>
// CORE imports
#include "../../core/utils/Symbols.h"
#include "../../core/FEM/PDE.h"
using fdaPDE::core::FEM::PDEBase;
#include "../space_time/SplineBasis.h"
#include "../../core/NLA/SparseBlockMatrix.h"
using fdaPDE::core::NLA::SparseBlockMatrix;
#include "../../core/NLA/SMW.h"
using fdaPDE::core::NLA::SMW;
#include "../../core/NLA/KroneckerProduct.h"
using fdaPDE::core::NLA::SparseKroneckerProduct;
using fdaPDE::core::NLA::Kronecker;
// calibration module imports
#include "../../calibration/iGCV.h"
using fdaPDE::calibration::iGCV;
// regression module imports
#include "RegressionBase.h"
using fdaPDE::models::RegressionBase;
#include "../ModelTraits.h"

namespace fdaPDE{
namespace models{

  // base class for STRPDE model
  template <typename PDE, typename TimeRegularization, Sampling SamplingDesign> class STRPDE;

  // implementation of STRPDE for separable space-time regularization
  template <typename PDE, Sampling SamplingDesign>
  class STRPDE<PDE, SpaceTimeSeparableTag, SamplingDesign>
    : public RegressionBase<STRPDE<PDE, SpaceTimeSeparableTag, SamplingDesign>>/*, public iGCV*/ {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
  private:
    typedef SpaceTimeSeparableTag TimeRegularization;
    typedef RegressionBase<STRPDE<PDE, TimeRegularization, SamplingDesign>> Base;
    SpMatrix<double> A_{}; // system matrix of non-parametric problem (2N x 2N matrix)
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> invA_; // factorization of matrix A
    DVector<double> b_{};  // right hand side of problem's linear system (1 x 2N vector)
  public:
    // import commonly defined symbols from base
    IMPORT_REGRESSION_SYMBOLS;
    using Base::lambdaS; // smoothing parameter in space
    using Base::lambdaT; // smoothing parameter in time
    using Base::Pt;      // time penalization matrix: [Pt_]_{ij} = \int_{[0,T]} (\phi_i)_tt*(\phi_j)_tt
    using Base::Rt;      // time mass matrix: [Rt_]_{ij} = \int_{[0,T]} \phi_i*\phi_j    
    // constructor
    STRPDE() = default;
    template <typename... SamplingData>
    STRPDE(const PDE& pde, const DMatrix<double>& time, const SamplingData&... s)
      : RegressionBase<STRPDE<PDE, SpaceTimeSeparableTag, SamplingDesign>>(pde, time, s...) {};
    
    // ModelBase interface implementation
    virtual void solve(); // finds a solution to the smoothing problem

    // iGCV interface implementation
    // virtual const DMatrix<double>& T(); // T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
    // virtual const DMatrix<double>& Q(); // Q = W(I - H) = W - W*X*(X^T*W*X)^{-1}X^T*W
    // returns the euclidian norm of y - \hat y
    // virtual double norm(const DMatrix<double>& obs, const DMatrix<double>& fitted) const;
    
    virtual ~STRPDE() = default;
  };
  template <typename PDE_, Sampling SamplingDesign>
  struct model_traits<STRPDE<PDE_, SpaceTimeSeparableTag, SamplingDesign>> {
    typedef PDE_ PDE;
    typedef SpaceTimeSeparableTag RegularizationType;
    typedef SplineBasis<3> TimeBasis; // use cubic B-splines
    static constexpr Sampling sampling = SamplingDesign;
  };

  // implementation of STRPDE for parabolic space-time regularization
  template <typename PDE, Sampling SamplingDesign>
  class STRPDE<PDE, SpaceTimeParabolicTag, SamplingDesign>
    : public RegressionBase<STRPDE<PDE, SpaceTimeParabolicTag, SamplingDesign>>/*, public iGCV*/ {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
  private:
    typedef SpaceTimeParabolicTag TimeRegularization;
    typedef RegressionBase<STRPDE<PDE, SpaceTimeParabolicTag, SamplingDesign>> Base;
    SpMatrix<double> A_{}; // system matrix of non-parametric problem (2N x 2N matrix)
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> invA_; // factorization of matrix A
    DVector<double> b_{};  // right hand side of problem's linear system (1 x 2N vector)
  public:
    // import commonly defined symbols from base
    IMPORT_REGRESSION_SYMBOLS;
    using Base::lambdaS; // smoothing parameter in space
    using Base::lambdaT; // smoothing parameter in time
    using Base::L;       // [L]_{ii} = 1/DeltaT for i \in {1 ... m} and [L]_{i,i-1} = -1/DeltaT for i \in {1 ... m-1}
    using Base::n_time;  // number of time instants m defined over [0,T]
    using Base::s;       // initial condition
    // constructor
    STRPDE() = default;
    template <typename... SamplingData>
    STRPDE(const PDE& pde, const DMatrix<double>& time, const SamplingData&... s)
      : RegressionBase<STRPDE<PDE, SpaceTimeParabolicTag, SamplingDesign>>(pde, time, s...) {};
    
    // ModelBase interface implementation
    virtual void solve(); // finds a solution to the smoothing problem

    // iGCV interface implementation
    // virtual const DMatrix<double>& T(); // T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
    // virtual const DMatrix<double>& Q(); // Q = W(I - H) = W - W*X*(X^T*W*X)^{-1}X^T*W
    // returns the euclidian norm of y - \hat y
    // virtual double norm(const DMatrix<double>& obs, const DMatrix<double>& fitted) const;
    
    virtual ~STRPDE() = default;
  };
  template <typename PDE_, Sampling SamplingDesign>
  struct model_traits<STRPDE<PDE_, SpaceTimeParabolicTag, SamplingDesign>> {
    typedef PDE_ PDE;
    typedef SpaceTimeParabolicTag RegularizationType;
    static constexpr Sampling sampling = SamplingDesign;
  };
  
#include "STRPDE.tpp"
}}
#endif // __STRPDE_H__

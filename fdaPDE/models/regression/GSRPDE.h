#ifndef __GSRPDE_H__
#define __GSRPDE_H__

#include <memory>
#include <type_traits>
// CORE imports
#include "../../core/utils/Symbols.h"
#include "../../core/FEM/PDE.h"
using fdaPDE::core::FEM::PDEBase;
#include "../ModelBase.h"
#include "../ModelTraits.h"
// calibration module imports
#include "../../calibration/iGCV.h"
using fdaPDE::calibration::iGCV;
#include "FPIRLS.h"
using fdaPDE::models::FPIRLS;

namespace fdaPDE{
namespace models{

  // base class for GSRPDE model
  template <typename PDE, typename RegularizationType, Sampling SamplingDesign,
	    SolverType Solver, typename Distribution>
  class GSRPDE : public RegressionBase<GSRPDE<PDE, RegularizationType, SamplingDesign, Solver, Distribution>>/*,
    public iGCV*/ {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
  private:
    typedef RegressionBase<GSRPDE<PDE, RegularizationType, SamplingDesign, Solver, Distribution>> Base;   
    // weight matrix obtained at FPIRLS convergence
    DiagMatrix<double> W_;

    // FPIRLS parameters (set to default)
    std::size_t max_iter_ = 15;
    double tol_ = 0.0002020;
  public:
    IMPORT_REGRESSION_SYMBOLS;
    using Base::lambda; // smoothing parameter in space
    // constructor
    GSRPDE() = default;
    template <typename... SamplingData>
    GSRPDE(const PDE& pde, const SamplingData&... s) : Base(pde, s...) {};

    // setter
    void setFPIRLSTolerance(double tol) { tol_ = tol; }
    void setFPIRLSMaxIterations(std::size_t max_iter) { max_iter_ = max_iter; }
    
    // ModelBase implementation
    virtual void solve(); // finds a solution to the smoothing problem
    
    // // iGCV interface implementation
    // virtual const DMatrix<double>& T(); // T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
    // virtual const DMatrix<double>& Q();
    // // returns the total deviance of the model as \sum dev(y - \hat \mu)
    // virtual double norm(const DMatrix<double>& obs, const DMatrix<double>& fitted) const;

    virtual ~GSRPDE() = default;
  };
  template <typename PDE_, typename RegularizationType_, Sampling SamplingDesign,
	    SolverType Solver, typename Distribution>
  struct model_traits<GSRPDE<PDE_, RegularizationType_, SamplingDesign, Solver, Distribution>> {
    typedef PDE_ PDE;
    typedef RegularizationType_ RegularizationType;
    static constexpr Sampling sampling = SamplingDesign;
    static constexpr SolverType solver = Solver;
  };
  
  #include "GSRPDE.tpp"
}}

#endif // __GSRPDE_H__

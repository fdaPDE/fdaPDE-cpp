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
  class GSRPDE : public RegressionBase<GSRPDE<PDE, RegularizationType, SamplingDesign, Solver, Distribution>>, public iGCV {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
  private:
    typedef RegressionBase<GSRPDE<PDE, RegularizationType, SamplingDesign, Solver, Distribution>> Base;   
    DiagMatrix<double> W_; // weight matrix obtained at FPIRLS convergence

    // FPIRLS parameters (set to default)
    std::size_t max_iter_ = 15;
    double tol_ = 0.0002020;
  public:
    IMPORT_REGRESSION_SYMBOLS;
    using Base::lambdaS; // smoothing parameter in space
    // constructor
    GSRPDE() = default;
    // space-only constructor
    template <typename U = RegularizationType,
	      typename std::enable_if< std::is_same<U, SpaceOnlyTag>::value, int>::type = 0> 
    GSRPDE(const PDE& pde) : Base(pde) {};
    // space-time constructor
    template <typename U = RegularizationType,
	      typename std::enable_if<!std::is_same<U, SpaceOnlyTag>::value, int>::type = 0> 
    GSRPDE(const PDE& pde, const DVector<double>& time) : Base(pde, time) {};

    // setter
    void setFPIRLSTolerance(double tol) { tol_ = tol; }
    void setFPIRLSMaxIterations(std::size_t max_iter) { max_iter_ = max_iter; }
    
    // ModelBase implementation
    void init_model() { return; }
    virtual void solve(); // finds a solution to the smoothing problem
    
    // iGCV interface implementation
    virtual const DMatrix<double>& T(); // T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
    virtual const DMatrix<double>& Q();
    // returns the total deviance of the model as \sum dev(y - \hat \mu)
    virtual double norm(const DMatrix<double>& obs, const DMatrix<double>& fitted) const;

    virtual ~GSRPDE() = default;
  };
  template <typename PDE_, typename RegularizationType_, Sampling SamplingDesign,
	    SolverType Solver, typename DistributionType_>
  struct model_traits<GSRPDE<PDE_, RegularizationType_, SamplingDesign, Solver, DistributionType_>> {
    typedef PDE_ PDE;
    typedef RegularizationType_ RegularizationType;
    static constexpr Sampling sampling = SamplingDesign;
    static constexpr SolverType solver = Solver;
    static constexpr int n_lambda = n_smoothing_parameters<RegularizationType_>::value;
    typedef DistributionType_ DistributionType;
  };

  // specialization for separable regularization
  template <typename PDE_, Sampling SamplingDesign, SolverType Solver, typename DistributionType_>
  struct model_traits<GSRPDE<PDE_, fdaPDE::models::SpaceTimeSeparableTag, SamplingDesign, Solver, DistributionType_>> {
    typedef PDE_ PDE;
    typedef fdaPDE::models::SpaceTimeSeparableTag RegularizationType;
    typedef SplineBasis<3> TimeBasis; // use cubic B-splines
    static constexpr Sampling sampling = SamplingDesign;
    static constexpr SolverType solver = Solver;
    static constexpr int n_lambda = 2;
    typedef DistributionType_ DistributionType;
  };
  
  #include "GSRPDE.tpp"
}}

#endif // __GSRPDE_H__

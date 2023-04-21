#ifndef __FPCA_H__
#define __FPCA_H__

#include <Eigen/SVD>
#include "../../core/utils/Symbols.h"
#include "../../core/OPT/optimizers/GridOptimizer.h"
using fdaPDE::core::OPT::GridOptimizer;
#include "../../calibration/GCV.h"
using fdaPDE::calibration::GCV;
#include "FunctionalBase.h"
using fdaPDE::models::FunctionalBase;
#include "ProfilingEstimation.h"

namespace fdaPDE {
namespace models {

  struct fixed_lambda {};
  struct gcv_lambda_selection {};
  struct kcv_lambda_selection {};
  
  // base class for any FPCA model
  template <typename PDE, typename RegularizationType, typename SamplingDesign, typename lambda_selection_strategy>
  class FPCA :
    public FunctionalBase< FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>> {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
  private:
    typedef FunctionalBase<FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>> Base;
    std::size_t n_pc_ = 3; // default number of principal components
    // ProfilingEstimation parameters
    double tol_ = 1e-6;         // relative tolerance between Jnew and Jold, used as stopping criterion
    std::size_t max_iter_ = 20; // maximum number of allowed iterations
    
    // problem solution
    DMatrix<double> loadings_;
    DMatrix<double> scores_;

    // tag dispatched private methods for computation of PCs, ** to be removed **
    void solve_(fixed_lambda);
    void solve_(gcv_lambda_selection);
    void solve_(kcv_lambda_selection);
  public:
    IMPORT_MODEL_SYMBOLS;
    using Base::lambda;
    using Base::lambdas;
    // constructor
    FPCA() = default;
    // space-only constructor
    template <typename U = RegularizationType,
	      typename std::enable_if< std::is_same<U, SpaceOnly>::value, int>::type = 0> 
    FPCA(const PDE& pde) : Base(pde) {};
    // space-time constructor
    template <typename U = RegularizationType,
	      typename std::enable_if<!std::is_same<U, SpaceOnly>::value, int>::type = 0> 
    FPCA(const PDE& pde, const DVector<double>& time) : Base(pde, time) {};
    
    void init_model() { return; };
    virtual void solve(); // compute principal components
    
    // getters
    const DMatrix<double>& loadings() const { return loadings_; }
    const DMatrix<double>& scores() const { return scores_; }

    // setters
    void set_tolerance(double tol) { tol_ = tol; }
    void set_max_iterations(std::size_t max_iter) { max_iter_ = max_iter; }
    void set_npc(std::size_t n_pc) { n_pc_ = n_pc; }
  };
  template <typename PDE_, typename SamplingDesign_, typename lambda_selection_strategy>
  struct model_traits<FPCA<PDE_, fdaPDE::models::SpaceOnly, SamplingDesign_, lambda_selection_strategy>> {
    typedef PDE_ PDE;
    typedef fdaPDE::models::SpaceOnly regularization;
    typedef SamplingDesign_ sampling;
    typedef fdaPDE::models::MonolithicSolver solver;
    static constexpr int n_lambda = 1;
  };
  // specialization for separable regularization
  template <typename PDE_, typename SamplingDesign_, typename lambda_selection_strategy>
  struct model_traits<FPCA<PDE_, fdaPDE::models::SpaceTimeSeparable, SamplingDesign_, lambda_selection_strategy>> {
    typedef PDE_ PDE;
    typedef fdaPDE::models::SpaceTimeSeparable regularization;
    typedef SplineBasis<3> TimeBasis; // use cubic B-splines
    typedef SamplingDesign_ sampling;
    typedef fdaPDE::models::MonolithicSolver solver;
    static constexpr int n_lambda = 2;
  };
  
  #include "fPCA.tpp"
  
}}

#endif // __FPCA_H__

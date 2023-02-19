#ifndef __FPCA_H__
#define __FPCA_H__

#include <Eigen/SVD>
#include "../../core/utils/Symbols.h"
#include "FunctionalBase.h"
using fdaPDE::models::FunctionalBase;
#include "../../calibration/GCV.h"
#include "../../calibration/KFoldCV.h"
using fdaPDE::calibration::GCV;
using fdaPDE::calibration::KFoldCV;
#include "FPIREM.h"
#include "PCScoreCV.h"
using fdaPDE::models::PCScoreCV;

namespace fdaPDE {
namespace models {

  struct fixed_lambda {};
  struct gcv_lambda_selection {};
  struct kcv_lambda_selection {};
  
  // base class for any FPCA model
  template <typename PDE, typename RegularizationType, Sampling
	    SamplingDesign, typename lambda_selection_strategy>
  class FPCA : public FunctionalBase<
    FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>> {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
  private:
    typedef FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy> ModelType;
    typedef FunctionalBase<ModelType> Base;
    typedef typename FPIREM<ModelType>::SmootherType SmootherType;
    // parameters used as stopping criterion by FPIREM algorithm
    std::size_t max_iter_ = 20; // maximum number of iterations before forced stop
    double tol_ = 1e-2;         // tolerance on |Jnew - Jold| used as convergence criterion
    std::size_t n_pc_ = 3;      // default number of principal components
    
    // problem solution
    DMatrix<double> loadings_;
    DMatrix<double> scores_;

    // tag dispatched private methods for computation of PCs
    void solve_(fixed_lambda);
    void solve_(gcv_lambda_selection);
    void solve_(kcv_lambda_selection);
    // required to support \lambda parameter selection
    std::vector<SVector<model_traits<SmootherType>::n_lambda>> lambda_vect_;
  public:
    IMPORT_MODEL_SYMBOLS;
    using Base::lambda;

    // constructor
    FPCA() = default;
    template <typename... SamplingData>
    FPCA(const PDE& pde, const SamplingData&... s) : Base(pde, s...) {};

    void init() { return; } // overload Base initialization since there is nothing to init
    virtual void solve();   // compute principal components
    
    // getters
    const DMatrix<double>& loadings() const { return loadings_; }
    const DMatrix<double>& scores() const { return scores_; }
    // setters
    void setTolerance(double tol) { tol_ = tol; }
    void setMaxIterations(std::size_t max_iter) { max_iter_ = max_iter; }
    void setNPC(std::size_t n_pc) { n_pc_ = n_pc; }
    // accepts a collection of \lambda parameters if a not fixed_lambda method is selected
    typename std::enable_if<
      !std::is_same<lambda_selection_strategy, fixed_lambda>::value,
      void>::type setLambda(const std::vector<SVector<model_traits<SmootherType>::n_lambda>>& lambda_vect) {
      lambda_vect_ = lambda_vect; }
  };
  template <typename PDE_, typename RegularizationType_,
	    Sampling SamplingDesign_, typename lambda_selection_strategy>
  struct model_traits<FPCA<PDE_, RegularizationType_, SamplingDesign_, lambda_selection_strategy>> {
    typedef PDE_ PDE;
    typedef RegularizationType_ RegularizationType;
    static constexpr Sampling sampling = SamplingDesign_;
    static constexpr SolverType solver = SolverType::Monolithic;
    static constexpr int n_lambda = n_smoothing_parameters<RegularizationType>::value;
  };

  #include "fPCA.tpp"
  
}}

#endif // __FPCA_H__

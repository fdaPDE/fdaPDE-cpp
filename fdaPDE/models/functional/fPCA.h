#ifndef __FPCA_H__
#define __FPCA_H__

#include <Eigen/SVD>
#include "../../core/utils/Symbols.h"
#include "FunctionalBase.h"
using fdaPDE::models::FunctionalBase;
#include "../ModelTraits.h"
using fdaPDE::models::SpaceOnlyTag;
using fdaPDE::models::SpaceTimeSeparableTag;
#include "../regression/SRPDE.h"
#include "../regression/STRPDE.h"

namespace fdaPDE {
namespace models {

  // trait to select internal solver
  template <typename Model>
  struct FPCA_internal_solver {
    typedef typename std::decay<Model>::type Model_;
    using type = typename std::conditional<
      !is_space_time<Model_>::value,
      // space-only problem
      SRPDE <typename model_traits<Model_>::PDE, model_traits<Model_>::sampling>,
      // space-time problem
      STRPDE<typename model_traits<Model_>::PDE, SpaceTimeSeparableTag, model_traits<Model_>::sampling,
	     SolverType::Monolithic>
      >::type;
  };
  
  // base class for any FPCA model
  template <typename PDE, typename RegularizationType, Sampling SamplingDesign>
  class FPCA : public FunctionalBase<FPCA<PDE, RegularizationType, SamplingDesign>> {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
  private:
    typedef FunctionalBase<FPCA<PDE, RegularizationType, SamplingDesign>> Base;

    std::size_t n_pc_ = 3; // number of principal components to compute, default to 3
    // parameters used as stopping criterion
    std::size_t max_iter_ = 20;
    double tol_ = 1e-3;

    // computed principal components
    DMatrix<double> loadings_;
    DMatrix<double> scores_;

  public:
    IMPORT_MODEL_SYMBOLS;
    using Base::locs;    // locations where data are sampled
    using Base::lambda;  // smoothing parameter vector
    using Base::pen;     // discretization of penalizing term
    
    // constructor
    FPCA() = default;
    template <typename... SamplingData>
    FPCA(const PDE& pde, const SamplingData&... s) : Base(pde, s...) {};
    
    // ModelBase implementation
    virtual void solve() { // finds a solution to the FPCA problem
      // define internal problem solver required for computation of the loadings
      typename FPCA_internal_solver<decltype(*this)>::type solver;
      if constexpr(std::is_same<RegularizationType, SpaceOnlyTag>::value) // space-only
	solver = typename FPCA_internal_solver<decltype(*this)>::type(pde(), locs());
      else // space-time
	solver = typename FPCA_internal_solver<decltype(*this)>::type(pde(), Base::time_domain(), locs());
      solver.setLambda(lambda());
      solver.init();

      // data are organized as a matrix n_obs() x n_subjects()
      DMatrix<double> Y = y().transpose(); // copy initial data
      loadings_.resize(Y.cols(), n_pc_);
      scores_.resize(Y.rows(), n_pc_);
      
      // compute principal components one at a time
      for(std::size_t i = 0; i < n_pc_; ++i){
	// perform Singular Value Decomposition of Y to initialize PC loadings and scores
	Eigen::JacobiSVD<DMatrix<double>> svd(Y, Eigen::ComputeThinU | Eigen::ComputeThinV); // Y = USV^T
	loadings_.col(i) = svd.matrixV().col(0);
	scores_.col(i) = svd.matrixU().col(0);
	
	// prepare for inner loop
	std::size_t iter_ = 0;
	double Jold = std::numeric_limits<double>::max();
	double Jnew = 0;

	while(std::abs(Jnew - Jold) > tol_ && iter_ < max_iter_){
	  // compute score vector s
	  scores_.col(i) = Y*loadings_.col(i);
	  scores_.col(i) = scores_.col(i)/scores_.col(i).norm();

	  // compute loadings by solving either an SRPDE or STRPDE (separable regularization) problem
	  BlockFrame<double, int> df;
	  df.insert<double>(OBSERVATIONS_BLK, Y.transpose()*scores_.col(i)); // Y^T*s
	  solver.setData(df);
	  solver.solve();

	  loadings_.col(i) = solver.f();
	  // prepare for next iteration
	  iter_++;
	  Jold = Jnew;
	  Jnew = (Y - scores_.col(i)*loadings_.col(i).transpose()).squaredNorm() + solver.g().squaredNorm();
	}
	// subtract computed PC from data
	Y = Y - scores_.col(i)*loadings_.col(i).transpose();
	// normalize loadings and unnormalize scores
	double norm = std::sqrt(loadings_.col(i).transpose()*R0()*loadings_.col(i));
	loadings_.col(i) = loadings_.col(i)/norm; 
	scores_.col(i) = scores_.col(i)*norm;
      }
    }

    // getters
    const DMatrix<double>& loadings() const { return loadings_; }
    const DMatrix<double>& scores() const { return scores_; }
    
    // setters
    void setNPC(std::size_t n_pc) { n_pc_ = n_pc; }
    void setTolerance(double tol) { tol_ = tol; }
    void setMaxIterations(std::size_t max_iter) { max_iter_ = max_iter; }
    
  };
  template <typename PDE_, typename RegularizationType_, Sampling SamplingDesign_>
  struct model_traits<FPCA<PDE_, RegularizationType_, SamplingDesign_>> {
    typedef PDE_ PDE;
    typedef RegularizationType_ RegularizationType;
    static constexpr Sampling sampling = SamplingDesign_;
    static constexpr SolverType solver = SolverType::Monolithic;
  };

  
}}

#endif // __FPCA_H__

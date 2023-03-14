#ifndef __FPIREM_H__
#define __FPIREM_H__

#include <Eigen/SVD>
#include "../../core/utils/Symbols.h"
#include "../ModelTraits.h"
using fdaPDE::models::SpaceOnly;
using fdaPDE::models::SpaceTimeSeparable;
#include "../regression/SRPDE.h"
#include "../regression/STRPDE.h"

namespace fdaPDE{
namespace models{

  // trait to select internal solver
  template <typename Model>
  struct FPIREM_internal_solver {
    typedef typename std::decay<Model>::type Model_;
    using type = typename std::conditional<
      !is_space_time<Model_>::value,
      // space-only problem
      SRPDE <typename model_traits<Model_>::PDE, model_traits<Model_>::sampling>,
      // space-time problem
      STRPDE<typename model_traits<Model_>::PDE, SpaceTimeSeparable, model_traits<Model_>::sampling,
	     SolverType::Monolithic>
      >::type;
  };
  
  // general implementation of Functional Penalized Iterative Reconstruction Error Minimization algorithm
  // finds minimum to \norm_F{Y - s^T*f}^2 + (s^T*s)*P(f) where P(f) is the penalization type of Model
  // and \norm_F{} denotes the Frobenius norm.
  template <typename Model>
  class FPIREM {
  private:
    typedef typename std::decay<Model>::type Model_;
    typename FPIREM_internal_solver<Model_>::type solver_; // internal solver
    BlockFrame<double,int> df_;  // BlockFrame of data
    
    // algorithm's parameter
    std::size_t max_iter_; // maximum number of iterations before forced stop
    double tol_;           // tolerance on |Jnew - Jold| used as convergence criterion
    // results at convergence
    DVector<double> loadings_;  // computed vector of loadings f
    DVector<double> scores_;    // computed vector of scores s
  public:
    // constructor
    FPIREM() = default;
    FPIREM(const Model& m) {
      // define internal problem solver required for computation of loadings
      if constexpr(!is_space_time<Model_>::value) // space-only
	solver_ = typename FPIREM_internal_solver<Model_>::type(m.pde());
      else // space-time
	solver_ = typename FPIREM_internal_solver<Model_>::type(m.pde(), m.time_domain());
      // partially initialize internal solver
      solver_.init_pde();
      std::cout << "initializzato" << std::endl;
    }
    
    // setters
    void setLambda(const SVector<model_traits<Model_>::n_lambda>& lambda) {
      std::cout << "set_lambda" << std::endl;

      solver_.setLambda(lambda); // set lambda
      solver_.init_model();      // update solver's matrices to reflect the update on smoothing parameters

      std::cout << "set_lambda" << std::endl;

    }
    void setData(const BlockFrame<double,int>& df) {
      std::cout << "set_data" << std::endl;

      df_ = df;
      // allocate space for solver's data.
      solver_.setData(df);
      solver_.init_regularization();
      solver_.init_sampling(); // init sampling informations (now is known how many data point we expect)
      // reserve space for solution
      loadings_.resize(df_.get<double>(OBSERVATIONS_BLK).rows());
      scores_.resize  (df_.get<double>(OBSERVATIONS_BLK).cols());
      // perform Singular Value Decomposition of Y to initialize PC loadings once
      Eigen::JacobiSVD<DMatrix<double>> svd(df_.get<double>(OBSERVATIONS_BLK).transpose(), Eigen::ComputeThinU|Eigen::ComputeThinV);
      loadings_ = svd.matrixV().col(0);

      std::cout << "set_data" << std::endl;
    }
    // control iterative algorithm parameters
    void setTolerance(double tol) { tol_ = tol; }
    void setMaxIterations(std::size_t max_iter) { max_iter_ = max_iter; }

    void analyze_nan() { solver_.analyze_nan(); }
    void set_nan() { solver_.set_nan(); }
    
    // solves the minimization problem \norm_F{Y - s^T*f}^2 + (s^T*s)*P(f)
    void solve() {
      std::cout << "solve" << std::endl;
      DMatrix<double> Y_ = df_.get<double>(OBSERVATIONS_BLK).transpose(); // n_subject() x n_obs() matrix
      // algorithm initialization
      std::size_t iter_ = 0;
      double Jold = std::numeric_limits<double>::max();
      double Jnew = 0;
      
      // Principal Component computation
      while(std::abs(Jnew - Jold) > tol_ && iter_ < max_iter_){
	// compute score vector s as Y*f/\norm(Y*f)
	scores_ = Y_*loadings_;
	scores_ = scores_/scores_.norm();

	std::cout << "scores: " << scores_.topRows(10) << std::endl;
	
	// compute loadings by solving either an SRPDE or STRPDE (separable regularization) problem
	solver_.data().template insert<double>(OBSERVATIONS_BLK, Y_.transpose()*scores_); // Y^T*s
	//solver_.analyze_nan(); solver_.set_nan();
	solver_.solve();
	// prepare for next iteration
	iter_++;
	Jold = Jnew;
	// update value of discretized functional
	loadings_ = solver_.fitted(); // \Psi*f

	std::cout << "loadings: " << loadings_.topRows(10) << std::endl;
	
	Jnew = (Y_ - scores_*loadings_.transpose()).squaredNorm() + solver_.g().squaredNorm();
      }
      // normalize loadings and unnormalize scores
      double norm = std::sqrt((solver_.f().transpose()*solver_.R0()*solver_.f()).coeff(0,0));
      loadings_ = loadings_/norm; scores_ = scores_*norm;
      return;
    }

    // getters
    const DVector<double>& scores() const { return scores_; }     // vector of computed scores
    const DVector<double>& loadings() const { return loadings_; } // vector of computed loadings
    const BlockFrame<double, int>& data() const { return df_; }   // BlockFrame of data (required by KFoldCV)
    // expose internal solver informations
    const DMatrix<double>& f()   const { return solver_.f(); } // estimated spatial/spatio-temporal field
    const DMatrix<double>& g()   const { return solver_.g(); } // estimated PDE misfit
    const SpMatrix<double>& R0() const { return solver_.R0();} // mass matrix
    typedef typename FPIREM_internal_solver<Model_>::type SmootherType;
    SmootherType& smoother() { return solver_; }
  };
  
}}

#endif // __FPIREM_H__

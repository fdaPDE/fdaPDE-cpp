#ifndef __PROFILING_ESTIMATION_H__
#define __PROFILING_ESTIMATION_H__

#include <Eigen/SVD>
#include "../../core/utils/Symbols.h"
#include "../ModelTraits.h"
using fdaPDE::models::SpaceOnly;
using fdaPDE::models::SpaceTimeSeparable;
#include "../regression/SRPDE.h"
#include "../regression/STRPDE.h"
#include "../../calibration/GCV.h"
using fdaPDE::calibration::GCV;
#include "../../calibration/StochasticEDF.h"
using fdaPDE::calibration::StochasticEDF;

namespace fdaPDE {
namespace models {

  // trait to select model type to use in the internal loop of ProfilingEstimation
  template <typename Model>
  class PE_internal_solver {
  private:
    typedef typename std::decay<Model>::type Model_;
    typedef typename model_traits<Model_>::PDE      PDE;
    typedef typename model_traits<Model_>::sampling sampling;
  public:
    using type = typename std::conditional<
      !is_space_time<Model_>::value,
      SRPDE <PDE, sampling>,
      STRPDE<PDE, fdaPDE::models::SpaceTimeSeparable, sampling, fdaPDE::models::MonolithicSolver>
      >::type;
  };
  
  // finds vectors s,f_n minimizing \norm{X - s*f_n^T}_F^2 + s^T*s P_{\lambda_{\mathcal{D}}, \lambda_T}(f)
  // being P_{\lambda_{\mathcal{D}}, \lambda_T}(f) the penalty term and \norm{}_F the Frobenius norm
  template <typename Model>
  class ProfilingEstimation {
  private:
    typedef typename std::decay<Model>::type Model_;
    typedef typename PE_internal_solver<Model_>::type SolverType;
    SolverType solver_;
    GCV<SolverType, StochasticEDF<SolverType>> gcv_; // gcv index associated to internal solver
    
    // algorithm's parameter
    double tol_ = 1e-6;         // relative tolerance between Jnew and Jold, used as stopping criterion
    std::size_t max_iter_ = 20; // maximum number of allowed iterations
    std::size_t k_ = 0;         // iteration index
    
    // parameters at convergence
    DVector<double> s_;   // estimate of vector s
    DVector<double> f_n_; // estimate of vector f_n (spatial field evaluated at observations' locations)
    DVector<double> f_;   // estimated spatial field
  public:
    // constructor
    ProfilingEstimation(const Model& m, double tol, std::size_t max_iter)
      : tol_(tol), max_iter_(max_iter), gcv_(solver_, 100) {
      // define internal problem solver required for smoothing step
      if constexpr(!is_space_time<Model_>::value) // space-only
	solver_ = typename PE_internal_solver<Model_>::type(m.pde());
      else // space-time
	solver_ = typename PE_internal_solver<Model_>::type(m.pde(), m.time_domain());
      solver_.init_pde();
      solver_.set_spatial_locations(m.locs());
      // set temporal locations here...
      solver_.setData(BlockFrame<double,int>(m.n_locs()));
      // initialize solver
      solver_.init_regularization();
      solver_.init_sampling();
    };

    // somewhere we need to provide a specilization for the case of missing data and a proper dispatcher

    // executes the ProfilingEstimation algorithm given data X and smoothing parameter \lambda
    void compute(const DMatrix<double>& X, const SVector<model_traits<Model_>::n_lambda>& lambda) {
      // solver initialization
      solver_.setLambda(lambda);
      solver_.init_model();
      // reserve space for solution
      f_n_.resize(X.cols()); s_.resize(X.rows());
      
      // initialization of f_ using SVD
      Eigen::JacobiSVD<DMatrix<double>> svd(X, Eigen::ComputeThinU|Eigen::ComputeThinV);
      f_n_ = svd.matrixV().col(0);
      // start iterative procedure
      DMatrix<double> X_ = X; // copy data to avoid side effects on caller state
      double Jold = std::numeric_limits<double>::max(); double Jnew = 1;
      while(!almost_equal(Jnew, Jold, tol_) && k_ < max_iter_){
	// compute score vector s as Y*f/\norm(Y*f)
	s_ = X_*f_n_;
	s_ = s_/s_.norm();
	// compute loadings by solving a proper smoothing problem
	solver_.data().template insert<double>(OBSERVATIONS_BLK, X_.transpose()*s_); // X^T*s
	solver_.solve();
	// prepare for next iteration
	k_++;
	Jold = Jnew;
	// update value of discretized functional
	f_n_ = solver_.fitted(); // \Psi*f
	Jnew = (X_ - s_*f_n_.transpose()).squaredNorm(); // Frobenius norm of reconstruction error
	if constexpr(is_space_only<Model>::value)
	  // for a space only problem we can leverage the following identity
	  // \int_D (Lf-u)^2 = g^\top*R_0*g = f^\top*P*f, being P = R_1^\top*(R_0)^{-1}*R_1
	  Jnew += solver_.g().dot(solver_.R0()*solver_.g());
	else
	  // space-time separable regularization requires to compute the penalty matrix
	  Jnew += solver_.f().dot(solver_.pen()*solver_.f());
      }
      // normalize loadings with respect to L^2 norm
      double norm = std::sqrt(solver_.f().dot(solver_.R0()*solver_.f()));
      f_n_ = f_n_/norm; s_ = s_*norm;
      // store estimated spatial field
      f_ = solver_.f();
      return;
    }
    // getters
    const DVector<double>& f_n() const { return f_n_; }  // vector f_n at convergence
    const DVector<double>& f() const { return f_; }      // estimated spatial field at convergence
    const DVector<double>& s() const { return s_; }      // vector s at convergence
    double gcv() { return gcv_.eval(); }                 // GCV index at convergence
    std::size_t n_iter() const { return k_ - 1; }        // number of iterations

    // setters
    void set_tolerance(double tol) { tol_ = tol; }
    void set_max_iter(std::size_t max_iter) { max_iter_ = max_iter; }
  };

}}

#endif // __PROFILING_ESTIMATION_H__

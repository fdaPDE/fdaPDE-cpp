#ifndef __PROFILING_ESTIMATION_H__
#define __PROFILING_ESTIMATION_H__

#include <Eigen/SVD>
#include "../../core/utils/Symbols.h"
#include "../ModelTraits.h"
using fdaPDE::models::SpaceOnly;
using fdaPDE::models::SpaceTimeSeparable;
#include "../SamplingDesign.h"
using fdaPDE::models::not_nan_corrected;
#include "../regression/SRPDE.h"
#include "../regression/STRPDE.h"
#include "../../calibration/GCV.h"
using fdaPDE::calibration::GCV;
#include "../../calibration/StochasticEDF.h"
using fdaPDE::calibration::StochasticEDF;

namespace fdaPDE {
namespace models {

  // forward declarations
  template <typename Model> class ProfilingEstimationStrategy;
  template <typename Model, typename ImplType> class ProfilingEstimationImpl;
  // tags for possible resolution strategies
  struct complete_data {}; 
  struct missing_data  {}; // specialized implementation for the missing data setting

  // finds vectors s,f_n minimizing \norm{X - s*f_n^T}_F^2 + s^T*s P_{\lambda_{\mathcal{D}}, \lambda_T}(f)
  // being P_{\lambda_{\mathcal{D}}, \lambda_T}(f) the penalty term and \norm{}_F the Frobenius norm
  template <typename Model>
  class ProfilingEstimation { // uses strategy pattern
  private:
    typedef typename std::decay<Model>::type Model_;
    std::unique_ptr<ProfilingEstimationStrategy<Model_>> pe_; // pointer to resolution strategy
  public:
    ProfilingEstimation(const Model& m, double tol, std::size_t max_iter) {
      if(m.has_nan()) // missing data
	pe_ = std::make_unique<
	  ProfilingEstimationImpl<Model_, missing_data>
	  >(m, tol, max_iter);
      else // fallback to complete data setting
	pe_ = std::make_unique<
	  ProfilingEstimationImpl<Model_, complete_data>
	  >(m, tol, max_iter);
    }
    // dynamically dispatch calls to instantiated strategy
    const DVector<double>& f_n() const { return pe_->f_n(); }  // vector f_n at convergence
    const DVector<double>& f() const { return pe_->f(); }      // estimated spatial field at convergence
    const DVector<double>& s() const { return pe_->s(); }      // vector s at convergence
    std::size_t n_iter() const { return pe_->n_iter(); }       // number of iterations
    // setters
    void set_tolerance(double tol) { pe_->set_tolerance(tol); }
    void set_max_iter(std::size_t max_iter) { pe_->set_max_iter(max_iter); }

    // apply profiling estimation algorithm on data matrix X and smoothing vector \lambda
    void compute(const DMatrix<double>& X, const SVector<model_traits<Model_>::n_lambda>& lambda) {
      pe_->compute(X, lambda); };
    double gcv() { return pe_->gcv(); }; // return gcv index at convergence
  };
  
  // base class for profiling estimation resolution strategy
  template <typename Model>
  class ProfilingEstimationStrategy {
  protected:
    typedef typename std::decay<Model>::type Model_;
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
    ProfilingEstimationStrategy() = default;
    ProfilingEstimationStrategy(double tol, std::size_t max_iter) :
      tol_(tol), max_iter_(max_iter) {};
    
    // getters
    const DVector<double>& f_n() const { return f_n_; }  // vector f_n at convergence
    const DVector<double>& f() const { return f_; }      // estimated spatial field at convergence
    const DVector<double>& s() const { return s_; }      // vector s at convergence
    std::size_t n_iter() const { return k_ - 1; }        // number of iterations
    // setters
    void set_tolerance(double tol) { tol_ = tol; }
    void set_max_iter(std::size_t max_iter) { max_iter_ = max_iter; }

    // methods discharged on actual resolution strategies
    virtual void compute(const DMatrix<double>& X, const SVector<model_traits<Model_>::n_lambda>& lambda) = 0;
    virtual double gcv() = 0;
  };

  // complete data setting  
  // trait to select model type to use in the internal loop of ProfilingEstimation
  template <typename Model>
  class PE_internal_solver {
  private:
    typedef typename std::decay<Model>::type Model_;
    typedef typename model_traits<Model_>::PDE PDE;
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
  class ProfilingEstimationImpl<Model, complete_data> : public ProfilingEstimationStrategy<Model> {
  private:
    typedef typename std::decay<Model>::type Model_;
    typedef typename PE_internal_solver<Model_>::type SolverType;
    typedef ProfilingEstimationStrategy<Model> Base;
    SolverType solver_;
    GCV<SolverType, StochasticEDF<SolverType>> gcv_; // gcv index associated to internal solver    
  public:
    using Base::f_n_; // spatial (spatio-temporal) field fitted values
    using Base::s_;   // scores vector
    using Base::f_;   // estimated spatial (spatio-temporal) field
    // constructor
    ProfilingEstimationImpl(const Model& m, double tol, std::size_t max_iter)
      : Base(tol, max_iter), gcv_(solver_, 100) {
      // define internal problem solver required for smoothing step
      if constexpr(!is_space_time<Model_>::value) // space-only
	solver_ = typename PE_internal_solver<Model_>::type(m.pde());
      else{ // space-time
	solver_ = typename PE_internal_solver<Model_>::type(m.pde(), m.time_domain());
	solver_.set_temporal_locations(m.time_locs());
      }
      solver_.set_spatial_locations(m.locs());
      solver_.init_pde();
      solver_.setData(BlockFrame<double,int>(m.n_locs()));
      // initialize solver
      solver_.init_regularization();
      solver_.init_sampling();
    };

    // executes the ProfilingEstimation algorithm given data X and smoothing parameter \lambda, assuming no missing data
    virtual void compute(const DMatrix<double>& X, const SVector<model_traits<Model_>::n_lambda>& lambda) {
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
      while(!almost_equal(Jnew, Jold, this->tol_) && this->k_ < this->max_iter_){
	// compute score vector s as Y*f/\norm(Y*f)
	s_ = X_*f_n_;
	s_ = s_/s_.norm();
	// compute loadings by solving a proper smoothing problem
	solver_.data().template insert<double>(OBSERVATIONS_BLK, X_.transpose()*s_); // X^T*s
	solver_.solve();
	// prepare for next iteration
	this->k_++;
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
      // compute L^2 norm of spatial field
      double L2norm = std::sqrt(solver_.f().dot(solver_.R0()*solver_.f()));;
      f_ = solver_.f()/L2norm;
      if constexpr(is_sampling_pointwise_at_mesh<Model_>::value){
	// store normalized f_n with respcet to L^2 norm
	f_n_ = f_; s_ = s_*L2norm;
      }else{ // use euclidean norm if L^2 norm of f_n vector cannot be computed 
	L2norm = f_n_.norm();
	f_n_ = f_n_/L2norm; s_ = s_*L2norm;
      } 
      return;
    }
    
    // getters
    virtual double gcv() { return gcv_.eval(); } // GCV index at convergence
  };

  // missing data setting
  // functional to minimize is different
  template <typename Model>
  class ProfilingEstimationImpl<Model, missing_data> : public ProfilingEstimationStrategy<Model> {
  private:
    typedef typename std::decay<Model>::type Model_;
    typedef ProfilingEstimationStrategy<Model> Base;
    const Model_& m_;
    std::vector<SpMatrix<double>> psi_vect_; // store \Psi matrix of i-th subject
  public:
    using Base::f_n_; // spatial (spatio-temporal) field fitted values
    using Base::s_;   // scores vector
    using Base::f_;   // estimated spatial (spatio-temporal) field
    // constructor
    ProfilingEstimationImpl(const Model& m, double tol, std::size_t max_iter)
      : Base(tol, max_iter), m_(m) {};
    
    // executes the ProfilingEstimation algorithm given data X and smoothing parameter \lambda
    virtual void compute(const DMatrix<double>& X, const SVector<model_traits<Model_>::n_lambda>& lambda) {
      if(psi_vect_.size() == 0){
	// compute \Psi matrix for each subject
	psi_vect_.resize(X.rows());
	for(std::size_t i = 0; i < X.rows(); ++i){
	  psi_vect_[i].resize(m_.n_basis(), m_.n_basis());
	  psi_vect_[i] = m_.Psi();
	  // zero \Psi_i rows depending on missingness pattern of i-th subject
	  for(std::size_t j = 0; j < X.cols(); ++j){
	    if(std::isnan(X(i,j))) psi_vect_[i].row(j) *= 0;
	  }
	}
      }

      DMatrix<double> X_ = X; // copy data to avoid side effects on caller state      
      // zero NaN entries in data matrix
      for(std::size_t i = 0; i < X_.rows(); ++i){
	for(std::size_t j = 0; j < X_.cols(); ++j)
	  if(std::isnan(X(i,j))) X_(i,j) = 0;
      }
      
      // solver initialization
      std::size_t N = m_.n_basis();
      SparseBlockMatrix<double,2,2>
	A(SpMatrix<double>(N,N),  lambda[0]*m_.R1().transpose(),
	  lambda[0]*m_.R1(),      lambda[0]*m_.R0()            );

      DVector<double> b_{};  // right hand side of problem's linear system (1 x 2N vector)
      b_.resize(2*N);
      b_.block(N,0, N,1) = DMatrix<double>::Zero(N, 1);
      
      // reserve space for solution
      f_n_.resize(X.cols()); s_.resize(X.rows());
      
      // initialization of f_ using SVD
      Eigen::JacobiSVD<DMatrix<double>> svd(X_, Eigen::ComputeThinU|Eigen::ComputeThinV);
      f_n_ = svd.matrixV().col(0);
      // start iterative procedure
      double Jold = std::numeric_limits<double>::max(); double Jnew = 1;
      while(!almost_equal(Jnew, Jold, this->tol_) && this->k_ < this->max_iter_){
	// compute score vector s as Y*f/\norm(Y*f)
	s_ = X_*f_n_;
	s_ = s_/s_.norm();

	// update north-west block of matrix A
	SpMatrix<double> L;
	L.resize(N,N);
	for(std::size_t i = 0; i < X.rows(); ++i)
	  L += (s_[i]*s_[i])*(psi_vect_[i].transpose()*psi_vect_[i]);
	
	SparseBlockMatrix<double,2,2>
	  B(-L,      SpMatrix<double>(N,N),
	    SpMatrix<double>(N,N), SpMatrix<double>(N,N));

	SpMatrix<double> AA = L.transpose()*L;
	
	SpMatrix<double> C = A.derived() + B.derived();
	// update rhs of linear system
	b_.block(0,0, N,1) = -m_.PsiTD()*X_.transpose()*s_;
	fdaPDE::SparseLU<SpMatrix<double>> invC_; // factorization of matrix C
	invC_.compute(C);
	// solve smoothing problem
	DVector<double> solution = invC_.solve(b_);
	f_ = solution.topRows(N);

	// prepare for next iteration
	this->k_++;
	Jold = Jnew;
	// update value of discretized functional
	f_n_ = m_.Psi(not_nan_corrected())*f_; // \Psi*f
	Jnew = (X_ - s_*f_n_.transpose()).squaredNorm(); // Frobenius norm of reconstruction error
	if constexpr(is_space_only<Model>::value)
	  // for a space only problem we can leverage the following identity
	  // \int_D (Lf-u)^2 = g^\top*R_0*g = f^\top*P*f, being P = R_1^\top*(R_0)^{-1}*R_1
	  Jnew += solution.bottomRows(N).dot(m_.R0()*solution.bottomRows(N));
	else
	  // space-time separable regularization requires to compute the penalty matrix
	  Jnew += f_.dot(m_.pen()*f_);
      }
      // normalize loadings with respect to L^2 norm
      double norm = std::sqrt(f_.dot(m_.R0()*f_));
      f_n_ = f_n_/norm; s_ = s_*norm;
      return;
    }
    // getters
    virtual double gcv() { return 0; } // GCV index at convergence
  };
  
}}

#endif // __PROFILING_ESTIMATION_H__

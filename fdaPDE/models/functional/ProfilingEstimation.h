#ifndef __PROFILING_ESTIMATION_H__
#define __PROFILING_ESTIMATION_H__

#include <Eigen/SVD>
#include "../../core/utils/Symbols.h"
#include "../ModelTraits.h"
using fdaPDE::models::SpaceOnly;
using fdaPDE::models::SpaceTimeSeparable;
#include "../SamplingDesign.h"
using fdaPDE::models::not_nan;
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
    const DVector<double>& f_n() const { return pe_->f_n(); } // vector f_n at convergence
    const DVector<double>& f() const { return pe_->f(); }     // estimated field at convergence
    const DVector<double>& g() const { return pe_->g(); }     // PDE misfit at convergence
    const DVector<double>& s() const { return pe_->s(); }     // vector s at convergence
    std::size_t n_iter() const { return pe_->n_iter(); }      // number of iterations
    double f_norm() const { return pe_->f_norm(); }           // L^2 norm of esimated field
    double f_n_norm() const { return pe_->f_n_norm(); }       // f_n vector norm
    // setters
    void set_tolerance(double tol) { pe_->set_tolerance(tol); }
    void set_max_iter(std::size_t max_iter) { pe_->set_max_iter(max_iter); }

    // apply profiling estimation algorithm on data matrix X and smoothing vector \lambda
    void compute(const BlockFrame<double,int>& df, const SVector<model_traits<Model_>::n_lambda>& lambda) { pe_->compute(df, lambda); };
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
    DVector<double> s_;   
    DVector<double> f_n_; 
    DVector<double> f_; // estimated spatial(spatio-temporal) field at convergence
    DVector<double> g_; // PDE misfit at convergence
    double f_norm_;     // L^2 norm of estimated field at converegence
    double f_n_norm_;   // f_n norm (use euclidean norm if is_sampling_pointwise_at_mesh<Model> evaluates false)
  public:
    // constructor
    ProfilingEstimationStrategy() = default;
    ProfilingEstimationStrategy(double tol, std::size_t max_iter) :
      tol_(tol), max_iter_(max_iter) {};
    
    // getters
    const DVector<double>& f_n() const { return f_n_; } 
    const DVector<double>& f() const { return f_; } 
    const DVector<double>& g() const { return g_; } 
    const DVector<double>& s() const { return s_; } 
    std::size_t n_iter() const { return k_ - 1; }   
    double f_norm() const { return f_norm_; }
    double f_n_norm() const { return f_n_norm_; }
    // setters
    void set_tolerance(double tol) { tol_ = tol; }
    void set_max_iter(std::size_t max_iter) { max_iter_ = max_iter; }

    // methods discharged on actual resolution strategies
    virtual void compute(const BlockFrame<double,int>& df, const SVector<model_traits<Model_>::n_lambda>& lambda) = 0;
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
      solver_.setData(BlockFrame<double,int>(m.n_locs()));
      solver_.init_pde();
      // initialize solver
      solver_.init_regularization();
      solver_.init_sampling();
    };

    // executes the ProfilingEstimation algorithm given data X and smoothing parameter \lambda, assuming no missing data
    virtual void compute(const BlockFrame<double, int>& df, const SVector<model_traits<Model_>::n_lambda>& lambda) {
      // solver initialization
      solver_.setLambda(lambda);
      solver_.init_model();
      const DMatrix<double>& X_ = df.template get<double>(OBSERVATIONS_BLK);
      // reserve space for solution
      f_n_.resize(X_.cols()); s_.resize(X_.rows());

      // initialization of f_ using SVD
      Eigen::JacobiSVD<DMatrix<double>> svd(X_, Eigen::ComputeThinU|Eigen::ComputeThinV);
      f_n_ = svd.matrixV().col(0);
      // start iterative procedure
	
      double Jold = std::numeric_limits<double>::max(); double Jnew = 1;
      this->k_ = 0; // reset iteration counter
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
	  Jnew += lambda[0]*solver_.g().dot(solver_.R0()*solver_.g());
	else
	  // space-time separable regularization requires to compute the penalty matrix
	  Jnew += solver_.f().dot(solver_.pen()*solver_.f());
      }
      // store results
      f_ = solver_.f(); // estimated field at convergence
      this->g_ = solver_.g(); // PDE misfit at convergence
      this->f_norm_ = std::sqrt(f_.dot(solver_.R0()*f_)); // L^2 norm of estimated field
      // compute norm of loadings vector
      if constexpr(is_sampling_pointwise_at_mesh<Model_>::value) this->f_n_norm_ = this->f_norm_;
      else this->f_n_norm_ = f_n_.norm(); // use euclidean norm if L^2 norm of f_n vector cannot be computed 
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
    Model_& m_;

    SpMatrix<double> P_; // Pt \kron R0, [Pt_]_{ij} = \int_{[0,T]} (\phi_i)_tt*(\phi_j)_tt, R0 FEM mass matrix
  public:
    using Base::f_n_; // spatial (spatio-temporal) field fitted values
    using Base::s_;   // scores vector
    using Base::f_;   // estimated spatial (spatio-temporal) field
    using Base::g_;   // PDE misfit
    // constructor
    ProfilingEstimationImpl(Model& m, double tol, std::size_t max_iter)
      : Base(tol, max_iter), m_(m) {
      // initialize P_ if regularization is space-time separable
      if constexpr(is_space_time<Model>::value) P_ = Kronecker(m_.Pt(), m_.pde().R0());
    };
    
    // executes the ProfilingEstimation algorithm given data X and smoothing parameter \lambda
    virtual void compute(const BlockFrame<double, int>& df, const SVector<model_traits<Model_>::n_lambda>& lambda) {
      // extract missingness pattern
      const DMatrix<double>& X_nan = df.template get<double>(OBSERVATIONS_BLK);
      auto nan_pattern = X_nan.array().isNaN(); 
      DMatrix<double> X_ = nan_pattern.select(0, X_nan); // set NaN to zero
      
      // solver initialization
      std::size_t N = m_.Psi(not_nan()).cols(); // number of basis, but for space-time (separable) problems should be something else....
      SparseBlockMatrix<double,2,2>
	A_(SpMatrix<double>(N,N),  lambda[0]*m_.R1().transpose(),
	   lambda[0]*m_.R1(),      lambda[0]*m_.R0()            );
      Eigen::SparseLU<SpMatrix<double>> invA_;

      DVector<double> b_;  // right hand side of problem's linear system (1 x 2N vector)
      b_.resize(2*N);
      b_.block(N,0, N,1) = DMatrix<double>::Zero(N, 1);

      // initialization of f_ using SVD
      f_n_.resize(X_.cols()); s_.resize(X_.rows());      
      Eigen::JacobiSVD<DMatrix<double>> svd(X_, Eigen::ComputeThinU|Eigen::ComputeThinV);
      f_n_ = svd.matrixV().col(0);
      
      // start iterative procedure
      double Jold = std::numeric_limits<double>::max(); double Jnew = 1;
      this->k_ = 0; // reset iteration counter
      while(!almost_equal(Jnew, Jold, this->tol_) && this->k_ < this->max_iter_){
	// compute score vector s as X*f/\norm(X*f)
	s_ = X_*f_n_;
	s_ = s_/s_.norm();
	// Assembly of matrix [L]_{ij} = \sum_{m=1}^M \sum_{n \in O_m}(s_m^2*\psi_i(p_n)*\psi_j(p_n))
	// being M : number of considered statistical units, O_l : observation locations' indexes for l-th unit
	SpMatrix<double> L; L.resize(N,N);
	for(std::size_t i = 0; i < s_.rows(); ++i){
	  // data might be possible shuffled (see for instance KFoldCV), recover correct stat unit index
	  int stat_unit = df.template get<int>(INDEXES_BLK)(i,0);
	  L += (s_[i]*s_[i])*(m_.PsiTD(stat_unit)*m_.Psi(stat_unit));
	}

	// set north-west block of system matrix A_ and factorize
	if constexpr(is_space_only<Model>::value) A_.block(0,0) = -L;
	else A_.block(0,0) = -L - lambda[1]*P_;
	invA_.compute(A_);
	// update rhs of linear system
	b_.block(0,0, N,1) = -m_.PsiTD(not_nan())*X_.transpose()*s_;
	// solve linear system and store results
	DVector<double> solution = invA_.solve(b_);
	f_ = solution.topRows(N); g_ = solution.bottomRows(N); // estimated field and PDE misfit
	f_n_ = m_.Psi(not_nan())*f_; // \Psi*f
	
	// prepare for next iteration
	this->k_++;
	Jold = Jnew;
	// Frobenius norm of reconstruction error
	double Jnew = nan_pattern.select(0, X_nan - s_*f_n_.transpose()).squaredNorm();
	if constexpr(is_space_only<Model>::value)
	  // for a space only problem we can leverage the following identity
	  // \int_D (Lf-u)^2 = g^\top*R_0*g = f^\top*P*f, being P = R_1^\top*(R_0)^{-1}*R_1
	  Jnew += lambda[0]*g_.dot(m_.R0()*g_);
	else
	  // space-time separable regularization requires to compute the penalty matrix
	  Jnew += f_.dot(m_.pen()*f_);
      }
      // store results
      this->f_norm_ = std::sqrt(f_.dot(m_.R0()*f_)); // L^2 norm of estimated field
      // compute norm of loadings vector
      if constexpr(is_sampling_pointwise_at_mesh<Model_>::value) this->f_n_norm_ = this->f_norm_;
      else this->f_n_norm_ = f_n_.norm(); // use euclidean norm if L^2 norm of f_n vector cannot be computed 
    }
    
    // **currently not working**
    // GCV is manually implemented here, should provide a general implementation of GCV model-independent (at least not of GCV
    // but of trace computation algorithms sure... )
    virtual double gcv() { return 0.0; };
    //     // compute trace of smoothing matrix using stochastic approximation
    //     // std::size_t n = m_.n_basis(); // number of basis functions
    //     // if(!init_){
    //     // 	// compute sample from Rademacher distribution
    //     // 	std::size_t seed = std::random_device()();
    //     // 	std::default_random_engine rng(seed);
    //     // 	std::bernoulli_distribution Be(0.5); // bernulli distribution with parameter p = 0.5
    //     // 	Us_.resize(m_.n_locs(), r_); // preallocate memory for matrix Us
    //     // 	// fill matrix
    //     // 	for(std::size_t i = 0; i < m_.n_locs(); ++i){
    //     // 	  for(std::size_t j = 0; j < r_; ++j){
    //     // 	    if(Be(rng)) Us_(i,j) =  1.0;
    //     // 	    else        Us_(i,j) = -1.0;
    //     // 	  }
    //     // 	}
    //     // 	// prepare matrix Bs_
    //     // 	Bs_ = DMatrix<double>::Zero(2*n, r_);
    //     // 	Bs_.topRows(n) = -m_.Psi(not_nan_corrected()).transpose()*Us_;
    //     // 	// prepare matrix Y
    //     // 	Y_ = Us_.transpose()*m_.Psi(not_nan_corrected());
    //     // 	init_ = true; // never reinitialize again
    //     // }
      
    //     // DMatrix<double> sol; // room for problem solution
    //     // sol = invC_.solve(Bs_);
    //     // // compute approximated Tr[S] using monte carlo mean
    //     // double MCmean = 0;
    //     // for(std::size_t i = 0; i < r_; ++i)
    //     // 	MCmean += Y_.row(i).dot(sol.col(i).head(n));
      
    //     // double trS = MCmean/r_;

    //     fdaPDE::SparseLU<SpMatrix<double>> invR0_;
    //     invR0_.compute(m_.R0());      
    //     DMatrix<double> T = L_ + lambda_[0]*( m_.R1().transpose() * invR0_.solve(m_.R1()) );

    //     Eigen::PartialPivLU<DMatrix<double>> invT_;
    //     invT_ = T.partialPivLu();

    //     // if locations are mesh nodes S = T^{-1}
    //     DMatrix<double> E = m_.Psi(not_nan()).transpose();
    //     DMatrix<double> S = m_.Psi(not_nan())*invT_.solve(E);

    //     double trS = S.trace();
      
    //     std::size_t nn = y_.rows();
    //     double dor = nn - trS;     // residual degrees of freedom
      
    //     // return gcv at point
    //     return (nn/std::pow(dor, 2))*( (f_n_ - y_).squaredNorm() ) ;
    //   } // GCV index at convergence
  };
  
}}

#endif // __PROFILING_ESTIMATION_H__

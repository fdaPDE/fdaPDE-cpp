// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __PROFILING_ESTIMATION_H__
#define __PROFILING_ESTIMATION_H__

#include <Eigen/SVD>
#include <fdaPDE/utils.h>

#include "../regression/gcv.h"
#include "../regression/stochastic_edf.h"
#include "../model_traits.h"
#include "../regression/srpde.h"
#include "../regression/strpde.h"
#include "../sampling_design.h"
using fdapde::models::GCV;
using fdapde::models::StochasticEDF;

namespace fdapde {
namespace models {

// forward declarations
template <typename Model> class ProfilingEstimationStrategy;
template <typename Model, typename ImplType> class ProfilingEstimationImpl;
// tags for possible resolution strategies
struct complete_data { };
struct missing_data { };   // specialized implementation for the missing data setting

// finds vectors s,f_n minimizing \norm{X - s*f_n^T}_F^2 + s^T*s P_{\lambda_{\mathcal{D}}, \lambda_T}(f)
// being P_{\lambda_{\mathcal{D}}, \lambda_T}(f) the penalty term and \norm{}_F the Frobenius norm
template <typename Model> class ProfilingEstimation {   // uses strategy pattern
   private:
    typedef typename std::decay<Model>::type Model_;
    std::unique_ptr<ProfilingEstimationStrategy<Model_>> pe_;   // pointer to resolution strategy
   public:
    // constructor
    ProfilingEstimation() = default;
    ProfilingEstimation(Model& m, double tol, std::size_t max_iter, int seed) {
        if (m.has_nan())   // missing data
            pe_ = std::make_unique<ProfilingEstimationImpl<Model_, missing_data>> (m, tol, max_iter, seed);
        else   // fallback to complete data setting
            pe_ = std::make_unique<ProfilingEstimationImpl<Model_, complete_data>>(m, tol, max_iter, seed);
    }
    // dynamically dispatch calls to instantiated strategy
    const DVector<double>& f_n() const { return pe_->f_n(); }   // vector f_n at convergence
    const DVector<double>& f() const { return pe_->f(); }       // estimated field at convergence
    const DVector<double>& g() const { return pe_->g(); }       // PDE misfit at convergence
    const DVector<double>& s() const { return pe_->s(); }       // vector s at convergence
    std::size_t n_iter() const { return pe_->n_iter(); }        // number of iterations
    double f_norm() const { return pe_->f_norm(); }             // L^2 norm of esimated field
    double f_n_norm() const { return pe_->f_n_norm(); }         // f_n vector norm
    // setters
    void set_tolerance(double tol) { pe_->set_tolerance(tol); }
    void set_max_iter(std::size_t max_iter) { pe_->set_max_iter(max_iter); }
  
    // apply profiling estimation algorithm on data matrix X and smoothing vector \lambda
    void compute(const BlockFrame<double, int>& df, const SVector<model_traits<Model_>::n_lambda>& lambda) {
        pe_->compute(df, lambda);
    };
    double gcv() { return pe_->gcv(); };   // return gcv index at convergence
};

// trait to select model type to use in the internal loop of ProfilingEstimation
template <typename Model> class PE_internal_solver {
   private:
    typedef typename std::decay<Model>::type Model_;
    typedef typename model_traits<Model_>::PDE PDE;
    typedef typename model_traits<Model_>::sampling sampling;
   public:
    using type = typename std::conditional<
      !is_space_time<Model_>::value, SRPDE<PDE, sampling>,
      STRPDE<PDE, SpaceTimeSeparable, sampling, MonolithicSolver>>::type;
};

// base class for profiling estimation resolution strategy
template <typename Model> class ProfilingEstimationStrategy {
   protected:
    typedef typename std::decay<Model>::type Model_;
    typedef typename PE_internal_solver<Model_>::type SolverType;
    typedef ProfilingEstimationStrategy<Model> Base;
    SolverType solver_;
    GCV<SolverType, StochasticEDF<SolverType>> gcv_;   // gcv index associated to internal solver

    // algorithm's parameter
    double tol_ = 1e-6;           // relative tolerance between Jnew and Jold, used as stopping criterion
    std::size_t max_iter_ = 20;   // maximum number of allowed iterations
    std::size_t k_ = 0;           // iteration index

    // parameters at convergence
    DVector<double> s_;
    DVector<double> f_n_;
    DVector<double> f_;   // estimated spatial(spatio-temporal) field at convergence
    DVector<double> g_;   // PDE misfit at convergence
    double f_norm_;       // L^2 norm of estimated field at converegence
    double f_n_norm_;     // f_n norm (use euclidean norm if is_sampling_pointwise_at_mesh<Model> evaluates false)
   public:
    // constructor
    ProfilingEstimationStrategy() = default;
    ProfilingEstimationStrategy(const Model& m, double tol, std::size_t max_iter, int seed) :
        tol_(tol), max_iter_(max_iter), gcv_(solver_, 100, seed) {
        // define internal problem solver required for smoothing step
        if constexpr (!is_space_time<Model_>::value)   // space-only
            solver_ = typename PE_internal_solver<Model_>::type(m.pde());
        else {   // space-time
            solver_ = typename PE_internal_solver<Model_>::type(m.pde(), m.time_domain());
            solver_.set_temporal_locations(m.time_locs());
        }
        solver_.set_spatial_locations(m.locs());
        solver_.set_data(BlockFrame<double, int>(m.n_locs()));
        solver_.init_pde();
        // initialize solver
        solver_.init_regularization();
        solver_.init_sampling();
    };

    // getters
    const DVector<double>& f_n() const { return f_n_; }
    const DVector<double>& f() const { return f_; }
    const DVector<double>& g() const { return g_; }
    const DVector<double>& s() const { return s_; }
    std::size_t n_iter() const { return k_ - 1; }
    double f_norm() const { return f_norm_; }
    double f_n_norm() const { return f_n_norm_; }
    double gcv() { return gcv_.eval(); }   // GCV index at convergence
    // setters
    void set_tolerance(double tol) { tol_ = tol; }
    void set_max_iter(std::size_t max_iter) { max_iter_ = max_iter; }
    // methods implemented by resolution schemes
    virtual void compute(const BlockFrame<double, int>& df, const SVector<model_traits<Model_>::n_lambda>& lambda) = 0;
};

// finds vectors s,f_n minimizing \norm{X - s*f_n^T}_F^2 + s^T*s P_{\lambda_{\mathcal{D}}, \lambda_T}(f)
// being P_{\lambda_{\mathcal{D}}, \lambda_T}(f) the penalty term and \norm{}_F the Frobenius norm
template <typename Model>
struct ProfilingEstimationImpl<Model, complete_data> : public ProfilingEstimationStrategy<Model> {
    typedef typename std::decay<Model>::type Model_;
    typedef ProfilingEstimationStrategy<Model> Base;
    using Base::f_;        // estimated spatial (spatio-temporal) field
    using Base::f_n_;      // spatial (spatio-temporal) field fitted values
    using Base::s_;        // scores vector
    using Base::solver_;   // internal solver used in smoothing step
    // constructor
    ProfilingEstimationImpl(Model& m, double tol, std::size_t max_iter, int seed) : Base(m, tol, max_iter, seed) {};

    // executes the ProfilingEstimation algorithm given data X and smoothing parameter \lambda, assuming no missing data
    virtual void compute(const BlockFrame<double, int>& df, const SVector<model_traits<Model_>::n_lambda>& lambda) {
        // solver initialization
        solver_.set_lambda(lambda);
        solver_.init_model();
        const DMatrix<double>& X_ = df.template get<double>(OBSERVATIONS_BLK);
        // reserve space for solution
        f_n_.resize(X_.cols());
        s_.resize(X_.rows());

        // initialization of f_ using SVD
        Eigen::JacobiSVD<DMatrix<double>> svd(X_, Eigen::ComputeThinU | Eigen::ComputeThinV);
        f_n_ = svd.matrixV().col(0);

        // start iterative procedure
        double Jold = std::numeric_limits<double>::max();
        double Jnew = 1;
        this->k_ = 0;   // reset iteration counter
        while (!almost_equal(Jnew, Jold, this->tol_) && this->k_ < this->max_iter_) {
            // compute score vector s as Y*f/\norm(Y*f)
            s_ = X_ * f_n_;
            s_ = s_ / s_.norm();
            // compute loadings by solving a proper smoothing problem
            solver_.data().template insert<double>(OBSERVATIONS_BLK, X_.transpose() * s_);   // X^T*s
            solver_.solve();
            // prepare for next iteration
            this->k_++;
            Jold = Jnew;
            // update value of discretized functional
            f_n_ = solver_.fitted();                             // \Psi*f
            Jnew = (X_ - s_ * f_n_.transpose()).squaredNorm();   // Frobenius norm of reconstruction error
            if constexpr (is_space_only<Model>::value)
                // for a space only problem we can leverage the following identity
                // \int_D (Lf-u)^2 = g^\top*R_0*g = f^\top*P*f, being P = R_1^\top*(R_0)^{-1}*R_1
                Jnew += lambda[0] * solver_.g().dot(solver_.R0() * solver_.g());
            else
                // space-time separable regularization requires to compute the penalty matrix
                Jnew += solver_.f().dot(solver_.pen() * solver_.f());
        }
        // store results
        f_ = solver_.f();                                       // estimated field at convergence
        this->g_ = solver_.g();                                 // PDE misfit at convergence
        this->f_norm_ = std::sqrt(f_.dot(solver_.R0() * f_));   // L^2 norm of estimated field
        // compute norm of loadings vector
        if constexpr (is_sampling_pointwise_at_mesh<Model_>::value)
            this->f_n_norm_ = this->f_norm_;
        else
            this->f_n_norm_ = f_n_.norm();   // use euclidean norm if L^2 norm of f_n vector cannot be computed
    }
};

// finds vectors s,f_n minimizing \norm{X - s*f_n^T}_F^2 + s^T*s P_{\lambda_{\mathcal{D}}, \lambda_T}(f)
// being P_{\lambda_{\mathcal{D}}, \lambda_T}(f) the penalty term and \norm{}_F the Frobenius norm
// X can have NaN. In this case a weighted least square problem has to be solved at each iteration
template <typename Model>
struct ProfilingEstimationImpl<Model, missing_data> : public ProfilingEstimationStrategy<Model> {
    typedef typename std::decay<Model>::type Model_;
    typedef ProfilingEstimationStrategy<Model> Base;
    using Base::f_;        // estimated spatial (spatio-temporal) field
    using Base::f_n_;      // spatial (spatio-temporal) field fitted values
    using Base::s_;        // scores vector
    using Base::solver_;   // internal solver used in smoothing step
    // constructor
    ProfilingEstimationImpl(Model& m, double tol, std::size_t max_iter, int seed) : Base(m, tol, max_iter, seed) {};

    // executes the ProfilingEstimation algorithm given data X and smoothing parameter \lambda, X can have NaN
    virtual void compute(const BlockFrame<double, int>& df, const SVector<model_traits<Model_>::n_lambda>& lambda) {
        // extract missingness pattern
        const DMatrix<double>& X_nan = df.template get<double>(OBSERVATIONS_BLK);
        auto nan_pattern = X_nan.array().isNaN();
        DMatrix<double> X_ = nan_pattern.select(0, X_nan);   // set NaN to zero
        solver_.set_lambda(lambda);
        solver_.init_model();
        // reserve space for solution
        f_n_.resize(X_.cols());
        s_.resize(X_.rows());
        // initialization of f_ using SVD
        Eigen::JacobiSVD<DMatrix<double>> svd(X_, Eigen::ComputeThinU | Eigen::ComputeThinV);
        f_n_ = svd.matrixV().col(0);

        // start iterative procedure
        double Jold = std::numeric_limits<double>::max();
        double Jnew = 1;
        this->k_ = 0;   // reset iteration counter
        while (!almost_equal(Jnew, Jold, this->tol_) && this->k_ < this->max_iter_) {
            // compute score vector s as Y*f/\norm(Y*f)
            s_ = X_ * f_n_;
            s_ = s_ / s_.norm();
            // compute weights matrix W_ = diag(\sum_k s_k^2 o_1^k, ..., \sum_k s_k^2 o_n^k).
            // o_j^k = 1 \iff data is observed at location j for unit k
            DMatrix<double> W_(X_.cols(), 1);
            for (std::size_t i = 0; i < W_.rows(); ++i) {
                W_(i, 0) = nan_pattern.col(i).select(0, s_.array().square()).sum();   // \sum_k s_k^2 o_i^k
            }
            // adjust observation vector as y_i = (X^T*s)_i/w_i. if w_i = 0, set y_i = 0
            DMatrix<double> y_ = (W_.array() == 0).select(0, (X_.transpose() * s_).array() / W_.array());
            // compute loadings by solving a proper **weighted** smoothing problem
            solver_.data().template insert<double>(OBSERVATIONS_BLK, y_);
            solver_.data().template insert<double>(WEIGHTS_BLK, W_);
            solver_.update_data();
            solver_.update_to_weights();   // update non-parametric matrix to change in the weights matrix
            solver_.solve();
            // prepare for next iteration
            this->k_++;
            Jold = Jnew;
            // update value of discretized functional
            f_n_ = solver_.fitted();                                                     // \Psi*f
            Jnew = nan_pattern.select(0, X_nan - s_ * f_n_.transpose()).squaredNorm();   // norm of reconstruction error
            if constexpr (is_space_only<Model>::value)
                // for a space only problem we can leverage the following identity
                // \int_D (Lf-u)^2 = g^\top*R_0*g = f^\top*P*f, being P = R_1^\top*(R_0)^{-1}*R_1
                Jnew += lambda[0] * solver_.g().dot(solver_.R0() * solver_.g());
            else
                // space-time separable regularization requires to compute the penalty matrix
                Jnew += solver_.f().dot(solver_.pen() * solver_.f());
        }
        f_ = solver_.f();                                       // estimated field at convergence
        this->g_ = solver_.g();                                 // PDE misfit at convergence
        this->f_norm_ = std::sqrt(f_.dot(solver_.R0() * f_));   // L^2 norm of estimated field
        // compute norm of loadings vector
        if constexpr (is_sampling_pointwise_at_mesh<Model_>::value)
            this->f_n_norm_ = this->f_norm_;
        else
            this->f_n_norm_ = f_n_.norm();   // use euclidean norm if L^2 norm of f_n vector cannot be computed
    }
};

}   // namespace models
}   // namespace fdapde

#endif   // __PROFILING_ESTIMATION_H__

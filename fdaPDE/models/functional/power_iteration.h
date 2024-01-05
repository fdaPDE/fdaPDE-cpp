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

#ifndef __POWER_ITERATION_H__
#define __POWER_ITERATION_H__

#include <fdaPDE/utils.h>
#include "../model_traits.h"
#include "../regression/gcv.h"
#include "../regression/srpde.h"
#include "../regression/stochastic_edf.h"
#include "../regression/strpde.h"
using fdapde::models::GCV;

namespace fdapde {
namespace models {

// finds vectors s,f_n minimizing \norm{X - s*f_n^T}_F^2 + s^T*s P_{\lambda_{\mathcal{D}}, \lambda_T}(f)
// being P_{\lambda_{\mathcal{D}}, \lambda_T}(f) the penalty term and \norm{}_F the Frobenius norm
template <typename Model_> class PowerIteration {
   private:
    using Model = typename std::decay<Model_>::type;
    using SolverType = typename std::conditional<
      is_space_only<Model>::value, SRPDE, STRPDE<typename Model::RegularizationType, fdapde::monolithic>>::type;
    static constexpr int n_lambda = Model::n_lambda;
    Model* m_;
    GCV gcv_;   // GCV functor
    // algorithm's parameters
    double tolerance_ = 1e-6;     // treshold on |Jnew - Jold| used as stopping criterion
    std::size_t max_iter_ = 20;   // maximum number of iterations before forced stop
    std::size_t k_ = 0;           // iteration index
    SolverType solver_;           // internal solver
    int seed_;

    DVector<double> s_;    // estimated score vector
    DVector<double> fn_;   // field evaluation at data location (\Psi*f)
    DVector<double> f_;    // field basis expansion at convergence
    DVector<double> g_;    // PDE misfit at convergence
    double f_norm_;        // L^2 norm of estimated field at converegence
    double fn_norm_;       // euclidean norm of field evaluation at data locations
   public:
    // constructors
    PowerIteration() = default;
    PowerIteration(Model* m, double tolerance, std::size_t max_iter, int seed) :
        m_(m), tolerance_(tolerance), max_iter_(max_iter),
        seed_((seed == fdapde::random_seed) ? std::random_device()() : seed) {};
    PowerIteration(Model* m, double tolerance, std::size_t max_iter) :
        PowerIteration(m, tolerance, max_iter, fdapde::random_seed) {};

    // initializes internal smoothing solver
    void init() {
        // define internal problem solver
        if constexpr (!is_space_time<Model>::value)   // space-only
            solver_ = SolverType(m_->pde(), m_->sampling());
        else {   // space-time
            if constexpr (is_space_time_parabolic<Model>::value) {
                solver_ = SolverType(m_->pde(), m_->sampling(), m_->time_domain());
                solver_.set_initial_condition(m_->s(), false);
            }
            if constexpr (is_space_time_separable<Model>::value) {
                solver_ = SolverType(m_->pde(), m_->time_pde(), m_->sampling());
                solver_.set_temporal_locations(m_->time_locs());
            }
        }
        solver_.set_spatial_locations(m_->locs());
        // initialize GCV functor
        gcv_ = solver_.template gcv<StochasticEDF>(100, seed_);
    }

    // executes the power itration algorithm on data X and smoothing parameter \lambda, starting from f0
    void compute(const DMatrix<double>& X, const SVector<n_lambda>& lambda, const DVector<double>& f0) {
        // initialization
        fn_.resize(X.cols()); s_.resize(X.rows());
	k_ = 0;   // reset iteration counter
        solver_.set_lambda(lambda);
        solver_.init();
        double Jold = std::numeric_limits<double>::max();
        double Jnew = 1;
        fn_ = f0;   // set starting point
        while (!almost_equal(Jnew, Jold, tolerance_) && k_ < max_iter_) {
            // compute score vector s as \frac{X*fn}{\norm(X*fn)}
            s_ = X * fn_;
            s_ = s_ / s_.norm();
            // compute loadings by solving a proper smoothing problem
            solver_.data().template insert<double>(OBSERVATIONS_BLK, X.transpose() * s_);   // X^\top*s
            solver_.solve();
            // prepare for next iteration
            k_++;
            Jold = Jnew;
            // update value of discretized functional
            fn_ = solver_.fitted();                            // \Psi*f
            Jnew = (X - s_ * fn_.transpose()).squaredNorm();   // frobenius norm of reconstruction error
            if constexpr (is_space_only<Model>::value) {
                // \int_D (Lf-u)^2 = g^\top*R_0*g = f^\top*P*f, being P = R_1^\top*(R_0)^{-1}*R_1
                Jnew += lambda[0] * solver_.g().dot(solver_.R0() * solver_.g());
	    } else {
	        // f^\top*P*f = g^\top*R0*g + f^\top*P_T*f
                Jnew += lambda[0] * solver_.g().dot(solver_.R0() * solver_.g()) +   // g^\top*R0*g
                        lambda[1] * solver_.f().dot(solver_.PT() * solver_.f());    // f^\top*P_T*f
	    }
        }
        // store results
        f_ = solver_.f();                                 // estimated field at convergence
        g_ = solver_.g();                                 // PDE misfit at convergence
        f_norm_ = std::sqrt(f_.dot(solver_.R0() * f_));   // L^2 norm of estimated field
        //fn_norm_ = fn_.norm();                            // euclidean norm of estimated field at data locations

	// TOOD: check this, should we always return the euclidean norm of the field evalutation?
	if (m_->sampling() == Sampling::mesh_nodes)
            fn_norm_ = f_norm_;
        else
            fn_norm_ = fn_.norm();
	
	return;
    }

    // getters
    const DVector<double>& f() const { return f_; }
    const DVector<double>& g() const { return g_; }
    const DVector<double>& s() const { return s_; }     // scores vector
    const DVector<double>& fn() const { return fn_; }   // loadings vector
    std::size_t n_iter() const { return k_; }
    double f_norm() const { return f_norm_; }
    double fn_norm() const { return fn_norm_; }
    double gcv() { return gcv_.eval(); }   // GCV index at convergence
    // setters
    void set_tolerance(double tolerance) { tolerance_ = tolerance; }
    void set_max_iter(std::size_t max_iter) { max_iter_ = max_iter; }
    void set_seed(std::size_t seed) { seed_ = seed; }
};

}   // namespace models
}   // namespace fdapde

#endif   // __POWER_ITERATION_H__

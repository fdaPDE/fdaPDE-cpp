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
    GCV gcv_;   // GCV functor
    // algorithm's parameters
    double tolerance_ = 1e-6;     // treshold on |Jnew - Jold| used as stopping criterion
    std::size_t max_iter_ = 20;   // maximum number of iterations before forced stop
    std::size_t k_ = 0;           // iteration index
    SolverType solver_;           // internal solver
    int seed_ = fdapde::random_seed;
    int n_mc_samples_ = 100;

    DVector<double> s_;    // estimated score vector
    DVector<double> fn_;   // field evaluation at data location \frac{\Psi*f}{\norm{f}_{L^2}}
    DVector<double> f_;    // field basis expansion at convergence
    DVector<double> g_;    // PDE misfit at convergence
    double f_norm_;        // L^2 norm of estimated field at converegence
   public:
    // constructors
    PowerIteration() = default;
    PowerIteration(const Model& m, double tolerance, std::size_t max_iter, int seed) :
        tolerance_(tolerance), max_iter_(max_iter),
        seed_((seed == fdapde::random_seed) ? std::random_device()() : seed) {
        // initialize internal smoothing solver
        if constexpr (is_space_only<SolverType>::value) { solver_ = SolverType(m.pde(), m.sampling()); }
        else {
            solver_ = SolverType(m.pde(), m.time_pde(), m.sampling());
            solver_.set_temporal_locations(m.time_locs());
        }
        solver_.set_spatial_locations(m.locs());
        return;
    };
    template <typename ModelType> PowerIteration(const ModelType& m, double tolerance, std::size_t max_iter) :
        PowerIteration(m, tolerance, max_iter, fdapde::random_seed) {};

    void init() { gcv_ = solver_.template gcv<StochasticEDF>(n_mc_samples_, seed_); }
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
            // s = \frac{X*fn}{\norm(X*fn)}
            s_ = X * fn_;
            s_ = s_ / s_.norm();
            // compute loadings (solve smoothing problem)
            solver_.data().template insert<double>(OBSERVATIONS_BLK, X.transpose() * s_);   // X^\top*s
            solver_.solve();
            // prepare for next iteration
            k_++;
            fn_ = solver_.fitted();   // \Psi*f
            // update value of discretized functional: \norm{X - s*f_n^\top}_F + f^\top*P(\lambda)*f
            Jold = Jnew;
            Jnew = (X - s_ * fn_.transpose()).squaredNorm() + solver_.ftPf(lambda);
        }
        // store results
        f_norm_ = std::sqrt(solver_.f().dot(solver_.R0() * solver_.f()));   // L^2 norm of estimated field
        f_ = solver_.f() / f_norm_;                                         // estimated field (L^2 normalized)
        g_ = solver_.g();                                                   // PDE misfit at convergence
        fn_ = solver_.Psi(not_nan()) * f_;   // evaluation of (L^2 unitary norm) estimated field at data locations
        return;
    }

    // getters
    const DVector<double>& f() const { return f_; }   // loadings vector
    const DVector<double>& s() const { return s_; }   // scores vector
    const DVector<double>& g() const { return g_; }
    const DVector<double>& fn() const { return fn_; }   
    double ftPf(const DVector<double>& lambda) const { return solver_.ftPf(lambda); }   // block f^\top*P(\lambda)*f
    std::size_t n_iter() const { return k_; }
    double f_norm() const { return f_norm_; }
    inline double f_squaredNorm() const { return std::pow(f_norm_, 2); }
    double gcv() { return gcv_.eval(); }   // GCV index at convergence
    // setters
    void set_tolerance(double tolerance) { tolerance_ = tolerance; }
    void set_max_iter(std::size_t max_iter) { max_iter_ = max_iter; }
    void set_seed(std::size_t seed) { seed_ = seed; }
};

}   // namespace models
}   // namespace fdapde

#endif   // __POWER_ITERATION_H__

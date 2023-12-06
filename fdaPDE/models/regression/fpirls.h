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

#ifndef __FPIRLS_H__
#define __FPIRLS_H__

#include <fdaPDE/linear_algebra.h>
#include <fdaPDE/utils.h>
using fdapde::core::SMW;
using fdapde::core::SparseBlockMatrix;

#include "../model_macros.h"
#include "../model_traits.h"
#include "distributions.h"
#include "srpde.h"
#include "strpde.h"

namespace fdapde {
namespace models {

// trait to select model type to use in the internal loop of FPIRLS
template <typename Model> class FPIRLS_internal_solver {
   private:
    typedef typename std::decay<Model>::type Model_;
    typedef typename model_traits<Model_>::PDE PDE;
    typedef typename model_traits<Model_>::sampling sampling;
    typedef typename model_traits<Model_>::solver solver;
    typedef typename model_traits<Model_>::regularization regularization;
   public:
    using type = typename std::conditional<
      !is_space_time<Model_>::value, SRPDE<PDE, sampling>,   // space-only problem
      STRPDE<PDE, regularization, sampling, solver>          // space-time problem
      >::type;
};

// a general implementation of the Functional Penalized Iterative Reweighted Least Square (FPIRLS) algorithm
template <typename Model> class FPIRLS {
   private:
    typedef typename std::decay<Model>::type Model_;
    typedef typename FPIRLS_internal_solver<Model_>::type SolverType;
    static_assert(is_regression_model<Model_>::value);
    Model& m_;
    // algorithm's parameters
    double tolerance_;       // treshold on objective functional J to convergence
    std::size_t max_iter_;   // maximum number of iterations before forced stop
    std::size_t k_ = 0;      // FPIRLS iteration index
    SolverType solver_;      // internal solver
   public:
    // constructor
    FPIRLS(const Model& m, double tolerance, std::size_t max_iter) :
        m_(m), tolerance_(tolerance), max_iter_(max_iter) {
        // define internal problem solver
        if constexpr (!is_space_time<Model_>::value)   // space-only
            solver_ = SolverType(m_.pde());
        else {   // space-time
            solver_ = SolverType(m_.pde(), m_.time_domain());
            if constexpr (is_space_time_parabolic<Model_>::value) { solver_.set_initial_condition(m_.s(), false); }
            if constexpr (is_space_time_separable<Model_>::value) { solver_.set_temporal_locations(m_.time_locs()); }
        }
        // solver initialization
        solver_.set_lambda(m_.lambda());
        solver_.set_spatial_locations(m_.locs());
        solver_.set_data(m_.data());   // possible covariates are passed from here

        // solver_.init()  ---> rimosso 
        solver_.init_pde();                      // init differential regularization
        solver_.init_regularization();   // init regularization term
        solver_.init_sampling(true);     // init \Psi matrix, always force recomputation
        solver_.init_nan(m_.nan_idxs());              // analyze and set missingness pattern
        solver_.init_model();
    };

    // executes the FPIRLS algorithm
    void compute() {
        // algorithm initialization
        //mu_ = m_.initialize_mu();
        //distr_.preprocess(mu_); 
	m_.fpirls_init(); 
        // objective functional value at consecutive iterations
        double J_old = tolerance_ + 1;
        double J_new = 0;
        while (k_ < max_iter_ && std::abs(J_new - J_old) > tolerance_) {
            //std::cout << "FPIRLS k = " << k_ << std::endl; 
            // compute weight matrix W and pseudo-observations \tilde{y}
            m_.fpirls_pre_solve_step();
            //auto pair = m_.compute(mu_);
            // solve weighted least square problem
            // \argmin_{\beta, f} [ \norm(W^{1/2}(y - X\beta - f_n))^2 + \lambda \int_D (Lf - u)^2 ]
            solver_.data().template insert<double>(OBSERVATIONS_BLK, m_.py()); //std::get<1>(pair));
            solver_.data().template insert<double>(WEIGHTS_BLK, m_.pW()); //std::get<0>(pair));
            // update solver to change in the weight matrix
            solver_.update_data();
            solver_.update_to_weights();
            solver_.solve();

            // \mu update
	    m_.fpirls_post_solve_step(solver_.fitted());
            //DVector<double> fitted = solver_.fitted();
            //mu_ = distr_.inv_link(fitted);

            // compute value of functional J for this pair (\beta, f)
            double J = m_.data_loss() ;     
            if constexpr (is_space_time_separable<Model>::value) 
                // space-time separable regularization requires to compute the penalty matrix
                // J += solver_.f().dot(m_.P()*solver_.f());  -> equivalente
                J += m_.lambda_D()*solver_.g().dot(m_.R0()*solver_.g()) + m_.lambda_T()*solver_.f().dot(Kronecker(solver_.P1(), solver_.pde().R0())*solver_.f());  
            
            else
                // for a space only problem and space time parabolic we can leverage the following identity
                // \int_D (Lf-u)^2 = g^\top*R_0*g = f^\top*P*f, being P = R_1^\top*(R_0)^{-1}*R_1
                J += m_.lambda_D()*solver_.g().dot(m_.R0() * solver_.g());


            // prepare for next iteration
            k_++; J_old = J_new; J_new = J;
        }
        std::cout << "niter fpirls = " << k_ << std::endl; 
        return;
    }

    // getters
    std::size_t n_iter() const { return k_; }
    const SolverType& solver() const { return solver_; }
};

}   // namespace models
}   // namespace fdapde

#endif   // __FPIRLS_H__

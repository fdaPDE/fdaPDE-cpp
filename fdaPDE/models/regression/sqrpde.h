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

#ifndef __SQRPDE_H__
#define __SQRPDE_H__

#include <memory>
#include <type_traits>

#include <fdaPDE/pde.h>
#include <fdaPDE/utils.h>
using fdapde::core::PDEBase;

#include "../model_base.h"
#include "../model_macros.h"
#include "../model_traits.h"
#include "../sampling_design.h"
#include "regression_base.h"
#include "fpirls.h"
#include "distributions.h"

namespace fdapde {
namespace models {

template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver>
class SQRPDE : public RegressionBase<SQRPDE<PDE, RegularizationType, SamplingDesign, Solver>> {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
   private:
    typedef RegressionBase<SQRPDE<PDE, RegularizationType, SamplingDesign, Solver>> Base;
    double alpha_;            // quantile order
    DVector<double> py_ {};   // y - (1-2*alpha)|y - X*beta - f|
    DVector<double> pW_ {};   // diagonal of W^k = 1/(2*n*|y - X*beta - f|)
    DVector<double> mu_;      // \mu^k = [ \mu^k_1, ..., \mu^k_n ] : mean vector at step k
    DMatrix<double> T_;       // T = \Psi^T*Q*\Psi + P

    double pinball_loss(double x) const { return 0.5 * std::abs(x) + (alpha_ - 0.5) * x; };  // quantile check function
   
    fdapde::SparseLU<SpMatrix<double>> invA_; // factorization of matrix A

    // FPIRLS parameters (set to default)
    std::size_t max_iter_ = 200;  
    double tol_weights_ = 1e-6;  
    double tol_ = 1e-6; 
  public:
   IMPORT_REGRESSION_SYMBOLS;
   using Base::lambda_D;   // smoothing parameter in space
   using Base::n_basis;    // number of spatial basis
   using Base::W_;         // weight matrix
   // constructor
   SQRPDE() = default;
   // space-only constructor
   SQRPDE(const PDE& pde, double alpha = 0.5) : Base(pde), alpha_(alpha) {};

   // setter
   void set_fpirls_tolerance(double tol) { tol_ = tol; }
   void set_fpirls_max_iter(std::size_t max_iter) { max_iter_ = max_iter; }
   void set_alpha(double alpha) { alpha_ = alpha; }

   // ModelBase implementation
   void init_model() { return; }
   void update_to_weights() { return; };   // update model object in case of changes in the weights matrix
   virtual void solve();   // finds a solution to the smoothing problem

   // required by FPIRLS (model_loss computes the unpenalized loss)

   // required by FPIRLS (initialize \mu for the first FPIRLS iteration)
   // non-parametric and semi-parametric cases coincide here, since beta^(0) = 0
   void fpirls_init() {
       // not consistent, fix needed

       // assemble system matrix
       SparseBlockMatrix<double, 2, 2> A(
         PsiTD() * Psi() / n_obs(), 2 * lambda_D() * R1().transpose(),
	 lambda_D() * R1(),         -lambda_D() * R0()               );
       // factorize
       fdapde::SparseLU<SpMatrix<double>> invA;
       invA.compute(A);
       DVector<double> b;
       b.resize(A.rows());
       b.block(n_basis(), 0, n_basis(), 1) = lambda_D() * u();
       b.block(0, 0, n_basis(), 1) = PsiTD() * y() / n_obs();

       mu_ = Psi(not_nan()) * (invA.solve(b)).head(n_basis());
   }

   // required by FPIRLS (computes weight matrix and vector of pseudo-observations)
   // returns a pair of references to W^k and \tilde y^k
   void fpirls_pre_solve_step() {
       DVector<double> abs_res(y().size());   // y().size() -> n_obs() ??

       // abs_res = (y() - mu_).array().abs(); // vectorized

       for (int i = 0; i < y().size(); ++i) abs_res(i) = std::abs(y()(i,0) - mu_[i]);

       pW_.resize(n_obs());

       // pW_ = (abs_res.array() < tol_weights_).select() ??

       for (int i = 0; i < y().size(); ++i) {
           if (abs_res(i) < tol_weights_) {   // avoid zero weight
               pW_[i] = (1. / (abs_res[i] + tol_weights_)) / (2. * n_obs());
           } else {
               pW_[i] = (1. / abs_res[i]) / (2. * n_obs());
           }
       }
       py_ = y() - (1 - 2. * alpha_) * abs_res;
   }
    void fpirls_post_solve_step(const DMatrix<double>& hat_f, const DMatrix<double>& hat_beta) {
        mu_ = hat_f; 
    }
    // returns the data loss ***
    double data_loss() const { return (pW_.cwiseSqrt().matrix().asDiagonal() * (py_ - mu_)).squaredNorm(); }
    const DVector<double>& py() const { return py_; }
    const DVector<double>& pW() const { return pW_; }

    // GCV support
    const DMatrix<double>& T();
    double norm(const DMatrix<double>& op1, const DMatrix<double>& op2) const {
        double result = 0;
        for (std::size_t i = 0; i < op2.rows(); ++i) { result += pinball_loss(op2.coeff(i, 0) - op1.coeff(i, 0)); }
        return result * result / n_obs();
    }

   // getters
   const fdapde::SparseLU<SpMatrix<double>>& invA() const { return invA_; }
   const DMatrix<double>& U() const { return U_; }
   const DMatrix<double>& V() const { return V_; }

   virtual ~SQRPDE() = default;
  };
 
template <
   typename PDE_, typename RegularizationType_, typename SamplingDesign_, typename Solver_>
struct model_traits<SQRPDE<PDE_, RegularizationType_, SamplingDesign_, Solver_>> {
    typedef PDE_ PDE;
    typedef SpaceOnly regularization;
    typedef SamplingDesign_ sampling;
    typedef MonolithicSolver solver;
    enum { N = PDE::N, M = PDE::M, n_lambda = 1 };
  };

  // finds a solution to the SQR-PDE smoothing problem
template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver>
void SQRPDE<PDE, RegularizationType, SamplingDesign, Solver>::solve() {
    // execute FPIRLS for minimization of the functional
    // \norm{V^{-1/2}(y - \mu)}^2 + \lambda \int_D (Lf - u)^2

    FPIRLS<decltype(*this)> fpirls(*this, tol_, max_iter_);   // FPIRLS engine
    fpirls.compute();

    // fpirls converged: extract matrix P and solution estimates
    // W_ = fpirls.weights().asDiagonal();
    W_ = fpirls.solver().W();

    // if (has_covariates()) {
    //     XtWX_ = X().transpose() * W_ * X(); // get from solver ?? 
    //     invXtWX_ = XtWX_.partialPivLu();
    // }

    invA_ = fpirls.solver().invA();

    // if (has_covariates()) {
    //     U_ = fpirls.solver().U();
    //     V_ = fpirls.solver().V();
    // }

    f_ = fpirls.solver().f();
    g_ = fpirls.solver().g();

    if (has_covariates()) beta_ = fpirls.solver().beta();
    return;
  }

// Non-parametric and semi-parametric cases coincide here, since beta^(0) = 0
// template <typename PDE, typename SamplingDesign>
// DVector<double> 
// SQRPDE<PDE, SamplingDesign>::initialize_mu() const{

//     // assemble system matrix 
//     SparseBlockMatrix<double,2,2>
//       A_init(PsiTD()*Psi()/n_obs(), 2*lambda_D()*R1().transpose(),
//         lambda_D()*R1(),     -lambda_D()*R0()            );

//     // cache non-parametric matrix and its factorization for reuse 
//     fdapde::SparseLU<SpMatrix<double>> invA_init;
//     invA_init.compute(A_init);
//     DVector<double> b_init; 
//     b_init.resize(A_init.rows());
//     b_init.block(n_basis(),0, n_basis(),1) = lambda_D()*u();  
//     b_init.block(0,0, n_basis(),1) = PsiTD()*y()/n_obs(); 
//     BLOCK_FRAME_SANITY_CHECKS;
//     DVector<double> f = (invA_init.solve(b_init)).head(n_basis());
//     DVector<double> fn =  Psi(not_nan())*f;  

//     return fn;
  
// }

// template <typename PDE, typename SamplingDesign>
// std::tuple<DVector<double>&, DVector<double>&>
// SQRPDE<PDE, SamplingDesign>::compute(const DVector<double>& mu){
//   // compute weight matrix and pseudo-observation vector
//   DVector<double> abs_res{};
//   abs_res.resize(y().size()); 

//   for(int i = 0; i < y().size(); ++i)
//     abs_res(i) = std::abs(y()(i) - mu(i));   

//   pW_.resize(n_obs());
//   for(int i = 0; i < y().size(); ++i) {
//     if (abs_res(i) < tol_weights_){
//       pW_(i) = ( 1./(abs_res(i) + tol_weights_) )/(2.*n_obs());

//     }    
//     else
//       pW_(i) = (1./abs_res(i))/(2.*n_obs()); 
//   }
 
//   py_ = y() - (1 - 2.*alpha_)*abs_res;

//   return std::tie(pW_, py_);
// }



// template <typename PDE, typename SamplingDesign>
// double
// SQRPDE<PDE, SamplingDesign>::model_loss(const DVector<double>& mu) const{
  
//   // compute value of functional J given mu: /(2*n) 
//     return (pW_.cwiseSqrt().matrix().asDiagonal()*(py_ - mu)).squaredNorm();     
// }

// required to support GCV based smoothing parameter selection
// in case of an SRPDE model we have T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
// template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver>
// const DMatrix<double>& SQRPDE<PDE, RegularizationType, SamplingDesign, Solver>::T() {
//   // compute value of R = R1^T*R0^{-1}*R1, cache for possible reuse
//   if(R_.size() == 0){
//       invR0_.compute(R0());
//       R_ = R1().transpose()*invR0_.solve(R1());
//   }
//   // compute and store matrix T for possible reuse
//   if(!has_covariates()) // case without covariates, Q is the identity matrix
//     T_ = PsiTD()*W()*Psi()   + lambda_D()*R_;
//   else // general case with covariates
//     T_ = PsiTD()*lmbQ(Psi()) + lambda_D()*R_;
//   return T_;
// }

// Q is computed on demand only when it is needed by GCV and cached for fast reacess (in general operations
// involving Q can be substituted with the more efficient routine lmbQ(), which is part of iRegressionModel interface)
// template <typename PDE, typename SamplingDesign>
// const DMatrix<double>& SQRPDE<PDE, SamplingDesign>::Q() {
//   // compute Q = W(I - H) = W ( I - X*(X^T*W*X)^{-1}*X^T*W ) 
//   Q_ = W()*(DMatrix<double>::Identity(n_obs(), n_obs()) - X()*invXtWX().solve(X().transpose()*W()));

//   return Q_;
// }


// returns the numerator of the GCV score
// template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver>
// double SQRPDE<PDE, SamplingDesign>::norm(const DMatrix<double>& op1, const DMatrix<double>& op2) const {
//   double result = 0;
//   for (std::size_t i = 0; i < op2.rows(); ++i) { result += pinball_loss(op2.coeff(i, 0) - op1.coeff(i, 0)); }
//   return result * result / n_op2();
// }

// returns the pinball loss at a specific x 
// template <typename PDE, typename SamplingDesign>
// double SQRPDE<PDE, SamplingDesign>::rho_alpha(const double& x) const{ 
//   return 0.5*std::abs(x) + (alpha_ - 0.5)*x; 
// }





  
// base class for SQRPDE model
// template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver, typename Distribution>
// class SQRPDE : public RegressionBase<SQRPDE<PDE, RegularizationType, SamplingDesign, Solver, Distribution>> {
//     // compile time checks
//     static_assert(std::is_base_of<PDEBase, PDE>::value);
//    private:
//     typedef RegressionBase<SQRPDE<PDE, RegularizationType, SamplingDesign, Solver, Distribution>> Base;
//     Distribution distr_ {};
//     DVector<double> py_;   // \tilde y^k = G^k(y-u^k) + \theta^k
//     DVector<double> pW_;   // diagonal of W^k = ((G^k)^{-2})*((V^k)^{-1})
//     DVector<double> mu_;   // \mu^k = [ \mu^k_1, ..., \mu^k_n ] : mean vector at step k
//     DMatrix<double> T_;    // T = \Psi^T*Q*\Psi + P

//     // FPIRLS parameters (set to default)
//     std::size_t max_iter_ = 200;
//     double tol_ = 1e-4;
//    public:
//     IMPORT_REGRESSION_SYMBOLS;
//     using Base::lambda_D;   // smoothing parameter in space
//     using Base::P;          // discretized penalty
//     using Base::W_;         // weight matrix
//     // constructor
//     SQRPDE() = default;
//     // space-only constructor
//     template <
//       typename U = RegularizationType, typename std::enable_if<std::is_same<U, SpaceOnly>::value, int>::type = 0>
//     SQRPDE(const PDE& pde) : Base(pde) {};
//     // space-time constructor
//     template <
//       typename U = RegularizationType, typename std::enable_if<!std::is_same<U, SpaceOnly>::value, int>::type = 0>
//     SQRPDE(const PDE& pde, const DVector<double>& time) : Base(pde, time) {};

//     // setter
//     void set_fpirls_tolerance(double tol) { tol_ = tol; }
//     void set_fpirls_max_iter(std::size_t max_iter) { max_iter_ = max_iter; }

//     void init_model() { return; };          // update model object in case of **structural** changes in its definition
//     void update_to_weights() { return; };   // update model object in case of changes in the weights matrix
//     virtual void solve();                   // finds a solution to the smoothing problem

//     // required by FPIRLS (see fpirls.h for details)
//     // initalizes mean vector \mu
//     void fpirls_init() {
//         mu_ = y();
//         distr_.preprocess(mu_);
//     };
//     // computes W^k = ((G^k)^{-2})*((V^k)^{-1}) and \tilde y^k = G^k(y-u^k) + \theta^k
//     void fpirls_pre_solve_step() {
//         DVector<double> theta_ = distr_.link(mu_);   // \theta^k = (g(\mu^k_1), ..., g(\mu^k_n))
//         DVector<double> G_ = distr_.der_link(mu_);   // G^k = diag(g'(\mu^k_1), ..., g'(\mu^k_n))
//         DVector<double> V_ = distr_.variance(mu_);   // V^k = diag(v(\mu^k_1), ..., v(\mu^k_n))
//         pW_ = ((G_.array().pow(2) * V_.array()).inverse()).matrix();
//         py_ = G_.asDiagonal() * (y() - mu_) + theta_;
//     }
//     // updates mean vector \mu after WLS solution
//     void fpirls_post_solve_step(const DMatrix<double>& hat_f, const DMatrix<double>& hat_beta) {
//         mu_ = distr_.inv_link(hat_f);
//     }
//     // returns the data loss \norm{V^{-1/2}(y - \mu)}^2
//     double data_loss() const {
//         DVector<double> V = distr_.variance(mu_).array().sqrt().inverse().matrix();
//         return (V.asDiagonal() * (y() - mu_)).squaredNorm();
//     }
//     const DVector<double>& py() const { return py_; }
//     const DVector<double>& pW() const { return pW_; }

//     // GCV support
//     const DMatrix<double>& T();                                                  // T = \Psi^T*Q*\Psi + P
//     double norm(const DMatrix<double>& op1, const DMatrix<double>& op2) const;   // total deviance \sum dev(op1 - op2)

//     virtual ~SQRPDE() = default;
// };

// // implementative details
  
// // finds a solution to the SQRPDE smoothing problem
// template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver, typename Distribution>
// void SQRPDE<PDE, RegularizationType, SamplingDesign, Solver, Distribution>::solve() {
//     // execute FPIRLS for minimization of functional \norm{V^{-1/2}(y - \mu)}^2 + \lambda \int_D (Lf - u)^2
//     FPIRLS<decltype(*this)> fpirls(*this, tol_, max_iter_);   // FPIRLS engine
//     fpirls.compute();

//     // fpirls converged: extract matrix W and solution estimates
//     W_ = fpirls.solver().W();
//     f_ = fpirls.solver().f();
//     if (has_covariates()) { beta_ = fpirls.solver().beta(); }
//     return;
// }

// template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver, typename Distribution>
// const DMatrix<double>& SQRPDE<PDE, RegularizationType, SamplingDesign, Solver, Distribution>::T() {
//     if (!has_covariates())   // case without covariates, Q is the identity matrix
//         T_ = PsiTD() * W() * Psi() + P();
//     else   // general case with covariates
//         T_ = PsiTD() * lmbQ(Psi()) + P();
//     return T_;
// }

// // returns the deviance of y - \hat y induced by the considered distribution
// template <typename PDE, typename RegularizationType, typename SamplingDesign, typename Solver, typename Distribution>
// double SQRPDE<PDE, RegularizationType, SamplingDesign, Solver, Distribution>::norm(
//   const DMatrix<double>& obs, const DMatrix<double>& fitted) const {
//     DMatrix<double> mu = distr_.inv_link(fitted);
//     double result = 0;
//     for (std::size_t i = 0; i < obs.rows(); ++i) { result += distr_.deviance(mu.coeff(i, 0), obs.coeff(i, 0)); }
//     return result;
// }

// template <
//   typename PDE_, typename RegularizationType_, typename SamplingDesign_, typename Solver_, typename DistributionType_>
// struct model_traits<SQRPDE<PDE_, RegularizationType_, SamplingDesign_, Solver_, DistributionType_>> {
//     typedef PDE_ PDE;
//     typedef RegularizationType_ regularization;
//     typedef SamplingDesign_ sampling;
//     typedef Solver_ solver;
//     typedef DistributionType_ DistributionType;
//     enum { N = PDE::N, M = PDE::M, n_lambda = n_smoothing_parameters<RegularizationType_>::value };
// };
// // specialization for separable regularization
// template <typename PDE_, typename SamplingDesign_, typename Solver_, typename DistributionType_>
// struct model_traits<SQRPDE<PDE_, fdapde::models::SpaceTimeSeparable, SamplingDesign_, Solver_, DistributionType_>> {
//     typedef PDE_ PDE;
//     typedef fdapde::models::SpaceTimeSeparable regularization;
//     typedef SplineBasis<3> TimeBasis;   // use cubic B-splines
//     typedef SamplingDesign_ sampling;
//     typedef Solver_ solver;
//     typedef DistributionType_ DistributionType;
//     enum { N = PDE::N, M = PDE::M, n_lambda = 2 };
// };

// // gsrpde trait
// template <typename Model> struct is_gsrpde {
//     static constexpr bool value = is_instance_of<Model, SQRPDE>::value;
// };
  
}   // namespace models
}   // namespace fdapde

#endif   // __SQRPDE_H__

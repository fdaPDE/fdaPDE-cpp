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

#ifndef __FPCA_H__
#define __FPCA_H__

#include <fdaPDE/optimization.h>
#include <fdaPDE/utils.h>
#include <Eigen/SVD>

#include "functional_base.h"
#include "power_iteration.h"
#include "rsvd.h"

namespace fdapde {
namespace models {

// FPCA (Functional Principal Components Analysis) model signature
template <typename RegularizationType, typename SolutionPolicy> class FPCA;
  
// implementation of FPCA for sequential approach, see e.g.
// Lila, E., Aston, J.A.D., Sangalli, L.M. (2016), Smooth Principal Component Analysis over two-dimensional manifolds
// with an application to Neuroimaging, Annals of Applied Statistics, 10 (4), 1854-1879.
template <typename RegularizationType_>
class FPCA<RegularizationType_, sequential> :
    public FunctionalBase<FPCA<RegularizationType_, sequential>, RegularizationType_> {
   public:
    using RegularizationType = std::decay_t<RegularizationType_>;
    using This = FPCA<RegularizationType, sequential>;
    using Base = FunctionalBase<This, RegularizationType>;
    IMPORT_MODEL_SYMBOLS;
    using Base::lambda;
    using Base::X;

    // constructors
    FPCA() = default;
    fdapde_enable_constructor_if(is_space_only, This) FPCA(const pde_ptr& pde, Sampling s) : Base(pde, s) {};
    fdapde_enable_constructor_if(is_space_time_separable, This)
      FPCA(const pde_ptr& space_penalty, const pde_ptr& time_penalty, Sampling s) :
        Base(space_penalty, time_penalty, s) {};

    void init_model() { return; };
    void solve() {
        loadings_.resize(X().cols(), n_pc_);
        scores_.resize(X().rows(), n_pc_);
        DMatrix<double> X_ = X();   // copy original data to avoid side effects

        // first guess of PCs set to a multivariate PCs (SVD)
        Eigen::JacobiSVD<DMatrix<double>> svd(X_, Eigen::ComputeThinU | Eigen::ComputeThinV);
        PowerIteration<This> solver(this, tolerance_, max_iter_);   // power iteration solver
        solver.init();
        // sequential extraction of principal components
        for (std::size_t i = 0; i < n_pc_; i++) {
            // find vectors s,f minimizing \norm_F{Y - s^T*f}^2 + (s^T*s)*P(f) fixed \lambda, using the multivariate
            // estimate of f as starting point
            solver.compute(X_, lambda(), svd.matrixV().col(i));
	    // store results and deflate
            loadings_.col(i) = solver.fn() / solver.fn_norm();
            scores_.col(i) = solver.s() * solver.fn_norm();
            X_ -= scores_.col(i) * loadings_.col(i).transpose();
        }
        return;
    }

    // getters
    const DMatrix<double>& loadings() const { return loadings_; }
    const DMatrix<double>& scores() const { return scores_; }
    // setters
    void set_npc(std::size_t n_pc) { n_pc_ = n_pc; }
    void set_tolerance(double tolerance) { tolerance_ = tolerance; }
    void set_max_iter(std::size_t max_iter) { max_iter_ = max_iter; }
    void set_seed(std::size_t seed) { seed_ = seed; }
   private:
    std::size_t n_pc_ = 3;   // number of principal components
    // power iteration parameters
    double tolerance_ = 1e-6;     // relative tolerance between Jnew and Jold, used as stopping criterion
    std::size_t max_iter_ = 20;   // maximum number of allowed iterations
    int seed_ = fdapde::random_seed;

    // problem solution
    DMatrix<double> loadings_;   // evaluation of the PC functions at data locations
    DMatrix<double> scores_;
};

// implementation of FPCA for monolithic approach, see e.g.
template <typename RegularizationType_>
class FPCA<RegularizationType_, monolithic> :
    public FunctionalBase<FPCA<RegularizationType_, monolithic>, RegularizationType_> {
   public:
    using RegularizationType = std::decay_t<RegularizationType_>;
    using This = FPCA<RegularizationType, monolithic>;
    using Base = FunctionalBase<This, RegularizationType>;
    IMPORT_MODEL_SYMBOLS;
    using Base::lambda;
    using Base::X;

    // constructors
    FPCA() = default;
    fdapde_enable_constructor_if(is_space_only, This) FPCA(const pde_ptr& pde, Sampling s) : Base(pde, s) {};
    fdapde_enable_constructor_if(is_space_time_separable, This)
      FPCA(const pde_ptr& space_penalty, const pde_ptr& time_penalty, Sampling s) :
        Base(space_penalty, time_penalty, s) {};

    void init_model() { return; };
    void solve() {
        // monolithic approach via Regularized SVD
	RSVD<This> solver(this);
	solver.compute(X(), lambda(), n_pc_);
	loadings_ = Psi() * solver.W();
	scores_ = solver.H();
	// normalization
	for(std::size_t i = 0; i < n_pc_; ++i) {
	  double fn_norm_ = std::sqrt(solver.W().col(i).dot(R0() * solver.W().col(i)));
	  loadings_.col(i) /= fn_norm_;
	  scores_.col(i) *= fn_norm_;
	}
        return;
    }

    // getters
    const DMatrix<double>& loadings() const { return loadings_; }
    const DMatrix<double>& scores() const { return scores_; }
    void set_npc(std::size_t n_pc) { n_pc_ = n_pc; }
   private:
    std::size_t n_pc_ = 3;   // number of principal components
    // problem solution
    DMatrix<double> loadings_;   // evaluation of the PC functions at data locations
    DMatrix<double> scores_;
};


}   // namespace models
}   // namespace fdapde

#endif   // __FPCA_H__

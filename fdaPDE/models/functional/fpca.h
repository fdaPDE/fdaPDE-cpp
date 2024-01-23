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
#include <fdaPDE/linear_algebra.h>
#include <Eigen/SVD>

#include "../../calibration/kfold_cv.h"
#include "../../calibration/symbols.h"
using fdapde::calibration::Calibration;
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
    using Base::lambda;         // vector of smoothing parameters
    using Base::n_basis;        // number of basis functions
    using Base::n_stat_units;   // number of statistical units
    using Base::X;              // n_stat_units \times n_locs data matrix

    // constructors
    FPCA() = default;
    fdapde_enable_constructor_if(is_space_only, This) FPCA(const pde_ptr& pde, Sampling s, Calibration c) :
        Base(pde, s), calibration_(c) {};
    fdapde_enable_constructor_if(is_space_time_separable, This)
      FPCA(const pde_ptr& space_penalty, const pde_ptr& time_penalty, Sampling s, Calibration c) :
        Base(space_penalty, time_penalty, s), calibration_(c) {};

    void init_model() { return; };
    void solve() {
        // preallocate space
        loadings_.resize(n_basis(), n_pc_);
        fitted_loadings_.resize(n_locs(), n_pc_);
        scores_.resize(n_stat_units(), n_pc_);
        DMatrix<double> X_ = X();   // copy original data to avoid side effects

        // first guess of PCs set to a multivariate PCA (SVD)
        Eigen::JacobiSVD<DMatrix<double>> svd(X_, Eigen::ComputeThinU | Eigen::ComputeThinV);
        PowerIteration<This> solver(this, tolerance_, max_iter_);   // power iteration solver
        solver.set_seed(seed_);
        solver.init();

        // sequential extraction of principal components
        for (std::size_t i = 0; i < n_pc_; i++) {
            DVector<double> f0 = svd.matrixV().col(i);
            switch (calibration_) {
            case Calibration::off: {
                // find vectors s,f minimizing \norm_F{Y - s^T*f}^2 + (s^T*s)*P(f) fixed \lambda
                solver.compute(X_, lambda(), f0);
            } break;
            case Calibration::gcv: {
                // select \lambda minimizing the GCV index
                ScalarField<Dynamic> gcv([&solver, &X_, &f0](const DVector<double>& lambda) -> double {
                    solver.compute(X_, lambda, f0);
                    return solver.gcv();   // return GCV index at convergence
                });
                solver.compute(X_, core::Grid<Dynamic> {}.optimize(gcv, lambda_grid_), f0);
            } break;
            case Calibration::kcv: {
                auto cv_score = [&solver, &X_, &f0](
                                  const DVector<double>& lambda, const core::BinaryVector<Dynamic>& train_set,
                                  const core::BinaryVector<Dynamic>& test_set) -> double {
                    solver.compute(train_set.blk_repeat(1, X_.cols()).select(X_), lambda, f0);   // fit on train set
                    // reconstruction error on test set: \norm{X_test * (I - fn*fn^\top/J)}_F/n_test, with
                    // J = \norm{f_n}_2^2 + f^\top*P(\lambda)*f (PS: the division of f^\top*P(\lambda)*f by
                    // \norm{f}_{L^2} is necessary to obtain J as expected)
                    return (test_set.blk_repeat(1, X_.cols()).select(X_) *
                            (DMatrix<double>::Identity(X_.cols(), X_.cols()) -
                             solver.fn() * solver.fn().transpose() /
                               (solver.fn().squaredNorm() + solver.ftPf(lambda) / solver.f_squaredNorm())))
                             .squaredNorm() /
                           test_set.count() * X_.cols();
                };
                solver.compute(X_, calibration::KCV {n_folds_, seed_}.fit(*this, lambda_grid_, cv_score), f0);
            } break;
            }
            // store results
            loadings_.col(i) = solver.f();
            fitted_loadings_.col(i) = solver.fn();   // \frac{\Psi * f}{\norm{f}_{L^2}}
            scores_.col(i) = solver.s() * solver.f_norm();
            X_ -= scores_.col(i) * fitted_loadings_.col(i).transpose();   // X <- X - s*f_n^\top (deflation step)
        }
        return;
    }

    // getters
    const DMatrix<double>& fitted_loadings() const { return fitted_loadings_; }
    const DMatrix<double>& loadings() const { return loadings_; }
    const DMatrix<double>& scores() const { return scores_; }
    // setters
    void set_npc(std::size_t n_pc) { n_pc_ = n_pc; }
    void set_tolerance(double tolerance) { tolerance_ = tolerance; }
    void set_max_iter(std::size_t max_iter) { max_iter_ = max_iter; }
    void set_seed(std::size_t seed) { seed_ = seed; }
    void set_lambda(const std::vector<DVector<double>>& lambda_grid) {
        fdapde_assert(calibration_ != Calibration::off);
        lambda_grid_ = lambda_grid;
    }
    void set_nfolds(std::size_t n_folds) {
        fdapde_assert(calibration_ == Calibration::kcv);
        n_folds_ = n_folds;
    }
   private:
    std::size_t n_pc_ = 3;       // number of principal components
    Calibration calibration_;    // PC function's smoothing parameter selection strategy
    std::size_t n_folds_ = 10;   // for a kcv calibration strategy, the number of folds
    std::vector<DVector<double>> lambda_grid_;
    // power iteration parameters
    double tolerance_ = 1e-6;     // relative tolerance between Jnew and Jold, used as stopping criterion
    std::size_t max_iter_ = 20;   // maximum number of allowed iterations
    int seed_ = fdapde::random_seed;

    // problem solution
    DMatrix<double> loadings_;          // PC functions' expansion coefficients
    DMatrix<double> fitted_loadings_;   // evaluation of the PC functions at data locations
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
    void solve() {   // monolithic approach via Regularized SVD
        RSVD<This> solver(this);
        solver.compute(X(), lambda(), n_pc_);
        // store results
	loadings_ = solver.W();
        fitted_loadings_ = Psi(not_nan()) * solver.W();   // evaluation of L^2 normalized PC functions at data locations
        scores_ = solver.H().array().rowwise() * solver.w_norm().transpose().array();
        return;
    }

    // getters
    const DMatrix<double>& fitted_loadings() const { return fitted_loadings_; }
    const DMatrix<double>& loadings() const { return loadings_; }
    const DMatrix<double>& scores() const { return scores_; }
    void set_npc(std::size_t n_pc) { n_pc_ = n_pc; }
   private:
    std::size_t n_pc_ = 3;   // number of principal components
    // problem solution
    DMatrix<double> loadings_;          // PC functions' expansion coefficients
    DMatrix<double> fitted_loadings_;   // evaluation of the PC functions at data locations
    DMatrix<double> scores_;
};


}   // namespace models
}   // namespace fdapde

#endif   // __FPCA_H__

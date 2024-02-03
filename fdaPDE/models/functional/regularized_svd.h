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

#ifndef __REGULARIZED_SVD_H__
#define __REGULARIZED_SVD_H__

#include <fdaPDE/utils.h>
#include <Eigen/SVD>

#include "../../calibration/kfold_cv.h"
#include "../../calibration/symbols.h"
using fdapde::calibration::Calibration;
#include "../model_traits.h"
#include "power_iteration.h"


namespace fdapde {
namespace models {

// Let X be a data matrix made of noisy and discrete measurements of smooth functions sampled from a random field
// \mathcal{X}. RegularizedSVD implements the computation of a low-rank approximation of X using some regularizing term
template <typename SolutionPolicy_> class RegularizedSVD;

// Finds a low-rank approximation of X while penalizing for the eigenfunctions of \mathcal{X} by sequentially solving
// \argmin_{s,f} \norm_F{X - s^\top*f}^2 + (s^\top*s)*P_{\lambda}(f), up to a desired rank
template <> class RegularizedSVD<sequential> {
   private:
    Calibration calibration_;    // PC function's smoothing parameter selection strategy
    std::size_t n_folds_ = 10;   // for a kcv calibration strategy, the number of folds
    std::vector<DVector<double>> lambda_grid_;
    // power iteration parameters
    double tolerance_ = 1e-6;     // relative tolerance between Jnew and Jold, used as stopping criterion
    std::size_t max_iter_ = 20;   // maximum number of allowed iterations
    int seed_ = fdapde::random_seed;

    // problem solution
    DMatrix<double> scores_;
    DMatrix<double> loadings_;        // PC functions' expansion coefficients
    DVector<double> loadings_norm_;   // L^2 norm of estimated fields

    // basic rsvd iterator (should be used only for the end() iterator)
    struct rsvd_iterator_base {
        friend RegularizedSVD;
        RegularizedSVD* rsvd_;
        int index_;   // current rank
        rsvd_iterator_base(RegularizedSVD* rsvd, int index) : rsvd_(rsvd), index_(index) {};
        friend bool operator!=(const rsvd_iterator_base& lhs, const rsvd_iterator_base& rhs) {
            return lhs.index_ != rhs.index_;
        }
    };

    // this iterator allows to range over each component of X, while triggering the computation only when required
    template <typename ModelType> struct rsvd_iterator : public rsvd_iterator_base {
       private:
        PowerIteration<ModelType> solver_;        // rank-one step solver
        DMatrix<double> X_;                       // deflated data
        Eigen::JacobiSVD<DMatrix<double>> svd_;   // thin monolithic Singular Value Decomposition (not regularized)
        ModelType& model_;

        // rank-one step: \argmin_{s,f} \norm_F{X - s^\top*f}^2 + (s^\top*s)*P_{\lambda}(f). calibration of \lambda
        // dispatched to desired strategy
        void rank_one_step() {
            DVector<double> f0 = svd_.matrixV().col(index_);
            switch (rsvd_->calibration_) {
            case Calibration::off: {
                // find vectors s,f minimizing \norm_F{Y - s^T*f}^2 + (s^T*s)*P(f) fixed \lambda
                solver_.compute(X_, model_.lambda(), f0);
            } break;
            case Calibration::gcv: {
                // select \lambda minimizing the GCV index
                ScalarField<Dynamic> gcv([&](const DVector<double>& lambda) -> double {
                    solver_.compute(X_, lambda, f0);
                    return solver_.gcv();   // return GCV index at convergence
                });
                solver_.compute(X_, core::Grid<Dynamic> {}.optimize(gcv, rsvd_->lambda_grid_), f0);
            } break;
            case Calibration::kcv: {
                // select \lambda minimizing the reconstruction error in cross-validation
                auto cv_score = [&](
                                  const DVector<double>& lambda, const core::BinaryVector<Dynamic>& train_set,
                                  const core::BinaryVector<Dynamic>& test_set) -> double {
                    solver_.compute(train_set.blk_repeat(1, X_.cols()).select(X_), lambda, f0);   // fit on train set
                    // reconstruction error on test set: \norm{X_test * (I - fn*fn^\top/J)}_F/n_test, with
                    // J = \norm{f_n}_2^2 + f^\top*P(\lambda)*f (PS: the division of f^\top*P(\lambda)*f by
                    // \norm{f}_{L^2} is necessary to obtain J as expected)
                    return (test_set.blk_repeat(1, X_.cols()).select(X_) *
                            (DMatrix<double>::Identity(X_.cols(), X_.cols()) -
                             solver_.fn() * solver_.fn().transpose() /
                               (solver_.fn().squaredNorm() + solver_.ftPf(lambda) / solver_.f_squaredNorm())))
                             .squaredNorm() /
                           test_set.count() * X_.cols();
                };
                solver_.compute(
                  X_, calibration::KCV {rsvd_->n_folds_, rsvd_->seed_}.fit(model_, rsvd_->lambda_grid_, cv_score), f0);
            } break;
            }
            X_ -= solver_.s() * solver_.fn().transpose() * solver_.f_norm();   // X <- X - s*f_n^\top (deflation step)
            return;
        }
       public:
        // constructor
        rsvd_iterator(RegularizedSVD* rsvd, int index, const DMatrix<double>& X, ModelType& model) :
            rsvd_iterator_base(rsvd, index),
            X_(X),
            solver_(model, rsvd->tolerance_, rsvd->max_iter_, rsvd->seed_),
            model_(model) {
            // first guess of PCs set to a multivariate PCA (SVD)
            svd_ = Eigen::JacobiSVD<DMatrix<double>>(X_, Eigen::ComputeThinU | Eigen::ComputeThinV);
            solver_.init();   // initialize power iteration solver
        };
        rsvd_iterator& operator++() {
            ++index_;
            return *this;
        }
        bool operator!=(const rsvd_iterator_base& rhs) {
            if (index_ != rhs.index_) { rank_one_step(); }   // trigger computation of next component only if not ended
            return index_ != rhs.index_;
        }
        // getters
        const DVector<double>& scores() const { return solver_.s(); }
        const DVector<double>& loading() const { return solver_.f(); }
        double norm() const { return solver_.f_norm(); }
    };
   public:
    // constructors
    RegularizedSVD(Calibration c) : calibration_(c) { }
    RegularizedSVD() : RegularizedSVD(Calibration::off) {};

    // sequentially solves \argmin_{s,f} \norm_F{X - s^\top*f}^2 + (s^\top*s)*P_{\lambda}(f), up to the specified rank,
    // selecting the level of smoothing of the component according to the desired strategy
    template <typename ModelType> void compute(const DMatrix<double>& X, std::size_t rank, ModelType& model) {
        // preallocate space
        loadings_.resize(model.n_basis(), rank);
        scores_.resize(X.rows(), rank);
        loadings_norm_.resize(rank);
        int i = 0;
        for (auto it = begin(X, model); it != end(rank); ++it, ++i) {
            loadings_.col(i) = it.loading();
            scores_.col(i) = it.scores() * it.norm();
            loadings_norm_[i] = it.norm();
        }
    }
    // iterator support
    template <typename ModelType> rsvd_iterator<ModelType> begin(const DMatrix<double>& X, ModelType&& model) {
        return rsvd_iterator<ModelType>(this, 0, X, model);
    }
    rsvd_iterator_base end(int rank) { return rsvd_iterator_base(this, rank); }
    // getters
    const DMatrix<double>& scores() const { return scores_; }
    const DMatrix<double>& loadings() const { return loadings_; }
    const DVector<double>& loadings_norm() const { return loadings_norm_; }
    // setters
    void set_tolerance(double tolerance) { tolerance_ = tolerance; }
    void set_max_iter(std::size_t max_iter) { max_iter_ = max_iter; }
    void set_seed(std::size_t seed) { seed_ = seed; }
    RegularizedSVD& set_lambda(const std::vector<DVector<double>>& lambda_grid) {
        fdapde_assert(calibration_ != Calibration::off);
        lambda_grid_ = lambda_grid;
        return *this;
    }
    RegularizedSVD& set_nfolds(std::size_t n_folds) {
        fdapde_assert(calibration_ == Calibration::kcv);
        n_folds_ = n_folds;
        return *this;
    }
};

// finds a rank r matrix U minimizing \norm{X - U*\Psi^\top}_F^2 + Tr[U*P_{\lambda}(f)*U^\top]
template <> class RegularizedSVD<monolithic> {
   public:
    // solves \norm{X - U*\Psi^\top}_F^2 + Tr[U*P_{\lambda}(f)*U^\top] retaining the first rank components
    template <typename ModelType> void compute(const DMatrix<double>& X, std::size_t rank, ModelType& model) {
        // compute matrix C = \Psi^\top*\Psi + P(\lambda)
        DMatrix<double> C = model.Psi().transpose() * model.Psi() + model.P();
        // compute the inverse of the cholesky factor of C, D^{-1}
        DMatrix<double> D = C.llt().matrixL();
        DMatrix<double> invD = D.inverse();

        // compute SVD of X*\Psi*(D^{-1})^\top
        Eigen::JacobiSVD<DMatrix<double>> svd(
          X * model.Psi() * invD.transpose(), Eigen::ComputeThinU | Eigen::ComputeThinV);
        // store results
        scores_ = svd.matrixU().leftCols(rank);
        loadings_ =
          (svd.singularValues().head(rank).asDiagonal() * svd.matrixV().leftCols(rank).transpose() * invD).transpose();
        loadings_norm_.resize(rank);
        for (std::size_t i = 0; i < rank; ++i) {
            loadings_norm_[i] = std::sqrt(loadings_.col(i).dot(model.R0() * loadings_.col(i)));   // L^2 norm
            loadings_.col(i) = loadings_.col(i) / loadings_norm_[i];
        }
        scores_ = scores_.array().rowwise() * loadings_norm_.transpose().array();
        return;
    }
    // getters
    const DMatrix<double>& scores() const { return scores_; }
    const DMatrix<double>& loadings() const { return loadings_; }
    const DVector<double>& loadings_norm() const { return loadings_norm_; }
   private:
    // let E*\Sigma*F^\top the reduced (rank r) SVD of X*\Psi*(D^{1})^\top, with D^{-1} the inverse of the cholesky
    // factor of \Psi^\top * \Psi + P(\lambda), then
    DMatrix<double> scores_;          // matrix E in the reduced SVD of X*\Psi*(D^{-1})^\top
    DMatrix<double> loadings_;        // \Sigma*F^\top*D^{-1} (PC functions expansion coefficients, L^2 normalized)
    DVector<double> loadings_norm_;   // L^2 norm of estimated fields
};

}   // namespace models
}   // namespace fdapde

#endif   // __REGULARIZED_SVD_H__

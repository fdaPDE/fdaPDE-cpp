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
using fdapde::core::Grid;

#include "../../calibration/gcv.h"
#include "../../calibration/kfold_cv.h"
#include "functional_base.h"
#include "profiling_estimation.h"
using fdapde::calibration::GCV;
using fdapde::calibration::KFoldCV;

namespace fdapde {
namespace models {

// Functional Principal Components Analysis
template <typename PDE, typename RegularizationType, typename SamplingDesign, typename lambda_selection_strategy>
class FPCA : public FunctionalBase<FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>> {
    static_assert(std::is_base_of<PDEBase, PDE>::value);
   private:
    typedef FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy> Model;
    typedef FunctionalBase<Model> Base;
    std::size_t n_pc_ = 3;   // number of principal components
    // profiling_estimation parameters
    double tol_ = 1e-6;           // relative tolerance between Jnew and Jold, used as stopping criterion
    std::size_t max_iter_ = 20;   // maximum number of allowed iterations

    // problem solution
    DMatrix<double> loadings_;
    DMatrix<double> scores_;

    // tag dispatched private methods for computation of PCs, ** to be removed **
    void solve_(fixed_lambda);
    void solve_(gcv_lambda_selection);
    void solve_(kcv_lambda_selection);
   public:
    IMPORT_MODEL_SYMBOLS;
    using Base::lambda;
    using Base::lambdas;
    using Base::nan_idxs;
    using Base::X;
    // constructor
    FPCA() = default;
    // space-only constructor
    template <
      typename U = RegularizationType, typename std::enable_if<std::is_same<U, SpaceOnly>::value, int>::type = 0>
    FPCA(const PDE& pde) : Base(pde) {};
    // space-time constructor
    template <
      typename U = RegularizationType, typename std::enable_if<!std::is_same<U, SpaceOnly>::value, int>::type = 0>
    FPCA(const PDE& pde, const DVector<double>& time) : Base(pde, time) {};

    void init_model() { return; };
    virtual void solve();   // compute principal components

    // getters
    const DMatrix<double>& loadings() const { return loadings_; }
    const DMatrix<double>& scores() const { return scores_; }

    // setters
    void set_tolerance(double tol) { tol_ = tol; }
    void set_max_iter(std::size_t max_iter) { max_iter_ = max_iter; }
    void set_npc(std::size_t n_pc) { n_pc_ = n_pc; }
};

// solution in case of fixed \lambda
template <typename PDE, typename RegularizationType, typename SamplingDesign, typename lambda_selection_strategy>
void FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>::solve_(fixed_lambda) {
    ProfilingEstimation<decltype(*this)> pe(*this, tol_, max_iter_);
    BlockFrame<double, int> data_ = data();
    // Principal Components computation
    for (std::size_t i = 0; i < n_pc_; i++) {
        // find vectors s,f minimizing \norm_F{Y - s^T*f}^2 + (s^T*s)*P(f) fixed \lambda
        pe.compute(data_, lambda());
        loadings_.col(i) = pe.f_n() / pe.f_n_norm();
        scores_.col(i) = pe.s() * pe.f_n_norm();
        // subtract computed PC from data
        data_.get<double>(OBSERVATIONS_BLK) -= scores_.col(i) * loadings_.col(i).transpose();
    }
    return;
}

// best \lambda for PC choosen according to GCV index
template <typename PDE, typename RegularizationType, typename SamplingDesign, typename lambda_selection_strategy>
void FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>::solve_(gcv_lambda_selection) {
    ProfilingEstimation<decltype(*this)> pe(*this, tol_, max_iter_);
    BlockFrame<double, int> data_ = data();
    // wrap objective into a ScalarField accepted by OPT module
    const std::size_t n_lambda = n_smoothing_parameters<RegularizationType>::value;
    ScalarField<n_lambda> f;
    f = [&pe, &data_](const SVector<n_lambda>& p) -> double {
      // find vectors s,f minimizing \norm_F{Y - s^T*f}^2 + (s^T*s)*P(f) fixed \lambda = p
      pe.compute(data_, p);
      return pe.gcv(); // return GCV at convergence
    };

    Grid<n_lambda> opt;   // optimization algorithm
    // Principal Components computation
    for (std::size_t i = 0; i < n_pc_; i++) {
        opt.optimize(f, lambdas());   // select optimal \lambda for i-th PC
        // compute and store results given estimated optimal \lambda
        pe.compute(data_, opt.optimum());
        loadings_.col(i) = pe.f_n() / pe.f_n_norm();
        scores_.col(i) = pe.s() * pe.f_n_norm();
        // subtract computed PC from data
        data_.get<double>(OBSERVATIONS_BLK) -= scores_.col(i) * loadings_.col(i).transpose();
    }
    return;
}

// best \lambda for PC choosen according to K-fold CV strategy, uses the reconstruction error on test set as CV score
template <typename PDE, typename RegularizationType, typename SamplingDesign, typename lambda_selection_strategy>
void FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>::solve_(kcv_lambda_selection) {
    ProfilingEstimation<decltype(*this)> pe(*this, tol_, max_iter_);
    // number of smoothing parameters
    const std::size_t n_lambda = n_smoothing_parameters<RegularizationType>::value;

    // routine executed by the CV-engine to produce the model score
    std::function<double(DVector<double>, BlockFrame<double, int>, BlockFrame<double, int>)> cv_score =
      [this](const DVector<double>& lambda,
	     const BlockFrame<double, int>& train_df,
	     const BlockFrame<double, int>& test_df) -> double {
        ProfilingEstimation<decltype(*this)> pe_(*this, tol_, max_iter_);
        // get references to train and test sets
        const DMatrix<double>& X_test = test_df.get<double>(OBSERVATIONS_BLK);
        SVector<n_lambda> p(lambda.data());
        // find vectors s,f minimizing \norm_F{Y - s^T*f}^2 + (s^T*s)*P(f)
        pe_.compute(train_df, p);
        // compute reconstruction error and scores (X_test * f_n)/(\norm{f_n} + \lambda*P(f)) on test set
        if (this->has_nan()) {
            auto nan_pattern = X_test.array().isNaN();   // missingness pattern
            std::size_t n = nan_pattern.count();         // number of not-NaN points in test set
            // scores vector
            DVector<double> s = (nan_pattern.select(0, X_test)) * pe_.f_n();
            double pen = p[0] * (pe_.g().dot(R0() * pe_.g()));   // \lambda*P(f)
            for (std::size_t i = 0; i < s.rows(); ++i) {
                // compute \norm{f_n} for i-th subject
                double sum_f = (X_test.row(i)).array().isNaN().select(0, pe_.f_n().transpose()).squaredNorm();
                s[i] /= (sum_f + pen);   // normalize i-th score
            }
            // evaluate reconstruction error on test set
            double err = nan_pattern.select(0, X_test - s * pe_.f_n().transpose()).squaredNorm();
            return err / n;
        } else {
            // scores vector
            DVector<double> s = (X_test * pe_.f_n()) / (pe_.f_n().squaredNorm() + p[0] * (pe_.g().dot(R0() * pe_.g())));
            // evaluate reconstruction error on test set
            return (X_test - s * pe_.f_n().transpose()).squaredNorm() / X_test.size();
        }
    };

    // define K-fold algorithm
    KFoldCV cv(10);   // allow user-defined number of folds!
    std::vector<DVector<double>> lambdas_;
    lambdas_.reserve(lambdas().size());
    for (const auto& l : lambdas()) { lambdas_.emplace_back(Eigen::Map<const DVector<double>>(l.data(), n_lambda, 1)); }
    // Principal Components computation
    BlockFrame<double, int> data_ = data();
    for (std::size_t i = 0; i < n_pc_; i++) {
        cv.compute(lambdas_, data_, cv_score);   // select optimal smoothing level
        pe.compute(data_, cv.optimum());         // execute profiling estimation given estimated optimal \lambda
        loadings_.col(i) = pe.f_n() / pe.f_n_norm();
        scores_.col(i) = pe.s() * pe.f_n_norm();
        // subtract computed PC from data
        data_.get<double>(OBSERVATIONS_BLK) -= pe.s() * pe.f_n().transpose();
    }
    return;
}

// finds solution to fPCA problem, dispatch to solver depending on \lambda selection criterion
template <typename PDE, typename RegularizationType, typename SamplingDesign, typename lambda_selection_strategy>
void FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>::solve() {
    loadings_.resize(X().cols(), n_pc_);
    scores_.resize(X().rows(), n_pc_);
    // dispatch to desired solution strategy
    solve_(lambda_selection_strategy());
    return;
}

template <typename PDE_, typename SamplingDesign_, typename lambda_selection_strategy>
struct model_traits<FPCA<PDE_, SpaceOnly, SamplingDesign_, lambda_selection_strategy>> {
    typedef PDE_ PDE;
    typedef SpaceOnly regularization;
    typedef SamplingDesign_ sampling;
    typedef MonolithicSolver solver;
    enum { N = PDE::N, M = PDE::M, R = PDE::R, n_lambda = 1 };
};
// specialization for separable regularization
template <typename PDE_, typename SamplingDesign_, typename lambda_selection_strategy>
struct model_traits<FPCA<PDE_, SpaceTimeSeparable, SamplingDesign_, lambda_selection_strategy>> {
    typedef PDE_ PDE;
    typedef SpaceTimeSeparable regularization;
    typedef SplineBasis<3> TimeBasis;   // use cubic B-splines
    typedef SamplingDesign_ sampling;
    typedef MonolithicSolver solver;
    enum { N = PDE::N, M = PDE::M, R = PDE::R, n_lambda = 2 };
};

}   // namespace models
}   // namespace fdapde

#endif   // __FPCA_H__

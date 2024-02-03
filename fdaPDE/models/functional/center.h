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

#ifndef __CENTER_H__
#define __CENTER_H__

#include <fdaPDE/utils.h>

namespace fdapde {
namespace models {

// computes the smooth weighted mean field from a set of functional data (stored rowwise in a data matrix X)
// solves \argmin_{f} \| X*\frac{w}{\norm{w}_2^2} - f \|_2^2 + P_{\lambda}(f) (using a linear smoother)
template <typename SmootherType_, typename CalibrationType_>
DMatrix<double> smooth_mean(
  const DMatrix<double>& X, const DVector<double>& w, SmootherType_&& smoother, CalibrationType_&& calibration) {
    fdapde_assert(X.rows() == w.rows());
    BlockFrame<double, int> df;
    // let O_{p_i} the set of index where x_j is observed at location p_i, compute smoother data {y_i}_i
    // y_i = \sum_{j \in O_{p_i}} x_j(p_i)*w_j / \sum_{j \in O_{p_i}} w_j
    DMatrix<double> X_ = X.array().isNaN().select(0, X).transpose() * w;
    for (std::size_t i = 0; i < X.cols(); ++i) { X_(i, 0) /= X.col(i).array().isNaN().select(0, w).squaredNorm(); }
    df.insert<double>(OBSERVATIONS_BLK, X_);
    smoother.set_data(df);
    smoother.set_lambda(calibration.fit(smoother));   // find optimal smoothing parameter
    smoother.init();
    smoother.solve();
    return smoother.f();
}

// computes the smooth mean field from a set of functional data
template <typename SmootherType_, typename CalibrationType_>
DMatrix<double> smooth_mean(const DMatrix<double>& X, SmootherType_&& smoother, CalibrationType_&& calibration) {
    return smooth_mean(X, DVector<double>::Ones(X.rows()), smoother, calibration);
}

// functional centering of a data matrix X
struct CenterReturnType {
    DMatrix<double> fitted;   // centred data, X - \mu
    DMatrix<double> mean;     // mean field expansion coefficients
};
template <typename SmootherType_, typename CalibrationType_>
CenterReturnType center(const DMatrix<double>& X, SmootherType_&& smoother, CalibrationType_&& calibration) {
    DMatrix<double> mean_field = smooth_mean(X, smoother, calibration);
    // compute mean matrix and return
    return {X - smoother.fitted().replicate(1, X.rows()).transpose(), smoother.f()};
}

}   // namespace models
}   // namespace fdapde

#endif   // __CENTER_H__

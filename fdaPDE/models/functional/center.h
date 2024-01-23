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

struct CenterReturnType {
    DMatrix<double> fitted;   // the centred data, X - \mu
    DMatrix<double> mean;     // the mean field \mu
};

// functional centering of the data matrix X
template <typename SmootherType_, typename CalibrationType_>
CenterReturnType center(const DMatrix<double>& X, SmootherType_&& smoother, CalibrationType_&& calibration) {
    BlockFrame<double, int> df;
    df.insert<double>(OBSERVATIONS_BLK, X.colwise().sum().transpose() / X.rows());
    smoother.set_data(df);
    // set optimal lambda according to choosen calibration strategy
    smoother.set_lambda(calibration.fit(smoother));
    smoother.init();
    smoother.solve();
    // compute mean matrix and return
    DMatrix<double> mean = smoother.fitted().replicate(1, X.rows()).transpose();
    return {X - mean, mean};
}

// pointwise centering of the data matrix X
CenterReturnType center(const DMatrix<double>& X) {
    DMatrix<double> mean = (X.colwise().sum().transpose() / X.rows()).replicate(1, X.rows()).transpose();
    return {X - mean, mean};
}

}   // namespace models
}   // namespace fdapde

#endif   // __CENTER_H__

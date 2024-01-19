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

template <typename CalibrationType_, typename SmootherType_>
std::pair<DMatrix<double>, DMatrix<double>> center(
  const DMatrix<double>& data, SmootherType_&& smoother, CalibrationType_&& calibration,
  const std::vector<DVector<double>>& lambdas) { // lambdas must be removed and put inside calibration strategy (use expression-templates and lazy eval)
    BlockFrame<double, int> df;
    df.insert<double>(OBSERVATIONS_BLK, data.colwise().sum().transpose() / data.rows());
    smoother.set_data(df);
    // set optimal lambda according to choosen calibration strategy (currently only GCV based supported)
    smoother.set_lambda(calibration.fit(smoother, lambdas));
    smoother.init();
    smoother.solve();   // should we include init() in solve() ?? cleaner interface

    std::cout << smoother.f() << std::endl;
    
    DMatrix<double> mean = smoother.fitted().replicate(1, data.rows()).transpose();
    return std::make_pair(data - mean, mean);
}

}   // namespace models
}   // namespace fdapde

#endif   // __CENTER_H__

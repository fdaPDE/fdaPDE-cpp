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

#ifndef __CALIBRATION_OFF_H__
#define __CALIBRATION_OFF_H__

#include "../models/regression/gcv.h"
#include "calibration_base.h"

namespace fdapde {
namespace calibration {

// a class representing a fixed lambda calibration strategy (no calibration)
class Off {
   private:
    DVector<double> lambda_;
   public:
    template <typename LambdaType> Off(LambdaType&& lambda) : lambda_(lambda) {
      static_assert(std::is_base_of<Eigen::MatrixBase<LambdaType>, LambdaType>::value);
      fdapde_assert(lambda.cols() == 1 && (lambda.rows() == 1 || lambda.rows() == 2));
    }
    template <typename ModelType_> DVector<double> fit(ModelType_& model) { return lambda_; }
};

}   // namespace calibration
}   // namespace fdapde

#endif   // __CALIBRATION_OFF_H__

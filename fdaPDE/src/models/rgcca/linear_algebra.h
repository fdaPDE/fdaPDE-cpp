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

#ifndef __FDAPDE_RGCCA_LINEAR_ALGEBRA_H__
#define __FDAPDE_RGCCA_LINEAR_ALGEBRA_H__

// TODO: consider moving these generic linear algebra helpers into core internals.

#include "../header_check.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SVD>

namespace fdapde {
namespace internals {

inline void ginv(const Eigen::MatrixXd& X, Eigen::MatrixXd& ginvX, double tol = std::sqrt(std::numeric_limits<double>::epsilon())){
    Eigen::BDCSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const auto& d = svd.singularValues();
    const double d1 = d(0);
    const double thresh = std::max(tol * d1, 0.0);

    std::vector<int> idx;
    idx.reserve(d.size());
    for (int i = 0; i < d.size(); ++i) {
        if (d(i) > thresh) idx.push_back(i);
    }

    Eigen::MatrixXd Upos(X.rows(), static_cast<int>(idx.size()));
    Eigen::MatrixXd Vpos(X.cols(), static_cast<int>(idx.size()));
    Eigen::VectorXd invd(static_cast<int>(idx.size()));

    for (int k = 0; k < static_cast<int>(idx.size()); ++k) {
        Upos.col(k) = svd.matrixU().col(idx[k]);
        Vpos.col(k) = svd.matrixV().col(idx[k]);
        invd(k) = 1.0 / d(idx[k]);
    }

    ginvX = Vpos * (invd.asDiagonal() * Upos.transpose());
}

} // namespace internals
} // namespace fdapde

#endif // __FDAPDE_RGCCA_LINEAR_ALGEBRA_H__

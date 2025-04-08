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

#ifndef __INDEXES_H__
#define __INDEXES_H__

namespace fdapde {

// collection of performace indexes

template <typename Lhs, typename Rhs>
    requires(internals::is_vector_like_v<Lhs> && internals::is_vector_like_v<Rhs>)
double RMSE(const Lhs& lhs, const Rhs& rhs) {
    fdapde_assert(lhs.size() == rhs.size());
    int n = lhs.size();
    double sse = 0;
    for (int i = 0; i < n; ++i) { sse += std::pow(lhs[i] - rhs[i], 2); }
    return std::sqrt((1. / n) * sse);
}

}   // namespace fdapde

#endif   // __INDEXES_H__

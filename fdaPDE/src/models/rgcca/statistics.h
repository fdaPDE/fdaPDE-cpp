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

#ifndef __FDAPDE_RGCCA_STATISTICS_H__
#define __FDAPDE_RGCCA_STATISTICS_H__

// TODO: consider moving these generic statistics helpers into core internals.

#include "../header_check.h"
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <limits>
#include <vector>

namespace fdapde {
namespace internals {

inline double median(std::vector<double>& x) {
    if (x.empty())
        return 0.0;

    const std::size_t n = x.size();
    const std::size_t mid = n / 2;

    std::nth_element(x.begin(), x.begin() + mid, x.end());

    if (n % 2 == 1)
        return x[mid];

    const double upper = x[mid];

    std::nth_element(x.begin(), x.begin() + mid - 1, x.end());
    const double lower = x[mid - 1];

    return 0.5 * (lower + upper);
}

inline double empirical_quantile(std::vector<double>& x, const double p) {
    x.erase(
        std::remove_if(x.begin(), x.end(), [](const double v) { return !std::isfinite(v); }),
        x.end()
    );

    if (x.empty())
        return std::numeric_limits<double>::quiet_NaN();

    std::sort(x.begin(), x.end());

    if (x.size() == 1)
        return x.front();

    const double pos = std::clamp(p, 0.0, 1.0) * static_cast<double>(x.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
    const std::size_t hi = static_cast<std::size_t>(std::ceil(pos));
    const double frac = pos - static_cast<double>(lo);

    return (1.0 - frac) * x[lo] + frac * x[hi];
}

} // namespace internals
} // namespace fdapde

#endif // __FDAPDE_RGCCA_STATISTICS_H__

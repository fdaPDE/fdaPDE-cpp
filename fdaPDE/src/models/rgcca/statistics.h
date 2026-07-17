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

inline double standard_normal_quantile(const double p) {
    if (!std::isfinite(p) || p <= 0.0 || p >= 1.0)
        return std::numeric_limits<double>::quiet_NaN();

    double lower = -8.0;
    double upper = 8.0;
    constexpr int iterations = 64;
    const double inv_sqrt_two = std::sqrt(0.5);

    for (int i = 0; i < iterations; ++i) {
        const double mid = 0.5 * (lower + upper);
        const double cdf = 0.5 * std::erfc(-mid * inv_sqrt_two);
        if (cdf < p)
            lower = mid;
        else
            upper = mid;
    }

    return 0.5 * (lower + upper);
}

inline double wilson_score_upper_bound(
    const int successes,
    const int trials,
    const double z
) {
    if (trials <= 0 || successes < 0 || successes > trials || !std::isfinite(z) || z < 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    const double n = static_cast<double>(trials);
    const double p = static_cast<double>(successes) / n;
    const double z2 = z * z;
    const double denominator = 1.0 + z2 / n;
    const double center = (p + z2 / (2.0 * n)) / denominator;
    const double margin = z / denominator *
        std::sqrt(p * (1.0 - p) / n + z2 / (4.0 * n * n));

    return std::clamp(center + margin, 0.0, 1.0);
}

inline double median_confidence_upper_bound(std::vector<double>& x, const double z) {
    x.erase(
        std::remove_if(x.begin(), x.end(), [](const double v) { return !std::isfinite(v); }),
        x.end()
    );

    if (x.empty() || !std::isfinite(z) || z < 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    // Normal approximation to the binomial order-statistic interval for a median.
    const double upper_probability = 0.5 + z / (2.0 * std::sqrt(static_cast<double>(x.size())));
    return empirical_quantile(x, std::clamp(upper_probability, 0.5, 1.0));
}

} // namespace internals
} // namespace fdapde

#endif // __FDAPDE_RGCCA_STATISTICS_H__

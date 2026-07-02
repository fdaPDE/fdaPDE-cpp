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

#ifndef __FDAPDE_RGCCA_GCV_H__
#define __FDAPDE_RGCCA_GCV_H__

#include "../header_check.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace fdapde {
namespace internals {

// TODO: this GCV path is generic enough for other smoothers; move it to the shared model/solver layer.
template<class Fun> inline std::pair<double,double> argmin_over_log_grid(Fun&& f, const double log10_min, const double log10_max, int n_grid) {
    if(n_grid<2) n_grid=2;
    double best_log=log10_min;
    double best_val=std::numeric_limits<double>::infinity();
    const double step=(log10_max-log10_min)/(n_grid-1);
    for(int i=0;i<n_grid;++i) {
        const double lg = log10_min+i*step;
        double lam = std::pow(10.0,lg);
        if (const double val = f(lam); val < best_val) {
            best_val=val;
            best_log=lg;
        }
    }
    return{std::pow(10.0,best_log),best_val};
}

struct GCVConfig {
    double log10_min = -9.0;
    double log10_max = 0.0;
    int grid = 20;
    int edf_r = 100;
    int edf_seed = 12345;
    double eps_dof = 1e-12;
};

template <class Smoother> struct GCVEval {
    Smoother* s;
    GCVConfig cfg;

    double operator()(double lambda) {
        s->fit(lambda);

        const int n = s->n_obs();
        const int q = s->n_covs();
        const double trS = s->edf(cfg.edf_r, cfg.edf_seed);

        const auto& y  = s->response();
        const auto yhat = s->fn();
        const double rss = (yhat - y).squaredNorm();

        const double dor = std::max(cfg.eps_dof, static_cast<double>(n) - (static_cast<double>(q) + trS));
        return (static_cast<double>(n) / (dor * dor)) * rss;
    }
};

template <typename SolverType> std::pair<bool, double> select_lambda_with_gcv(SolverType& solver, const GCVConfig& gcv_cfg) {
    GCVEval<SolverType> gcv{ &solver, gcv_cfg };
    auto [lambda_opt, gcv_opt] = argmin_over_log_grid(
        [&](double lam){ return gcv(lam); },
        gcv_cfg.log10_min, gcv_cfg.log10_max, gcv_cfg.grid
    );
    return {lambda_opt < std::pow(10.0, gcv_cfg.log10_max), lambda_opt};
}

} // namespace internals
} // namespace fdapde

#endif // __FDAPDE_RGCCA_GCV_H__

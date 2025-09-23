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

#ifndef __DENSITY_ESTIMATION_H__
#define __DENSITY_ESTIMATION_H__

#include "header_check.h"

namespace fdapde {

template <typename VariationalSolver>
    requires(std::is_same_v<typename VariationalSolver::solver_category, de_solver>)
class DEPDE {
   private:
    using solver_t = std::decay_t<VariationalSolver>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
   public:
    static constexpr int n_lambda = solver_t::n_lambda;

    DEPDE() noexcept = default;
    template <typename GeoFrame, typename Penalty> DEPDE(const GeoFrame& gf, Penalty&& penalty) noexcept : solver_() {
        discretize(penalty.get());
        analyze_data(gf);
    }
    template <typename... Args> const vector_t& fit(Args&&... args) { return solver_.fit(std::forward<Args>(args)...); }
    template <typename... Args> void discretize(Args&&... args) { solver_.discretize(std::forward<Args>(args)...); }
    template <typename GeoFrame> void analyze_data(const GeoFrame& gf) {
        fdapde_assert(gf.n_layers() == 1);
        n_obs_ = gf[0].rows();
        solver_.analyze_data(gf);
    }
    // observers
    const vector_t& log_density() const { return solver_.log_density(); }
    vector_t density() const { return solver_.density(); }
    vector_t gn() const { return solver_.gn(); }
    vector_t fn() const { return solver_.fn(); }
    int n_obs() const { return n_obs_; }
    // modifiers
    void set_llik_tolerance(double tol) { solver_.set_llik_tolerance(tol); }
   private:
    int n_obs_ = 0;
    solver_t solver_;
};

// deduction guide
template <typename GeoFrame, typename Penalty>
DEPDE(const GeoFrame& gf, Penalty&& solver) -> DEPDE<typename std::decay_t<Penalty>::solver_t>;

}   // namespace fdapde

#endif   // __DENSITY_ESTIMATION_H__

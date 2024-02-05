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

#ifndef __FPCA_H__
#define __FPCA_H__

#include <fdaPDE/optimization.h>
#include <fdaPDE/utils.h>
#include <fdaPDE/linear_algebra.h>
#include <Eigen/SVD>

#include "../../calibration/kfold_cv.h"
#include "../../calibration/symbols.h"
using fdapde::calibration::Calibration;
#include "functional_base.h"
#include "power_iteration.h"
#include "regularized_svd.h"

namespace fdapde {
namespace models {

// Functional Principal Components Analysis
template <typename RegularizationType_>
class FPCA : public FunctionalBase<FPCA<RegularizationType_>, RegularizationType_> {
   private:
    std::size_t n_pc_ = 3;   // number of principal components

    struct SolverType__ {   // type erased solver strategy
        using This_ = FPCA<RegularizationType_>;
        template <typename T> using fn_ptrs = mem_fn_ptrs<&T::template compute<This_>, &T::loadings, &T::scores>;
        void compute(const DMatrix<double>& X, This_& model, std::size_t rank) {
            invoke<void, 0>(*this, X, model, rank);
        }
        decltype(auto) loadings() const { return invoke<const DMatrix<double>&, 1>(*this); }
        decltype(auto) scores()   const { return invoke<const DMatrix<double>&, 2>(*this); }
    };
    using SolverType = fdapde::erase<heap_storage, SolverType__>;
    SolverType solver_;
   public:
    using RegularizationType = std::decay_t<RegularizationType_>;
    using This = FPCA<RegularizationType>;
    using Base = FunctionalBase<This, RegularizationType>;
    IMPORT_MODEL_SYMBOLS;
    using Base::X;              // n_stat_units \times n_locs data matrix

    // constructors
    FPCA() = default;
    fdapde_enable_constructor_if(is_space_only, This) FPCA(const pde_ptr& pde, Sampling s, SolverType solver) :
        Base(pde, s), solver_(solver) {};
    fdapde_enable_constructor_if(is_space_time_separable, This)
      FPCA(const pde_ptr& space_penalty, const pde_ptr& time_penalty, Sampling s, SolverType solver) :
        Base(space_penalty, time_penalty, s), solver_(solver) {};

    void init_model() { fdapde_assert(bool(solver_) == true); };   // check solver is assigned
    void solve() { solver_.compute(X(), *this, n_pc_); }
    // getters
    const DMatrix<double>& loadings() const { return solver_.loadings(); }
    const DMatrix<double>& scores() const { return solver_.scores(); }
    // setters
    void set_npc(std::size_t n_pc) { n_pc_ = n_pc; }
    void set_solver(SolverType solver) { solver_ = solver; }
};

}   // namespace models
}   // namespace fdapde

#endif   // __FPCA_H__

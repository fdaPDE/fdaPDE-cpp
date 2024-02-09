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

#ifndef __REGRESSION_TYPE_ERASURE_H__
#define __REGRESSION_TYPE_ERASURE_H__

#include <fdaPDE/utils.h>
#include <fdaPDE/linear_algebra.h>
#include "../model_traits.h"
#include "../model_type_erasure.h"
#include "../sampling_design.h"
using fdapde::models::is_space_only;
using fdapde::models::Sampling;
using fdapde::models::SpaceOnly;
using fdapde::core::BinaryVector;
using fdapde::Dynamic;

namespace fdapde {
namespace models {

// type erased wrapper for regression models
struct RegressionModel__ {
    template <typename M>
    using fn_ptrs = fdapde::mem_fn_ptrs<
      &M::f, &M::beta, &M::g, &M::fitted, &M::W, &M::XtWX, &M::U, &M::V, &M::invXtWX, &M::invA, &M::q, &M::n_obs,
      &M::norm, &M::y, &M::T, &M::lmbQ, &M::has_covariates, &M::nan_mask, &M::set_mask, &M::X>;
    // interface implementation
    decltype(auto) f()       const { return invoke<const DVector<double>&   , 0>(*this); }
    decltype(auto) beta()    const { return invoke<const DVector<double>&   , 1>(*this); }
    decltype(auto) g()       const { return invoke<const DVector<double>&   , 2>(*this); }
    decltype(auto) fitted()  const { return invoke<DMatrix<double>          , 3>(*this); }
    decltype(auto) W()       const { return invoke<const DiagMatrix<double>&, 4>(*this); }
    decltype(auto) XtWX()    const { return invoke<const DMatrix<double>&   , 5>(*this); }
    decltype(auto) U()       const { return invoke<const DMatrix<double>&   , 6>(*this); }
    decltype(auto) V()       const { return invoke<const DMatrix<double>&   , 7>(*this); }
    decltype(auto) invXtWX() const { return invoke<const Eigen::PartialPivLU<DMatrix<double>>&, 8>(*this); }
    decltype(auto) invA()    const { return invoke<const fdapde::SparseLU<SpMatrix<double>>&  , 9>(*this); }
    decltype(auto) q()       const { return invoke<std::size_t, 10>(*this); }
    decltype(auto) n_obs()   const { return invoke<std::size_t, 11>(*this); }
    decltype(auto) norm(const DMatrix<double>& op1, const DMatrix<double>& op2) const {
        return invoke<double, 12> (*this, op1, op2);
    }
    decltype(auto) y() const { return invoke<const DMatrix<double>&, 13>(*this); }
    decltype(auto) T() { return invoke<const DMatrix<double>&, 14>(*this); }
    decltype(auto) lmbQ(const DMatrix<double>& x) const { return invoke<DMatrix<double>, 15>(*this, x); }
    decltype(auto) has_covariates() const { return invoke<bool, 16>(*this); }
    decltype(auto) nan_mask() const { return invoke<const BinaryVector<Dynamic>&, 17>(*this); }
    decltype(auto) set_mask(const BinaryVector<Dynamic>& mask) { return invoke<void, 18>(*this, mask); }
    decltype(auto) X() const { return invoke<const DMatrix<double>&, 19>(*this); }
};

template <typename RegularizationType>
using RegressionModel = fdapde::erase<fdapde::heap_storage, StatisticalModel__<RegularizationType>, RegressionModel__>;
template <typename RegularizationType>
using RegressionView =
  fdapde::erase<fdapde::non_owning_storage, StatisticalModel__<RegularizationType>, RegressionModel__>;

}   // namespace models
}   // namespace fdapde

#endif   // __REGRESSION_TYPE_ERASURE_H__

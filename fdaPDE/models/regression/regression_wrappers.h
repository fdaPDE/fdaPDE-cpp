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

#ifndef __REGRESSION_WRAPPERS_H__
#define __REGRESSION_WRAPPERS_H__

#include <fdaPDE/utils.h>
#include "../model_traits.h"
#include "../model_wrappers.h"
#include "../sampling_design.h"
using fdapde::models::SpaceOnly;
using fdapde::models::is_space_only;
using fdapde::models::Sampling;

namespace fdapde {
namespace models {
  
// type erased wrapper for regression models
struct IRegression {
    template <typename M>
    using fn_ptrs = fdapde::mem_fn_ptrs<
      &M::f, &M::beta, &M::g, &M::fitted, &M::W, &M::XtWX, &M::U, &M::V, &M::invXtWX, &M::invA, &M::q, &M::n_obs,
      &M::norm, &M::y, &M::T, &M::lmbQ, &M::has_covariates, &M::P1, &M::R0, &M::nan_idxs>;
    // interface implementation
    decltype(auto) f()       const { return fdapde::invoke<const DVector<double>&   , 0>(*this); }
    decltype(auto) beta()    const { return fdapde::invoke<const DVector<double>&   , 1>(*this); }
    decltype(auto) g()       const { return fdapde::invoke<const DVector<double>&   , 2>(*this); }
    decltype(auto) fitted()  const { return fdapde::invoke<DMatrix<double>          , 3>(*this); }
    decltype(auto) W()       const { return fdapde::invoke<const DiagMatrix<double>&, 4>(*this); }
    decltype(auto) XtWX()    const { return fdapde::invoke<const DMatrix<double>&   , 5>(*this); }
    decltype(auto) U()       const { return fdapde::invoke<const DMatrix<double>&   , 6>(*this); }
    decltype(auto) V()       const { return fdapde::invoke<const DMatrix<double>&   , 7>(*this); }
    decltype(auto) invXtWX() const { return fdapde::invoke<const Eigen::PartialPivLU<DMatrix<double>>&, 8>(*this); }
    decltype(auto) invA()    const { return fdapde::invoke<const fdapde::SparseLU<SpMatrix<double>>&  , 9>(*this); }
    decltype(auto) q()       const { return fdapde::invoke<std::size_t, 10>(*this); }
    decltype(auto) n_obs()   const { return fdapde::invoke<std::size_t, 11>(*this); }
    decltype(auto) norm(const DMatrix<double>& op1, const DMatrix<double>& op2) const {
        return fdapde::invoke<double, 12> (*this, op1, op2);
    }
    decltype(auto) y() const { return fdapde::invoke<const DMatrix<double>&, 13>(*this); }
    decltype(auto) T() { return fdapde::invoke<const DMatrix<double>&, 14>(*this); }
    decltype(auto) lmbQ(const DMatrix<double>& x) const { return fdapde::invoke<DMatrix<double>, 15>(*this, x); }
    decltype(auto) has_covariates() const { return fdapde::invoke<bool, 16>(*this); }
    // M
    decltype(auto) P1() const {return fdapde::invoke<const SpMatrix<double>&, 17>(*this);  }
    decltype(auto) R0() const {return fdapde::invoke<const SpMatrix<double>&, 18>(*this);  }
    decltype(auto) nan_idxs() const {return fdapde::invoke<const SpMatrix<double>&, 19>(*this);  }
}; 

template <typename RegularizationType>
using RegressionModel = fdapde::erase<fdapde::heap_storage, IStatModel<RegularizationType>, IRegression>;
template <typename RegularizationType>
using RegressionView  = fdapde::erase<fdapde::non_owning_storage, IStatModel<RegularizationType>, IRegression>;

}   // namespace models

template <typename Model, typename PDE>
typename std::enable_if< is_space_only<Model>::value, models::RegressionModel<SpaceOnly>>::type
make_model(PDE&& args, Sampling s) {
    return models::RegressionModel<SpaceOnly>(Model(std::forward<PDE>(args), s));
}
template <typename Model, typename PDE>
typename std::enable_if<!is_space_only<Model>::value, models::RegressionModel<typename Model::RegularizationType>>::type
make_model(PDE&& args, Sampling s, const DVector<double>& time) {
    return models::RegressionModel<typename Model::RegularizationType>(Model(std::forward<PDE>(args), s, time));
}

}   // namespace fdapde

#endif   // __REGRESSION_WRAPPERS_H__

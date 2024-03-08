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

#ifndef __FUNCTIONAL_BASE_H__
#define __FUNCTIONAL_BASE_H__

#include <fdaPDE/linear_algebra.h>
#include <fdaPDE/utils.h>

#include "../model_macros.h"
#include "../model_traits.h"
#include "../sampling_design.h"
#include "../space_only_base.h"
#include "../space_time_base.h"
#include "../space_time_parabolic_base.h"
#include "../space_time_separable_base.h"

namespace fdapde {
namespace models {

// base class for any *functional* model
template <typename Model, typename RegularizationType>
class FunctionalBase : public select_regularization_base<Model, RegularizationType>::type, public SamplingBase<Model> {
   public:
    using Base = typename select_regularization_base<Model, RegularizationType>::type;
    using Base::df_;                    // BlockFrame for problem's data storage
    using SamplingBase<Model>::Psi;     // matrix of basis evaluations at locations p_1 ... p_n
    using SamplingBase<Model>::PsiTD;   // block \Psi^T*D

    FunctionalBase() = default;
    // space-only and space-time parabolic constructor (they require only one PDE)
    FunctionalBase(const Base::PDE& pde, Sampling s)
        requires(is_space_only<Model>::value || is_space_time_parabolic<Model>::value)
        : Base(pde), SamplingBase<Model>(s) {};
    // space-time separable constructor
    FunctionalBase(const Base::PDE& space_penalty, const Base::PDE& time_penalty, Sampling s)
        requires(is_space_time_separable<Model>::value)
        : Base(space_penalty, time_penalty), SamplingBase<Model>(s) {};

    // getters
    const DMatrix<double>& X() const { return df_.template get<double>(OBSERVATIONS_BLK); }   // observation matrix X
    std::size_t n_stat_units() const { return X().rows(); }
    std::size_t n_obs() const { return X().size(); };
    const SpMatrix<double>& Psi() const { return Psi(not_nan()); }
    const SpMatrix<double>& PsiTD() const { return PsiTD(not_nan()); }

    // initialization methods
    void analyze_data() { return; }
    void correct_psi() { return; }
};

}   // namespace models
}   // namespace fdapde

#endif   // __FUNCTIONAL_BASE_H__

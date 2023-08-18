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

#include <algorithm>

#include <fdaPDE/utils.h>
#include "../model_base.h"
#include "../model_traits.h"
#include "../sampling_design.h"

namespace fdapde {
namespace models {

// ** to be moved **
struct fixed_lambda { };
struct gcv_lambda_selection { };
struct kcv_lambda_selection { };

// base class for any *functional* fdaPDE model
template <typename Model>
class FunctionalBase :
    public select_regularization_type<Model>::type,
    public SamplingDesign<Model, typename model_traits<Model>::sampling> {
   protected:
    // vector of smoothing parameters
    std::vector<SVector<model_traits<Model>::n_lambda>> lambdas_;

    // quantites related to missing data
    std::vector<std::unordered_set<std::size_t>> nan_idxs_;   // missing observations indexes for each statistical unit
    std::vector<SpMatrix<double>> B_;                         // i-th statistical unit \Psi matrix
   public:
    typedef typename model_traits<Model>::PDE PDE;   // PDE used for regularization in space
    typedef typename select_regularization_type<Model>::type Base;
    typedef SamplingDesign<Model, typename model_traits<Model>::sampling> SamplingBase;
    using Base::df_;             // BlockFrame for problem's data storage
    using Base::n_locs;          // number of spatial (spatio-temporal) data locations
    using Base::pde_;            // differential operator L
    using SamplingBase::D;       // matrix of subdomains measures (for areal sampling)
    using SamplingBase::Psi;     // matrix of basis evaluations at locations p_1 ... p_n
    using SamplingBase::PsiTD;   // block \Psi^T*D

    FunctionalBase() = default;
    // space-only constructor
    template <
      typename U = Model,   // fake type to enable substitution in SFINAE
      typename std::enable_if<std::is_same<typename model_traits<U>::regularization, SpaceOnly>::value, int>::type = 0>
    FunctionalBase(const PDE& pde) :
        select_regularization_type<Model>::type(pde),
        SamplingDesign<Model, typename model_traits<Model>::sampling>() {};
    // space-time constructor
    template <
      typename U = Model,   // fake type to enable substitution in SFINAE
      typename std::enable_if<!std::is_same<typename model_traits<U>::regularization, SpaceOnly>::value, int>::type = 0>
    FunctionalBase(const PDE& pde, const DVector<double>& time) :
        select_regularization_type<Model>::type(pde, time),
        SamplingDesign<Model, typename model_traits<Model>::sampling>() {};
    // copy constructor, copy only pde object (as a consequence also the problem domain)
    FunctionalBase(const FunctionalBase& rhs) { pde_ = rhs.pde_; }

    // getters
    const DMatrix<double>& X() const { return df_.template get<double>(OBSERVATIONS_BLK); }   // observation matrix y
    std::size_t n_stat_units() const { return X().rows(); }
    const std::vector<std::unordered_set<std::size_t>>& nan_idxs() const { return nan_idxs_; }   // missing data indexes
    std::size_t n_nan() const {   // total number of missing data points
        return std::accumulate(
          nan_idxs_.begin(), nan_idxs_.end(), 0,
          [](std::size_t n, const std::unordered_set<std::size_t>& stat_unit_nan) { return n + stat_unit_nan.size(); });
    }
    std::size_t n_obs() const { return X().size() - n_nan(); };   // number of not missing observations
    // access to NaN corrected \Psi and \Psi^T*D matrices of i-th unit
    const SpMatrix<double>& Psi(std::size_t i) const { return has_nan(i) ? B_[i] : Psi(not_nan()); }
    auto PsiTD(std::size_t i) const { return has_nan(i) ? B_[i].transpose() * D() : Psi(not_nan()).transpose() * D(); }

    // accepts a collection of \lambda parameters if a not fixed_lambda method is selected ** not stable **
    void set_lambda(const std::vector<SVector<model_traits<Model>::n_lambda>>& lambdas) { lambdas_ = lambdas; }
    const std::vector<SVector<model_traits<Model>::n_lambda>>& lambdas() const { return lambdas_; }

    // utilities
    bool has_nan(std::size_t i) const { return !nan_idxs_[i].empty(); }   // true if the i-th unit has missing data
    bool has_nan() const;   // true if any of the statistical unit has missing data

    // initialization methods
    void update_data() { return; }
    void init_nan();   // functional models' missing data logic (called by SamplingBase::init_sampling())
};

// implementative details
  
// missing data logic
template <typename Model> void FunctionalBase<Model>::init_nan() {
    nan_idxs_.clear();   // clean previous missingness structure
    nan_idxs_.resize(n_stat_units());
    B_.resize(n_stat_units());
    // \Psi matrix dimensions
    std::size_t n = Psi(not_nan()).rows();
    std::size_t N = Psi(not_nan()).cols();
    // for i-th statistical unit, analyze missingness structure and set \Psi_i
    for (std::size_t i = 0; i < n_stat_units(); ++i) {
        // derive missingness pattern for i-th statistical unit
        for (std::size_t j = 0; j < n_locs(); ++j) {
            if (std::isnan(X()(i, j)))   // requires -ffast-math compiler flag to be disabled
                nan_idxs_[i].insert(j);
        }

        // NaN detected for this unit, start assembly
        if (!nan_idxs_[i].empty()) {
            for (std::size_t i = 0; i < n_stat_units(); ++i) {
                B_[i].resize(n, N);   // reserve space
                std::vector<fdapde::Triplet<double>> tripletList;
                tripletList.reserve(n * N);
                for (int k = 0; k < Psi(not_nan()).outerSize(); ++k)
                    for (SpMatrix<double>::InnerIterator it(Psi(not_nan()), k); it; ++it) {
                        if (nan_idxs_[i].find(it.row()) == nan_idxs_[i].end())
                            // no missing data at this location for i-th statistical unit
                            tripletList.emplace_back(it.row(), it.col(), it.value());
                    }
                // finalize construction
                B_[i].setFromTriplets(tripletList.begin(), tripletList.end());
                B_[i].makeCompressed();
            }
        }
        // otherwise no matrix is assembled, full \Psi is returned by Psi(std::size_t) getter
    }
    return;
}

// true if there are missing data in any of the statistical unit
template <typename Model> bool FunctionalBase<Model>::has_nan() const {
    for (auto s : nan_idxs_) {
        if (!s.empty()) return true;
    }
    return false;
}

}   // namespace models
}   // namespace fdapde

#endif   // __FUNCTIONAL_BASE_H__

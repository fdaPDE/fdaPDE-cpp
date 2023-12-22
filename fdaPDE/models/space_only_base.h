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

#ifndef __SPACE_ONLY_BASE_H__
#define __SPACE_ONLY_BASE_H__

#include <fdaPDE/utils.h>
#include "model_base.h"

namespace fdapde {
namespace models {

// abstract base interface for any *space-only* fdaPDE statistical model.
template <typename Model> class SpaceOnlyBase : public ModelBase<Model> {
   protected:
    typedef ModelBase<Model> Base;
    using Base::pde_;      // regularizing PDE
    using Base::model;     // underlying model object
    static constexpr int n_lambda = n_smoothing_parameters<SpaceOnly>::value;
  
    SpMatrix<double> P_;   // discretization of penalty term: R1^T*R0^{-1}*R1
    SVector<n_lambda> lambda_ = SVector<n_lambda>::Zero();
   public:
    using Base::lambda;       // dynamic sized smoothing parameter vector
    using Base::set_lambda;   // dynamic sized setter for \lambda
    // constructor
    SpaceOnlyBase() = default;
    SpaceOnlyBase(const pde_ptr& pde) : ModelBase<Model>(pde) {};
    void init_regularization() { return; }   // do nothing

    // setters
    void set_lambda(const SVector<n_lambda>& lambda) {
        if(lambda_ == lambda) return;
        model().runtime().set(runtime_status::is_lambda_changed);
        lambda_ = lambda;
    }
    void set_lambda_D(double lambda_D) { set_lambda(SVector<n_lambda>(lambda_D)); }
    // getters
    SVector<n_lambda> lambda() const { return lambda_; }
    double lambda_D() const { return lambda_[0]; }                 // smoothing parameter
    const SpMatrix<double>& R0() const { return pde_.mass(); }     // mass matrix in space
    const SpMatrix<double>& R1() const { return pde_.stiff(); }    // discretization of differential operator L
    const DMatrix<double>& u() const { return pde_.force(); }      // discretization of forcing term u
    inline std::size_t n_temporal_locs() const { return 1; }       // number of time instants
    std::size_t n_basis() const { return pde_.n_dofs(); };         // number of basis functions
    std::size_t n_spatial_basis() const { return n_basis(); }      // number of basis functions in space

    // computes and cache R1^T*R0^{-1}*R1. Returns the discretized penalty P = \lambda_D*(R1^T*R0^{-1}*R1)
    auto P() {
        if (is_empty(P_)) {
            fdapde::SparseLU<SpMatrix<double>> invR0_;
            invR0_.compute(pde_.mass());
            P_ = R1().transpose() * invR0_.solve(R1());   // R1^T*R0^{-1}*R1
        }
        return lambda_D() * P_;
    }
  
    // destructor
    virtual ~SpaceOnlyBase() = default;
};

}   // namespace models
}   // namespace fdapde

#endif   // __SPACE_ONLY_BASE_H__

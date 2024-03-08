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
#include <fdaPDE/pde.h>
#include <fdaPDE/linear_algebra.h>
#include "model_base.h"
using fdapde::core::lump;

namespace fdapde {
namespace models {

// abstract base interface for any *space-only* fdaPDE statistical model.
template <typename Model> class SpaceOnlyBase : public ModelBase<Model> {
   public:
    using PDE = erase<heap_storage, core::PDE__>;
    using Base = ModelBase<Model>;
    static constexpr int n_lambda = n_smoothing_parameters<SpaceOnly>::value;
    using Base::lambda;       // dynamic sized smoothing parameter vector
    using Base::model;        // underlying model object
    using Base::set_lambda;   // dynamic sized setter for \lambda
    // constructor
    SpaceOnlyBase() = default;
    SpaceOnlyBase(const PDE& space_penalty) : pde_(space_penalty) {};
    void init_regularization() {
        pde_.init();
        if (mass_lumping) { R0_lumped_ = lump(pde_.mass()); }   // lump mass matrix if requested
    }
    // public flags
    bool mass_lumping = false;
    // setters
    void set_lambda(const SVector<n_lambda>& lambda) {
        if(lambda_ == lambda) return;
        model().runtime().set(runtime_status::is_lambda_changed);
        lambda_ = lambda;
    }
    void set_lambda_D(double lambda_D) { set_lambda(SVector<n_lambda>(lambda_D)); }
    void set_penalty(const PDE& pde) {
        pde_ = pde;
        model().runtime().set(runtime_status::require_penalty_init);
    }
    // getters
    SVector<n_lambda> lambda() const { return lambda_; }
    double lambda_D() const { return lambda_[0]; }
    const SpMatrix<double>& R0() const { return mass_lumping ? R0_lumped_ : pde_.mass(); }   // mass matrix
    const SpMatrix<double>& R1() const { return pde_.stiff(); }    // discretization of differential operator L
    const DMatrix<double>& u() const { return pde_.force(); }      // discretization of forcing term u
    inline int n_temporal_locs() const { return 1; }               // number of time instants
    int n_basis() const { return pde_.n_dofs(); };                 // number of basis functions
    int n_spatial_basis() const { return n_basis(); }              // number of basis functions in space
    const PDE& pde() const { return pde_; }                        // regularizing term Lf - u
    const fdapde::SparseLU<SpMatrix<double>>& invR0() const {      // LU factorization of mass matrix R0
        if (!invR0_) { invR0_.compute(R0()); }
        return invR0_;
    }
    const SpMatrix<double>& PD() const {   // space-penalty component (R1^T*R0^{-1}*R1)
        if (is_empty(P_)) { P_ = R1().transpose() * invR0().solve(R1()); }
        return P_;
    }
    // evaluation of penalty term \lambda*(R1^\top*R0^{-1}*R1) at \lambda
    auto P(const SVector<n_lambda>& lambda) const { return lambda[0] * PD(); }
    auto P() const { return P(lambda()); }
    // destructor
    virtual ~SpaceOnlyBase() = default;
   protected:
    PDE pde_ {};                   // differential penalty in space Lf - u
    mutable SpMatrix<double> P_;   // discretization of penalty term: R1^T*R0^{-1}*R1
    mutable fdapde::SparseLU<SpMatrix<double>> invR0_;
    SpMatrix<double> R0_lumped_;   // lumped mass matrix, if mass_lumping == true, empty otherwise
    SVector<n_lambda> lambda_ = SVector<n_lambda>::Zero();
};

}   // namespace models
}   // namespace fdapde

#endif   // __SPACE_ONLY_BASE_H__

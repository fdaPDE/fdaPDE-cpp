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

#ifndef __SAMPLING_DESIGN_H__
#define __SAMPLING_DESIGN_H__

#include <fdaPDE/linear_algebra.h>
#include <fdaPDE/mesh.h>
#include <fdaPDE/utils.h>
using fdapde::core::ct_nnodes;
using fdapde::core::Kronecker;
using fdapde::core::PointLocationStrategy;
using fdapde::core::PointLocator;

#include "model_base.h"
#include "model_macros.h"
#include "model_traits.h"

namespace fdapde {
namespace models {

// base classes for the implemetation of the different sampling designs.
// Here is computed the matrix of spatial basis evaluations \Psi = [\Psi]_{ij} = \psi_i(p_j)
template <typename Model, typename S> class SamplingDesign { };

// tag to request the not-NaN corrected version of matrix \Psi
struct not_nan { };

// base class for all sampling strategies implementing common operations on \Psi matrix
template <typename Model> class SamplingBase {
   protected:
    FDAPDE_DEFINE_MODEL_GETTER;   // import model() method (const and non-const access)
    SpMatrix<double> Psi_;        // n x N matrix \Psi = [\psi_{ij}] = \psi_j(p_i)
   public:
    // if the model is space-time, perform a proper tensorization of matrix \Psi
    void tensorize() {
        if constexpr (is_solver_monolithic<Model>::value) {
            if constexpr (is_space_time_separable<Model>::value) Psi_ = Kronecker(model().Phi(), Psi_);
            if constexpr (is_space_time_parabolic<Model>::value) {
                SpMatrix<double> Im(model().n_temporal_locs(), model().n_temporal_locs());   // m x m identity matrix
                Im.setIdentity();
                Psi_ = Kronecker(Im, Psi_);
            }
        }
    }
    // getters to not NaN corrected \Psi and \Psi^T*D matrices (\Psi^T*D redefined for Areal sampling)
    const SpMatrix<double>& Psi(not_nan) const { return Psi_; }
    auto PsiTD(not_nan) const { return Psi_.transpose(); }
};

// data sampled at mesh nodes
template <typename Model> class SamplingDesign<Model, GeoStatMeshNodes> : public SamplingBase<Model> {
   private:
    FDAPDE_DEFINE_MODEL_GETTER;   // import model() method (const and non-const access)
    typedef SamplingBase<Model> Base;
    using Base::Psi_;
    using Base::tensorize;   // tensorize matrix \Psi for space-time problems
   public:
    // constructor
    SamplingDesign() = default;
    // init sampling data structures
    void init_sampling(bool forced = false) {
        // compute once if not forced to recompute
        if (!is_empty(Psi_) && forced == false) return;
        // preallocate space for Psi matrix
        std::size_t n = model().n_spatial_basis();
        std::size_t N = model().n_spatial_basis();
        Psi_.resize(n, N);
        std::vector<fdapde::Triplet<double>> triplet_list;
        triplet_list.reserve(n);

        // if data locations are equal to mesh nodes then \Psi is the identity matrix.
        // \psi_i(p_i) = 1 and \psi_i(p_j) = 0 \forall i \neq j
        for (std::size_t i = 0; i < n; ++i) triplet_list.emplace_back(i, i, 1.0);
        // finalize construction
        Psi_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        Psi_.makeCompressed();
        tensorize();   // tensorize \Psi for space-time problems
    }

    // getters
    std::size_t n_spatial_locs() const { return model().domain().dof(); }
    auto D() const { return DVector<double>::Ones(Psi_.rows()).asDiagonal(); }
    DMatrix<double> locs() const { return model().domain().dof_coords(); }
    // set locations (nothing to do, locations are implicitly set to mesh nodes)
    template <typename Derived> void set_spatial_locations(const DMatrix<Derived>& locs) { return; }
};

// data sampled at general locations p_1, p_2, ... p_n
template <typename Model> class SamplingDesign<Model, GeoStatLocations> : public SamplingBase<Model> {
   private:
    static constexpr std::size_t M = model_traits<Model>::M, N = model_traits<Model>::N, R = model_traits<Model>::R;
    FDAPDE_DEFINE_MODEL_GETTER;   // import model() method (const and non-const access)
    DMatrix<double> locs_;        // matrix of spatial locations p_1, p2_, ... p_n
    // point location strategy over triangulation
    PointLocationStrategy point_location_strategy_;
    std::shared_ptr<PointLocator<M, N, R>> point_locator_ = nullptr;
    typedef SamplingBase<Model> Base;
    using Base::Psi_;
    using Base::tensorize;   // tensorize matrix \Psi for space-time problems
   public:
    // constructor
    SamplingDesign() = default;
    // init sampling data structures
    void init_sampling(bool forced = false) {
        fdapde_assert(locs_.size() != 0);   // enable custom message
        // compute once if not forced to recompute
        if (!is_empty(Psi_) && forced == false) return;
        // set-up point location strategy
        switch (point_location_strategy_) {   // strategy pattern
        case PointLocationStrategy::naive_search:
            point_locator_ = std::make_shared<core::NaiveSearch<M, N, R>>(model().domain());
            break;
        case PointLocationStrategy::barycentric_walk:
            point_locator_ = std::make_shared<core::BarycentricWalk<M, N, R>>(model().domain());
            break;
        case PointLocationStrategy::tree_search:
            point_locator_ = std::make_shared<core::ADT<M, N, R>>(model().domain());
            break;
        default:   // default to tree search
            point_locator_ = std::make_shared<core::ADT<M, N, R>>(model().domain());
        }
        // preallocate space for Psi matrix
        std::size_t n_locs = locs_.rows();
        std::size_t n_basis = model().n_spatial_basis();
        Psi_.resize(n_locs, n_basis);
        std::vector<fdapde::Triplet<double>> triplet_list;
        triplet_list.reserve(n_locs * ct_nnodes(M, R));

        // cycle over all locations
        for (std::size_t i = 0; i < n_locs; ++i) {
            SVector<N> p_i(locs_.row(i));
            // search element containing the point
            auto e = point_locator_->locate(p_i);
            // update \Psi matrix
            for (std::size_t j = 0; j < ct_nnodes(M, R); ++j) {
                std::size_t h = e->node_ids()[j];   // column index of \Psi matrix
                // extract \psi_h from basis and evaluate at p_i
                auto psi_h = model().pde().basis()(*e, j);
                triplet_list.emplace_back(i, h, psi_h(p_i));
            }
        }
        // finalize construction
        Psi_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        Psi_.makeCompressed();
        tensorize();   // tensorize \Psi for space-time problems
    };

    // getters
    std::size_t n_spatial_locs() const { return locs_.rows(); }
    auto D() const { return DVector<double>::Ones(Psi_.rows()).asDiagonal(); }
    const DMatrix<double>& locs() const { return locs_; }
    // setter
    void set_spatial_locations(const DMatrix<double>& locs) { locs_ = locs; }
    void set_point_location_strategy(PointLocationStrategy strategy) { point_location_strategy_ = strategy; }
};

// data sampled at subdomains D_1, D_2, ... D_d
template <typename Model> class SamplingDesign<Model, Areal> : public SamplingBase<Model> {
   private:
    FDAPDE_DEFINE_MODEL_GETTER;   // import model() method (const and non-const access)
    DMatrix<int> subdomains_;     // incidence matrix D = [D]_{ij} = 1 \iff element j belongs to subdomain i.
    DiagMatrix<double> D_;        // diagonal matrix of subdomains' measures
    typedef SamplingBase<Model> Base;
    using Base::Psi_;
    using Base::tensorize;   // tensorize matrix \Psi for space-time problems
   public:
    // constructor
    SamplingDesign() = default;
    // init sampling data structures
    void init_sampling(bool forced = false) {
        fdapde_assert(subdomains_.size() != 0);   // allow for custom message
        // compute once if not forced to recompute
        if (!is_empty(Psi_) && forced == false) return;
        // preallocate space for Psi matrix
        std::size_t n = subdomains_.rows();
        std::size_t N = model().n_spatial_basis();
        Psi_.resize(n, N);
        std::vector<fdapde::Triplet<double>> triplet_list;
        triplet_list.reserve(n * N);   // n is small, should not cause any bad_alloc

        DVector<double> D;   // store measure of subdomains, this will be ported to a diagonal matrix at the end
        D.resize(subdomains_.rows());
        // start construction of \Psi matrix
        std::size_t tail = 0;
        for (std::size_t k = 0; k < n; ++k) {
            std::size_t head = 0;
            double Di = 0;   // measure of subdomain D_i
            for (std::size_t l = 0; l < subdomains_.cols(); ++l) {
                if (subdomains_(k, l) == 1) {   // element with ID l belongs to k-th subdomain
                    // get element with this ID
                    auto e = model().domain().element(l);
                    // compute \int_e \psi_h \forall \psi_h defined on e
		    for (std::size_t j = 0; j < ct_nnodes(Model::M, Model::K); ++j) {
		      std::size_t h = e.node_ids()[j];   // column index of \Psi matrix
		      auto psi_h = model().pde().basis()(e, j);
		      triplet_list.emplace_back(k, h, model().pde().integrator().integrate(e, psi_h));
		      head++;
                    }
                    Di += e.measure();   // update measure of subdomain D_i
                }
            }
            // divide each \int_{D_i} \psi_j by the measure of subdomain D_i
            for (std::size_t j = 0; j < head; ++j) { triplet_list[tail + j].value() /= Di; }
            D[k] = Di;   // store measure of subdomain
            tail += head;
        }
        // here we must be carefull of the type of model (space-only or space-time) we are handling
        if constexpr (is_space_time<Model>::value) {
            // store I_m \kron D
            std::size_t m = model().n_temporal_locs();
            std::size_t n = n_spatial_locs();
            DVector<double> IkronD(n * m);
            for (std::size_t i = 0; i < m; ++i) IkronD.segment(i * n, n) = D;
            // compute and store result
            D_ = IkronD.asDiagonal();
        } else {
            // for space-only problems store diagonal matrix D_ = diag(D_1, D_2, ... ,D_d) as it is
            D_ = D.asDiagonal();
        }
        // finalize construction
        Psi_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        Psi_.makeCompressed();
        tensorize();   // tensorize \Psi for space-time problems
    };

    // getters
    auto PsiTD(not_nan) const { return Psi_.transpose() * D_; }
    std::size_t n_spatial_locs() const { return subdomains_.rows(); }
    const DiagMatrix<double>& D() const { return D_; }
    const DMatrix<int>& locs() const { return subdomains_; }
    // setter
    void set_spatial_locations(const DMatrix<int>& subdomains) { subdomains_ = subdomains; }
};

}   // namespace models
}   // namespace fdapde

#endif   // __SAMPLING_DESIGN_H__

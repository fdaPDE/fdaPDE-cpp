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

#ifndef __MODEL_BASE_H__
#define __MODEL_BASE_H__

#include <memory>
#include <fdaPDE/mesh.h>
#include <fdaPDE/utils.h>
using fdapde::core::BlockFrame;
using fdapde::core::Mesh;

#include "model_macros.h"
#include "model_traits.h"

namespace fdapde {
namespace models {

// abstract base interface for any fdaPDE statistical model. Uses CRTP pattern
template <typename Model> class ModelBase {
   public:
    typedef typename model_traits<Model>::PDE PDE;   // differential regularization
    typedef typename PDE::DomainType DomainType;
    static constexpr std::size_t M = PDE::M;   // tangent space dimension
    static constexpr std::size_t N = PDE::N;   // embedding space dimension

    // constructors
    ModelBase() = default;
    ModelBase(const PDE& pde) : pde_(std::make_shared<PDE>(pde)) {};

    void init_pde() { pde_->init(); };   // entry point for pde initialization
    void init();                         // entry point for full model stack initialization
    virtual void solve() = 0;            // finds a solution to the problem, whatever the problem is.

    // setters
    void set_dirichlet_bc(SpMatrix<double>& A, DMatrix<double>& b);
    void set_data(const BlockFrame<double, int>& df, bool reindex = false);
    void set_lambda(const SVector<model_traits<Model>::n_lambda>& lambda) { lambda_ = lambda; }
    void set_pde(const PDE& pde) { pde_ = std::make_shared<PDE>(pde); }
  
    // getters
    const BlockFrame<double, int>& data() const { return df_; }
    BlockFrame<double, int>& data() { return df_; }   // direct write-access to model's internal data storage
    const DMatrix<int>& idx() const { return df_.get<int>(INDEXES_BLK); }   // data indices
    const PDE& pde() const { return *pde_; }   // regularizing term Lf - u (defined on some domain \Omega)
    const DomainType& domain() const { return pde_->domain(); }
    std::size_t n_locs() const { return model().n_spatial_locs() * model().n_temporal_locs(); }
    SVector<model_traits<Model>::n_lambda> lambda() const { return lambda_; }
  
    virtual ~ModelBase() = default;
   protected:
    std::shared_ptr<PDE> pde_ {};     // regularizing term Lf - u and domain definition D
    BlockFrame<double, int> df_ {};   // blockframe for data storage
    SVector<model_traits<Model>::n_lambda> lambda_ = SVector<model_traits<Model>::n_lambda>::Zero();

    // getter to underlying model object
    inline Model& model() { return static_cast<Model&>(*this); }
    inline const Model& model() const { return static_cast<const Model&>(*this); }
};

// implementative details

// perform initialization of model object, must be called before call to .solve()
template <typename Model> void ModelBase<Model>::init() {
    init_pde();                      // init differential regularization
    model().init_regularization();   // init regularization term
    model().init_sampling(true);     // init \Psi matrix, always force recomputation
    model().init_nan();              // analyze and set missingness patternl; 
    model().init_model();            // specific model's initialization
}

// set model's data from blockframe
template <typename Model> void ModelBase<Model>::set_data(const BlockFrame<double, int>& df, bool reindex) {
    df_ = df;
    // insert an index row (if not yet present or requested)
    if (!df_.has_block(INDEXES_BLK) || reindex) {
        std::size_t n = df_.rows();
        DMatrix<int> idx(n, 1);
        for (std::size_t i = 0; i < n; ++i) idx(i, 0) = i;
        df_.insert(INDEXES_BLK, idx);
    }
    model().update_data();   // specific data initialization requested by the model
}

// set boundary conditions on problem's linear system
// BUG: not working - fix needed due to SparseBlockMatrix interface
template <typename Model> void ModelBase<Model>::set_dirichlet_bc(SpMatrix<double>& A, DMatrix<double>& b) {
    std::size_t n = A.rows() / 2;

    for (std::size_t i = 0; i < n; ++i) {
        if (pde_->domain().is_on_boundary(i)) {
            A.row(i) *= 0;          // zero all entries of this row
            A.coeffRef(i, i) = 1;   // set diagonal element to 1 to impose equation u_j = b_j

            A.row(i + n) *= 0;
            A.coeffRef(i + n, i + n) = 1;

            // boundaryDatum is a pair (nodeID, boundary value)
            double boundaryDatum = pde_->boundaryData().empty() ? 0 : pde_->boundaryData().at(i)[0];
            b.coeffRef(i, 0) = boundaryDatum;   // impose boundary value
            b.coeffRef(i + n, 0) = 0;
        }
    }
    return;
}

}   // namespace models
}   // namespace fdapde

#endif   // __MODEL_BASE_H__

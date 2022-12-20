#ifndef __MODEL_BASE_H__
#define __MODEL_BASE_H__

#include <memory>
#include <type_traits>

#include "../core/utils/Symbols.h"
#include "../core/MESH/Mesh.h"
using fdaPDE::core::MESH::Mesh;
#include "../core/MESH/engines/AlternatingDigitalTree/ADT.h"
using fdaPDE::core::MESH::ADT;
#include "../core/utils/DataStructures/BlockFrame.h"
#include "ModelTraits.h"
using fdaPDE::models::model_traits;

namespace fdaPDE {
namespace models {
   
  // abstract base interface for any fdaPDE statistical model. Uses CRTP pattern
  template <typename Model>
  class ModelBase {
  public:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    static constexpr std::size_t M = PDE::local_dimension;
    static constexpr std::size_t N = PDE::embedding_dimension;
    static constexpr std::size_t K = PDE::basis_order;

    // constructor
    ModelBase() = default;
    ModelBase(const PDE& pde) : pde_(std::make_shared<PDE>(pde)) {};
    // copy constructor
    ModelBase(const ModelBase& rhs) { pde_ = rhs.pde_; }
    
    // setters
    void setDirichletBC(SpMatrix<double>& A, DMatrix<double>& b);
    void setData(const BlockFrame<double, int>& df); // initialize model's data
    
    // getters
    const BlockFrame<double, int>& data() const { return df_; }
    const DMatrix<double>& y() const { return df_.get<double>(OBSERVATIONS_BLK); } // observation vector zy
    const DMatrix<int>& idx() const { return df_.get<int>(INDEXES_BLK); } // data indices
    // informations related to discretization of regularization term
    const PDE& pde() const { return *pde_; } // regularizing term Lf - u (defined on some domain \Omega)
    const Mesh<M,N,K>& domain() const { return pde_->domain(); }
    std::size_t n_basis() const { return pde_->domain().dof(); }; // number of basis functions used in space discretization
    std::size_t n_obs() const { return df_.template get<double>(OBSERVATIONS_BLK).rows(); } // number of observations
    const ADT<M,N,K>& gse() { if(gse_ == nullptr){ gse_ = std::make_unique<ADT<M,N,K>>(domain()); } return *gse_; }
    
    // abstract part of the interface, must be implemented by concrete models
    virtual void solve() = 0; // finds a solution to the problem, whatever the problem is.
    // destructor
    virtual ~ModelBase() = default;
  protected:   
    std::shared_ptr<PDE> pde_; // regularizing term Lf - u and domain definition D
    std::unique_ptr<ADT<M,N,K>> gse_; // geometric search engine
    BlockFrame<double, int> df_; // blockframe for data storage

    // getter to underlying model object
    inline Model& model() { return static_cast<Model&>(*this); }
    inline const Model& model() const { return static_cast<const Model&>(*this); } // const version
  };

#include "ModelBase.tpp"
  
}}

#endif // __MODEL_BASE_H__

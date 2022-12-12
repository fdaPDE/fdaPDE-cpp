#ifndef __I_STAT_MODEL__
#define __I_STAT_MODEL__

#include <cstddef>
#include <memory>
#include <Eigen/LU>
#include <type_traits>

#include "../core/utils/Symbols.h"
#include "../core/utils/Traits.h"
using fdaPDE::is_base_of_template;
#include "../core/MESH/engines/AlternatingDigitalTree/ADT.h"
using fdaPDE::core::MESH::ADT;
#include "../core/utils/DataStructures/BlockFrame.h"

namespace fdaPDE {
namespace models {

  // enumeration of possible sampling strategies
  enum SamplingStrategy{ GeostatisticalAtNodes, GeostatisticalAtLocations, Areal };
  
  // standardized definitions for stat model BlockFrame. layers below will make heavy assumptions on the layout of the BlockFrame,
  // use these instead of manually typing the block name when accessing df_
#define OBSERVATIONS_BLK "y" // matrix of observations
#define INDEXES_BLK "i"      // vector of observation indices
#define LOCATIONS_BLK "P"    // matrix of spatial locations coordinates
#define AREAL_BLK "D"        // incidence matrix for areal observations

  // type of regularization
  struct SpaceOnly {};
  struct SpaceTimeSeparable {};
  struct SpaceTimeParabolic {};

  // base class for model traits
  template <typename B> struct model_traits;
  
  // abstract base interface for any fdaPDE statistical model. Uses CRTP pattern
  template <typename Model>
  class iStatModel {
  protected:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    
    // helper to check if data member contains valid data
    template <typename T>
    bool isAlloc(const T& t) const { return t.size() != 0; }
    
    // problem's data
    static constexpr std::size_t M = PDE::local_dimension;
    static constexpr std::size_t N = PDE::embedding_dimension;
    static constexpr std::size_t K = PDE::basis_order;
    std::shared_ptr<PDE> pde_; // regularizing term Lf - u and domain definition \Omega
    BlockFrame<double, int> df_; // blockframe for data storage

    // algorithm used by MESH to solve search queries
    std::shared_ptr<ADT<M,N,K>> searchEngine_;
    // n x N matrix \Psi = [\psi_{ij}] = \psi_j(p_i) whose (i,j)-th entry is the
    // evaluation of the j-th basis function at the i-th spatial location
    SpMatrix<double> Psi_{};
    const SpMatrix<double>& __Psi(); // pointer to n x N sparse matrix \Psi. This computes \Psi if not available
    // diagonal matrix of subdomains' measure, used only in case of areal sampling
    Eigen::DiagonalMatrix<double, Eigen::Dynamic, Eigen::Dynamic> D_;
    SpMatrix<double> PsiTD_{}; // matrix \Psi corrected for areal observations (stores \Psi^T*D if D \neq I)
  public:
    // constructor
    iStatModel() = default;
    iStatModel(const PDE& pde)
      : pde_(std::make_shared<PDE>(pde)) {};
    // copy constructor, copy only pde object (as a consequence also the problem domain)
    iStatModel(const iStatModel& rhs) { pde_ = rhs.pde_; }
    
    // setters
    void setData(const BlockFrame<double, int>& df); // initialize model's data directly from a BlockFrame object 
    void setDirichletBC(SpMatrix<double>& A, DMatrix<double>& b); // set dirichlet boundary conditions. boundary data are passed by pde_ object
  
    // getters
    const Model& get() const { return static_cast<const Model&>(*this); } // get underlying model object
    const BlockFrame<double, int>& data() const { return df_; }
    std::size_t obs() const { return df_.template get<double>(OBSERVATIONS_BLK).rows(); } // number of observations
    std::size_t locs() const; // number of geostatistical locations where data are observed
    std::size_t nbasis() const { return pde_->domain().dof(); }; // number of basis functions used in space discretization
    const PDE& pde() const { return *pde_; } // regularizing term Lf - u (defined on some domain \Omega)
    const DMatrix<double>& y() const { return df_.get<double>(OBSERVATIONS_BLK); } // observation vector zy
    const DMatrix<int>& idx() const { return df_.get<int>(INDEXES_BLK); } // data indices
    SamplingStrategy sampling() const; // how data are sampled
    const ADT<M,N,K>& searchEngine(); // algorithm used to search element over the mesh
    // available if data are sampled at general locations inside the domain
    DMatrix<double> locations() const {
      if constexpr(std::is_same<typename model_traits<Model>::RegularizationType, SpaceOnly>::value)
	return df_.get<double>(LOCATIONS_BLK);
      else return df_(0, df_.rows()/get().time_domain().rows()-1).template get<double>(LOCATIONS_BLK);
    }
    // available if data are sampled at subdomains (areal observations)
    const DMatrix<int>& subdomains() const { return df_.get<int>(AREAL_BLK); }

    // utilities
    bool dataAtNodes() const { return !df_.hasBlock(LOCATIONS_BLK); } // true if locations are a subset of mesh nodes
    
    // abstract part of the interface, must be implemented by concrete models
    virtual void solve() = 0; // finds a solution to the problem, whatever the problem is.
    // destructor
    virtual ~iStatModel() = default;  
  };

#include "iStatModel.tpp"
  
  // this macro is intended to import all **common** symbols a model can expect from its parent Base class
  // a type Base must be in the scope of the macro
#define IMPORT_STAT_MODEL_SYMBOLS	\
  /* direct accessible fields */	\
  using Base::pde_;			\
  using Base::df_;			\
  using Base::locs;			\
  using Base::obs;			\
  using Base::y;			\
  using Base::Psi;			\
  using Base::R0;			\
  using Base::R1;			\
  using Base::u;			\
  /* utilities */			\
  using Base::idx;			\
  using Base::isAlloc;			\
  using Base::sampling;			\
  using Base::dataAtNodes;		\
  
  // trait to detect if a type implements iStatModel
  template <typename T>
  struct is_stat_model {
    static constexpr bool value = fdaPDE::is_base_of_template<iStatModel, T>::value;
  };

  // forward declarations
  template <typename Model> class iSpaceOnlyModel;
  template <typename Model> class iSpaceTimeSeparableModel;
  template <typename Model> class iSpaceTimeParabolicModel;
  // trait for the selection of the type of regularization on the base of the property of a model
  template <typename Model>
  struct select_regularization_type {
    using type = typename std::conditional<
      std::is_same<typename model_traits<Model>::RegularizationType, SpaceOnly>::value,
      iSpaceOnlyModel<Model>,
      typename std::conditional<
        std::is_same<typename model_traits<Model>::RegularizationType, SpaceTimeSeparable>::value,
        iSpaceTimeSeparableModel<Model>,
        iSpaceTimeParabolicModel<Model>>::type
      >::type;
  };
  
  
}}

#endif // __I_STAT_MODEL__

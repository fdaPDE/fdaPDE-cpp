#ifndef __I_STAT_MODEL__
#define __I_STAT_MODEL__

#include <cstddef>
#include <memory>
#include <Eigen/LU>
#include <type_traits>

#include "../core/utils/Symbols.h"
#include "../core/MESH/engines/AlternatingDigitalTree/ADT.h"
using fdaPDE::core::MESH::ADT;
#include "../core/utils/DataStructures/BlockFrame.h"
#include "../core/FEM/PDE.h"
using fdaPDE::core::FEM::PDE;

namespace fdaPDE {
namespace models {

  // enumeration of possible sampling strategies
  enum SamplingStrategy{ GeostatisticalAtNodes, GeostatisticalAtLocations, Areal };
  
  // standardized definitions for stat model BlockFrame. layers below will make heavy assumptions on the layout of the BlockFrame,
  // use these instead of manually typing the block name when accessing df_
#define STAT_MODEL_Y_BLK "y" // matrix of observations
#define STAT_MODEL_I_BLK "i" // vector of observation indices
#define STAT_MODEL_P_BLK "P" // matrix of spatial locations coordinates
#define STAT_MODEL_D_BLK "D" // incidence matrix for areal observations
  
  // abstract base interface for any fdaPDE statistical model.
  //   * n: the number of observations
  //   * N: the number of locations where data are observed
  //   * q: the number of covariates
  template <typename PDE>
  class iStatModel {
  protected:
    // helper to check if data member contains valid data
    template <typename T>
    bool isAlloc(const T& t) const { return t.size() != 0; }
    
    // problem's data
    static constexpr std::size_t M = PDE::local_dimension;
    static constexpr std::size_t N = PDE::embedding_dimension;
    static constexpr std::size_t K = PDE::basis_order;
    std::shared_ptr<PDE> pde_; // regularizing term Lf - u and domain definition \Omega
    double lambda_; // smoothing parameter
    BlockFrame<double, int> df_; // blockframe for data storage

    // algorithm used by MESH to solve search queries
    std::shared_ptr<ADT<M,N,K>> searchEngine_;
    // n x N matrix \Psi = [\psi_{ij}] = \psi_j(p_i) whose (i,j)-th entry is the
    // evaluation of the j-th basis function at the i-th spatial location
    SpMatrix<double> Psi_{};
    // diagonal matrix of subdomains' measure, used only in case of areal sampling
    Eigen::DiagonalMatrix<double, Eigen::Dynamic, Eigen::Dynamic> D_;
    SpMatrix<double> PsiTD_{}; // matrix \Psi corrected for areal observations (stores \Psi^T*D if D \neq I)
  public:
    // constructor
    iStatModel() = default;
    iStatModel(const PDE& pde, double lambda)
      : pde_(std::make_shared<PDE>(pde)), lambda_(lambda) {};
    // copy constructor, copy only pde object (as a consequence also the problem domain)
    iStatModel(const iStatModel& rhs) { pde_ = rhs.pde_; }

    // setters
    void setData(const BlockFrame<double, int>& df); // initialize model's data directly from a BlockFrame object
    void setLambda(double lambda) { lambda_ = lambda; } 
    void setDirichletBC(SpMatrix<double>& A, DMatrix<double>& b); // set dirichlet boundary conditions. boundary data are passed by pde_ object
  
    // getters
    const BlockFrame<double, int>& data() const { return df_; } 
    std::size_t loc() const { return pde_->domain().nodes();} // ???
    std::size_t obs() const { return df_.get<double>(STAT_MODEL_Y_BLK).rows(); } // number of observations
    const PDE& pde() const { return *pde_; } // regularizing term Lf - u (defined on some domain \Omega)
    const DMatrix<double>& y() const { return df_.get<double>(STAT_MODEL_Y_BLK); } // observation vector z
    const DMatrix<int>& idx() const { return df_.get<int>(STAT_MODEL_I_BLK); } // data indices
    double lambda() const { return lambda_; } // smoothing parameter \lambda
    SamplingStrategy sampling() const; // sampling strategy adopted from the model.
    const ADT<M,N,K>& searchEngine(); // algorithm used to search element over the mesh
    // available if data are sampled at general locations inside the domain
    const DMatrix<double>& locations() const { return df_.get<double>(STAT_MODEL_P_BLK); }
    // available if data are sampled at subdomains (areal observations)
    const DMatrix<int>& subdomains() const { return df_.get<int>(STAT_MODEL_D_BLK); }

    // pointers to FEM related quantites
    const SpMatrix<double>& R0() const { return *(pde_->R0()); }
    const SpMatrix<double>& R1() const { return *(pde_->R1()); }
    const DMatrix<double>&  u()  const { return *(pde_->force()); }
    const SpMatrix<double>& Psi(); // pointer to n x N sparse matrix \Psi. This computes \Psi if not available
    
    // utilities
    bool dataAtNodes() const { return !df_.hasBlock(STAT_MODEL_P_BLK); } // true if locations are a subset of mesh nodes
    // an efficient implementation of left multiplication by \Psi
    DMatrix<double> lmbPsi(const DMatrix<double>& x) const;
    auto PsiTD() const { // returns the block \Psi^T*D as eigen expression, if D = I returns \Psi^T
      return sampling() == SamplingStrategy::Areal ? PsiTD_ : Psi_.transpose(); }; 
    
    // abstract part of the interface, must be implemented by concrete models
    virtual void solve() = 0; // finds a solution to the problem, whatever the problem is.
    // destructor
    virtual ~iStatModel() = default;  
  };

#include "iStatModel.tpp"
  
  // import all symbols from iStatModel interface in derived classes
#define IMPORT_STAT_MODEL_SYMBOLS(E)		 \
  /* direct accessible fields */		 \
  using iStatModel<E>::pde_;			 \
  using iStatModel<E>::df_;			 \
  /* those info should be accessed by methods */ \
  using iStatModel<E>::lambda;			 \
  using iStatModel<E>::Psi;			 \
  using iStatModel<E>::loc;			 \
  using iStatModel<E>::obs;			 \
  using iStatModel<E>::y;			 \
  using iStatModel<E>::idx;			 \
  using iStatModel<E>::R0;			 \
  using iStatModel<E>::R1;			 \
  using iStatModel<E>::u;			 \
  using iStatModel<E>::isAlloc;			 \
  using iStatModel<E>::PsiTD;			 \
  
  // trait to detect if a type implements iStatModel
  template <typename T>
  struct is_stat_model {
    static constexpr bool value = fdaPDE::is_base_of_template<iStatModel, T>::value;
  };
  
}}
#endif // __I_STAT_MODEL__

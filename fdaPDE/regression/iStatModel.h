#ifndef __I_STAT_MODEL__
#define __I_STAT_MODEL__

#include <cstddef>
#include <memory>
#include <Eigen/LU>

#include "../core/utils/Symbols.h"
#include "../core/FEM/PDE.h"
using fdaPDE::core::FEM::PDE;
#include "Internals.h"
using fdaPDE::regression::internal::psi;
#include "../core/MESH/engines/AlternatingDigitalTree/ADT.h"
using fdaPDE::core::MESH::ADT;

// trait to detect if a type is a specialization of iStatModel
template <typename M> struct is_stat_model {
  // query the type M for a flag which is owned only by iStatModel and its derived classes
  static constexpr bool value = M::stat_model ? true : false;
};

// abstract base class for any fdaPDE statistical model
template <unsigned int M, unsigned int N, unsigned int K, typename E, typename B>
class iStatModel {
protected:
  // helper to check if data member contains valid data
  template <typename T>
  bool isAlloc(const std::shared_ptr<T>& t) const { return t != nullptr && t->size() != 0; }

  // data of the problem
  std::shared_ptr<PDE<M,N,K,E,B>> pde_;   // regularizing term and domain information
  double lambda_;                         // smoothing parameter
  std::shared_ptr<DVector<double>> z_{};  // vector of observations
  std::vector<std::size_t> z_idx_{};      // subset of locations where data are observed
  std::shared_ptr<DMatrix<double>> W_{};  // design matrix

  // matrix of locations. If this matrix contains any data it is assumed that data are observed at general locations in space.
  // Otherwise this information is not provided and data are assumed to be observed at mesh nodes.
  std::shared_ptr<DMatrix<double>> locations_{};

  // notation:
  //   * n: number of observations
  //   * N: number of locations where data are observed
  //   * q: number of regressors

  // algorithm used by MESH to solve search queries
  std::shared_ptr<ADT<M,N,K>> searchEngine_;
  
  // data coming from FEM module
  std::shared_ptr<SpMatrix<double>> R0_{};  // mass matrix, result of the discretization of the identity operator (N x N matrix)
  std::shared_ptr<SpMatrix<double>> R1_{};  // discretization matrix of the differential operator L (N x N matrix)
  std::shared_ptr<DMatrix<double>>  u_{};   // discretization of forcing term (1 x N vector)
  
  // n x N matrix \Psi = [\psi_{ij}] = \psi_j(p_i) whose ij-th entry is the evaluation of the j-th basis function at the i-th spatial location
  std::shared_ptr<SpMatrix<double>> Psi_{};
  
  // system matrix of non-parametric problem (2N x 2N matrix)
  //     | -\Psi^T*\Psi  \lambda*R1^T |
  // A = |                            |
  //     | \lambda*R1    \lambda*R0   |
  std::shared_ptr<SpMatrix<double>> A_{};
  // right hand side of problem's linear system (1 x 2N vector)
  //     | -\Psi^T*Q*z |
  // b = |             |
  //     |  \lambda*u  |
  std::shared_ptr<DVector<double>>  b_{};

  std::shared_ptr<DMatrix<double>>  H_{};   // hat matrix of the problem: H_ = W*(W*T*W)^{-1}*W^T (n x n matrix)
  std::shared_ptr<DMatrix<double>>  WTW_{}; // q x q dense matrix W^T*W
  std::shared_ptr<DMatrix<double>>  Q_{};   // n x n projection matrix onto the orthogonal space of Im(W), Q_ = I - H_
  
  // partial LU (with pivoting) factorization of the dense (square invertible) q x q matrix W^T*W.
  Eigen::PartialPivLU<DMatrix<double>> invWTW_{};

  // problem solution
  std::shared_ptr<DVector<double>> f_;      // estimate of the spatial field (1 x N vector)
  std::shared_ptr<DVector<double>> beta_;   // estimate of the coefficient vector (1 x q vector)
  
public:
  // constructor
  iStatModel() = default;
  iStatModel(const PDE<M,N,K,E,B>& pde, double lambda)
    : pde_(std::make_shared<PDE<M,N,K,E,B>>(pde)), lambda_(lambda) {};

  // copy constructor, copy only pde object (as a consequence also the problem domain)
  iStatModel(const iStatModel& rhs) {
    pde_ = rhs.pde_;
  }
  
  // set problem design matrix and precomputes all related quantites
  void setCovariates(const DMatrix<double>& W) {
    // store design matrix
    W_ = std::make_shared<DMatrix<double>>(W);
    // compute q x q dense matrix
    WTW_ = std::make_shared<DMatrix<double>>(W.transpose()*W);
    // compute the factorization of W^T*W
    invWTW_ = WTW_->partialPivLu();
    return;
  }
  // set observation vector (assumes there is an observation for each mesh node)
  void setObservations(const DVector<double>& z) {
    z_ = std::make_shared<DVector<double>>(z);
    return;
  }
  // set observation vector in case data are observed at general locations
  void setObservations(const DVector<double>& z, const DMatrix<double>& locations){
    z_ = std::make_shared<DVector<double>>(z);
    locations_ = std::make_shared<DMatrix<double>>(locations);
    return;
  }
  
  // set observation indexes: the i-th element in z_idx denotes the spatial point p_i where the i-th datum in z is observed
  // this allows to set observation and locations and then filter out some of them (used e.g. in KFoldCV)
  void setObsIndexes(const std::vector<std::size_t>& z_idx){
    z_idx_ = z_idx;
    return;
  }
  
  // set smoothing parameter
  void setLambda(double lambda) { lambda_ = lambda; }

  // getters
  std::size_t q() const { return isAlloc(W_) ? W_->cols() : 0; } // q, the number of covariates
  std::size_t loc() const { return pde_->domain().nodes(); }     // N, the number of locations
  std::size_t obs() const { return z_->rows(); }                 // n, the number of observations
  // pointers to problem data
  std::shared_ptr<DVector<double>> z() const { return z_; } // observation vector
  const std::vector<std::size_t>& z_idx() const { return z_idx_; } // data locations 
  std::shared_ptr<DMatrix<double>> W() const { return W_; } // design matrix
  double lambda() const { return lambda_; } // smoothing parameter
  std::shared_ptr<PDE<M,N,K,E,B>> pde() const { return pde_; }
  std::shared_ptr<ADT<M,N,K>> searchEngine() {
    // if still not available, initialize the search engine
    if(searchEngine_ == nullptr){
      searchEngine_ = std::make_shared<ADT<M,N,K>>(pde_->domain());
    }
    return searchEngine_;
  }
  std::shared_ptr<DMatrix<double>> locations() const { return locations_; }
  
  // pointer to projection matrix. Q is computed on demand only when it is needed (in general operations involving Q can be substituted
  // with the more efficient routine lmbQ())
  std::shared_ptr<DMatrix<double>> Q() {
    if(!isAlloc(Q_)){ // Q is computed on request since not needed in general
      // compute transpose of W once here
      DMatrix<double> Wt = W_->transpose();

      // compute hat matrix H = W*(W*W^T)^{-1}*W^T
      H_ = std::make_shared<DMatrix<double>>((*W_)*invWTW_.solve(Wt));
      // compute Q = I - H_
      Q_ = std::make_shared<DMatrix<double>>
	(DMatrix<double>::Identity(H_->rows(), H_->cols()) - *H_);
    }
    return Q_;
  }

  // pointer to q x q dense matrix W^T*W and its inverse
  std::shared_ptr<DMatrix<double>>  WTW() const { return WTW_; }
  Eigen::PartialPivLU<DMatrix<double>> invWTW() const { return invWTW_; }

  // pointers to FEM related quantites
  std::shared_ptr<SpMatrix<double>> R0() const { return pde_->R0(); }
  std::shared_ptr<SpMatrix<double>> R1() const { return pde_->R1(); }
  std::shared_ptr<DMatrix<double>>  u()  const { return pde_->force(); }
  
  std::shared_ptr<SpMatrix<double>> A() const { return A_; }   // pointer to non-parametric part of the problem
  std::shared_ptr<DVector<double>>  b() const { return b_; }   // pointer to rhs of the problem

  // pointer to n x N sparse matrix \Psi.
  std::shared_ptr<SpMatrix<double>> Psi() {
    if(!isAlloc(Psi_)){ // compute \Psi if not already available
      // call to Internals::psi already takes into account of the sampling strategy of the model
      Psi_ = psi(*this);
    }
    return Psi_;
  } 

  // pointers to problem solution
  std::shared_ptr<DVector<double>> f_hat() const { return f_; }
  std::shared_ptr<DVector<double>> beta_hat() const { return beta_; }

  // compile time informations
  static constexpr unsigned int local_dimension = M;
  static constexpr unsigned int embedding_dimension = N;
  static constexpr unsigned int mesh_order = K;
  static constexpr bool stat_model = true; // used by trait is_stat_model to assert if a type is a statistical model
  
  // methods
  bool hasCovariates() const { return q() != 0; } // true if the model has a parametric part
  bool dataAtNodes() const { return !isAlloc(locations_); } // true if locations are subset of mesh nodes

  // an efficient way to perform a left multiplication by Q. The following method is based on the following strategy
  //  given the design matrix W and x
  //    compute v = W^T*x
  //    solve Yz = v
  //    return x - Wz = Qx
  // it is required to having assigned a design matrix W to the model before calling this method
  DMatrix<double> lmbQ(const DMatrix<double>& x){
    DMatrix<double> v = W_->transpose()*x; // W^T*x
    DMatrix<double> z = invWTW_.solve(v);  // (W^T*W)^{-1}*W^T*x
    // compute x - W*z = x - (W*(W^T*W)^{-1}*W^T)*x = (I - H)*x = Q*x
    return x - (*W_)*z;
  }

  // an efficient implementation of left multiplication by \Psi
  std::shared_ptr<DMatrix<double>> lmbPsi(const DMatrix<double>& x){
    // compute dimensions of resulting matrix
    std::size_t n = Psi_->rows();
    std::size_t m = x.cols();
    // preallocate space for n x m result matrix
    std::shared_ptr<DMatrix<double>> result = std::make_shared<DMatrix<double>>(n,m);

    // if data are sampled at mesh nodes (or a subset of them) then \Psi is a permutation matrix
    if(dataAtNodes()){
      // just permute input matrix columns
      for(std::size_t k = 0; k < Psi_->outerSize(); ++k){
	for (SpMatrix<double>::InnerIterator it(*Psi_,k); it; ++it){
	  result->row(it.row()) = x.row(it.col());
	}
      }
    }else{
      // in the general case no optimization can be put in place
      *result = (*Psi_)*x;
    }
    return result;
  }
  
  // abstract part of the interface, must be implemented by concrete models

  // finds a solution to the smoothing problem.
  // After a call to smooth() all quantites related to the solution of the problem must contain valid data
  virtual void smooth() = 0;
  // computes \hat z, the fitted values at the observations' locations
  virtual DVector<double> fitted() const = 0;
  // compute prediction at new location
  virtual double predict(const DVector<double>& covs, const std::size_t loc) const = 0;
  
  virtual ~iStatModel() = default;
};

// import all symbols from iStatModel interface in derived classes
#define IMPORT_STAT_MODEL_SYMBOLS(M,N,K,E,B)		\
  using iStatModel<M,N,K,E,B>::pde_;			\
  using iStatModel<M,N,K,E,B>::lambda_;			\
  using iStatModel<M,N,K,E,B>::z_;			\
  using iStatModel<M,N,K,E,B>::W_;			\
  using iStatModel<M,N,K,E,B>::R0_;			\
  using iStatModel<M,N,K,E,B>::R1_;			\
  using iStatModel<M,N,K,E,B>::u_;			\
  using iStatModel<M,N,K,E,B>::Psi_;			\
  using iStatModel<M,N,K,E,B>::A_;			\
  using iStatModel<M,N,K,E,B>::b_;			\
  using iStatModel<M,N,K,E,B>::H_;			\
  using iStatModel<M,N,K,E,B>::WTW_;			\
  using iStatModel<M,N,K,E,B>::Q_;			\
  using iStatModel<M,N,K,E,B>::invWTW_;			\
  using iStatModel<M,N,K,E,B>::f_;			\
  using iStatModel<M,N,K,E,B>::beta_;			\

#endif // __I_STAT_MODEL__

#ifndef __SRPDE_H__
#define __SRPDE_H__

#include <memory>
#include "../core/utils/Symbols.h"
#include "../core/FEM/PDE.h"
using fdaPDE::core::FEM::PDE;
#include "../core/NLA/SparseBlockMatrix.h"
using fdaPDE::core::NLA::SparseBlockMatrix;
#include "../core/NLA/SMW.h"
using fdaPDE::core::NLA::SMW;
#include "Internals.h"
using fdaPDE::regression::internal::psi;
using fdaPDE::regression::internal::lmbQ;

template <unsigned int M, unsigned int N, unsigned int R, typename E>
class SRPDE{
private:
  // data of the problem
  const PDE<M,N,R,E>& pde_;
  double lambda_;

  // internal information required to solve an SR-PDE problem. Not all calls to .smooth() will fill all the
  // structures here listed. Anyway there is no additional overhead in declaring them here, the default constructor of a dynamic eigen
  // matrix "never performs any dynamic memory allocation, and never initializes the matrix coefficients." (from eigen documentation).
  // We assume that anything which has a size larger than 0 contains valid informations. Helper .isAlloc() perform exactly this check
  template <typename T>
  bool isAlloc(const std::shared_ptr<T>& t) const { return t != nullptr && t->size() != 0; }

  // notation:
  //   * n: number of observations
  //   * N: number of locations where data are observed (typically the number of mesh nodes)
  //   * q: number of regressors
  // In the following matrices R1 and R0 denotes respectively the N x N discretization matrix of the differential operator and the N x N
  // mass matrix computed by the PDE solver. You can access them from any pde object using .R1() and .R0() methods (see FEM module for details)

  // n x N matrix \Psi = [\psi_{ij}] = \psi_j(p_i) whose ij-th entry is the evaluation of the j-th basis function at the i-th spatial location 
  // (typically \psi_j denotes a Lagrangian basis function centered in the j-th node of the mesh)
  std::shared_ptr<SpMatrix<double>> Psi_;

  // nonparametric part of the SR-PDE linear system (2N x 2N matrix)
  //     | -\Psi^T*\Psi  \lambda*R1^T |
  // A = |                            |  
  //     | \lambda*R1    \lambda*R0   |
  std::shared_ptr<SpMatrix<double>> A_;

  // design matrix
  std::shared_ptr<DMatrix<double>> W_;
  // hat matrix of the regression.
  // Letting W the n x q design matrix of the problem, we have H_ = W*(W*T*W)^{-1}*W^T (n x n matrix)
  std::shared_ptr<DMatrix<double>> H_;

  // q x q dense matrix W^T*W
  std::shared_ptr<DMatrix<double>> WTW_;
  // partial LU (with pivoting) factorization of the dense (square invertible) q x q matrix W^T*W. Cached to compute expressions
  // involving the (W^T*W)^{-1} factor without the need to recompute it
  Eigen::PartialPivLU<DMatrix<double>> invWTW_;
  
  // n x n projection matrix onto the orthogonal space of Im(W), Q_ = I - H_
  std::shared_ptr<DMatrix<double>> Q_;

  // results of the regression problem
  std::shared_ptr<DVector<double>> f_;    // estimate of the spatial field (nonparametric part)
  std::shared_ptr<DVector<double>> beta_; // estimate of the coefficient vector (parametric part)

  DVector<double> z_; // vector of observations
public:
  // constructor
  SRPDE() = default;
  SRPDE(const PDE<M,N,R,E>& pde, double lambda) : pde_(pde), lambda_(lambda) {
    Psi_ = psi(pde_);
  };

  // finds a solution to the regression problem
  void smooth();
  // computes the fitted values \hat z
  DVector<double> fitted() const;
  
  // set problem design matrix and precomputes all related quantites
  void setCovariates(const DMatrix<double>& W) {
    // store design matrix
    W_ = std::make_shared<DMatrix<double>>(W);

    // compute transpose of W once here
    DMatrix<double> Wt = W.transpose();
    // compute q x q dense matrix
    WTW_ = std::make_unique<DMatrix<double>>(Wt*W);
    // compute the factorization of the dense q x q W^T*W matrix
    invWTW_ = WTW_->partialPivLu();
    // compute hat matrix H = W*(W*W^T)^{-1}*W^T
    H_ = std::make_unique<DMatrix<double>>(W*invWTW_.solve(Wt));
    // compute Q = I - H_
    Q_ = std::make_unique<DMatrix<double>>
      (DMatrix<double>::Identity(H_->rows(), H_->cols()) - *H_);
    return;
  }

  void setObservations(const DVector<double>& z) { z_ = z; }
  void setLambda(double lambda) { lambda_ = lambda; }
  
  // getters
  std::size_t q() const { return isAlloc(W_) ? W_->size() : 0; } // number of covariates
  std::size_t n() const { return pde_.domain().nodes(); } // number of observation locations (need to change if locations != nodes)
  std::shared_ptr<DMatrix<double>> W() const { return W_; }
  std::shared_ptr<DMatrix<double>> Q() const { return Q_; }
  Eigen::PartialPivLU<DMatrix<double>> invWTW() const { return invWTW_; }
  std::shared_ptr<SpMatrix<double>> Psi() const { return Psi_; }

  std::shared_ptr<DMatrix<double>> T(double lambda) const;
  DVector<double> z() const { return z_; }
  
  // solution of smoothing problem
  DVector<double> f() const { return *f_; }
  DVector<double> beta() const { return *beta_; }
};

// template argument deduction guide
template <unsigned int M, unsigned int N, unsigned int R, typename E>
SRPDE(const PDE<M,N,R,E>& pde_, double lambda_) -> SRPDE<M,N,R,E>;

// Let z the vector of observations passed as argument to the smooth function
// For a nonparametric SR-PDE model we have to solve the linear system A*x = b with
//
//         | -\Psi^T*\Psi  \lambda*R1^T |        | -Psi^T*z  |
//     A = |                            |    b = |           |
//         | \lambda*R1    \lambda*R0   |        | \lambda*u |
//
// The system is solved in an efficient way exploiting the sparsity of matrix A.
// For a parametric SR-PDE problem we need instead to solve a slightly different linear system M*x = b where
// 
//         | -\Psi^T*Q*\Psi  \lambda*R1^T |      | -Psi^T*Q*z |
//     M = |                              |  b = |            |
//         | \lambda*R1    \lambda*R0     |      | \lambda*u  |
//
// Due to the north-west dense block of matrix M an efficient solution of the system cannot be obtained using a sparse solver, 
// instead letting A the system matrix of a nonparametric SR-PDE problem, we can show that M = A + U*C*V where
//
//         | Psi^T*W |
//     U = |         |  C = (W^T*W)^{-1}  V = | W^T*Psi  O_N |
//         |   O_N   |
//
// using the above decomposition an efficient solution to the system M*x = b is obtain by SMW decomposition (from core/NLA module)
template <unsigned int M, unsigned int N, unsigned int R, typename E>
void SRPDE<M, N, R, E>::smooth() {
  // assemble system matrix for the nonparameteric part of the model
  SparseBlockMatrix<double,2,2>
    A(-Psi_->transpose()*(*Psi_), lambda_ * pde_.R1().transpose(),
      lambda_ * pde_.R1(),        lambda_ * pde_.R0()            );
  // cache system matrix for reuse
  A_ = std::make_shared<SpMatrix<double>>(A.derived());
  DVector<double> solution; // where the system solution will be stored
  
  if(!isAlloc(W_)){ // nonparametric case
    // rhs of SR-PDE linear system
    DVector<double> b;
    b.resize(A.rows());
    b << -Psi_->transpose()*z_,
      lambda_ * pde_.force();

    // define system solver. Use a sparse solver
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver{};
    solver.analyzePattern(*A_);
    solver.factorize(*A_);
  
    // solve linear system A_*x = b
    solution = solver.solve(b);
  }else{ // parametric case
    // rhs of SR-PDE linear system
    DVector<double> b;
    b.resize(A_->rows());
    b << -Psi_->transpose()*lmbQ(*this, z_),
      lambda_ * pde_.force();

    std::size_t q_ = W_->cols(); // number of covariates
    // definition of matrices U and V 
    DMatrix<double> U = DMatrix<double>::Zero(A.rows(), q_);
    U.block(0,0, A.rows()/2, q_) = Psi_->transpose()*(*W_);
  
    DMatrix<double> V = DMatrix<double>::Zero(q_, A.rows());
    V.block(0,0, q_, A.rows()/2) = W_->transpose()*(*Psi_);

    // Define system solver. Use SMW solver from NLA module
    SMW<> solver{};
    solver.compute(*A_);
    // solve system Mx = b
    solution = solver.solve(U, *WTW_, V, b);

    // store estimation of coefficient vector for the parametric part
    beta_ = std::make_shared<DVector<double>>( invWTW_.solve(W_->transpose())*(z_ - *f_) );
  }
  // store estimation of coefficient vector for the nonparametric part
  f_ = std::make_shared<DVector<double>>( (*Psi_)*solution.head(A_->rows()/2) );
  return;
}

// it is asssumed that smooth has already been called on the model object
// in general fitted values \hat z are equal to f_ + W * beta_, in case the parametric part is absent W * beta_ is omitted
template <unsigned int M, unsigned int N, unsigned int R, typename E>
DVector<double> SRPDE<M, N, R, E>::fitted() const {
  DVector<double> hat_z = *f_;
  // if the model has a parametric part, we need to sum its contribute
  if(isAlloc(W_))
    hat_z += (*W_)*(*beta_);
  
  return hat_z;
}

// required to support GCV based smoothing parameter selection
// in case of an SRPDE model we have T = \Psi^T*Q*\Psi + \lambda*R1_^T*R0_^{-1}*R1_
template <unsigned int M, unsigned int N, unsigned int R, typename E>
std::shared_ptr<DMatrix<double>> SRPDE<M, N, R, E>::T(double lambda) const{
  // compute value of R_ = R1_^T*R0_^{-1}*R1_
  Eigen::SparseLU<SpMatrix<double>> invR0{};
  invR0.analyzePattern(pde_.R0());
  invR0.factorize(pde_.R0());
  DMatrix<double> R_ = pde_.R1().transpose()*invR0.solve(pde_.R1());

  if(!isAlloc(W_)) // case without covariates, Q is the identity matrix
    return std::make_shared<DMatrix<double>>
      (Psi_->transpose()*(*Psi_) + lambda*R_);
  else // general case with covariates
    return std::make_shared<DMatrix<double>>
      (Psi_->transpose()*lmbQ(*this, *Psi_) + lambda*R_);
}

#endif // __SRPDE_H__

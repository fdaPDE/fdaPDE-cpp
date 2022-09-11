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
  bool isAlloc(const std::unique_ptr<T>& t) { return t->size() != 0; }

  // notation:
  //   * n: number of observations
  //   * N: number of locations where data are observed (typically the number of mesh nodes)
  //   * q: number of regressors
  // In the following matrices R1 and R0 denotes respectively the N x N discretization matrix of the differential operator and the N x N
  // mass matrix computed by the PDE solver. You can access them from any pde object using .R1() and .R0() methods (see FEM module for details)

  // n x N matrix \Psi = [\psi_{ij}] = \psi_j(p_i) whose ij-th entry is the evaluation of the j-th basis function at the i-th spatial location 
  // (typically \psi_j denotes a Lagrangian basis function centered in the j-th node of the mesh)
  std::unique_ptr<SpMatrix<double>> Psi_;

  // nonparametric part of the SR-PDE linear system (2N x 2N matrix)
  //     | -\Psi^T*\Psi  \lambda*R1^T |
  // A = |                            |  
  //     | \lambda*R1    \lambda*R0   |
  std::unique_ptr<SpMatrix<double>> A_;

  // hat matrix of the regression.
  // Letting W the n x q design matrix of the problem, we have H_ = W*(W*T*W)^{-1}*W^T (n x n matrix)
  std::unique_ptr<DMatrix<double>> H_;

  // q x q dense matrix W^T*W
  std::unique_ptr<DMatrix<double>> WTW_;
  // partial LU (with pivoting) factorization of the dense (square invertible) q x q matrix W^T*W. Used to compute expressions
  // involving the (W^T*W)^{-1} factor
  Eigen::PartialPivLU<DMatrix<double>> invWTW_;
  
  // n x n projection matrix onto the orthogonal space of Im(W), Q_ = I - H_
  std::unique_ptr<DMatrix<double>> Q_;

  // results of the regression problem
  std::unique_ptr<DVector<double>> f_;    // estimate of the spatial field (nonparametric part)
  std::unique_ptr<DVector<double>> beta_; // estimate of the coefficient vector (parametric part)
  
public:
  // constructor
  SRPDE() = default;
  SRPDE(const PDE<M,N,R,E>& pde, double lambda) : pde_(pde), lambda_(lambda) {};

  // finds a solution to the regression problem
  void smooth(const DVector<double>& data);
  void smooth(const DVector<double>& data, const DMatrix<double>& covariates);

  // getters
  DVector<double> f() const { return *f_; }
  DVector<double> beta() const { return *beta_; }
};

template <unsigned int M, unsigned int N, unsigned int R, typename E>
SRPDE(const PDE<M,N,R,E>& pde_, double lambda_) -> SRPDE<M,N,R,E>;

// z is the vector of observations
template <unsigned int M, unsigned int N, unsigned int R, typename E>
void SRPDE<M, N, R, E>::smooth(const DVector<double>& z) {
  // for an SR-PDE model we have to solve A*x = b with
  //
  //         | -\Psi^T*\Psi  \lambda*R1^T |      | -Psi^T*z  |
  //     A = |                            |  b = |           |
  //         | \lambda*R1    \lambda*R0   |      | \lambda*u |
  //
  // The system is solved in an efficient way exploiting the sparsity of matrix A
  
  // assemble system matrix
  Psi_ = psi(pde_);
  SparseBlockMatrix<double,2,2>
    A(-Psi_->transpose()*(*Psi_), lambda_ * pde_.R1().transpose(),
      lambda_ * pde_.R1(),        lambda_ * pde_.R0()            );
  // develop SparseBlockMatrix here once (expensive operation), cache for reuse
  A_ = std::make_unique<SpMatrix<double>>(A.derived());
  
  // rhs of SR-PDE linear system
  DVector<double> b;
  b.resize(A.rows());
  b << -Psi_->transpose()*z,
       lambda_ * pde_.force();

  // define system solver. Matrix A is sparse!
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver{};
  solver.analyzePattern(*A_);
  solver.factorize(*A_);
  
  // solve linear system A_*x = b
  DVector<double> sol = solver.solve(b);
  std::cout << sol << std::endl;
  f_ = std::make_unique<DVector<double>>( (*Psi_)*sol.head(A_->rows()/2) );
}

// smoothing in case covariates are supplied to the model. In this case the efficient solution of the linear system requires the application
// of the SMW decomposition. See NLA module for more details.
template <unsigned int M, unsigned int N, unsigned int R, typename E>
void SRPDE<M, N, R, E>::smooth(const DVector<double>& z, const DMatrix<double>& W) {
  // for an SR-PDE problem with covariates we need to solve M*x = b where
  // 
  //         | -\Psi^T*Q*\Psi  \lambda*R1^T |      | -Psi^T*Q*z |
  //     M = |                              |  b = |            |
  //         | \lambda*R1    \lambda*R0     |      | \lambda*u  |
  //
  // Let A_ the system matrix of a nonparametric SR-PDE problem, we can show that M = A + U*C*V where
  //
  //      | Psi^T*W |
  // U  = |         |  C = (W^T*W)^{-1}  V = | W^T*Psi  O_N |
  //      |   O_N   |
  //
  // using the above decomposition an efficient solution to the system M*x = b is obtain using SMW decomposition (from NLA module)

  std::size_t q = W.cols(); // number of covariates
  // assemble system matrix for the nonparametric part
  Psi_ = psi(pde_);
  SparseBlockMatrix<double,2,2>
    A(-Psi_->transpose()*(*Psi_), lambda_ * pde_.R1().transpose(),
      lambda_ * pde_.R1(),        lambda_ * pde_.R0()            );
  // develop SparseBlockMatrix here once (expensive operation), cache for reuse
  A_ = std::make_unique<SpMatrix<double>>(A.derived());
  
  // Define SMW solver
  SMW<> solver{};
  solver.compute(*A_);

  // compute transpose of W once here
  DMatrix<double> Wt = W.transpose();

  // compute q x q dense matrix
  WTW_ = std::make_unique<DMatrix<double>>(Wt*W);
  // compute the factorization of the dense q x q W^T*W matrix
  invWTW_ = WTW_->partialPivLu();
  // compute hat matrix H = W*(W*W^T)^{-1}*W^T
  H_ = std::make_unique<DMatrix<double>>(W*invWTW_.solve(Wt));
  // compute Q = I - H_
  Q_ = std::make_unique<DMatrix<double>>(DMatrix<double>::Identity(H_->rows(), H_->cols()) - *H_);

  // compute rhs of SR-PDE linear system
  DVector<double> b;
  b.resize(A_->rows());
  b << -Psi_->transpose()*(*Q_)*z,
       lambda_ * pde_.force();

  // definition of matrices U and V required for SMW decomposition
  DMatrix<double> U = DMatrix<double>::Zero(A.rows(), q);
  U.block(0,0, A.rows()/2, q) = Psi_->transpose()*W;
  
  DMatrix<double> V = DMatrix<double>::Zero(q, A.rows());
  V.block(0,0, q, A.rows()/2) = W.transpose()*(*Psi_);
  // SMW implementation requires directly the inversion of C, computed before as WTW_ = W^T*W
  DVector<double> sol = solver.solve(U, *WTW_, V, b);

  // store result of regression problem
  f_    = std::make_unique<DVector<double>>( (*Psi_)*sol.head(A.rows()/2) );
  beta_ = std::make_unique<DVector<double>>( invWTW_.solve(W.transpose())*(z - *f_) );  
  return;
}

#endif // __SRPDE_H__

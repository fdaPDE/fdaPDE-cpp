// finds a solution to the SR-PDE smoothing problem
template <typename PDE, Sampling SamplingDesign>
void SRPDE<PDE, SamplingDesign>::solve() {
  // assemble system matrix for the nonparameteric part of the model
  SparseBlockMatrix<double,2,2>
    A(-PsiTD()*W()*Psi(), lambda()*R1().transpose(),
      lambda()*R1(),      lambda()*R0()            );
  // cache system matrix for reuse
  A_ = A.derived();
  b_.resize(A_.rows());
  DVector<double> sol; // room for problem' solution
  
  if(!hasCovariates()){ // nonparametric case
    // rhs of SR-PDE linear system
    b_ << -PsiTD()*W()*y(),
          lambda()*u();
    
    // define system solver. Use a sparse solver
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver{};
    solver.compute(A_);
    // solve linear system A_*x = b_
    sol = solver.solve(b_);
    
    // store result of smoothing
    f_ = sol.head(A_.rows()/2);
  }else{ // parametric case
    // rhs of SR-PDE linear system
    b_ << -PsiTD()*lmbQ(y()), // -\Psi^T*D*Q*z
          lambda()*u();
    
    std::size_t q_ = X().cols(); // number of covariates
    // definition of matrices U and V  for application of woodbury formula
    DMatrix<double> U = DMatrix<double>::Zero(A_.rows(), q_);
    U.block(0,0, A_.rows()/2, q_) = PsiTD()*W()*X();
    DMatrix<double> V = DMatrix<double>::Zero(q_, A_.rows());
    V.block(0,0, q_, A_.rows()/2) = X().transpose()*W()*Psi();

    // Define system solver. Use SMW solver from NLA module
    SMW<> solver{};
    solver.compute(A_);
    // solve system Mx = b
    sol = solver.solve(U, XtWX(), V, b_);
    // store result of smoothing 
    f_    = sol.head(A_.rows()/2);
    beta_ = invXtWX().solve(X().transpose()*W())*(y() - Psi()*f_);
  }
  // store PDE misfit
  g_ = sol.tail(A_.rows()/2);
  return;
}

// required to support GCV based smoothing parameter selection
// in case of an SRPDE model we have T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
template <typename PDE, Sampling SamplingDesign>
const DMatrix<double>& SRPDE<PDE, SamplingDesign>::T() {
  if(T_.size() == 0){ // compute only at first time
    // compute value of R = R1^T*R0^{-1}*R1, cache for possible reuse
    invR0_->compute(R0());
    R_ = R1().transpose()*invR0_->solve(R1());

    // compute and store matrix T for possible reuse
    if(!hasCovariates()) // case without covariates, Q is the identity matrix
      T_ = Psi().transpose()*W()*Psi() + lambda()*R_;
    else // general case with covariates
      T_ = Psi().transpose()*lmbQ(Psi()) + lambda()*R_;
  }
  return T_;
}

// Q is computed on demand only when it is needed by GCV and cached for fast reacess (in general operations
// involving Q can be substituted with the more efficient routine lmbQ(), which is part of iRegressionModel interface)
template <typename PDE, Sampling SamplingDesign>
const DMatrix<double>& SRPDE<PDE, SamplingDesign>::Q() {
  if(Q_.size() == 0){ // Q is computed on request since not needed in general
    // compute Q = W(I - H) = W - W*X*(X*W*X^T)^{-1}*X^T*W
    Q_ = W()*(DMatrix<double>::Identity(n_obs(), n_obs()) - X()*invXtWX().solve(X().transpose()*W()));
  }
  return Q_;
}

// returns the euclidean norm of y - \hat y
template <typename PDE, Sampling SamplingDesign>
double SRPDE<PDE, SamplingDesign>::norm(const DMatrix<double>& obs, const DMatrix<double>& fitted) const {
  return (obs - fitted).squaredNorm();
}

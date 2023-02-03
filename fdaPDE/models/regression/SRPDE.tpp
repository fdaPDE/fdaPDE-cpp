// finds a solution to the SR-PDE smoothing problem
template <typename PDE, Sampling SamplingDesign>
void SRPDE<PDE, SamplingDesign>::solve() {
  // assemble system matrix for the nonparameteric part of the model
  SparseBlockMatrix<double,2,2>
    A(-PsiTD()*W()*Psi(), lambdaS()*R1().transpose(),
      lambdaS()*R1(),     lambdaS()*R0()            );
  // cache non-parametric matrix and its factorization for reuse
  A_ = A.derived();
  invA_.compute(A_);
  b_.resize(A_.rows());
  DVector<double> sol; // room for problem' solution
  
  if(!hasCovariates()){ // nonparametric case
    // rhs of SR-PDE linear system
    b_ << -PsiTD()*W()*y(),
          lambdaS()*u();
    
    // solve linear system A_*x = b_
    sol = invA_.solve(b_);
    // store result of smoothing
    f_ = sol.head(A_.rows()/2);
  }else{ // parametric case
    // rhs of SR-PDE linear system
    b_ << -PsiTD()*lmbQ(y()), // -\Psi^T*D*Q*z
          lambdaS()*u();
    
    // definition of matrices U and V  for application of woodbury formula
    U_ = DMatrix<double>::Zero(A_.rows(), q());
    U_.block(0,0, A_.rows()/2, q()) = PsiTD()*W()*X();
    V_ = DMatrix<double>::Zero(q(), A_.rows());
    V_.block(0,0, q(), A_.rows()/2) = X().transpose()*W()*Psi();
    // solve system (A_ + U_*(X^T*W_*X)*V_)x = b using woodbury formula from NLA module
    sol = SMW<>().solve(invA_, U_, XtWX(), V_, b_);
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
  // compute value of R = R1^T*R0^{-1}*R1, cache for possible reuse
  if(R_.size() == 0){
    invR0_->compute(R0());
    R_ = R1().transpose()*invR0_->solve(R1());
  }
  // compute and store matrix T for possible reuse
  if(!hasCovariates()) // case without covariates, Q is the identity matrix
    T_ = PsiTD()*W()*Psi()   + lambdaS()*R_;
  else // general case with covariates
    T_ = PsiTD()*lmbQ(Psi()) + lambdaS()*R_;
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

// returns the euclidean norm of op1 - op2
template <typename PDE, Sampling SamplingDesign>
double SRPDE<PDE, SamplingDesign>::norm(const DMatrix<double>& op1, const DMatrix<double>& op2) const {
  return (op1 - op2).squaredNorm();
}

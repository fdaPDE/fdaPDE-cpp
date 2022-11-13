// an efficient way to perform a left multiplication by Q implementing the following
//  given the design matrix W and x
//    compute v = W^T*x
//    solve Yz = v
//    return x - Wz = Qx
// it is required to having assigned a design matrix W to the model before calling this method
template <typename PDE>
DMatrix<double> SRPDE<PDE>::lmbQ(const DMatrix<double>& x){
  DMatrix<double> v = W().transpose()*x; // W^T*x
  DMatrix<double> z = invWTW_.solve(v);  // (W^T*W)^{-1}*W^T*x
  // compute x - W*z = x - (W*(W^T*W)^{-1}*W^T)*x = (I - H)*x = Q*x
  return x - W()*z;
}

// finds a solution to the SR-PDE smoothing problem
template <typename PDE>
void SRPDE<PDE>::solve() {
  // ask to iStatModel to return \Psi matrix
  this->Psi();
  // assemble system matrix for the nonparameteric part of the model
  SparseBlockMatrix<double,2,2>
    A(-PsiTD()*Psi_,  lambda_ * R1().transpose(),
      lambda_ * R1(), lambda_ * R0()            );
  // cache system matrix for reuse
  A_ = A.derived();
 
  b_.resize(A_.rows());
  DVector<double> sol; // room for problem' solution
  
  if(!hasCovariates()){ // nonparametric case
    // rhs of SR-PDE linear system
    b_ << -PsiTD()*z(),
      lambda_*u();
    
    // define system solver. Use a sparse solver
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver{};
    solver.compute(A_);
    // solve linear system A_*x = b_
    sol = solver.solve(b_);
    
    // store result of smoothing
    f_ = sol.head(A_.rows()/2);
  }else{ // parametric case
    if(!isAlloc(WTW_)){
      // compute q x q dense matrix W^T*W and its factorization
      WTW_ = W().transpose()*W();
      invWTW_ = WTW_.partialPivLu();
    }
    // rhs of SR-PDE linear system
    b_ << -PsiTD()*lmbQ(z()), // -\Psi^T*D*Q*z
      lambda_*u();
    
    std::size_t q_ = W().cols(); // number of covariates
    // definition of matrices U and V  for application of woodbury formula
    DMatrix<double> U = DMatrix<double>::Zero(A_.rows(), q_);
    U.block(0,0, A_.rows()/2, q_) = PsiTD()*W();
    DMatrix<double> V = DMatrix<double>::Zero(q_, A_.rows());
    V.block(0,0, q_, A_.rows()/2) = W().transpose()*Psi_;

    // Define system solver. Use SMW solver from NLA module
    SMW<> solver{};
    solver.compute(A_);
    // solve system Mx = b
    sol = solver.solve(U, WTW_, V, b_);
    // store result of smoothing 
    f_    = sol.head(A_.rows()/2);
    beta_ = invWTW_.solve(W().transpose())*(z() - Psi_*f_);
  }
  return;
}

// it is asssumed that smooth has already been called on the model object
// computes fitted values \hat z = \Psi*f_ + W*beta_
template <typename PDE>
DMatrix<double> SRPDE<PDE>::fitted() const {
  DMatrix<double> hat_z = Psi_*f_;
  // if the model has a parametric part, we need to sum its contribute
  if(hasCovariates()) hat_z += W()*beta_;
  return hat_z;
}

// compute prediction of model at new unseen data location (location equal to mesh node)
// W_{n+1}^T * \beta + f_*\psi(p_{n+1})
template <typename PDE>
double SRPDE<PDE>::predict(const DVector<double>& covs, std::size_t loc) const {
  double prediction = f_.coeff(loc,0);
  // parametetric contribute of the model, if any is present
  if(hasCovariates()) prediction += (covs.transpose()*beta_).coeff(0,0);
  return prediction;
}

// required to support GCV based smoothing parameter selection
// in case of an SRPDE model we have T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
template <typename PDE>
const DMatrix<double>& SRPDE<PDE>::T() {
  if(!isAlloc(T_)){ // compute only at first time
    // compute value of R = R1^T*R0^{-1}*R1, cache for possible reuse
    invR0_->compute(R0());
    R_ = R1().transpose()*invR0_->solve(R1());

    // compute and store matrix T for possible reuse
    if(!hasCovariates()) // case without covariates, Q is the identity matrix
      T_ = Psi_.transpose()*Psi_ + lambda_*R_;
    else // general case with covariates
      T_ = Psi_.transpose()*lmbQ(Psi_) + lambda_*R_;
  }
  return T_;
}

// Q is computed on demand only when it is needed by GCV and cached for fast reacess (in general operations
// involving Q can be substituted with the more efficient routine lmbQ(), which is part of iRegressionModel interface)
template <typename PDE>
const DMatrix<double>& SRPDE<PDE>::Q() {
  if(!isAlloc(Q_)){ // Q is computed on request since not needed in general
    // compute Q = I - H = I - W*(W*W^T)^{-1}*W^T
    Q_ = DMatrix<double>::Identity(obs(), obs()) - W()*invWTW_.solve(W().transpose());
  }
  return Q_;
}

// returns the euclidean norm of y - \hat y
template <typename PDE>
double SRPDE<PDE>::norm(const DMatrix<double>& obs, const DMatrix<double>& fitted) const {
  return (obs - fitted).squaredNorm();
}

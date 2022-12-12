// perform proper initialization of model
template <typename PDE>
void SRPDE<PDE>::init() {
  Psi(); // ask to iStatModel to compute \Psi matrix

  // if heteroscedastic observations are provided as datum, prepare weights matrix
  if(hasWeights()) W_ = W().asDiagonal();
  else W_ = DVector<double>::Ones(obs()).asDiagonal(); // homoscedastic observations
}

// finds a solution to the SR-PDE smoothing problem
template <typename PDE>
void SRPDE<PDE>::solve() {
  init(); // initialize the model, analyze input data
  
  // assemble system matrix for the nonparameteric part of the model
  SparseBlockMatrix<double,2,2>
    A(-PsiTD() * W_ * Psi(), lambda() * R1().transpose(),
      lambda() * R1(),       lambda() * R0()            );
  // cache system matrix for reuse
  A_ = A.derived();
  b_.resize(A_.rows());
  DVector<double> sol; // room for problem' solution
  
  if(!hasCovariates()){ // nonparametric case
    // rhs of SR-PDE linear system
    b_ << -PsiTD()*W_*y(),
          lambda()*u();
    
    // define system solver. Use a sparse solver
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver{};
    solver.compute(A_);
    // solve linear system A_*x = b_
    sol = solver.solve(b_);
    
    // store result of smoothing
    f_ = sol.head(A_.rows()/2);
  }else{ // parametric case
    if(!isAlloc(XTX_)){
      // compute q x q dense matrix W^T*W and its factorization
      XTX_ = X().transpose()*W_*X();
      invXTX_ = XTX_.partialPivLu();
    }
    // rhs of SR-PDE linear system
    b_ << -PsiTD()*lmbQ(y()), // -\Psi^T*D*Q*z
          lambda()*u();
    
    std::size_t q_ = X().cols(); // number of covariates
    // definition of matrices U and V  for application of woodbury formula
    DMatrix<double> U = DMatrix<double>::Zero(A_.rows(), q_);
    U.block(0,0, A_.rows()/2, q_) = PsiTD()*W_*X();
    DMatrix<double> V = DMatrix<double>::Zero(q_, A_.rows());
    V.block(0,0, q_, A_.rows()/2) = X().transpose()*W_*Psi();

    // Define system solver. Use SMW solver from NLA module
    SMW<> solver{};
    solver.compute(A_);
    // solve system Mx = b
    sol = solver.solve(U, XTX_, V, b_);
    // store result of smoothing 
    f_    = sol.head(A_.rows()/2);
    beta_ = invXTX_.solve(X().transpose()*W_)*(y() - Psi()*f_);
  }
  // store PDE misfit
  g_ = sol.tail(A_.rows()/2);
  return;
}

// it is asssumed that smooth has already been called on the model object
// computes fitted values \hat y = \Psi*f_ + X*beta_
template <typename PDE>
DMatrix<double> SRPDE<PDE>::fitted() {
  DMatrix<double> hat_y = Psi()*f_;
  // if the model has a parametric part, we need to sum its contribute
  if(hasCovariates()) hat_y += X()*beta_;
  return hat_y;
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
      T_ = Psi().transpose()*W_*Psi() + lambda()*R_;
    else // general case with covariates
      T_ = Psi().transpose()*lmbQ(Psi()) + lambda()*R_;
  }
  return T_;
}

// Q is computed on demand only when it is needed by GCV and cached for fast reacess (in general operations
// involving Q can be substituted with the more efficient routine lmbQ(), which is part of iRegressionModel interface)
template <typename PDE>
const DMatrix<double>& SRPDE<PDE>::Q() {
  if(!isAlloc(Q_)){ // Q is computed on request since not needed in general
    // compute Q = W(I - H) = W - W*X*(X*W*X^T)^{-1}*X^T*W
    Q_ = W_*(DMatrix<double>::Identity(obs(), obs()) - X()*invXTX_.solve(X().transpose()*W_));
  }
  return Q_;
}

// returns the euclidean norm of y - \hat y
template <typename PDE>
double SRPDE<PDE>::norm(const DMatrix<double>& obs, const DMatrix<double>& fitted) const {
  return (obs - fitted).squaredNorm();
}

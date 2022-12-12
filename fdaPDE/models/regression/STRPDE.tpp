template<typename PDE>
DMatrix<double> STRPDE<PDE, SpaceTimeSeparable>::fitted() {
  return DMatrix<double>::Zero(1,1);
}

template<typename PDE>
double STRPDE<PDE, SpaceTimeSeparable>::predict(const DVector<double>& covs, const std::size_t loc) const {
  return 0;
}

// finds a solution to the STR-PDE smoothing problem (separable penalization)
template <typename PDE>
void STRPDE<PDE, SpaceTimeSeparable>::solve() {
  // if heteroscedastic observations are provided as datum, prepare weights matrix
  if(hasWeights()) W_ = W().asDiagonal();
  else W_ = DVector<double>::Ones(locs()*this->time_domain().rows()).asDiagonal(); // homoscedastic observations
    
  // assemble system matrix for the nonparameteric part of the model
  SparseKroneckerProduct<> P = Kronecker(Pt(), pde_->R0());
  SparseBlockMatrix<double,2,2>
    A(-PsiTD()*Psi()-lambdaT()*P, lambdaS()*R1().transpose(),
      lambdaS()*R1(),             lambdaS()*R0()            );
  // cache system matrix for reuse
  A_ = A.derived();
  b_.resize(A_.rows());
  DVector<double> sol; // room for problem' solution
   
  if(!Base::hasCovariates()){ // nonparametric case
    // rhs of STR-PDE linear system
    b_ << -PsiTD()*y(),
      lambdaS()*u();
      
    // define system solver. Use a sparse solver
    Eigen::SparseLU<SpMatrix<double>, Eigen::COLAMDOrdering<int>> solver{};
    solver.compute(A_);
    // solve linear system A_*x = b_
    sol = solver.solve(b_);

    // store result of smoothing
    f_ = sol.head(A_.rows()/2);
  }else{ // parametric case
    if(!isAlloc(XTX_)){
      // compute q x q dense matrix X^T*W*X and its factorization
      XTX_ = X().transpose()*X();
      invXTX_ = XTX_.partialPivLu();
    }
    // rhs of STR-PDE linear system
    b_ << -PsiTD()*lmbQ(y()), // -\Psi^T*D*Q*z
      lambdaS()*u();

    std::size_t q_ = X().cols(); // number of covariates
    // definition of matrices U and V  for application of woodbury formula
    DMatrix<double> U = DMatrix<double>::Zero(A_.rows(), q_);
    U.block(0,0, A_.rows()/2, q_) = PsiTD()*X();
    DMatrix<double> V = DMatrix<double>::Zero(q_, A_.rows());
    V.block(0,0, q_, A_.rows()/2) = X().transpose()*Psi();

    // Define system solver. Use SMW solver from NLA module
    SMW<> solver{};
    solver.compute(A_);
    // solve system Mx = b
    sol = solver.solve(U, XTX_, V, b_);
    // store result of smoothing
    f_    = sol.head(A_.rows()/2);
    beta_ = invXTX_.solve(X().transpose())*(y() - Psi()*f_);
  }
  // store PDE misfit
  g_ = sol.tail(A_.rows()/2);
  return;
}

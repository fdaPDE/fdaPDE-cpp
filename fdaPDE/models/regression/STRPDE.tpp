// perform proper initialization and update of model. Computes quantites which can be reused
// across many calls to solve() and are **not affected by a change in the data**.
// It is implicitly called by ModelBase::init() as part of the initialization process.
// NB: a change in the smoothing parameter must trigger a re-initialization of the model
template <typename PDE, Sampling SamplingDesign>
void STRPDE<PDE, SpaceTimeSeparableTag, SamplingDesign, SolverType::Monolithic>::init_model() {
  // assemble system matrix for the nonparameteric part of the model
  SparseKroneckerProduct<> P = Kronecker(Pt(), pde().R0());  
  SparseBlockMatrix<double,2,2>
    A(-PsiTD()*W()*Psi()-lambdaT()*P, lambdaS()*R1().transpose(),
      lambdaS()*R1(),                 lambdaS()*R0()            );
  // cache system matrix for reuse
  A_ = A.derived();
  invA_.compute(A_);
  // prepare rhs of linear system
  b_.resize(A_.rows());
  b_.block(A_.rows()/2,0, A_.rows()/2,1) = lambdaS()*u();
  return;
}

// finds a solution to the STR-PDE smoothing problem (separable penalization)
template <typename PDE, Sampling SamplingDesign>
void STRPDE<PDE, SpaceTimeSeparableTag, SamplingDesign, SolverType::Monolithic>::solve() {
  BLOCK_FRAME_SANITY_CHECKS;
  DVector<double> sol; // room for problem' solution
   
  if(!Base::hasCovariates()){ // nonparametric case
    // update rhs of STR-PDE linear system
    b_.block(0,0, A_.rows()/2,1) = -PsiTD()*W()*y();
    // solve linear system A_*x = b_
    sol = invA_.solve(b_);
    // store result of smoothing
    f_ = sol.head(A_.rows()/2);
  }else{ // parametric case
    // update rhs of STR-PDE linear system
    b_.block(0,0, A_.rows()/2,1) = -PsiTD()*lmbQ(y()); // -\Psi^T*D*Q*z

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

// perform proper initialization and update of model. Computes quantites which can be reused
// across many calls to solve() and are **not affected by a change in the data**.
// It is implicitly called by ModelBase::init() as part of the initialization process.
// NB: a change in the smoothing parameter must trigger a re-initialization of the model
template <typename PDE, Sampling SamplingDesign>
void STRPDE<PDE, SpaceTimeParabolicTag, SamplingDesign, SolverType::Monolithic>::init_model() {
  // assemble system matrix for the nonparameteric part of the model
  SparseKroneckerProduct<> L_ = Kronecker(L(), pde().R0());
  SparseBlockMatrix<double,2,2>
    A(-PsiTD()*W()*Psi(),              lambdaS()*(R1() + lambdaT()*L_).transpose(),
      lambdaS()*(R1() + lambdaT()*L_), lambdaS()*R0()                             );
  // cache system matrix for reuse
  A_ = A.derived();
  invA_.compute(A_);
  // prepare rhs of linear system  
  b_.resize(A_.rows());
  b_.block(A_.rows()/2,0, A_.rows()/2,1) = lambdaS()*u();
  return;
}

// finds a solution to the STR-PDE smoothing problem (parabolic penalization, monolithic solution)
template <typename PDE, Sampling SamplingDesign>
void STRPDE<PDE, SpaceTimeParabolicTag, SamplingDesign, SolverType::Monolithic>::solve() {
  BLOCK_FRAME_SANITY_CHECKS;
  DVector<double> sol; // room for problem' solution
  
  if(!Base::hasCovariates()){ // nonparametric case
    // update rhs of STR-PDE linear system
    b_.block(0,0, A_.rows()/2,1) = -PsiTD()*W()*y();
    // solve linear system A_*x = b_
    sol = invA_.solve(b_);
    // store result of smoothing
    f_ = sol.head(A_.rows()/2);
  }else{ // parametric case
    // rhs of STR-PDE linear system
    b_.block(0,0, A_.rows()/2,1) = -PsiTD()*lmbQ(y()); // -\Psi^T*D*Q*z

    // definition of matrices U and V  for application of woodbury formula
    U_ = DMatrix<double>::Zero(A_.rows(), q());
    U_.block(0,0, A_.rows()/2, q()) = PsiTD()*W()*X();
    V_ = DMatrix<double>::Zero(q(), A_.rows());
    V_.block(0,0, q(), A_.rows()/2) = X().transpose()*W()*Psi();
    // solve system (A_ + U_*(X^T*W_*X)*V_)x = b using woodbury formula from NLA module
    sol = SMW<>().solve(invA_, U_, XtWX(), V_, b_);
    // store result of smoothing
    f_ = sol.head(A_.rows()/2);
    beta_ = invXtWX().solve(X().transpose()*W())*(y() - Psi()*f_);
  }
  // store PDE misfit
  g_ = sol.tail(A_.rows()/2);
  return;
}

// J(f,g) = \sum_{k=1}^m (z^k - \Psi*f^k)^T*(z^k - \Psi*f^k) + \lambda_S*(g^k)^T*(g^k)
template <typename PDE, Sampling SamplingDesign>
double STRPDE<PDE, SpaceTimeParabolicTag, SamplingDesign, SolverType::Iterative>::
J(const DMatrix<double>& f, const DMatrix<double>& g) const {
  double SSE = 0;
  // SSE = \sum_{k=1}^m (z^k - \Psi*f^k)^T*(z^k - \Psi*f^k)
  for(std::size_t t = 0; t < n_time(); ++t){
    SSE += (y(t) - Psi()*f.block(n_basis()*t,0, n_basis(),1)).squaredNorm();
  }
  return SSE + lambdaS()*g.squaredNorm();
}

// internal solve routine used by the iterative method
template <typename PDE, Sampling SamplingDesign>
void STRPDE<PDE, SpaceTimeParabolicTag, SamplingDesign, SolverType::Iterative>::solve
(std::size_t t, BlockVector<double>& f_new, BlockVector<double>& g_new) const {
  DVector<double> x = invA_.solve(b_);
  f_new(t) = x.topRows(n_basis()); g_new(t) = x.bottomRows(n_basis());
  return;
}

// finds a solution to the STR-PDE smoothing problem (parabolic penalization, iterative solution)
template <typename PDE, Sampling SamplingDesign>
void STRPDE<PDE, SpaceTimeParabolicTag, SamplingDesign, SolverType::Iterative>::solve() {  
  // compute starting point (f^(k,0), g^(k,0)) k = 1 ... m for iterative minimization of functional J(f,g)
  SparseBlockMatrix<double,2,2>
    A(PsiTD()*Psi(),   lambdaS()*R1().transpose(), 
      lambdaS()*R1(), -lambdaS()*R0()           );
  // cache system matrix and its factorization
  A_ = A.derived();
  invA_.compute(A_);
  b_.resize(A_.rows());
  
  // compute f^(k,0), k = 1 ... m as solution of Ax = b_(k)
  BlockVector<double> f_old(n_time(), n_basis());
  // solve n_time() space only linear systems
  for(std::size_t t = 0; t < n_time(); ++t){
    // right hand side at time step t
    b_ << PsiTD()*y(t), // should put W()
      lambdaS()*lambdaT()*u(t);
    // solve linear system Ax = b_(t) and store estimate of spatial field
    f_old(t) = invA_.solve(b_).head(A_.rows()/2);
  }
  
  // compute g^(k,0), k = 1 ... m as solution of the system
  //    G0 = [(\lambda_S*\lambda_T)/DeltaT * R_0 + \lambda_S*R_1^T]
  //    G0*g^(k,0) = \Psi^T*y^k + (\lambda_S*\lambda_T/DeltaT*R_0)*g^(k+1,0) - \Psi^T*\Psi*f^(k,0)
  SpMatrix<double> G0 = (lambdaS()*lambdaT()/DeltaT())*R0() + SpMatrix<double>((lambdaS()*R1()).transpose());
  Eigen::SparseLU<SpMatrix<double>, Eigen::COLAMDOrdering<int>> invG0;
  invG0.compute(G0); // compute factorization of matrix G0
  
  BlockVector<double> g_old(n_time(), n_basis());
  // solve n_time() distinct problems (in backward order)
  // at last step g^(t+1,0) is zero
  b_ = PsiTD()*(y(n_time()-1) - Psi()*f_old(n_time()-1)); 
  g_old(n_time()-1) = invG0.solve(b_);
  // general step
  for(int t = n_time()-2; t >= 0; --t){
    // compute rhs at time t: \Psi^T*y^t + (\lambda_S*\lambda_T/DeltaT*R_0)*g^(t+1,0) - \Psi^T*\Psi*f^(t,0)
    b_ = PsiTD()*(y(t) - Psi()*f_old(t)) + (lambdaS()*lambdaT()/DeltaT())*R0()*g_old(t+1);
    // solve linear system G0*g^(t,1) = b_t and store estimate of PDE misfit
    g_old(t) = invG0.solve(b_);
  }
  
  // initialize value of functional J to minimize
  double Jold = std::numeric_limits<double>::max();
  double Jnew = J(f_old.get(), g_old.get());
  std::size_t i = 1; // iteration number

  // build system matrix for the iterative scheme
  std::size_t N = n_basis();
  SparseBlockMatrix<double,2,2>
    M(SpMatrix<double>(N, N),            lambdaS()*lambdaT()/DeltaT()*R0(),
      lambdaS()*lambdaT()/DeltaT()*R0(), SpMatrix<double>(N, N)           );
  A_ = A_ + M.derived();
  invA_.compute(A_);
  b_.resize(A_.rows());
  
  // internal iteration variables
  BlockVector<double> f_new(n_time(), n_basis()), g_new(n_time(), n_basis());
  // iterative scheme for minimization of functional J
  while(i < max_iter_ && std::abs((Jnew-Jold)/Jnew) > tol_){
    // at step 0 f^(k-1,i-1) is zero    
    b_ << PsiTD()*y(0) + (lambdaS()*lambdaT()/DeltaT())*R0()*g_old(1),
      lambdaS()*u(0);
    // solve linear system
    solve(0, f_new, g_new);
    
    // general step
    for(std::size_t t = 1; t < n_time()-1; ++t){
      // \Psi^T*y^k   + (\lambdaS*\lambdaT/DeltaT)*R_0*g^(k+1,i-1),
      // \lambdaS*u^k + (\lambdaS*\lambdaT/DeltaT)*R_0*f^(k-1,i-1)
      b_ << PsiTD()*y(t) + (lambdaS()*lambdaT()/DeltaT())*R0()*g_old(t+1),
	lambdaS()*(lambdaT()/DeltaT()*R0()*f_old(t-1) + u(t));
      // solve linear system
      solve(t, f_new, g_new);
    }
    
    // at last step g^(k+1,i-1) is zero
    b_ << PsiTD()*y(n_time()-1),
      lambdaS()*(lambdaT()/DeltaT()*R0()*f_old(n_time()-2) + u(n_time()-1));
    // solve linear system
    solve(n_time()-1, f_new, g_new);

    // prepare for next iteration
    Jold = Jnew;
    f_old = f_new; g_old = g_new;
    Jnew = J(f_old.get(), g_old.get()); 
    i++;
  }
  
  // store solution
  f_ = f_old.get(); g_ = g_old.get();
  return;
}


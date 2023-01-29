// finds a solution to the STR-PDE smoothing problem (separable penalization)
template <typename PDE, Sampling SamplingDesign>
void STRPDE<PDE, SpaceTimeSeparableTag, SamplingDesign, SolverType::Monolithic>::solve() {
  // assemble system matrix for the nonparameteric part of the model
  SparseKroneckerProduct<> P = Kronecker(Pt(), pde().R0());
  SparseBlockMatrix<double,2,2>
    A(-PsiTD()*Psi()-lambdaT()*P, lambdaS()*R1().transpose(),
      lambdaS()*R1(),             lambdaS()*R0()            );
  // cache system matrix for reuse
  A_ = A.derived();
  invA_.compute(A_);
  b_.resize(A_.rows());
  DVector<double> sol; // room for problem' solution
   
  if(!Base::hasCovariates()){ // nonparametric case
    // rhs of STR-PDE linear system
    b_ << -PsiTD()*y(),
          lambdaS()*u();
      
    // solve linear system A_*x = b_
    sol = invA_.solve(b_);

    // store result of smoothing
    f_ = sol.head(A_.rows()/2);
  }else{ // parametric case
    // rhs of STR-PDE linear system
    b_ << -PsiTD()*lmbQ(y()), // -\Psi^T*D*Q*z
          lambdaS()*u();

    // definition of matrices U and V  for application of woodbury formula
    DMatrix<double> U = DMatrix<double>::Zero(A_.rows(), q());
    U.block(0,0, A_.rows()/2, q()) = PsiTD()*X();
    DMatrix<double> V = DMatrix<double>::Zero(q(), A_.rows());
    V.block(0,0, q(), A_.rows()/2) = X().transpose()*Psi();

    // Define system solver. Use SMW solver from NLA module
    SMW<> solver{};
    // solve system Mx = b
    sol = solver.solve(invA_, U, XtWX(), V, b_);
    // store result of smoothing
    f_    = sol.head(A_.rows()/2);
    beta_ = invXtWX().solve(X().transpose())*(y() - Psi()*f_);
  }
  // store PDE misfit
  g_ = sol.tail(A_.rows()/2);
  return;
}

// finds a solution to the STR-PDE smoothing problem (parabolic penalization, monolithic solution)
template <typename PDE, Sampling SamplingDesign>
void STRPDE<PDE, SpaceTimeParabolicTag, SamplingDesign, SolverType::Monolithic>::solve() {
  // set first n points of the estimation problem to the initial condition
  f_.resize((n_time()+1)*n_basis(), 1);
  f_.block(0,0, n_basis(),1) = s();
  
  // assemble system matrix for the nonparameteric part of the model
  SparseKroneckerProduct<> L_ = Kronecker(L(), pde().R0());
  SparseBlockMatrix<double,2,2>
    A(-PsiTD()*Psi(),                  lambdaS()*(R1() + lambdaT()*L_).transpose(),
      lambdaS()*(R1() + lambdaT()*L_), lambdaS()*R0()                             );
  // cache system matrix for reuse
  A_ = A.derived();
  invA_.compute(A_);
  b_.resize(A_.rows());
  DVector<double> sol; // room for problem' solution

  if(!Base::hasCovariates()){ // nonparametric case
    // rhs of STR-PDE linear system
    b_ << -PsiTD()*y(),
          lambdaS()*u();
      
    // solve linear system A_*x = b_
    sol = invA_.solve(b_);
    // store result of smoothing
    f_.block(n_basis(),0, n_time()*n_basis(),1) = sol.head(A_.rows()/2);
  }else{ // parametric case
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
    // solve system Mx = b
    sol = solver.solve(invA_, U, XtWX(), V, b_);
    // store result of smoothing
    f_.block(n_basis(),0, n_time()*n_basis(),1) = sol.head(A_.rows()/2);
    beta_ = invXtWX().solve(X().transpose())*(y() - Psi()*f_);
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


// finds a solution to the STR-PDE smoothing problem (parabolic penalization, iterative solution)
template <typename PDE, Sampling SamplingDesign>
void STRPDE<PDE, SpaceTimeParabolicTag, SamplingDesign, SolverType::Iterative>::solve() {  
  // compute starting point (f^(k,0), g^(k,0)) k = 1 ... m for iterative minimization of functional J(f,g)
  // assemble system matrix for the nonparameteric part of the model
  SparseBlockMatrix<double,2,2>
    A(-PsiTD()*Psi(), lambdaS()*R1().transpose(), 
      lambdaS()*R1(), lambdaS()*R0()            );
  // cache system matrix for reuse
  A_ = A.derived();
  invA_.compute(A_);
  b_.resize(A_.rows());
  
  // compute f^(k,0), k = 1 ... m as solution of Ax = b_(k)
  DMatrix<double> f_hat;
  f_hat.resize((n_time()+1)*n_basis(), 1);
  f_hat.block(0,0, n_basis(),1) = DMatrix<double>::Zero(n_basis(),1); // ????
  
  // solve n_time() space only linear systems
  for(std::size_t t = 0; t < n_time(); ++t){
    // right hand side at time step t
    b_ << -PsiTD()*y(t), // should put W()
      lambdaS()*lambdaT()*u(t);

    // solve linear system Ax = b_(t) and store estimate of spatial field
    f_hat.block(n_basis()*(t+1),0, n_basis(),1) = invA_.solve(b_).head(A_.rows()/2);
  }
  
  // compute g^(k,0), k = 1 ... m as solution of the system
  // [(\lambda_S*\lambda_T)/DeltaT * R_0 + \lambda_S*R_1^T]g^(k,0) = \Psi^T*z^k + (\lambda_S*\lambda_T/DeltaT*R_0)*g^(k+1,0) - \Psi^T*\Psi*f^(k,0)
  // need to evaluate R1().transpose() in a temporay to have compatible storage (see eigen docs)
  SpMatrix<double> G0 = (lambdaS()*lambdaT()/DeltaT())*R0() + SpMatrix<double>((lambdaS()*R1()).transpose());
  Eigen::SparseLU<SpMatrix<double>, Eigen::COLAMDOrdering<int>> invG0;
  invG0.compute(G0); // compute factorization of matrix G0
  
  DMatrix<double> g_hat;
  g_hat.resize((n_time()+1)*n_basis(), 1);
  g_hat.block(n_basis()*n_time(),0, n_basis(),1) = DMatrix<double>::Zero(n_basis(),1); // final condition for g is g = 0
  
  // solve n_time() distinct problems (in backward order, g follows a backward PDE)
  for(std::size_t t = n_time(); t > 0; --t){
    // compute rhs at time t: \Psi^T*z^t + (\lambda_S*\lambda_T/DeltaT*R_0)*g^(t+1,0) - \Psi^T*\Psi*f^(t,0)
    b_ = PsiTD()*(y(t-1) - Psi()*f_hat.block(n_basis()*t,0, n_basis(),1)) - // meno qui davanti a lambdaT() per problema R0_lambda
      (lambdaS()*lambdaT()/DeltaT())*R0()*g_hat.block(n_basis()*t,0, n_basis(),1);
    // solve linear system G0*g^(t,1) = b_t and store estimate of PDE misfit
    g_hat.block(n_basis()*(t-1),0, n_basis(),1) = invG0.solve(b_);
  }

  // initialize value of functional J to minimize
  double Jold = std::numeric_limits<double>::max();
  double Jnew = J(f_hat.bottomRows(n_time()*n_basis()), g_hat.topRows(n_time()*n_basis()));
  std::size_t i = 1; // iteration number

  // build system matrix for the iterative scheme
  SparseBlockMatrix<double,2,2>
    M(PsiTD()*Psi(), SpMatrix<double>(lambdaS()*R1().transpose()) + (lambdaS()*lambdaT()/DeltaT())*R0()   ,
      lambdaS()*R1()+(lambdaS()*lambdaT()/DeltaT())*R0(),   -lambdaS()*R0());
  SpMatrix<double> M_ = M.derived();
  Eigen::SparseLU<SpMatrix<double>> invM;
  invM.compute(M_);
  b_.resize(M_.rows());

  DMatrix<double> f_new, g_new;
  f_new = DMatrix<double>::Zero((n_time()+1)*n_basis(), 1);
  g_new = DMatrix<double>::Zero((n_time()+1)*n_basis(), 1);
  
  // iterative scheme for minimization of functional J
  while(i < max_iter_ && std::abs((Jnew-Jold)/Jnew) > tol_){

    for(std::size_t t = 1; t < n_time()+1; ++t){
      // update rhs of linear system as
      //    \Psi^T*y^k + (\lambdaS*\lambdaT/DeltaT)*R_0*g^(k+1)
      //    (\lambdaS*\lambdaT/DeltaT)*R_0*f^(k-1) + \lambdaS*u^k
      b_ << PsiTD()*y(t-1) + (lambdaS()*lambdaT()/DeltaT())*(-R0())*g_hat.block(n_basis()*t,0, n_basis(),1), // meno qui davanti a R0() per problema R0_lambda
	lambdaS()*(-lambdaT()/DeltaT()*R0()*f_hat.block(n_basis()*(t-1),0, n_basis(),1) - u(t-1)); // meno qui davanti a lambdaT() per problema R0_lambda

      // solve linear system
      DVector<double> sol = invM.solve(b_);
      f_new.block(n_basis()*t,0, n_basis(),1) = sol.topRows(n_basis());
      g_new.block(n_basis()*(t-1),0, n_basis(),1) = sol.bottomRows(n_basis());
    }
    // update value of functional J
    Jold = Jnew;
    f_hat = f_new; g_hat = g_new;
    Jnew = J(f_hat.bottomRows(n_time()*n_basis()), g_hat.topRows(n_time()*n_basis()));
    i++; // increase iteration count
  }
  // store solution
  f_ = f_hat.bottomRows(n_time()*n_basis());
  g_ = g_hat.topRows(n_time()*n_basis());
  return;
}

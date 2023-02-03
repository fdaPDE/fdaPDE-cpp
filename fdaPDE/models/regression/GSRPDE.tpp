// finds a solution to the GSR-PDE smoothing problem
template <typename PDE, typename RegularizationType, Sampling SamplingDesign,
	  SolverType Solver, typename Distribution>
void GSRPDE<PDE, RegularizationType, SamplingDesign, Solver, Distribution>::solve() {
  FPIRLS<decltype(*this), Distribution> fpirls(*this, tol_, max_iter_); // FPIRLS engine
  
  // call FPIRLS for minimization of the functional
  // \norm{V^{-1/2}(y - \mu)}^2 + \lambda \int_D (Lf - u)^2
  fpirls.compute();
  
  // fpirls converged: extract matrix P and solution estimates
  W_ = fpirls.weights().asDiagonal();
  f_ = fpirls.f();
  if(hasCovariates()) beta_ = fpirls.beta();
  return;
}

// required to support GCV based smoothing parameter selection
// in case of a GSRPDE model we have T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1), with Q = W*(I-H)
/*template <typename PDE, typename Distribution>
const DMatrix<double>& GSRPDE<PDE, Distribution>::T() {
  if(!isAlloc(T_)){ // compute only at first time
    // compute value of R = R1^T*R0^{-1}*R1, cache for possible reuse
    invR0_->compute(R0());
    R_ = R1().transpose()*invR0_->solve(R1());

    // compute and store matrix T for possible reuse
    if(!hasCovariates()) // case without covariates, Q is the identity matrix
      T_ = Psi().transpose()*Psi() + lambda()*R_;
    else // general case with covariates
      T_ = Psi().transpose()*lmbQ(Psi()) + lambda()*R_; // lmbQ will take care of the weight matrix W
  }
  return T_;
}

// Q is computed on demand only when it is needed by GCV and cached for fast reacess (in general operations
// involving Q can be substituted with the more efficient routine lmbQ(), which is part of iRegressionModel interface)
template <typename PDE, typename Distribution>
const DMatrix<double>& GSRPDE<PDE, Distribution>::Q() {
  if(!isAlloc(Q_)){ // Q is computed on request since not needed in general
    // compute Q = W(I - H) = W(I - (X*(X^T*W*X)^{-1}*X^T*W))
    Q_ = W_*(DMatrix<double>::Identity(obs(), obs()) - X()*invXTX_.solve(X().transpose()*W_));
  }
  return Q_;
}

// returns the deviance of y - \hat y induced by the specific distribution considered.
template <typename PDE, typename Distribution>
double GSRPDE<PDE, Distribution>::norm(const DMatrix<double>& obs, const DMatrix<double>& fitted) const {
  Distribution distr_{}; // define distribution object
  // total deviance computation
  double result = 0;
  for(std::size_t i = 0; i < obs.rows(); ++i)
    result += distr_.deviance(obs.coeff(i,0), fitted.coeff(i,0));
  return result;
  }*/

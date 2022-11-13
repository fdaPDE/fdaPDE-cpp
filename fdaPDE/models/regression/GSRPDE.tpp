// an efficient way to perform a left multiplication by Q implementing the following
//  given the design matrix W, the P matrix obtained from FPIRLS at convergence and x
//    compute v = W^T*P*x
//    solve Yz = v
//    return Px - PWz = P(I-H)x = Qx
// it is required to having assigned a design matrix W to the model before calling this method
template <typename PDE, typename Distribution>
DMatrix<double> GSRPDE<PDE, Distribution>::lmbQ(const DMatrix<double>& x){
  DMatrix<double> v = W().transpose()*P_*x; // W^T*P*x
  DMatrix<double> z = invWTW_.solve(v);  // (W^T*P*W)^{-1}*W^T*P*x
  // compute P*x - P*W*z = P*x - (P*W*(W^T*P*W)^{-1}*W^T*P)*x = P(I - H)*x = Q*x
  return P_*x - P_*W()*z;
}

// finds a solution to the GSR-PDE smoothing problem
template <typename PDE, typename Distribution>
void GSRPDE<PDE, Distribution>::solve() {
  // require iStatModel to compute Psi matrix
  this->Psi();
  
  // call FPIRLS for minimization of the functional
  // \norm{V^{-1/2}(y - \mu)}^2 + \lambda \int_D (Lf - u)^2
  fpirls.compute(*this);
  
  // fpirls converged: extract matrix P and solution estimates
  P_ = fpirls.weights().asDiagonal();
  f_ = fpirls.f();
  if(hasCovariates()) beta_ = fpirls.beta();
  return;
}

// it is asssumed that smooth has already been called on the model object
// computes fitted values \hat z = \Psi*f_ + W*beta_
template <typename PDE, typename Distribution>
DMatrix<double> GSRPDE<PDE, Distribution>::fitted() const {
  DMatrix<double> hat_z = Psi_*f_;
  // if the model has a parametric part, we need to sum its contribute
  if(hasCovariates()) hat_z += W()*beta_;
  return hat_z;
}

// compute prediction of model at new unseen data location (location equal to mesh node)
// W_{n+1}^T * \beta + f_*\psi(p_{n+1})
template <typename PDE, typename Distribution>
double GSRPDE<PDE, Distribution>::predict(const DVector<double>& covs, std::size_t loc) const {
  double prediction = f_.coeff(loc,0);
  // parametetric contribute of the model, if any is present
  if(hasCovariates()) prediction += (covs.transpose()*beta_).coeff(0,0);
  return prediction;
}

// required to support GCV based smoothing parameter selection
// in case of a GSRPDE model we have T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1), with Q = W*(I-H)
template <typename PDE, typename Distribution>
const DMatrix<double>& GSRPDE<PDE, Distribution>::T() {
  if(!isAlloc(T_)){ // compute only at first time
    // compute value of R = R1^T*R0^{-1}*R1, cache for possible reuse
    invR0_->compute(R0());
    R_ = R1().transpose()*invR0_->solve(R1());

    // compute and store matrix T for possible reuse
    if(!hasCovariates()) // case without covariates, Q is the identity matrix
      T_ = Psi_.transpose()*Psi_ + lambda_*R_;
    else // general case with covariates
      T_ = Psi_.transpose()*lmbQ(Psi_) + lambda_*R_; // lmbQ will take care of the weight matrix W
  }
  return T_;
}

// Q is computed on demand only when it is needed by GCV and cached for fast reacess (in general operations
// involving Q can be substituted with the more efficient routine lmbQ(), which is part of iRegressionModel interface)
template <typename PDE, typename Distribution>
const DMatrix<double>& GSRPDE<PDE, Distribution>::Q() {
  if(!isAlloc(Q_)){ // Q is computed on request since not needed in general
    // compute Q = P(I - H) = P(I - (W*(W^T*P*W)^{-1}*W^T*P))
    Q_ = P_*(DMatrix<double>::Identity(obs(), obs()) - W()*invWTW_.solve(W().transpose()*P));
  }
  return Q_;
}

// returns the deviance of y - \hat \mu induced by the specific distribution considered.
template <typename PDE, typename Distribution>
double GSRPDE<PDE, Distribution>::norm(const DMatrix<double>& obs, const DMatrix<double>& fitted) const {
  Distribution distr_{}; // define distribution object
  // total deviance computation
  double result = 0;
  for(std::size_t i = 0; i < obs.rows(); ++i)
    result += distr_.deviance(obs.coeff(i,0), fitted.coeff(i,0));
  return result;
}

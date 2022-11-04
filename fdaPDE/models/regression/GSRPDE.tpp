// an efficient way to perform a left multiplication by Q implementing the following
//  given the design matrix W, the P matrix obtained from FPIRLS at convergence and x
//    compute v = W^T*P*x
//    solve Yz = v
//    return Px - PWz = P(I-H)x = Qx
// it is required to having assigned a design matrix W to the model before calling this method
template <typename PDE, typename Dist>
DMatrix<double> GSRPDE<PDE, Dist>::lmbQ(const DMatrix<double>& x){
  DMatrix<double> v = W().transpose()*P_*x; // W^T*P*x
  DMatrix<double> z = invWTW_.solve(v);  // (W^T*P*W)^{-1}*W^T*P*x
  // compute P*x - P*W*z = P*x - (P*W*(W^T*P*W)^{-1}*W^T*P)*x = P(I - H)*x = Q*x
  return P_*x - P_*W()*z;
}

// finds a solution to the GSR-PDE smoothing problem
template <typename PDE, typename Dist>
void GSRPDE<PDE, Dist>::solve() {
  // require iStatModel to compute Psi matrix
  this->Psi();
  
  // call to FPIRLS for the minimization of
  // \norm{V^{-1/2}(y - \mu)}^2 + \lambda \int_D (Lf - u)^2
  fpirls.compute(*this);

  // at fpirls convergence we can get matrix P and solution estimates
  P_ = fpirls.weights().asDiagonal();
  f_ = fpirls.f();
  if(hasCovariates()) beta_ = fpirls.beta();

  return;
}

// it is asssumed that smooth has already been called on the model object
// computes fitted values \hat z = \Psi*f_ + W*beta_
template <typename PDE, typename Dist>
DMatrix<double> GSRPDE<PDE, Dist>::fitted() const {
  DMatrix<double> hat_z = Psi_*f_;
  // if the model has a parametric part, we need to sum its contribute
  if(hasCovariates())
    hat_z += W()*beta_;
  return hat_z;
}

// compute prediction of model at new unseen data location (location equal to mesh node)
// W_{n+1}^T * \beta + f_*\psi(p_{n+1})
template <typename PDE, typename Dist>
double GSRPDE<PDE, Dist>::predict(const DVector<double>& covs, std::size_t loc) const {
  double prediction = f_.coeff(loc,0);
  // parametetric contribute of the model, if any is present
  if(hasCovariates())
    prediction += (covs.transpose()*beta_).coeff(0,0);
  return prediction;
}

/*

// required to support GCV based smoothing parameter selection
// in case of an SRPDE model we have T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
template <typename PDE>
std::shared_ptr<DMatrix<double>> SRPDE<PDE>::T() {
  // compute value of R = R1^T*R0^{-1}*R1, cache for possible reuse
  invR0_->compute(R0());
  R_ = std::make_shared<DMatrix<double>>( R1().transpose()*invR0_->solve(R1()) );

  // compute and store matrix T for possible reuse
  if(!hasCovariates()) // case without covariates, Q is the identity matrix
    T_ = std::make_shared<DMatrix<double>>
      ( Psi_.transpose()*Psi_ + lambda_*(*R_) );
  else // general case with covariates
    T_ = std::make_shared<DMatrix<double>>
      ( Psi_.transpose()*lmbQ(Psi_) + lambda_*(*R_) );

  // return pointer to T
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
*/

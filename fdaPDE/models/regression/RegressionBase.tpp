// an efficient way to perform a left multiplication by Q implementing the following
//  given the design matrix X, the weight matrix W and x
//    compute v = X^T*W*x
//    solve Yz = v
//    return Wx - WXz = W(I-H)x = Qx
// it is required to having assigned a design matrix X to the model before calling this method
template <typename Model>
DMatrix<double> RegressionBase<Model>::lmbQ(const DMatrix<double>& x) const {
  if(!hasCovariates()) return W_*x;
  DMatrix<double> v = X().transpose()*W_*x; // X^T*W*x
  DMatrix<double> z = invXtWX_.solve(v);  // (X^T*W*X)^{-1}*X^T*W*x
  // compute W*x - W*X*z = W*x - (W*X*(X^T*W*X)^{-1}*X^T*W)*x = W(I - H)*x = Q*x
  return W_*x - W_*X()*z;
}

// initialization stuffs depending on the data
template <typename Model>
void RegressionBase<Model>::init_data() {
  // default to homoscedastic observations
  DVector<double> W = DVector<double>::Ones(Base::n_locs());
  if(hasWeights()) // update observations' weights if provided
    W = df_.template get<double>(WEIGHTS_BLK).col(0);
  W_ = W.asDiagonal();
  // model is semi-parametric
  if(hasCovariates()){
    // compute q x q dense matrix X^T*W*X and its factorization
    XtWX_ = X().transpose()*W_*X();
    invXtWX_ = XtWX_.partialPivLu();
  }
}

// computes fitted values \hat y = \Psi*f_ + X*beta_
template <typename Model>
DMatrix<double> RegressionBase<Model>::fitted() const {
  DMatrix<double> hat_y = Psi(not_nan())*f_;
  if(hasCovariates()) hat_y += X()*beta_;
  return hat_y;
}

// missing data logic
template <typename Model>
void RegressionBase<Model>::init_nan() {
  // derive missingness pattern
  nan_idxs_.clear(); // empty nan indexes set
  for(std::size_t i = 0; i < n_obs(); ++i){
    if(std::isnan(y()(i,0))){ // requires -ffast-math compiler flag to be disabled
      nan_idxs_.insert(i);
      df_.template get<double>(OBSERVATIONS_BLK)(i,0) = 0.0; // zero out NaN
    }
  }
  // matrix B assembly logic (set to zero rows corresponding to missing observations)
  if(has_nan()){
    // reserve space
    std::size_t n = Psi(not_nan()).rows();
    std::size_t N = Psi(not_nan()).cols();
    B_.resize(n, N);
    // triplet list to fill sparse matrix
    std::vector<fdaPDE::Triplet<double>> tripletList;
    tripletList.reserve(n*N);
    for (int k = 0; k < Psi(not_nan()).outerSize(); ++k)
      for (SpMatrix<double>::InnerIterator it(Psi(not_nan()),k); it; ++it){
	if(nan_idxs_.find(it.row()) == nan_idxs_.end()){
	  // no missing data at this location
	  tripletList.emplace_back(it.row(), it.col(), it.value());
	}
      }
    // finalize construction
    B_.setFromTriplets(tripletList.begin(), tripletList.end());
    B_.makeCompressed();
  }
  return;
}

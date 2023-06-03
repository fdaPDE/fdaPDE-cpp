// solution in case of fixed \lambda
template <typename PDE, typename RegularizationType,
	  typename SamplingDesign, typename lambda_selection_strategy>
void FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>::solve_(fixed_lambda) {
  ProfilingEstimation<decltype(*this)> pe(*this, tol_, max_iter_);
  DMatrix<double> X = data().template get<double>(OBSERVATIONS_BLK);
  // Principal Components computation
  for(std::size_t i = 0; i < n_pc_; i++){
    // find vectors s,f minimizing \norm_F{Y - s^T*f}^2 + (s^T*s)*P(f) fixed \lambda
    pe.compute(X, lambda()); 
    loadings_.col(i) = pe.f_n(); scores_.col(i) = pe.s();
    // subtract computed PC from data	
    X -= scores_.col(i)*loadings_.col(i).transpose();
  }
  return;
}

// best \lambda for PC choosen according to GCV index
template <typename PDE, typename RegularizationType,
	  typename SamplingDesign, typename lambda_selection_strategy>
void FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>::solve_(gcv_lambda_selection) {
  ProfilingEstimation<decltype(*this)> pe(*this, tol_, max_iter_);
  DMatrix<double> X = data().template get<double>(OBSERVATIONS_BLK);
  // wrap GCV into a ScalarField accepted by OPT module
  const std::size_t n_lambda = n_smoothing_parameters<RegularizationType>::value;
  ScalarField<n_lambda> f;
  f = [&pe, &X](const SVector<n_lambda>& p) -> double {
    // find vectors s,f minimizing \norm_F{Y - s^T*f}^2 + (s^T*s)*P(f) fixed \lambda = p
    pe.compute(X, p);
    return pe.gcv(); // return GCV at convergence
  };
  GridOptimizer<n_lambda> opt; // optimization algorithm
  // Principal Components computation
  for(std::size_t i = 0; i < n_pc_; i++){
    opt.optimize(f, lambdas()); // select optimal \lambda for i-th PC
    // compute and store results given estimated optimal \lambda
    pe.compute(X, opt.optimum());
    loadings_.col(i) = pe.f_n(); scores_.col(i) = pe.s();
    // subtract computed PC from data
    X -= scores_.col(i)*loadings_.col(i).transpose();
  }
  return;
}

// best \lambda for PC choosen according to K-fold CV strategy, uses the reconstruction error on test set as CV score
template <typename PDE, typename RegularizationType,
	  typename SamplingDesign, typename lambda_selection_strategy>
void FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>::solve_(kcv_lambda_selection) {
  ProfilingEstimation<decltype(*this)> pe(*this, tol_, max_iter_);
  // number of smoothing parameters
  const std::size_t n_lambda = n_smoothing_parameters<RegularizationType>::value;

  // routine executed by the CV-engine to produce the model score
  std::function<double(DVector<double>, BlockFrame<double, int>, BlockFrame<double, int>)> cv_score = 
    [&pe, this](const DVector<double>& lambda,
		const BlockFrame<double, int>& train_df,
		const BlockFrame<double, int>& test_df) -> double {
    // get references to train and test sets
    //const DMatrix<double>& X_train = train_df.get<double>(OBSERVATIONS_BLK);
    const DMatrix<double>& X_test  = test_df. get<double>(OBSERVATIONS_BLK);
    SVector<n_lambda> p(lambda.data());
    // find vectors s,f minimizing \norm_F{Y - s^T*f}^2 + (s^T*s)*P(f)
    pe.compute(train_df, p);
    // compute scores (X_test * f_n)/(\norm{f_n} + \lambda*P(f)) on test set
    DVector<double> s = X_test*pe.f_n();
    if(this->has_nan()){
      std::size_t n = 0; // number of not-NaN points in test set
      double pen = p[0]*(pe.g().dot(R0()*pe.g())); // \lambda*P(f)
      for(std::size_t i = 0; i < s.rows(); ++i){
	double sum_f = 0;
	// compute \norm{f_n} for i-th subject and update count of not-NaN points
	int stat_unit = test_df.get<int>(INDEXES_BLK)(i,0);
	for(std::size_t j = 0; j < n_locs(); ++j){
	  if(nan_idxs()[stat_unit].find(j) == nan_idxs()[stat_unit].end()) {
	    sum_f += std::pow(pe.f_n()[j], 2); n++;
	  }
	}
	s[i] /= (sum_f + pen);
      }
      // evaluate reconstruction error on test set
      double reconstruction_error = 0;
      for(std::size_t i = 0; i < X_test.rows(); ++i){
	int stat_unit = test_df.get<int>(INDEXES_BLK)(i,0);
	for(std::size_t j = 0; j < n_locs(); ++j)
	  if(nan_idxs()[stat_unit].find(j) == nan_idxs()[stat_unit].end()){
	    reconstruction_error += std::pow(X_test(i,j) - s[i]*pe.f_n()[j], 2);
	  }
      }
      return reconstruction_error/n;
    }else{
      // normalize scores with respect to \norm{f_n} + \lambda*P(f)
      s /= pe.f_n().squaredNorm() + p[0]*(pe.g().dot(R0()*pe.g())); 
      return (X_test - s*pe.f_n().transpose()).squaredNorm()/X_test.size(); // reconstruction error on test set
    }
  };
  
  // define K-fold algorithm
  KFoldCV cv(10); // allow user-defined number of folds!
  std::vector<DVector<double>> lambdas_;
  lambdas_.reserve(lambdas().size());
  for(const auto& l : lambdas()) lambdas_.emplace_back(Eigen::Map<const DVector<double>>(l.data(), n_lambda, 1));
  // Principal Components computation
  BlockFrame<double, int> data_ = data();//.shuffle(); // copy data to avoid side effects on caller state. ** can remove this **
  for(std::size_t i = 0; i < n_pc_; i++){
    cv.compute(lambdas_, data_, cv_score, false); // select optimal smoothing level
    // compute and store results given estimated optimal \lambda
    pe.compute(data_, cv.optimum());
    std::cout << cv.optimum() << std::endl;
    loadings_.col(i) = pe.f_n(); scores_.col(i) = pe.s();
    // subtract computed PC from data
    data_.get<double>(OBSERVATIONS_BLK) -= pe.s()*pe.f_n().transpose();
    // set NaN back
    for(std::size_t i = 0; i < data_.rows(); ++i){
      for(auto j : nan_idxs()[i])
	data_.get<double>(OBSERVATIONS_BLK)(i,j) = 0.0;
    }
  }
  return;
}

// finds solution to fPCA problem, dispatch to solver depending on \lambda selection criterion
template <typename PDE, typename RegularizationType,
	  typename SamplingDesign, typename lambda_selection_strategy>
void FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>::solve() {
  // pre-allocate space
  loadings_.resize(X().cols(), n_pc_);
  scores_.resize  (X().rows(), n_pc_);

  // dispatch to desired solution strategy
  solve_(lambda_selection_strategy());
  return;
}

// solution in case of fixed \lambda
template <typename PDE, typename RegularizationType,
	  Sampling SamplingDesign, typename lambda_selection_strategy>
void FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>::solve_(fixed_lambda) {
  // define internal solver
  FPIREM<decltype(*this)> solver(*this);
  solver.setTolerance(tol_);
  solver.setMaxIterations(max_iter_);

  DMatrix<double> Y = y().transpose(); // move original data to a n_subject() x n_obs() matrix
  // Principal Components computation
  for(std::size_t i = 0; i < n_pc_; i++){
    solver.setData(Y);
    solver.setLambda(lambda()); // take \lambda from FPCA object
    solver.solve(); // find minimum of \norm_F{Y - s^T*f}^2 + (s^T*s)*P(f) fixed \lambda
    // store result
    loadings_.col(i) = solver.loadings();
    scores_.col(i)   = solver.scores();
    // subtract computed PC from data	
    Y -= loadings_.col(i)*scores_.col(i).transpose();		
  }
  return;
}

// best \lambda for PC choosen according to GCV index
template <typename PDE, typename RegularizationType,
	  Sampling SamplingDesign, typename lambda_selection_strategy>
void FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>::solve_(gcv_lambda_selection) {
  // define internal solver
  typedef typename FPIREM<decltype(*this)>::SmootherType SmootherType; // solver used to smooth loadings
  FPIREM<decltype(*this)> solver(*this);
  solver.setTolerance(tol_);
  solver.setMaxIterations(max_iter_);

  BlockFrame<double, int> df;
  df.insert<double>(OBSERVATIONS_BLK, y().transpose()); // move original data to a n_subject() x n_obs() matrix
  // define GCV objective for internal smoothing model used by FPIREM
  GCV<SmootherType, fdaPDE::calibration::StochasticEDF<SmootherType>> GCV(solver.smoother(), 100);
  // Principal Components computation
  for(std::size_t i = 0; i < n_pc_; i++){
    solver.setData(df);
    // initialize GCV minimization
    double min_gcv = std::numeric_limits<double>::max();
    SVector<model_traits<SmootherType>::n_lambda> min_lambda;

    // explore vector of \lambda
    for(auto l : lambda_vect_){
      solver.setLambda(l);
      solver.solve(); // find minimum of \norm_F{Y - s^T*f}^2 + (s^T*s)*P(f) fixed \lambda
      // compute GCV of smoother for this value of \lambda
      double gcv = GCV.eval();
      if(gcv < min_gcv){ // update found minimium
	min_gcv = gcv;
	min_lambda = l;
	// store result
	loadings_.col(i) = solver.loadings();
	scores_.col(i)   = solver.scores();
      }
    }
    // subtract computed PC from data
    df.get<double>(OBSERVATIONS_BLK) -= scores_.col(i)*loadings_.col(i).transpose();
  }
  return;      
}

// best \lambda for PC choosen according to K fold CV error
template <typename PDE, typename RegularizationType,
	  Sampling SamplingDesign, typename lambda_selection_strategy>
void FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>::solve_(kcv_lambda_selection) {
  // define internal solver
  FPIREM<decltype(*this)> solver(*this);
  solver.setTolerance(tol_);
  solver.setMaxIterations(max_iter_);

  BlockFrame<double, int> df;
  df.insert<double>(OBSERVATIONS_BLK, y().transpose()); // n_subject() x n_obs() matrix
  KFoldCV<decltype(solver)> lambda_selector(solver, 4);
  
  for(std::size_t i = 0; i < n_pc_; i++){
    solver.setData(df);
    SVector<1> optimal_lambda = lambda_selector.compute(lambda_vect_, PCScoreCV());
    // compute result given estimated optimal \lambda
    solver.setLambda(optimal_lambda);
    solver.solve(); // find minimum of \norm_F{Y - s^T*f}^2 + (s^T*s)*P(f)
    // store result
    loadings_.col(i) = solver.loadings();
    scores_.col(i)   = solver.scores();
    // subtract computed PC from data
    df.get<double>(OBSERVATIONS_BLK) -= scores_.col(i)*loadings_.col(i).transpose();
  }

  std::cout << scores_ << std::endl;
  return;
}

// finds solution to fPCA problem, dispatch to solver depending on \lambda selection criterion
template <typename PDE, typename RegularizationType,
	  Sampling SamplingDesign, typename lambda_selection_strategy>
void FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>::solve() {
  // pre-allocate space
  loadings_.resize(y().rows(), n_pc_);
  scores_.resize(y().cols(),   n_pc_);
  // dispatch to desired solution strategy
  solve_(lambda_selection_strategy());
}

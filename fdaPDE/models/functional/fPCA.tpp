// solution in case of fixed \lambda
template <typename PDE, typename RegularizationType,
	  Sampling SamplingDesign, typename lambda_selection_strategy>
void FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>::solve_(fixed_lambda) {
  // define internal solver
  FPIREM<decltype(*this)> solver(*this);
  solver.setTolerance(tol_);
  solver.setMaxIterations(max_iter_);
  
  BlockFrame<double, int> df;
  df.insert<double>(OBSERVATIONS_BLK, y().transpose()); // move original data to a n_subject() x n_obs() matrix
  // Principal Components computation
  for(std::size_t i = 0; i < n_pc_; i++){
    solver.setData(df);
    solver.setLambda(lambda()); // take \lambda from FPCA object
    solver.solve(); // find minimum of \norm_F{Y - s^T*f}^2 + (s^T*s)*P(f) fixed \lambda

    // store result
    loadings_.col(i) = solver.loadings();
    scores_.col(i)   = solver.scores();
    // subtract computed PC from data	
    df.get<double>(OBSERVATIONS_BLK) -= scores_.col(i)*loadings_.col(i).transpose();
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
  GCV<SmootherType, StochasticEDF<SmootherType>> GCV(solver.smoother(), 100);
  // wrap GCV into a ScalarField accepted by the OPT module
  ScalarField<model_traits<SmootherType>::n_lambda> f;
  f = [&GCV, &solver](const SVector<model_traits<SmootherType>::n_lambda>& p) -> double {
    // set \lambda and solve smoothing problem on solver
    std::cout << "lambda: " << p << std::endl;
    solver.setLambda(p);
    solver.solve();
    // return evaluation of GCV at point
    return GCV.eval();
  };
    
  // define GCV optimization algorithm
  GridOptimizer<model_traits<SmootherType>::n_lambda> opt;
  // Principal Components computation
  for(std::size_t i = 0; i < n_pc_; i++){
    solver.setData(df);
    opt.optimize(f, lambda_vect_); // select optimal lambda for this PC
    // compute result given estimated optimal \lambda
    solver.setLambda(opt.optimum());
    solver.solve(); // find minimum of \norm_F{Y - s^T*f}^2 + (s^T*s)*P(f)
    // store result
    loadings_.col(i) = solver.loadings();
    scores_.col(i)   = solver.scores();    
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
  KFoldCV<decltype(solver)> lambda_selector(solver, 5);
  
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

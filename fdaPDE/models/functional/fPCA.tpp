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
    loadings_.col(i) = pe.f(); scores_.col(i) = pe.s();
    // subtract computed PC from data	
    X -= loadings_.col(i)*scores_.col(i).transpose();
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
    loadings_.col(i) = pe.f(); scores_.col(i) = pe.s();
    // subtract computed PC from data
    X -= loadings_.col(i)*scores_.col(i).transpose();
  }
  return;
}

// finds solution to fPCA problem, dispatch to solver depending on \lambda selection criterion
template <typename PDE, typename RegularizationType,
	  typename SamplingDesign, typename lambda_selection_strategy>
void FPCA<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>::solve() {
  // pre-allocate space
  loadings_.resize(y().rows(), n_pc_);
  scores_.resize(y().cols(),   n_pc_);

  // dispatch to desired solution strategy
  solve_(lambda_selection_strategy());
}

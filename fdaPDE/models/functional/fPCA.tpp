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

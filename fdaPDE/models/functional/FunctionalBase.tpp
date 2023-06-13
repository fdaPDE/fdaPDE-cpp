// missing data logic
template <typename Model>
void FunctionalBase<Model>::init_nan() {
  nan_idxs_.clear(); // clean previous missingness structure
  nan_idxs_.resize(n_stat_units()); B_.resize(n_stat_units());
  // \Psi matrix dimensions
  std::size_t n = Psi(not_nan()).rows();
  std::size_t N = Psi(not_nan()).cols();
  // for i-th statistical unit, analyze missingness structure and set \Psi_i
  for(std::size_t i = 0; i < n_stat_units(); ++i){
    // derive missingness pattern for i-th statistical unit
    for(std::size_t j = 0; j < n_locs(); ++j){
      if(std::isnan(X()(i,j))) // requires -ffast-math compiler flag to be disabled
	nan_idxs_[i].insert(j);
    }

    // NaN detected for this unit, start assembly
    if(!nan_idxs_[i].empty()){
      for(std::size_t i = 0; i < n_stat_units(); ++i){
	B_[i].resize(n, N); // reserve space
	std::vector<fdaPDE::Triplet<double>> tripletList;
	tripletList.reserve(n*N);
	for(int k = 0; k < Psi(not_nan()).outerSize(); ++k)
	  for(SpMatrix<double>::InnerIterator it(Psi(not_nan()),k); it; ++it){
	    if(nan_idxs_[i].find(it.row()) == nan_idxs_[i].end())
	      // no missing data at this location for i-th statistical unit
	      tripletList.emplace_back(it.row(), it.col(), it.value());
	  }
	// finalize construction
	B_[i].setFromTriplets(tripletList.begin(), tripletList.end());
	B_[i].makeCompressed();
      }	  
    }
    // otherwise no matrix is assembled, full \Psi is returned by Psi(std::size_t) getter
  }
  return;
}

// true if there are missing data in any of the statistical unit
template <typename Model>
bool FunctionalBase<Model>::has_nan() const { 
  for(auto s : nan_idxs_) {
    if(!s.empty()) return true;
  }
  return false;
}

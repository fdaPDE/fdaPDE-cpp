// set model's data from blockframe
template <typename E>
void iStatModel<E>::setData(const BlockFrame<double, int>& df) {
  df_ = df;
  // insert an index row (if not yet present)
  if(!df_.hasBlock(INDEXES_BLK)){
    std::size_t n = df_.rows();
    DMatrix<int> idx(n,1);
    for(std::size_t i = 0; i < n; ++i) idx(i,0) = i;
    df_.insert(INDEXES_BLK, idx);
  }
  return;
}

// return the sampling strategy adopted for the model.
template <typename E>
SamplingStrategy iStatModel<E>::sampling() const {
  if  (df_.hasBlock(AREAL_BLK)) // subdomains are given as datum
    return   SamplingStrategy::Areal;
  else{ // fallback to a Geostatistical sampling
    if(df_.hasBlock(LOCATIONS_BLK))
      return SamplingStrategy::GeostatisticalAtLocations;
    else // if no information is provided at all, assume data sampled at mesh nodes
      return SamplingStrategy::GeostatisticalAtNodes;
  }
}
// return a reference to the search engine over the problem's domain. Initializes it if still not available
template <typename E>
const ADT<iStatModel<E>::M, iStatModel<E>::N, iStatModel<E>::K>&
iStatModel<E>::searchEngine() {
  // if still not available, initialize the search engine
  if(searchEngine_ == nullptr){
    searchEngine_ = std::make_shared<ADT<M,N,K>>(pde_->domain());
  }
  return *searchEngine_;
}

// compute n x N \Psi matrix depending on sampling strategy
template <typename E>
const SpMatrix<double>& iStatModel<E>::Psi() {
  if(!isAlloc(Psi_)){ // compute \Psi if not already available
    // preallocate space for Psi matrix
    // detect number of observations and number of nodes from model object
    std::size_t n = obs();
    std::size_t N = loc();    
    Psi_.resize(n, N);    
    // triplet list to fill sparse matrix
    std::vector<fdaPDE::Triplet<double>> tripletList;
    tripletList.reserve(n*N);
    
    switch(sampling()){
    case SamplingStrategy::GeostatisticalAtNodes:
      // if data locations are equal to mesh nodes then \Psi is the identity matrix.
      //   (\psi_i(p_i) = 1 and \psi_i(p_j) = 0 \forall i \neq j)
      // if data are observed in a subset of the spatial locations then \Psi will have some of its columns set to zero
      //   (\psi_i(p_j) = 0 \forall j \in {1, ..., n} such that no data is observed at location i)
      for(std::size_t i = 0; i < n; ++i){
	tripletList.push_back(fdaPDE::Triplet<double>(i, idx()(i,0), 1.0));
      }
      break;
    case SamplingStrategy::GeostatisticalAtLocations:
      // general case in which locations are provided as plain coordinates to the model. In this case \Psi has no particular structure
      for(std::size_t i = 0; i < locations().rows(); ++i){ // cycle over all locations
	SVector<M> p_i(locations().row(i));
	// search element containing the point
	auto e = searchEngine().search(p_i);
	// update \Psi matrix
	for(std::size_t j = 0; j < pde_->basis()[e->ID()].size(); ++j){
	  std::size_t h = e->nodeIDs()[j]; // column index of \Psi matrix
	  // extract \phi_h from basis
	  auto psi_h = pde_->basis()[e->ID()][j];
	  // evaluate \phi_h(p_i) (value of the basis function centered in mesh node h and evaluated in point p_i)
	  tripletList.push_back(fdaPDE::Triplet<double>(i, h, psi_h(p_i)));
	}
      }
      break;
    case SamplingStrategy::Areal:
      DVector<double> D; // store measure of subdomains, this will be ported to a diagonal matrix at the end
      D.resize(subdomains().rows());
      
      // incidence matrix is a dxM sparse matrix (d number of subdomains, M number of elements) where D_{ij} = 1 \iff element j belongs
      // to subdomain i. In this case matrix \Psi is such that [\Psi]_{ij} = \int_{D_i} \psi_j = \sum_{e \in D_i} \int_{e} \psi_j
      std::size_t tail = 0;
      for(std::size_t k = 0; k < subdomains().rows(); ++k){
	std::size_t head = 0;
	double Di = 0; // measure of subdomain D_i
	for(std::size_t l = 0; l < subdomains().cols(); ++l){
	  if(subdomains()(k,l) == 1){ // element with ID l belongs to k-th subdomain
	    // get element with this ID
	    auto e = pde_->domain().element(l);
	    // compute \int_e \phi_h \forall \phi_h defined on e
	    std::size_t j = 0;
	    for(const auto& phi : pde_->basis()[e->ID()]){
	      std::size_t h = phi.node(); // if we write the finite element as \phi_h, this is h
	      // evaluate \int_e \phi_h and insert in tripletList. summation is implicitly resolved by Eigen::setFromTriplets
	      tripletList.push_back(fdaPDE::Triplet<double>(k, h, pde_->integrator().integrate(*e, phi)));
	      head++, j++; // increment counters
	    }
	    Di += e->measure(); // update measure of subdomain D_i
	  }
	}
	// divide each \int_{D_i} \psi_j by the measure of subdomain D_i
	for(std::size_t j = 0; j < head; ++j){
	  tripletList[tail + j].value() /= Di;
	}
	D[k] = Di; // store measure of subdomain
	tail += head;
      }
      // store diagonal matrix D_ = diag(D1, D2, ... ,Dd)
      D_ = D.asDiagonal();
      break;
    }
    // finalize construction
    Psi_.setFromTriplets(tripletList.begin(), tripletList.end());
    Psi_.makeCompressed();
    if(isAlloc(D_)) PsiTD_ = Psi_.transpose()*D_; // store \Psi^T*D in case of areal sampling
  }
  return Psi_;
}

template <typename E>
DMatrix<double> iStatModel<E>::lmbPsi(const DMatrix<double>& x) const {
  // compute dimensions of resulting matrix
  std::size_t n = Psi_.rows();
  std::size_t m = x.cols();
  // preallocate space for n x m result matrix
  DMatrix<double> result(n,m);
  // if data are sampled at mesh nodes (or a subset of them) then \Psi is a permutation matrix
  if(dataAtNodes()){
    // just permute input matrix columns
    for(std::size_t k = 0; k < Psi_.outerSize(); ++k){
      for (SpMatrix<double>::InnerIterator it(Psi_, k); it; ++it){
	result.row(it.row()) = x.row(it.col());
      }
    }
  }else{
    // in the general case no optimization can be put in place
    result = Psi_*x;
  }
  return result;
}


// set boundary conditions on problem's linear system
template <typename E>
void iStatModel<E>::setDirichletBC(SpMatrix<double>& A, DMatrix<double>& b){
  std::size_t n = A.rows()/2;

  for(std::size_t i = 0; i < n; ++i){
    if(pde_->domain().isOnBoundary(i)){
      A.row(i) *= 0;       // zero all entries of this row
      A.coeffRef(i,i) = 1; // set diagonal element to 1 to impose equation u_j = b_j

      A.row(i+n) *= 0;
      A.coeffRef(i+n,i+n) = 1;

      // boundaryDatum is a pair (nodeID, boundary value)
      double boundaryDatum = pde_->boundaryData().empty() ? 0 : pde_->boundaryData().at(i)[0];
      b.coeffRef(i,0) = boundaryDatum; // impose boundary value
      b.coeffRef(i+n, 0) = 0;
    }
  }
  return;
}

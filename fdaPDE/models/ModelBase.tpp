// perform initialization of model object, must be called before call to .solve()
template <typename Model>
void ModelBase<Model>::init(){
  pde_->init(); // init pde object
  model().init_regularization(); // init regularization
  model().init_sampling(); // init sampling design
}

// a trait to detect if a model requires a preprocessing step
template <typename Model, typename T = void>
struct requires_preprocess : std::false_type {};

template <typename Model> 
struct requires_preprocess<
  Model, std::void_t<decltype(std::declval<Model>().preprocess())>
  > : std::true_type {};

// set model's data from blockframe
template <typename Model>
void ModelBase<Model>::setData(const BlockFrame<double, int>& df) {
  // stop if incoming data has no observations
  if(!df.hasBlock(OBSERVATIONS_BLK))
    throw std::logic_error("model without observations is ill-formed");
  
  df_ = df;
  // insert an index row (if not yet present)
  if(!df_.hasBlock(INDEXES_BLK)){
    std::size_t n = df_.rows();
    DMatrix<int> idx(n,1);
    for(std::size_t i = 0; i < n; ++i) idx(i,0) = i;
    df_.insert(INDEXES_BLK, idx);
  }
  // perform preprocessing of data depending on model type, if model has a defined preprocess() step
  if constexpr(requires_preprocess<Model>::value){
    model().preprocess();
  }
  return;
}

// set boundary conditions on problem's linear system
template <typename Model>
void ModelBase<Model>::setDirichletBC(SpMatrix<double>& A, DMatrix<double>& b){
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

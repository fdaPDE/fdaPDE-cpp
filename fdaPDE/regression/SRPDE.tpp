// finds a solution to the SR-PDE smoothing problem
template <unsigned int M, unsigned int N, unsigned int K, typename E>
void SRPDE<M, N, K, E>::smooth() {
  // assemble system matrix for the nonparameteric part of the model
  SparseBlockMatrix<double,2,2>
    A(-Psi_->transpose()*(*Psi_), lambda_ * pde_.R1()->transpose(),
      lambda_ * (*pde_.R1()),     lambda_ * (*pde_.R0())          );
  // cache system matrix for reuse
  A_ = std::make_shared<SpMatrix<double>>(A.derived());
  DVector<double> solution; // where the system solution will be stored
  
  if(!this->isAlloc(W_)){ // nonparametric case
    // rhs of SR-PDE linear system
    DVector<double> b;
    b.resize(A.rows());
    b << -Psi_->transpose()*(*z_),
      lambda_ * (*pde_.force());
    b_ = std::make_shared<DVector<double>>(b);
    
    // define system solver. Use a sparse solver
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver{};
    solver.compute(*A_);
    // solve linear system A_*x = b
    solution = solver.solve(b);
    
    // store result of smoothing
    f_ = std::make_shared<DVector<double>>( (*Psi_)*solution.head(A_->rows()/2) );
  }else{ // parametric case
    // rhs of SR-PDE linear system
    DVector<double> b;
    b.resize(A_->rows());
    b << -Psi_->transpose()*lmbQ(*this, *z_),
      lambda_ * (*pde_.force());
    b_ = std::make_shared<DVector<double>>(b);
    
    std::size_t q_ = W_->cols(); // number of covariates
    // definition of matrices U and V 
    DMatrix<double> U = DMatrix<double>::Zero(A.rows(), q_);
    U.block(0,0, A.rows()/2, q_) = Psi_->transpose()*(*W_);
  
    DMatrix<double> V = DMatrix<double>::Zero(q_, A.rows());
    V.block(0,0, q_, A.rows()/2) = W_->transpose()*(*Psi_);

    // Define system solver. Use SMW solver from NLA module
    SMW<> solver{};
    solver.compute(*A_);
    // solve system Mx = b
    solution = solver.solve(U, *WTW_, V, b);

    // store result of smoothing 
    f_    = std::make_shared<DVector<double>>( (*Psi_)*solution.head(A_->rows()/2) );
    beta_ = std::make_shared<DVector<double>>( invWTW_.solve(W_->transpose())*(*z_ - *f_) );
  }
  return;
}

// it is asssumed that smooth has already been called on the model object
// in general fitted values \hat z are equal to f_ + W * beta_, in case the parametric part is absent W * beta_ is omitted
template <unsigned int M, unsigned int N, unsigned int K, typename E>
DVector<double> SRPDE<M, N, K, E>::fitted() const {
  DVector<double> hat_z = *f_;
  // if the model has a parametric part, we need to sum its contribute
  if(this->isAlloc(W_))
    hat_z += (*W_)*(*beta_);
  
  return hat_z;
}

// required to support GCV based smoothing parameter selection
// in case of an SRPDE model we have T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
template <unsigned int M, unsigned int N, unsigned int K, typename E>
std::shared_ptr<DMatrix<double>> SRPDE<M, N, K, E>::T() {
  // compute value of R = R1^T*R0^{-1}*R1, cache for possible reuse
  invR0_.analyzePattern(*pde_.R0());
  invR0_.factorize(*pde_.R0());
  R_ = std::make_shared<DMatrix<double>>( pde_.R1()->transpose()*invR0_.solve(*pde_.R1()) );

  // compute and store matrix T for possible reuse
  if(!this->hasCovariates()) // case without covariates, Q is the identity matrix
    T_ = std::make_shared<DMatrix<double>>
      ( Psi_->transpose()*(*Psi_) + lambda_* (*R_) );
  else // general case with covariates
    T_ = std::make_shared<DMatrix<double>>
      ( Psi_->transpose()*lmbQ(*this, *Psi_) + lambda_*(*R_) );

  // return pointer to T
  return T_;
}

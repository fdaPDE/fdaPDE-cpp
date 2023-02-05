// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/core/utils/Symbols.h>
#include <fdaPDE/models/regression/SRPDE.h>
using fdaPDE::models::SRPDE;
#include <fdaPDE/core/utils/DataStructures/BlockFrame.h>
#include <fdaPDE/models/ModelTraits.h>
using fdaPDE::models::Sampling;
#include <fdaPDE/core/FEM/PDE.h>
using fdaPDE::core::FEM::PDE;
using fdaPDE::core::FEM::DefaultOperator;
#include <fdaPDE/core/MESH/Mesh.h>
using fdaPDE::core::MESH::Mesh;

// this file contains the R wrapper for the SRPDE model

template <unsigned int M, unsigned int N, unsigned int R, typename F>
class RegularizingPDE {
private:
  typedef typename std::decay<F>::type BilinearFormType;
  // internal data
  Mesh<M,N,R> domain_;
  PDE<M,N,R, BilinearFormType, DMatrix<double>> pde_;
public:
  // constructor
  RegularizingPDE(const Rcpp::List& R_Mesh) :
    // initialize domain
    domain_(Rcpp::as<DMatrix<double>>(R_Mesh["nodes"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["edges"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["elements"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["neigh"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["boundary"])),
    pde_(domain_) { pde_.setBilinearForm(BilinearFormType()); };
  
  // setters
  void set_dirichlet_bc(const DMatrix<double>& data){ pde_.setDirichletBC(data); }
  void set_forcing_term(const DMatrix<double>& data){ pde_.setForcing(data); }
  // getters
  DMatrix<double> get_quadrature_nodes() const { return pde_.integrator().quadratureNodes(domain_); };
  DMatrix<double> get_dofs_coordinates() const { return domain_.dofCoords(); };
  const PDE<M,N,R, BilinearFormType, DMatrix<double>>& pde() { return pde_; }
  
  // compile time informations
  typedef PDE<M,N,R, BilinearFormType, DMatrix<double>> PDEType;
};
// define 2D simple Laplacian regularization.
typedef RegularizingPDE<2,2,1, decltype( std::declval<Laplacian<DefaultOperator>>() )> Laplacian_2D_Order1;
// expose RegularizingPDE as possible argument to other Rcpp modules
RCPP_EXPOSED_AS  (Laplacian_2D_Order1)
RCPP_EXPOSED_WRAP(Laplacian_2D_Order1)

RCPP_MODULE(Laplacian_2D_Order1) {
  Rcpp::class_<Laplacian_2D_Order1>("Laplacian_2D_Order1")
    .constructor<Rcpp::List>()
    // getters
    .method("get_quadrature_nodes", &Laplacian_2D_Order1::get_quadrature_nodes)
    .method("get_dofs_coordinates", &Laplacian_2D_Order1::get_dofs_coordinates)
    // setters
    .method("set_dirichlet_bc",     &Laplacian_2D_Order1::set_dirichlet_bc)
    .method("set_forcing_term",     &Laplacian_2D_Order1::set_forcing_term);
}

// wrapper for SRPDE module
template <typename RegularizingPDE_, Sampling S_>
class R_SRPDE {
private:
  RegularizingPDE_ regularization_;
  // the model this Rcpp module wraps
  SRPDE<typename RegularizingPDE_::PDEType, S_> model_;
  BlockFrame<double, int> df_;
public:
  // constructor
  R_SRPDE(const RegularizingPDE_& regularization) : regularization_(regularization) {
    model_.setPDE(regularization_.pde());
  };

  // setters
  void set_lambda_s(double lambdaS) { model_.setLambdaS(lambdaS); }
  void set_observations(const DMatrix<double>& data) { df_.template insert<double>(OBSERVATIONS_BLK, data); }
  void set_covariates(const DMatrix<double>& data) { df_.template insert<double>(DESIGN_MATRIX_BLK, data); }
  // getters
  DMatrix<double> f() const { return model_.f(); }
  DMatrix<double> beta() const { return model_.beta(); }
  
  // initialize model and solve smoothing problem
  void solve() {
    model_.setData(df_);
    model_.init(); model_.solve();
  }
};

// SRPDE with simple laplacian regularization, pointwise at mesh nodes
typedef R_SRPDE<Laplacian_2D_Order1, fdaPDE::models::Sampling::GeoStatMeshNodes> SRPDE_Laplacian_2D_GeoStatNodes;

RCPP_MODULE(SRPDE_Laplacian_2D_GeoStatNodes) {
  Rcpp::class_<SRPDE_Laplacian_2D_GeoStatNodes>("SRPDE_Laplacian_2D_GeoStatNodes")
    .constructor<Laplacian_2D_Order1>()
    // getters
    .method("f",                &SRPDE_Laplacian_2D_GeoStatNodes::f)
    .method("beta",             &SRPDE_Laplacian_2D_GeoStatNodes::beta)
    // setters
    .method("set_lambda_s",     &SRPDE_Laplacian_2D_GeoStatNodes::set_lambda_s)
    .method("set_observations", &SRPDE_Laplacian_2D_GeoStatNodes::set_observations)
    .method("set_covariates",   &SRPDE_Laplacian_2D_GeoStatNodes::set_covariates)
    // solve method
    .method("solve",            &SRPDE_Laplacian_2D_GeoStatNodes::solve);
}

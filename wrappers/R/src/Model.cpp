// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <fdaPDE/Core.h>
using fdaPDE::core::ScalarField;
#include <fdaPDE/core/MESH/Mesh.h>
using fdaPDE::core::MESH::Mesh;
using fdaPDE::core::MESH::neighboring_structure;
#include <fdaPDE/core/FEM/Evaluator.h>
using fdaPDE::core::FEM::Evaluator;
using fdaPDE::core::FEM::Raster;

#include <fdaPDE/FEM.h>
#include <fdaPDE/regression/SRPDE.h>

using namespace Rcpp;

// Model is the main interface of the problem

// mesh is stored as a set of matrices
// pde as a set of parameters
// data as a set of matrices/vectors
// lambda as a single scalar
template <unsigned int N, unsigned int M>
struct MeshWrapper{
  // internal mesh data structures
  DMatrix<double> points_{};
  DMatrix<int> edges_{};
  DMatrix<int> elements_{};
  typename neighboring_structure<M, N>::type neighbors_{};
  DMatrix<int> boundary_{};

  // constructor, take raw data directly from R layer
  MeshWrapper() = default;
  MeshWrapper(DMatrix<double> points, DMatrix<int> edges, DMatrix<int> elements,
	      DMatrix<int> neighbors, DMatrix<int> boundary) :
    points_(points), edges_(edges), elements_(elements), neighbors_(neighbors), boundary_(boundary) {}  
};

// define 2D mesh wrapper.
typedef MeshWrapper<2,2> MeshWrapper2D;
// Expose Rcpp:as to be accepted as parameter from other Rcpp modules
RCPP_EXPOSED_AS(MeshWrapper2D)
RCPP_EXPOSED_WRAP(MeshWrapper2D)

RCPP_MODULE(Mesh2D){
  Rcpp::class_<MeshWrapper2D>("Mesh2D")
    .constructor()
    .constructor<DMatrix<double>, DMatrix<int>, DMatrix<int>, DMatrix<int>, DMatrix<int>>()
    .field("points", &MeshWrapper2D::points_)
    .field("edges", &MeshWrapper2D::edges_)
    .field("elements", &MeshWrapper2D::elements_)
    .field("neighbors", &MeshWrapper2D::neighbors_)
    .field("boundary", &MeshWrapper2D::boundary_);
}

class Model{
private:
  DVector<double> z_{}; // vector of observations
  DVector<double> lambda_{};
  MeshWrapper<2,2> domain_;
public:
  Model() {};
  Model(MeshWrapper<2,2> domain, DVector<double> observations, DVector<double> lambda) :
    domain_(domain), z_(observations), lambda_(lambda) {};

  DVector<double> getObs() { return z_; }
  DMatrix<int> nodes() { return domain_.elements_; }
  MeshWrapper<2,2> domain() { return domain_; }
  List plotProva(DVector<double> coeff, double h) {
    // define mesh object
    Mesh<2,2,1> m(domain_.points_, domain_.edges_, domain_.elements_, domain_.neighbors_, domain_.boundary_);
    constexpr double pi = 3.14159265358979323846;
    // true spatial field
    std::function<double(SVector<2>)> f = [](SVector<2> x) -> double {
      double z = 1;
      std::function <double(double x, double y)> coe = [](double x, double y) -> double {
	return 1./2 * std::sin(5*pi*x)*std::exp(-std::pow(x,2)) + 1;
      };
      return std::sin(2 * pi * (coe(x[1],1) * x[0] * std::cos(z-2) - x[1] * std::sin(z-2))) *
	std::cos(2 * pi * (coe(x[1],1) * x[0] * std::cos(z-2+pi/2) + coe(x[0],1) * x[1] * sin((z-2) * pi/2)));
    };
    
    Evaluator<2,2,1> eval;
    Raster<2> img = eval.toRaster(m, f, h);
    // print stuffs required by R
    std::vector<double> x_coord{};
    x_coord.resize(img.size());
    std::vector<double> y_coord{};
    y_coord.resize(img.size());
    // prepare for R layer
    for(std::size_t i = 0; i < img.size(); ++i){
      x_coord[i] = img.coords_[i][0];
      y_coord[i] = img.coords_[i][1];
    }
    return List::create(Named("x") = x_coord, Named("y") = y_coord, Named("field") = img.data_);
  }

  List smooth(DVector<double> observations, Rcpp::Nullable<Rcpp::NumericMatrix> covariates = R_NilValue) {
    // define mesh object
    Mesh<2,2,1> m(domain_.points_, domain_.edges_, domain_.elements_, domain_.neighbors_, domain_.boundary_);

    // definition of PDE Lf = u used in the regularization term
    auto form = Laplacian(); // differential operator L
    DMatrix<double> forcingData; // forcing term u
    forcingData.resize(m.elements(), 1);
    forcingData.fill(4);
    // define differential problem
    PDE problem(m, form, forcingData);
    DMatrix<double> bc; // set dirichelt BC
    bc.resize(m.elements(), 1);
    bc.fill(0);
    problem.setDirichletBC(bc);
    // solve differential problem
    // define functional basis for FEM approximation of PDE solution
    LagrangianBasis<2, 2, 1> basis{};
    // define integrator for numerical approximations of integrals
    Integrator<2, 6> integrator{};
    problem.solve(basis, integrator);

    // define statistical model: simple spatial regression with laplacian regularization
    SRPDE model(problem, lambda_[0]);
    DVector<double> f;
    if(covariates.isNotNull()){
      DMatrix<double> W = Rcpp::as<DMatrix<double>>(covariates);
      f = model.smooth(observations, W);
    }else      
      f = model.smooth(observations); // solve the problem
    
    // solution evaluation
    Evaluator<2,2,1> eval;
    double h = 0.01;
    Raster<2> img = eval.toRaster(m, f, h);
    // print stuffs required by R
    std::vector<double> x_coord{};
    x_coord.resize(img.size());
    std::vector<double> y_coord{};
    y_coord.resize(img.size());
    // prepare for R layer
    for(std::size_t i = 0; i < img.size(); ++i){
      x_coord[i] = img.coords_[i][0];
      y_coord[i] = img.coords_[i][1];
    }
    return List::create(Named("x") = x_coord, Named("y") = y_coord, Named("field") = img.data_, Named("coeff") = f);
  }
};

// exposed interface to R layer
RCPP_MODULE(model) {
  Rcpp::class_<Model>("Model")
    .constructor()
    .constructor<MeshWrapper<2,2>, DVector<double>, DVector<double>>()
    .method("plotProva", &Model::plotProva)
    .method("domain", &Model::domain)
    .method("smooth", &Model::smooth)
    .method("nodes", &Model::nodes);
}

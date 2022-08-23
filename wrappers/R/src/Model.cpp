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
    
    NumericVector result;
    NumericVector x_coord;
    NumericVector y_coord;

    for(double x = 0; x <= 1.0; x+=0.05){
      for(double y = 0; y <= 1.0; y+=0.05){
	double eval = f(SVector<2>(x,y));
	result.push_back(eval);
	x_coord.push_back(x);
	y_coord.push_back(y);
      }
    }

    return List::create(Named("x") = x_coord, Named("y") = y_coord, Named("field") = result);

  }

  List smooth(DVector<double> data, double h) {
    // define mesh object
    Mesh<2,2,1> m(domain_.points_, domain_.edges_, domain_.elements_, domain_.neighbors_, domain_.boundary_);
    // define differential operator
    auto form = Laplacian();
    // forcing term
    DMatrix<double> forcingData;
    forcingData.resize(m.elements(), 1);
    forcingData.fill(4);

    PDE problem(m, form, forcingData);
    // set dirichelt BC
    DMatrix<double> bc;
    bc.resize(m.elements(), 1);
    bc.fill(0);
  
    problem.setDirichletBC(bc);

    LagrangianBasis<2, 2, 1> basis{};  // define functional basis
    Integrator<2, 6> integrator{};  // define integrator
    // solve problem
    problem.solve(basis, integrator);

    SRPDE model(problem, 0.001);
    DVector<double> f = model.smooth(data);
    // evaluate solution
    Evaluator<2,2,1> eval;
    Raster<2> img = eval.toRaster(m, f, h);

    // print stuffs required by R
    std::vector<double> x_coord{};
    x_coord.resize(img.size());
    std::vector<double> y_coord{};
    y_coord.resize(img.size());

    for(std::size_t i = 0; i < img.size(); ++i){
      x_coord[i] = img.coords_[i][0];
      y_coord[i] = img.coords_[i][1];
    }
    return List::create(Named("x") = x_coord, Named("y") = y_coord, Named("field") = img.data_, Named("coeff") = f);
  }

};

// void Model::plotProva(DVector<double> coeff, double h) {
//   // // define mesh object
//   Mesh<2,2,1> m(domain_.points_, domain_.edges_, domain_.elements_, domain_.neighbors_, domain_.boundary_);
//   // // evaluate field
//   // Evaluator<2,2,1> eval;
//   // Raster<2> img = eval.toRaster(m, coeff, h);
//   // Rcout << "ci sono" << std::endl;
  
//   // NumericVector result;
//   // NumericVector x_coord;
//   // NumericVector y_coord;

//   // Rcout << "sono proprio io :)" << std::endl;
  
//   // // for(std::size_t j = 0; j < img.data_.size(); ++j){
//   // //   result.push_back(img.data_[j]);
//   // //   x_coord.push_back(img.coords_[j][0]);
//   // //   y_coord.push_back(img.coords_[j][1]);
//   // // }
//   // return List::create(Named("x") = x_coord, Named("y") = y_coord, Named("field") = result);

//   // define differential operator
//   auto form = Laplacian();
//   // forcing term
//   DMatrix<double> forcingData;
//   forcingData.resize(m.elements(), 1);
//   forcingData.fill(4);

//   PDE problem(m, form, forcingData);
//   // set dirichelt BC
//   DMatrix<double> bc;
//   bc.resize(m.elements(), 1);
//   bc.fill(0);
  
//   problem.setDirichletBC(bc);

//   LagrangianBasis<2, 2, 1> basis{};  // define functional basis
//   Integrator<2, 6> integrator{};  // define integrator
//   // solve problem
//   problem.solve(basis, integrator);
//   // evaluate solution
//   //Evaluator<2,2,1> eval;
//   //Raster<2> img = eval.toRaster(problem, 0.1);

//   // print stuffs required by R
//   NumericVector result{};
//   NumericVector x_coord{};
//   NumericVector y_coord{};

//   // for(std::size_t i = 0; i < img.data_.size(); ++i){
//   //   result.push_back(img.data_[i]);
//   //   x_coord.push_back(img.coords_[i][0]);
//   //   y_coord.push_back(img.coords_[i][1]);
//   // }
  
//   return List::create(Named("x") = x_coord, Named("y") = y_coord, Named("field") = result);  
// }

RCPP_MODULE(model) {
  Rcpp::class_<Model>("Model")
    .constructor()
    .constructor<MeshWrapper<2,2>, DVector<double>, DVector<double>>()
    .method("plotProva", &Model::plotProva)
    .method("domain", &Model::domain)
    .method("smooth", &Model::smooth)
    .method("nodes", &Model::nodes);
}

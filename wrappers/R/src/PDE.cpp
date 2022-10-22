// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <fdaPDE/Core.h>
using fdaPDE::core::ScalarField;

#include <fdaPDE/core/MESH/Mesh.h>
using fdaPDE::core::MESH::Mesh;
using fdaPDE::core::MESH::neighboring_structure;

#include <fdaPDE/FEM.h>
using fdaPDE::core::FEM::PDE;

#include <fdaPDE/core/FEM/Evaluator.h>
using fdaPDE::core::FEM::Evaluator;
using fdaPDE::core::FEM::Raster;

#include "FieldParser.h"
#include "MeshWrapper.h"

class R_PDE{
private:
  Mesh<2,2,1> m_;
  MeshWrapper<2,2> domain_;
  FieldParser<2> forcing_;

  DMatrix<double> solution_;
public:
  R_PDE(MeshWrapper<2,2> domain, std::string forcing)
    : m_(domain.points_, domain.edges_, domain.elements_, domain.neighbors_, domain.boundary_),
      forcing_(forcing) {};

  // solve 2D simple laplacian
  void solve() {
    auto L = Laplacian();
    
    // define differential problem
    ScalarField<2> forcing = forcing_.get();
    
    PDE problem(m_, L, forcing);
    DMatrix<double> bc; // set dirichlet BC
    bc.resize(m_.elements(), 1);
    bc.fill(0);
    problem.setDirichletBC(bc);
    // compute R1_ and R0_ matrices
    problem.solve();

    // store solution
    solution_ = *problem.solution();
    return;
  }

  Rcpp::List plot(double resolution) const {
    // solution evaluation
    Evaluator<2,2,1> eval;
    Raster<2> img = eval.toRaster(m_, solution_, resolution);
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

    return Rcpp::List::create(Rcpp::Named("x") = x_coord, Rcpp::Named("y") = y_coord, Rcpp::Named("solution") = img.data_);
  }
  
};

// exposed interface to R layer
RCPP_MODULE(pde) {
  Rcpp::class_<R_PDE>("PDE.2D")
    .constructor<MeshWrapper<2,2>, std::string>()
    .method("solve", &R_PDE::solve)
    .method("plot", &R_PDE::plot);
}

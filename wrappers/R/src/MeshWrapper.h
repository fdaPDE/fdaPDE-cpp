#ifndef __MESH_WRAPPER_H__
#define __MESH_WRAPPER_H__

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <fdaPDE/Core.h>
#include <fdaPDE/core/MESH/Mesh.h>
using fdaPDE::core::MESH::Mesh;
using fdaPDE::core::MESH::neighboring_structure;

using namespace Rcpp;

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
    .constructor<DMatrix<double>, DMatrix<int>, DMatrix<int>, DMatrix<int>, DMatrix<int>>();
}

#endif

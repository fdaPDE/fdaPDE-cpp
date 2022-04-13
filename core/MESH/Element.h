#ifndef __ELEMENT_H__
#define __ELEMENT_H__

#include <iostream>
#include "../OPT/Utils.h"
#include <memory>
#include "MeshUtils.h"
#include <Eigen/LU>

// a mesh element
template <unsigned int N, unsigned int M>
class Element{

 private:
  unsigned int ID;
  std::array<SVector<N>, N_VERTICES(N,M)> coords;
  std::array<unsigned int, N+1> neighbors;

  // matrix defining the affine transformation from cartesian to barycentric coordinates
  Eigen::Matrix<double, M, M> baryMatrix;
  Eigen::Matrix<double, M, M> invBaryMatrix{};
  
 public:
  Element() = default;
  
  Element(int ID_, std::array<SVector<N>, N_VERTICES(N,M)> coords_, std::array<unsigned int, N+1> neighbors_) :
    ID(ID_), coords(coords_), neighbors(neighbors_) {

    // precompute barycentric coordinate matrix for fast access
    // use first point as reference
    SVector<N> ref = coords[0];
    for(size_t j = 0; j < M; ++j){
      baryMatrix.col(j) = coords[j+1] - ref;
    }
    
    // find barycentric coordinates is equivalent to solve a linear system
    // for efficiecy reasons caching the inverse of baryMatrix can be usefull
    // expetially if there is the need to continuously access to barycentric
    // coordinates. Moreover Eigen inverse() is fast for very small matrices
    // at most 4x4
    // see documentation to learn what system we are trying to solve
    invBaryMatrix = baryMatrix.inverse();
  };

  std::array<SVector<N>, N_VERTICES(N,M)> getCoords() const { return coords; }
  std::array<unsigned int, N+1> getNeighbors() const { return neighbors; }
  unsigned int getID() { return ID; }

  // computes the baricentric coordinates of the element
  SVector<N+1> computeBarycentricCoordinates(const SVector<N>& x) const;

  // check if a given point is inside the element
  bool contains(const SVector<M>& x) const;

  // compute bounding box of element
  std::pair<SVector<N>, SVector<N>> computeBoundingBox() const;
};

template <unsigned int N, unsigned int M>
SVector<N + 1> Element<N, M>::computeBarycentricCoordinates(const SVector<N>& x) const {
  // solve linear system baryMatrix*z = (x - ref) by using the precomputed inverse of baryMatrix
  SVector<N> z = invBaryMatrix*(x - coords[0]);
  // compute barycentric coordinate of reference element
  double z0 = 1 - z.sum();
  
  SVector<N+1> result;
  result << SVector<1>(z0), z;
  
  return result;
}

template <unsigned int N, unsigned int M>
bool Element<N, M>::contains(const SVector<M> &x) const {
  // you can prove that a point is inside the element if all its barycentric coordinates are positive
  
  // get barycentric coordinates of input point
  SVector<M+1> baryCoord = computeBarycentricCoordinates(x);

  // use Eigen visitor to check for positiveness of elements
  return (baryCoord.array() >= 0).all();
}

template <unsigned int N, unsigned int M>
std::pair<SVector<N>, SVector<N>> Element<N,M>::computeBoundingBox() const{

  // define lower-left and upper-right corner of bounding box
  SVector<N> ll, ur;

  // projection of each vertex coordinate on reference axis
  std::array<std::array<double, N_VERTICES(N,M)>, N> projCoords;

  for(size_t j = 0; j < N_VERTICES(N,M); ++j){
    for(size_t dim = 0; dim < N; ++dim){
      projCoords[dim][j] = coords[j][dim];
    }
  }

  // take minimum and maximum value along each dimension, those values define the lower-left and
  // upper-right corner of the bounding box
  for(size_t dim = 0; dim < N; ++dim){
    ll[dim] = *std::min_element(projCoords[dim].begin(), projCoords[dim].end());      
    ur[dim] = *std::max_element(projCoords[dim].begin(), projCoords[dim].end());
  }

  return std::make_pair(ll, ur);
}

#endif // __ELEMENT_H__

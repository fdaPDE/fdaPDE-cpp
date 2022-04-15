#ifndef __ELEMENT_H__
#define __ELEMENT_H__

#include "../OPT/Utils.h"
#include "MeshUtils.h"
#include "Geometry.h"

#include <Eigen/LU>
#include <limits>
#include <memory>
#include <iostream>
#include <type_traits>

template <unsigned int M, unsigned int N> class Element; // forward declaration

// using declarations for manifold specialization
using SurfaceElement = Element<2,3>;
using NetworkElement = Element<1,2>;

// trait to detect if we are dealing with a manifold mesh (N != M)
template <unsigned int N, unsigned int M> struct is_manifold {
public: 
  static constexpr bool value = (N != M);
};

// a mesh element
template <unsigned int M, unsigned int N>
class Element{

 private:
  unsigned int ID;
  std::array<SVector<N>, N_VERTICES(M,N)> coords;
  std::array<unsigned int, N+1> neighbors;

  // matrix defining the affine transformation from cartesian to barycentric coordinates
  Eigen::Matrix<double, N, M> baryMatrix;
  Eigen::Matrix<double, M, N> invBaryMatrix{};

  // compute the inverse of the barycentric matrix. This will be specialization by manifold elements
  Eigen::Matrix<double, M, N> computeInvBaryMatrix(const Eigen::Matrix<double, N, M>& baryMatrix) {
    return baryMatrix.inverse();
  };
  
 public:
  Element() = default;
  
  Element(int ID_, std::array<SVector<N>, N_VERTICES(M,N)> coords_, std::array<unsigned int, N+1> neighbors_) :
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
    // (at most 4x4)
    // see documentation to learn what system we are trying to solve
    invBaryMatrix = computeInvBaryMatrix(baryMatrix);
  };

  std::array<SVector<N>, N_VERTICES(M,N)> getCoords() const { return coords; }
  std::array<unsigned int, N+1> getNeighbors() const { return neighbors; }
  unsigned int getID() { return ID; }

  // computes the baricentric coordinates of the element
  SVector<M+1> computeBarycentricCoordinates(const SVector<N>& x) const;

  // check if a given point is inside the element
  bool contains(const SVector<N>& x) const;

  // compute bounding box of element
  std::pair<SVector<N>, SVector<N>> computeBoundingBox() const;
};

template <unsigned int M, unsigned int N>
SVector<M+1> Element<M, N>::computeBarycentricCoordinates(const SVector<N>& x) const {
  // solve linear system baryMatrix*z = (x - ref) by using the precomputed inverse of baryMatrix
  SVector<M> z = invBaryMatrix*(x - coords[0]);
  // compute barycentric coordinate of reference element
  double z0 = 1 - z.sum();
  
  SVector<M+1> result;
  result << SVector<1>(z0), z;

  return result;
}

template <unsigned int M, unsigned int N>
std::pair<SVector<N>, SVector<N>> Element<M,N>::computeBoundingBox() const{

  // define lower-left and upper-right corner of bounding box
  SVector<N> ll, ur;

  // projection of each vertex coordinate on reference axis
  std::array<std::array<double, N_VERTICES(M,N)>, N> projCoords;

  for(size_t j = 0; j < N_VERTICES(M,N); ++j){
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

template <unsigned int M, unsigned int N>
bool Element<M, N>::contains(const SVector<N> &x) const {
  // you can prove that a point is inside the element if all its barycentric coordinates are positive
  
  // get barycentric coordinates of input point
  SVector<N+1> baryCoord = computeBarycentricCoordinates(x);

  // use Eigen visitor to check for positiveness of elements
  return (baryCoord.array() >= 0).all();
}

// specialization for 2.5D domains (surfaces)
template <> 
Eigen::Matrix<double, 2, 3> SurfaceElement::computeInvBaryMatrix(const Eigen::Matrix<double, 3, 2>& baryMatrix) {
  // returns the generalized inverse of baryMatrix (which is a rectangular for surface elements)
  return (baryMatrix.transpose()*baryMatrix).inverse()*baryMatrix.transpose();
}
  
template <>
bool SurfaceElement::contains(const SVector<3>& x) const {
  // we start checking if the point is contained in the plane spanned by the mesh element

  // basis for the plane passing by the element, observe that the spase spanned by this set of basis is a
  // vector space in the proper sense (the plane passes throught zero). To cope with affine spaces getL2Distance
  // of mdule Geometry accepts an offset parameter representing the point throught which the space
  // spanned by this basis set has to pass
  SVector<3> a = (coords[1] - coords[0]);
  SVector<3> b = (coords[2] - coords[0]);
  
  // if the distance between the point projection into the plane and the point itself is larger than 0
  // return false, the point does not belong to the plane and therefore cannot belong to the surface element
  if(Geometry<3>::getL2Distance({a,b}, coords[0], x) > std::numeric_limits<double>::epsilon()){
    return false;
  }

  // if the point belongs to the spanned plane, check if its barycentric coordinates are all positive
  SVector<3> baryCoord = computeBarycentricCoordinates(x);

  // use Eigen visitor to check for positiveness of elements
  return (baryCoord.array() >= 0).all();
}

#endif // __ELEMENT_H__

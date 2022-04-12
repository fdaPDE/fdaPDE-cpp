#ifndef __MESH_H__
#define __MESH_H__

#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>
#include <memory>
#include <array>

#include "Element.h"
#include "CSVReader.h"
#include "../OPT/Utils.h"
#include "MeshUtils.h"

// at the moment used for barycentric walk

/* at the moment mesh data from R require the development version of RTriangle, download from
   https://github.com/davidcsterratt/RTriangle

   install from R

   install.packages("devtools")
   devtools::install_github("davidcsterratt/RTriangle", subdir="pkg")

   reason: CRAN version of RTriangle does not outputs neighboring information. Computing the neighboring
   matrix here could be expensive expetially for very large meshes.

   instruction to produce the triangulation
   library(fdaPDE)
   data(quasicircle2d)

   library(RTriangle)
   X = rbind(quasicircle2d$boundary_nodes, quasicircle2d$locations)

   X_pslg = pslg(P = X)

   pQ = RTriangle::triangulate(p = X_pslg)
 */

template <typename T> using DynamicMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

// this class offers a point of access to mesh information from an high level
// reasoning, i.e. without considering the internal representation of a mesh in
// memory.
// N is the order of the mesh, M is the dimension of the space where the mesh is embedded
template <unsigned int N, unsigned int M>
class Mesh{

 private:
  // raw mesh informations

  // coordinates of points constituting the vertices of mesh elements
  DynamicMatrix<double> points_;
  // matrix of edges. Each row of the matrix contains the row numbers in points_ matrix
  // of the points which form the edge
  DynamicMatrix<int> edges_;
  // matrix of triangles in the triangulation. Each row of the matrix contains the row
  // numbers in points_ matrix of the points which form the triangle
  DynamicMatrix<int> triangles_;
  unsigned int numElements = 0;
  // row i in this matrix contains the indexes as row number in triangles_ matrix of the
  // neighboring triangles to triangle i (all triangles in the triangulation which share an
  // edge with i)
  DynamicMatrix<int> neighbors_;

  // store min-max values for each dimension of the mesh
  std::array<std::pair<double, double>, N> meshRange;
  // is often required to access just to the minimum value along each dimension and to the quantity
  // 1/(max[dim] - min[dim]) = 1/(meshRange[dim].second - meshRange[dim].first). Compute here
  // once and cache results for efficiency
  std::array<double, N> minMeshRange;
  std::array<double, N> kMeshRange; // kMeshRange[dim] = 1/(meshRange[dim].second - meshRange[dim].first)
  
 public:

  // constructor from .csv files
  Mesh(const std::string& pointsFile,    const std::string& edgesFile,
       const std::string& trianglesFile, const std::string& neighborsFile){
    // open and parse CSV files
    CSVReader reader;
    CSVFile<double> points = reader.parseFile<double>(pointsFile);
    CSVFile<int> edges     = reader.parseFile<int>(edgesFile);
    CSVFile<int> triangles = reader.parseFile<int>(trianglesFile);
    CSVFile<int> neighbors = reader.parseFile<int>(neighborsFile);
    
    // set mesh internal representation to eigen matrix
    points_    = points.toEigen();

    // compute mesh limits
    for(size_t dim = 0; dim < N; ++dim){
      meshRange[dim].first  = points_.col(dim).minCoeff();
      meshRange[dim].second = points_.col(dim).maxCoeff();

      minMeshRange[dim] = meshRange[dim].first;
      kMeshRange[dim] = 1/(meshRange[dim].second - meshRange[dim].first);
    }
    
    // need to subtract 1 from all indexes since triangle indexes start from 1
    // C++ start counting from 0 instead
    edges_          = edges.toEigen();
    edges_          = (edges_.array() - 1).matrix();

    triangles_      = triangles.toEigen();
    triangles_      = (triangles_.array() -1).matrix();
    numElements     = triangles_.rows();
    
    // a negative value means no neighbor
    neighbors_      = neighbors.toEigen();
    neighbors_      = (neighbors_.array() - 1).matrix();
  }

  // get an element object given its ID (its row number in the triangles_ matrix)
  // this creates an element abstraction only when it is actually needed to avoid waste
  // of resources. Mesh class does not represent explicitly elements of the mesh,
  // instead directly manages raw representation built on matrices
  std::shared_ptr<Element<N,M>> requestElementById(unsigned int ID) const;

  // allow range-for loop over mesh elements
  // we can write something like
  // for(Element& e : mesh) { ... e is a Mesh element ... }
  struct iterator{
  private:
    friend Mesh;
    Mesh* meshContainer; // pointer to mesh object
    int index;           // keep track of current iteration during for-loop
    // constructor
    iterator(Mesh* container_, int index_) : meshContainer(container_), index(index_) {}; 
  public:
    // just increment the current iteration and return this iterator
    iterator& operator++() {
      ++index;
      return *this;
    }
    // dereference the iterator means to create Element object at current index
    std::shared_ptr<Element<N,M>> operator*() const {
      return meshContainer->requestElementById(index);
    }
    // two iterators are different when their indexes are different
    friend bool operator!=(const iterator& lhs, const iterator& rhs) {
      return lhs.index != rhs.index;
    }
  };

  // provide begin() and end() methods
  iterator begin() { return iterator(this, 0); }
  iterator end()   { return iterator(this, triangles_.rows()); }

  // getters
  unsigned int getNumberOfElements() const { return numElements; }
  std::array<std::pair<double, double>, N> getMeshRange() const { return meshRange; }
  std::array<double, N> getMinMeshRange() const { return minMeshRange; }
  std::array<double, N> getKMeshRange() const { return kMeshRange; }
  
};

// implementation for 2D case only...
template <unsigned int N, unsigned int M>
std::shared_ptr<Element<N,M>> Mesh<N,M>::requestElementById(unsigned int ID) const {
  // in the following use auto to take advantage of eigen acceleration
  
  // get the indexes of vertices from triangles_
  auto pointIndexes = triangles_.row(ID);
  // get neighbors information
  auto elementNeighbors = neighbors_.row(ID);
  
  // get vertices coordinates
  std::array<SVector<N>, N_VERTICES(N,M)> coords;
  std::array<unsigned int, N+1> neighbors;

  for(size_t i = 0; i < pointIndexes.size(); ++i){
    SVector<N> vertex(points_.row(pointIndexes[i]));
    coords[i] = vertex;

    // from triangle documentation: The first neighbor of triangle i is opposite the first corner of triangle i, and so on.
    // by storing neighboring informations as they come from triangle we have that neighbor[0] is the
    // triangle adjacent to the face opposite to coords[0]

    neighbors[i] = elementNeighbors[i];
  }
         
  return std::make_shared<Element<N,M>>(ID, coords, neighbors);
}

#endif // __MESH_H__

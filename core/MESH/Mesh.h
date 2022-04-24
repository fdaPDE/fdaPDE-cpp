#ifndef __MESH_H__
#define __MESH_H__

#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <unordered_map>
#include <vector>
#include <memory>
#include <array>

#include "Element.h"
#include "CSVReader.h"
#include "../OPT/Utils.h"
#include "MeshUtils.h"

namespace fdaPDE{
namespace core{
namespace MESH{

  template <typename T> using DynamicMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

  // this class offers a point of access to mesh information from an high level
  // reasoning, i.e. without considering the internal representation of a mesh in memory.
  // M is the order of the mesh (M = 1 linear network elements, M = 2 surface elements, M = 3 volumetric elemnts),
  // N is the dimension of the space where the mesh is embedded. If N != M the mesh is a manifold
  template <unsigned int M, unsigned int N>
  class Mesh{
  private:
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
	 const std::string& trianglesFile, const std::string& neighborsFile);

    // get an element object given its ID (its row number in the triangles_ matrix)
    // this creates an element abstraction only when it is actually needed to avoid waste
    // of resources. Mesh class does not represent explicitly elements of the mesh,
    // instead directly manages raw representation built on matrices
    const std::shared_ptr<Element<M,N>> requestElementById(unsigned int ID) const;

    // allow range-for loop over mesh elements
    struct iterator{
    private:
      friend Mesh;
      const Mesh* meshContainer; // pointer to mesh object
      int index;           // keep track of current iteration during for-loop
      // constructor
      iterator(const Mesh* container_, int index_) : meshContainer(container_), index(index_) {}; 
    public:
      // just increment the current iteration and return this iterator
      iterator& operator++() {
	++index;
	return *this;
      }
      // dereference the iterator means to create Element object at current index
      std::shared_ptr<Element<M,N>> operator*() {
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
    unsigned int getNumberOfElements()                      const { return numElements;  }
    std::array<std::pair<double, double>, N> getMeshRange() const { return meshRange;    }
    std::array<double, N> getMinMeshRange()                 const { return minMeshRange; }
    std::array<double, N> getKMeshRange()                   const { return kMeshRange;   }  
  };

  // export some aliases
  using Mesh2D             = Mesh<2,2>;
  using Mesh3D             = Mesh<3,3>;
  using SurfaceMesh        = Mesh<2,3>; // manifold cases
  using LinearNetworkMesh  = Mesh<1,2>;

#include "Mesh.tpp"
}}}
  
#endif // __MESH_H__

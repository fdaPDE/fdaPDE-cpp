#ifndef __MESH_H__
#define __MESH_H__

#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <unordered_map>
#include <vector>
#include <memory>

#include "Element.h"
#include "CSVReader.h"

/* at the moment mesh data from R require the development version of RTriangle, download from
   https://github.com/davidcsterratt/RTriangle

   install from R

   install.packages("devtools")
   devtools::install_github("davidcsterratt/RTriangle", subdir="pkg")

   reason: CRAN version of RTriangle does not outputs neighboring information. Computing the neighboring
   matrix here could be expensive expetially for very large meshes.
 */

template <typename T> using DynamicMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

// this class offers a point of access to mesh information from an high level reasoning, i.e.
// without considering the internal representation of a mesh in memory
class Mesh{

 private:
  // raw mesh information

  // coordinates of points constituting the vertices of mesh elements
  DynamicMatrix<double> points_;
  // matrix of edges. Each row of the matrix contains the row numbers in points_ matrix
  // of the points which form the edge
  DynamicMatrix<int> edges_;
  // matrix of triangles in the triangulation. Each row of the matrix contains the row
  // numbers in points_ matrix of the points which form the triangle
  DynamicMatrix<int> triangles_;
  // row i in this matrix contains the indexes as row number in triangles_ matrix of the
  // neighboring triangles to triangle i (all triangles in the triangulation which share an
  // edge with i)
  DynamicMatrix<int> neighbors_;
  
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
    edges_     = edges.toEigen();
    triangles_ = triangles.toEigen();
    neighbors_ = neighbors.toEigen();
  }

  // get an element object given its ID (its row number in the triangles_ matrix)
  // this creates an element abstraction only when it is actually needed to avoid waste
  // of resources. Mesh class does not represent explicitly elements of the mesh,
  // instead directly manages raw representation built on matrices
  std::shared_ptr<Element> getElementById(unsigned int ID) const;

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
    std::shared_ptr<Element> operator*() const {
      return meshContainer->getElementById(index);
    }
    // two iterators are different when their indexes are different
    friend bool operator!=(const iterator& lhs, const iterator& rhs) {
      return lhs.index != rhs.index;
    }
  };

  // provide begin() and end() methods
  iterator begin() { return iterator(this, 0); }
  iterator end()   { return iterator(this, triangles_.rows()); }
  
  // search for an element containg a point (big part here: naive, barycentric walking, ADT, structured meshes)
};

std::shared_ptr<Element> Mesh::getElementById(unsigned int ID) const {
  return std::make_shared<Element>(ID);
}


#endif // __MESH_H__

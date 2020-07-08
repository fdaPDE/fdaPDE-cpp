#ifndef __MESH_INPUT_HELPER_HPP__
#define __MESH_INPUT_HELPER_HPP__

#include <array>
#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
#include <utility>
#include <type_traits>


template<UInt mydim>
class simplex{
public:

  using nodeIndices = std::array<UInt, mydim>;
  using const_iterator = typename nodeIndices::const_iterator;
  using const_reverse_iterator = typename nodeIndices::const_reverse_iterator;

  simplex()=delete;
  simplex(UInt elementID_, UInt subelementID_, std::array<UInt, mydim> nodes_) :
    elementID(elementID_), subelementID(subelementID_), nodes(nodes_) {}

  UInt i() const {return elementID;}
  UInt j() const {return subelementID;}
  const UInt& operator[](UInt i) const {return nodes[i];}

  friend bool operator==(const simplex& lhs, const simplex& rhs) {return std::equal(lhs.rbegin(),lhs.rend(), rhs.rbegin());}
  friend bool operator!=(const simplex& lhs, const simplex& rhs) {return !(lhs==rhs);}

  const_iterator begin() const {return nodes.begin();}
  const_iterator end() const {return nodes.end();}
  const_reverse_iterator rbegin() const {return nodes.rbegin();}
  const_reverse_iterator rend() const {return nodes.rend();}

private:
  UInt elementID;
  UInt subelementID;
  nodeIndices nodes;
};

template<UInt mydim>
class simplex_container{
  static_assert(mydim==2 || mydim==3,
    "ERROR! TRYING TO INSTANTIATE SIMPLEX_CONTAINER IN DIMENSION OTHER THAN 2 OR 3! See mesh_input_helper.h");

public:

  using OutputType=std::tuple<std::vector<UInt>, std::vector<bool>, std::vector<bool>, std::vector<int> >;
  using simplex_t = simplex<mydim>;
  using simplex_container_t = std::vector<simplex_t>;
  using const_iterator = typename simplex_container_t::const_iterator;

  simplex_container()=delete;

  template<std::size_t SIZE>
  simplex_container(const UInt* const elements_, UInt  num_elements_, UInt num_points_, const std::array<UInt, SIZE>& ORDERING) :
      elements(elements_), num_elements(num_elements_), num_points(num_points_) {this->fill_container(elements_, ORDERING);}

  OutputType assemble_output() const;
  std::vector<UInt> get_simplexes() const {return this->assemble_subs();};

  const simplex_t& operator[](UInt i) const {return simplexes[i];}
  const_iterator begin() const {return simplexes.begin();}
  const_iterator end() const {return simplexes.end();}

  bool is_repeated(UInt i) const {return duplicates[i];}

  UInt size() const {return simplexes.size();}
  UInt get_num_points() const {return num_points;}
  UInt get_num_elements() const {return num_elements;}

private:
  simplex_container_t simplexes;
  std::vector<bool> duplicates;
  std::vector<UInt> distinct_indexes;
  const UInt num_elements;
  const UInt num_points;
  const UInt* const elements;

  template<std::size_t SIZE>
  void fill_container(const UInt* const, const std::array<UInt, SIZE>&);
  std::vector<UInt> compute_offsets(const UInt, std::vector<UInt>&);
  void bin_sort_(const UInt, std::vector<UInt>&);
  void bin_sort();
  void check_duplicates();
  void store_indexes();
  std::vector<bool> mark_boundary() const;
  std::vector<UInt> assemble_subs() const;
  std::vector<int> compute_neighbors() const;


};

std::vector<UInt> order2extend(const simplex_container<2> &edge_container){
  std::vector<UInt> edges_extended(edge_container.size());
  UInt offset{edge_container.get_num_points()};
  {
    UInt pos=0;
    for(auto const &curr : edge_container){
      offset += !edge_container.is_repeated(pos++);
      edges_extended[curr.i()+edge_container.get_num_elements()*curr.j()]=offset;
    }
  }
  return edges_extended;
}

std::vector<double> compute_midpoints(const double* const points, const std::vector<UInt>& edges, const UInt num_points){
  const UInt num_edges=edges.size()/2;
  std::vector<double> midpoints(3*num_edges);
  for (int i=0; i<num_edges; ++i)
    for (int j=0; j<3; ++j)
      midpoints[i+j*num_edges]=(points[edges[i]+j*num_points]+points[edges[i+num_edges]+j*num_points])/2;
  return midpoints;
}

std::vector<double> compute_midpoints2D(const double* const points, const std::vector<UInt>& edges, const UInt num_points){
  const UInt num_edges=edges.size()/2;
  std::vector<double> midpoints(2*num_edges);
  for (int i=0; i<num_edges; ++i)
    for (int j=0; j<2; ++j)
      midpoints[i+j*num_edges]=(points[edges[i]+j*num_points]+points[edges[i+num_edges]+j*num_points])/2;
  return midpoints;
}


std::vector<UInt> split(const std::vector<UInt>& extended_triangles, const int* const triangles, const UInt num_triangles){

  std::vector<UInt> splitted_triangles;
  splitted_triangles.reserve(12*num_triangles);

  for(int i=0; i<3*num_triangles; ++i)
    splitted_triangles.push_back(triangles[i]+1);

  for (auto const j : {0,2,0,1,1,1,2,0,2})
    for (int i=0; i<num_triangles; ++i)
      splitted_triangles.push_back(extended_triangles[i+j*num_triangles]);

  return splitted_triangles;
}

std::vector<UInt> split3D(const std::vector<UInt>& extended_tetrahedrons, const int* const tetrahedrons, const UInt num_tetrahedrons){

  std::vector<UInt> splitted_tetrahedrons;
  splitted_tetrahedrons.reserve(32*num_tetrahedrons);

  for(int i=0; i<num_tetrahedrons; ++i)
    splitted_tetrahedrons.push_back(tetrahedrons[i]+1);

  for (auto const j : {0,1,2,0,0,1,1,0})
    for (int i=0; i<num_tetrahedrons; ++i)
      splitted_tetrahedrons.push_back(extended_tetrahedrons[i+j*num_tetrahedrons]);

  for(int i=0; i<num_tetrahedrons; ++i)
    splitted_tetrahedrons.push_back(tetrahedrons[i+num_tetrahedrons]+1);

  for (auto const j : {3,5,1,1,2,3,1,3})
    for (int i=0; i<num_tetrahedrons; ++i)
      splitted_tetrahedrons.push_back(extended_tetrahedrons[i+j*num_tetrahedrons]);

  for(int i=0; i<num_tetrahedrons; ++i)
    splitted_tetrahedrons.push_back(tetrahedrons[i+2*num_tetrahedrons]+1);

  for (auto const j : {4,2,3,5,5,2,5,4})
    for (int i=0; i<num_tetrahedrons; ++i)
      splitted_tetrahedrons.push_back(extended_tetrahedrons[i+j*num_tetrahedrons]);

  for(int i=0; i<num_tetrahedrons; ++i)
    splitted_tetrahedrons.push_back(tetrahedrons[i+3*num_tetrahedrons]+1);

  for (auto const j : {5,5,4,4})
    for (int i=0; i<num_tetrahedrons; ++i)
      splitted_tetrahedrons.push_back(extended_tetrahedrons[i+j*num_tetrahedrons]);

  return splitted_tetrahedrons;
}



#include "mesh_input_helper_imp.h"

#endif

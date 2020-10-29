#ifndef __MESH_INPUT_HELPER_H__
#define __MESH_INPUT_HELPER_H__

#include <array>
#include <vector>
#include <iterator>
#include <algorithm>

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

  using simplex_t = simplex<mydim>;
  using simplex_container_t = std::vector<simplex_t>;
  using const_iterator = typename simplex_container_t::const_iterator;

  simplex_container()=delete;

  template<std::size_t SIZE>
  simplex_container(RIntegerMatrix elements_, RNumericMatrix nodes_, const std::array<UInt, SIZE>& ORDERING) :
      elements(elements_), nodes(nodes_), isTriangleContainer(SIZE==6) {this->fill_container(ORDERING);}

  template<std::size_t SIZE>
  simplex_container(SEXP Relements, SEXP Rnodes, const std::array<UInt, SIZE>& ORDERING) :
      elements(Relements), nodes(Rnodes), isTriangleContainer(SIZE==6) {this->fill_container(ORDERING);}

  const simplex_t& operator[](UInt i) const {return simplexes[i];}
  const UInt& distinct(UInt i, UInt j) const {return simplexes[distinct_indexes[i]][j];}
  const_iterator begin() const {return simplexes.begin();}
  const_iterator end() const {return simplexes.end();}

  bool is_repeated(UInt i) const {return duplicates[i];}

  UInt size() const {return simplexes.size();}
  UInt num_distinct() const {return distinct_indexes.size();}
  UInt get_num_points() const {return nodes.nrows();}
  UInt get_num_elements() const {return elements.nrows();}

  void mark_boundary(SEXP Routput, UInt index) const;
  void assemble_subs(SEXP Routput, UInt index) const;
  void compute_neighbors(SEXP Routput, UInt index) const;
  void order2extend(SEXP Routput, UInt index) const;


private:
  simplex_container_t simplexes;
  std::vector<bool> duplicates;
  std::vector<UInt> distinct_indexes;

  const RNumericMatrix nodes;
  const RIntegerMatrix elements;

  const bool isTriangleContainer;

  template<std::size_t SIZE>
  void fill_container(const std::array<UInt, SIZE>&);
  
  std::vector<UInt> compute_offsets(const UInt, std::vector<UInt>&);
  void bin_sort_(const UInt, std::vector<UInt>&);
  void bin_sort();
  
  void check_duplicates();
  void store_indexes();

};


void mark_boundary_nodes(SEXP Routput, SEXP Rnodes, UInt index, UInt index_subs, UInt index_markers) {
  
  const RNumericMatrix nodes(Rnodes);
  SET_VECTOR_ELT(Routput, index, Rf_allocMatrix(LGLSXP, nodes.nrows(), 1));
  const RIntegerMatrix subs(VECTOR_ELT(Routput, index_subs));
  const RIntegerMatrix submarkers(VECTOR_ELT(Routput, index_markers));
  RIntegerMatrix nodesmarkers(VECTOR_ELT(Routput, index));

  for (UInt i=0; i<nodes.nrows(); ++i)
    nodesmarkers[i]=0;

  for(UInt j=0; j<subs.ncols(); ++j)
    for(UInt i=0; i<subs.nrows(); ++i)
      if(nodesmarkers[subs(i,j)-1]==0)
        nodesmarkers[subs(i,j)-1] = submarkers[i];
    
}


void compute_midpoints(SEXP Routput, SEXP Rnodes, UInt index, UInt index_edges){
  
  const RNumericMatrix nodes(Rnodes);
  const RIntegerMatrix edges(VECTOR_ELT(Routput, index_edges));

  SET_VECTOR_ELT(Routput, index, Rf_allocMatrix(REALSXP, edges.nrows(), nodes.ncols()));
  RNumericMatrix midpoints(VECTOR_ELT(Routput, index));

  for (int i=0; i<midpoints.nrows(); ++i)
    for (int j=0; j<midpoints.ncols(); ++j)
      midpoints(i,j) = .5*(nodes(edges(i,0)-1, j)+nodes(edges(i,1)-1, j));
}

void compute_midpoints(SEXP Routput, SEXP Rnodes, UInt index, const simplex_container<2> &edge_container){
  
  const RNumericMatrix nodes(Rnodes);

  SET_VECTOR_ELT(Routput, index, Rf_allocMatrix(REALSXP, edge_container.num_distinct(), nodes.ncols()));
  RNumericMatrix midpoints(VECTOR_ELT(Routput, index));

  for (int i=0; i<midpoints.nrows(); ++i)
    for (int j=0; j<midpoints.ncols(); ++j)
      midpoints(i,j) = .5*(nodes(edge_container.distinct(i,0), j)+nodes(edge_container.distinct(i,1), j));
}


void split(SEXP Routput, SEXP Rtriangles, UInt index, const simplex_container<2> &edge_container){

  std::vector<UInt> extended_triangles(edge_container.size());
  {
    UInt offset{edge_container.get_num_points()};
    UInt pos=0;
    for(auto const &curr : edge_container){
      offset += !edge_container.is_repeated(pos++);
      extended_triangles[curr.i()+edge_container.get_num_elements()*curr.j()]=offset;
    }
  }

  const RIntegerMatrix triangles(Rtriangles);

  SET_VECTOR_ELT(Routput, index, Rf_allocMatrix(INTSXP, 4*triangles.nrows(), 3));
  RIntegerMatrix splitted_triangles(VECTOR_ELT(Routput, index));

  int i=0;
  for( ; i<3*triangles.nrows(); ++i)
    splitted_triangles[i] = triangles[i]+1;

  for (auto const j : {0,2,0,1,1,1,2,0,2})
    for (int k=0; k<triangles.nrows(); ++i, ++k)
      splitted_triangles[i] = extended_triangles[k+j*triangles.nrows()];

}

void split3D(SEXP Routput, SEXP Rtetrahedrons, UInt index, const simplex_container<2> &edge_container){


  std::vector<UInt> extended_tetrahedrons(edge_container.size());
  {
    UInt offset{edge_container.get_num_points()};
    UInt pos=0;
    for(auto const &curr : edge_container){
      offset += !edge_container.is_repeated(pos++);
      extended_tetrahedrons[curr.i()+edge_container.get_num_elements()*curr.j()]=offset;
    }
  }

  const RIntegerMatrix tetrahedrons(Rtetrahedrons);

  SET_VECTOR_ELT(Routput, index, Rf_allocMatrix(INTSXP, 8*tetrahedrons.nrows(), 4));
  RIntegerMatrix splitted_tetrahedrons(VECTOR_ELT(Routput, index));
  
  int i=0;
  for( ; i<tetrahedrons.nrows(); ++i)
    splitted_tetrahedrons[i] = tetrahedrons[i]+1;

  for (auto const j : {0,1,2,0,0,1,1,0})
    for (int k=0; k<tetrahedrons.nrows(); ++k, ++i)
      splitted_tetrahedrons[i]=extended_tetrahedrons[k+j*tetrahedrons.nrows()];

  for(int k=0; k<tetrahedrons.nrows(); ++k, ++i)
    splitted_tetrahedrons[i]=tetrahedrons[k+tetrahedrons.nrows()]+1;

  for (auto const j : {3,5,1,1,2,3,1,3})
    for (int k=0; k<tetrahedrons.nrows(); ++k, ++i)
      splitted_tetrahedrons[i]=extended_tetrahedrons[k+j*tetrahedrons.nrows()];

  for(int k=0; k<tetrahedrons.nrows(); ++k, ++i)
    splitted_tetrahedrons[i]=tetrahedrons[k+2*tetrahedrons.nrows()]+1;

  for (auto const j : {4,2,3,5,5,2,5,4})
    for (int k=0; k<tetrahedrons.nrows(); ++k, ++i)
      splitted_tetrahedrons[i]=extended_tetrahedrons[k+j*tetrahedrons.nrows()];

  for(int k=0; k<tetrahedrons.nrows(); ++k, ++i)
    splitted_tetrahedrons[i]=tetrahedrons[k+3*tetrahedrons.nrows()]+1;

  for (auto const j : {5,5,4,4})
    for (int k=0; k<tetrahedrons.nrows(); ++k, ++i)
      splitted_tetrahedrons[i]=extended_tetrahedrons[k+j*tetrahedrons.nrows()];

}

#include "Mesh_Input_Helper_imp.h"

#endif

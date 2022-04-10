#ifndef __TREE_SEARCH_H__
#define __TREE_SEARCH_H__

#include <cstddef>
#include <list>
#include <stack>
#include <utility>

#include "Tree.h"
#include "../OPT/Utils.h"

template <unsigned int N> using rectangle = std::pair<SVector<N>, SVector<N>>;

// a specific node data structure for easy the handling of an ADT structure
template <unsigned int N>
class MeshNode {
private:
  SVector<N>   point_;       // the point stored in this node 
  rectangle<N> node_range_;  // the range in the hypercube this node refers to
public:
  // constructor
  MeshNode(const SVector<N>& point, const rectangle<N>& node_range) : point_(point), node_range_(node_range) {}
  // getters
  SVector<N> getPoint()   const { return point_; }
  rectangle<N> getRange() const { return node_range_; }
};

// Alternating Decision Tree implementation for tree-based search of elements over an unstructured mesh
template <unsigned int N>
class TreeSearch{
private:
  Tree<MeshNode<N>> tree;

  // returns true if the two hyper-rectangles intersect
  bool intersect(const rectangle<N>& rect_1, const rectangle<N>& rect_2) const;
  bool intersect(const SVector<N>& point, const rectangle<N>& rectangle) const;  
public:
  // initialize an empty decision tree
  TreeSearch() = default;
  // initialize the ADT given a set of N-dimensional points
  TreeSearch(const std::vector<SVector<N>>& data) { init(data); }

  // build the Alternating Decision Tree for unstructured mesh given a set of N-dimensional points
  // recall that given a mesh element of dimension N, we can project it into a point of dimension 2N
  void init(const std::vector<SVector<N>>& data);

  std::list<SVector<N>> search(const rectangle<N>& query);

  Tree<MeshNode<N>> getTree() const { return tree; }
};

// data is already scaled in the range [0,1]
template <unsigned int N> void TreeSearch<N>::init(const std::vector<SVector<N>>& data) {

  // initialize here once
  SVector<N> left_lower_corner = SVector<N>::Zero(), right_upper_corner = SVector<N>::Ones();
  
  // initialize tree data structure
  tree = Tree(MeshNode<N>(data[0], std::make_pair(left_lower_corner, right_upper_corner)));

  // process all points inside data one by one and insert them in the correct position
  for(size_t j = 1; j < data.size(); ++j){
    SVector<N> nodeData = data[j];
    
    // at first we consider the whole hypercube
    rectangle<N> nodeRange = std::make_pair(left_lower_corner, right_upper_corner);

    // traverse the tree based on coordinates of data and insert the corresponding
    // node at right position in the tree
    node_ptr<MeshNode<N>> current = tree.getNode(0); // start from root

    // init visit
    bool inserted = false;         // stop iterating when an insertion point has been found
    unsigned int iteration = 1;    // defines the granularity of the split dimension
    std::array<double,N> offset{}; // used to keep track of the virtual splits of the domain at each iteration
  
    // search for the right insertion location in the tree
    while(!inserted){
      for(size_t dim = 0; dim < N; ++dim){       // cycle over dimensions
	double treshold = offset[dim] + std::pow(0.5, iteration);
	
	if(nodeData[dim] < treshold){
	  nodeRange.second[dim] = treshold;      // shrink node range on the left
	  if(tree.insert(MeshNode<N>(nodeData, nodeRange), current->getKey(), LinkDirection::LEFT )){ // O(1) operation
	    inserted = true;                     // stop searching for location
	    break;
	  }
	  else
	    current = current->getChildren()[0]; // move to left child
	}
	else{
	  nodeRange.first[dim] = treshold;       // shrink node range on the right
	  if(tree.insert(MeshNode<N>(nodeData, nodeRange), current->getKey(), LinkDirection::RIGHT)){ // O(1) operation
	    inserted = true;                     // stop searching for location
	    break;
	  }
	  else{
	    current = current->getChildren()[1]; // move to right child
	    offset[dim] += std::pow(0.5, iteration);
	  }
	}
      }
      // virtually perform an half split of the hyper-cube
      iteration++;
    }
  }
  // construction ended
  return;
}

template <unsigned int N>
bool TreeSearch<N>::intersect(const rectangle<N> &r1, const rectangle<N> &r2) const {
  // if we find a dimension along which the two rectangles have a common intersection return true

  // cycle over dimensions
  for(size_t dim = 0; dim < N; ++dim){
    // check intersection conditions
    if(r2.first[dim] < r1.second[dim] || r1.first [dim] < r2.second[dim]  || // partially overlapping sides
      (r1.first[dim] < r2.first [dim] && r2.second[dim] < r1.second[dim]) || // r2 side contained into r1 side
      (r2.first[dim] < r1.first [dim] && r1.second[dim] < r1.second[dim]))   // r1 side contained into r2 side
      return true;
  }
  return false;
}

template <unsigned int N>
bool TreeSearch<N>::intersect(const SVector<N> &p, const rectangle<N> &r) const {
  // cycle over dimensions
  for(size_t dim = 0; dim < N; ++dim){
    // return false if point lies outside 
    if(p[dim] > r.second[dim] || p[dim] < r.first[dim])
      return false;
  }
  return true;
}

// a searching range is supplied as a pair of points (a,b) where a is the
// lower-left corner and b the upper-right corner of the query rectangle
template <unsigned int N>
std::list<SVector<N>> TreeSearch<N>::search(const rectangle<N> &query) {
  // initialize result
  std::list<SVector<N>> searchResult;

  // use a stack to assist the searching process
  std::stack< node_ptr<MeshNode<N>> > stack;
  stack.push(tree.getNode(0)); // start from root
  
  unsigned int dim = 0; 
  while(!stack.empty()){
    node_ptr<MeshNode<N>> current = stack.top();
    stack.pop();

    // add to solution if point is contained in rectangle
    if(intersect(current->getData().getPoint(), query)){
      searchResult.push_back(current->getData().getPoint());
    }
    
    // get children at node
    std::array<node_ptr<MeshNode<N>>, 2> children = current->getChildren();
    
    bool left_child_test  = children[0] != nullptr ? intersect(children[0]->getData().getRange(), query) : false;
    bool right_child_test = children[1] != nullptr ? intersect(children[1]->getData().getRange(), query) : false;
    
    if(left_child_test)    {   // test if left child range intersects query range
      stack.push(children[0]); // push on stack for later inspection
    }
    if(right_child_test)   {   // test if right child range intersects query range
      stack.push(children[1]); // push on stack for later inspection
    }      
    //dim = (dim+1)%N;
  }
  
  return searchResult;
}

#endif // __TREE_SEARCH_H__

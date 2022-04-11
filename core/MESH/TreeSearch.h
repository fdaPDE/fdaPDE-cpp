#ifndef __TREE_SEARCH_H__
#define __TREE_SEARCH_H__

#include <cstddef>
#include <list>
#include <stack>
#include <utility>

#include "Tree.h"
#include "../OPT/Utils.h"

template <unsigned int N> using rectangle = std::pair<SVector<N>, SVector<N>>;

// in geometric search a query is an N-dimensional rectangle described by its left-lower corner
// and upper-right corner. This class wraps a rectangle providing some facilities to perform geometric
// operations between queries
template <unsigned int N> class Query {
private:
  rectangle<N> queryRange_;
public:
  // constructor
  Query(const rectangle<N>& queryRange) : queryRange_(queryRange) {}

  // returns true if a given point lies inside the query
  bool contains(const SVector<N>& point) const;
  // returns true if the query intersects a given rectangle along a given dimension
  bool intersect(const rectangle<N>& rectangle) const; 

  rectangle<N> getRange() const { return queryRange_; }
};

template <unsigned int N>
bool Query<N>::intersect(const rectangle<N> &rect) const {
  // keep track along which dimension query and rect intersects
  std::array<bool, N> boolVector{};
  
  for(size_t dim = 0; dim < N; ++dim){
    if(// partially overlapping sides
       (queryRange_.second[dim] > rect.first[dim] && rect.second[dim] > queryRange_.first[dim])  ||
       (rect.second[dim] > queryRange_.first[dim] && queryRange_.second[dim] > rect.first[dim])  ||
       // queryRange_ contained into rect 
       (rect.first[dim] < queryRange_.first [dim] && queryRange_.second[dim] < rect.second[dim]) ||
       // rect contained into queryRange_ 
       (queryRange_.first[dim] < rect.first [dim] && rect.second[dim] < rect.second[dim])){

      boolVector[dim] = true; // query intersect rectangle along this dimension
    }
  }

  // query and rect intersects if and only if they intersects along each dimension
  bool result = true;
  for(bool b : boolVector) result &= b;
  
  return result;
}

template <unsigned int N>
bool Query<N>::contains(const SVector<N> &point) const {
  // cycle over dimensions
  for(size_t dim = 0; dim < N; ++dim){
    // return false if point lies outside 
    if(point[dim] > queryRange_.second[dim] || point[dim] < queryRange_.first[dim])
      return false;
  }
  return true;
}

// a specific node data structure for easy management of an ADT structure
template <unsigned int N>
class MeshNode {
private:
  unsigned int elementID_;    // the element ID to which this node referes to
  SVector<N>   point_;       // the point stored in this node 
  rectangle<N> node_range_;  // the range in the hypercube this node refers to
public:
  // constructor
  MeshNode(unsigned int elementID, const SVector<N>& point, const rectangle<N>& node_range)
    : elementID_(elementID), point_(point), node_range_(node_range) {}
  // getters
  SVector<N>   getPoint()     const { return point_;      }
  rectangle<N> getRange()     const { return node_range_; }
  unsigned int getElementID() const { return elementID_;  }
};

// Alternating Decision Tree implementation for tree-based search of elements over an unstructured mesh
template <unsigned int N>
class TreeSearch{
private:
  Tree<MeshNode<N>> tree;

public:
  // initialize an empty decision tree
  TreeSearch() = default;
  // initialize the ADT given a set of N-dimensional points
  TreeSearch(const std::vector<std::pair<SVector<N>, unsigned int>>& data) { init(data); }

  // build an Alternating Decision Tree given a set of N-dimensional points
  void init(const std::vector<std::pair<SVector<N>, unsigned int>>& data);
  // perform a geometric search returning all points which lie in a given query
  std::list<unsigned int> search(const Query<N>& query);

  // getters
  Tree<MeshNode<N>> getTree() const { return tree; }
};

// data is already scaled in the range [0,1]
template <unsigned int N>
void TreeSearch<N>::init(const std::vector<std::pair<SVector<N>, unsigned int>>& data) {

  // initialize here once
  SVector<N> left_lower_corner = SVector<N>::Zero(), right_upper_corner = SVector<N>::Ones();
  
  // initialize tree data structure
  tree = Tree(MeshNode<N>(data[0].second, data[0].first, std::make_pair(left_lower_corner, right_upper_corner)));

  // process all points inside data one by one and insert them in the correct position
  for(size_t j = 1; j < data.size(); ++j){
    SVector<N>   nodeData  = data[j].first;
    unsigned int nodeID    = data[j].second;
    rectangle<N> nodeRange = std::make_pair(left_lower_corner, right_upper_corner);

    // traverse the tree based on coordinates of data and insert the corresponding node at right position in the tree
    node_ptr<MeshNode<N>> current = tree.getNode(0); // root node

    bool inserted = false;         // stop iterating when an insertion point has been found
    unsigned int iteration = 1;    // defines the granularity of the split dimension
    std::array<double,N> offset{}; // keep track of the virtual splits of the domain at each iteration
  
    // search for the right insertion location in the tree
    while(!inserted){
      for(size_t dim = 0; dim < N; ++dim){       // cycle over dimensions
	double split_point = offset[dim] + std::pow(0.5, iteration); // split point
	if(nodeData[dim] < split_point){
      	  nodeRange.second[dim] = split_point;      // shrink node range on the left
	  if(tree.insert(MeshNode<N>(nodeID, nodeData, nodeRange), current->getKey(), LinkDirection::LEFT )){ // O(1) operation
	    inserted = true;                     // stop searching for location
	    break;
	  }else
	    current = current->getChildren()[0]; // move to left child
	}
	else{
	  nodeRange.first[dim] = split_point;       // shrink node range on the right
	  if(tree.insert(MeshNode<N>(nodeID, nodeData, nodeRange), current->getKey(), LinkDirection::RIGHT)){ // O(1) operation
	    inserted = true;                     // stop searching for location
	    break;
	  }else{
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

// a searching range is supplied as a pair of points (a,b) where a is the
// lower-left corner and b the upper-right corner of the query rectangle
template <unsigned int N>
std::list<unsigned int> TreeSearch<N>::search(const Query<N> &query) {
  // initialize result
  std::list<unsigned int> searchResult;

  // use a stack to assist the searching process
  std::stack< node_ptr<MeshNode<N>> > stack;
  stack.push(tree.getNode(0)); // start from root
  
  while(!stack.empty()){
    node_ptr<MeshNode<N>> current = stack.top();
    stack.pop();
        
    // add to solution if point is contained in query range
    if(query.contains(current->getData().getPoint()))
      searchResult.push_back(current->getData().getElementID());
    
    // get children at node
    std::array<node_ptr<MeshNode<N>>, 2> children = current->getChildren();
    
    bool left_child_test  = children[0] != nullptr ? query.intersect(children[0]->getData().getRange()) : false;
    bool right_child_test = children[1] != nullptr ? query.intersect(children[1]->getData().getRange()) : false;
    
    if(left_child_test)        // test if left child range intersects query range
      stack.push(children[0]); 
    if(right_child_test)       // test if right child range intersects query range
      stack.push(children[1]); 
  }
  // search completed
  return searchResult;
}

#endif // __TREE_SEARCH_H__

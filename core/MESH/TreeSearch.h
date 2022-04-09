#ifndef __TREE_SEARCH_H__
#define __TREE_SEARCH_H__

#include "Tree.h"
#include "../OPT/Utils.h"

// Alternating Decision Tree implementation for tree-based search over an unstructured mesh
template <unsigned int N>
class TreeSearch{
private:

  Tree<SVector<N>> tree;
  
public:

  TreeSearch() = default;
  
  void init(std::vector<SVector<N>> data);
};


// data is already scaled in the range [0,1]
template <unsigned int N> void TreeSearch<N>::init(std::vector<SVector<N>> data) {

  // insert root
  tree.insert(data[0]);

  // keeps track of the dimension we have to split
  unsigned int dimension = 1;
  
  for(size_t j = 1; j < data.size(); ++j){

    SVector<N> nodeData = data[j];

    // traverse the tree based on coordinates of data and insert at right position
    node_ptr<SVector<N>> current = tree.getNode(0); // start from root

    // init visit
    bool inserted = false;
    unsigned int iteration = 1;
    std::array<double,N> anchor{}; // used to virtually split the domain at each iteration (spiega meglio)
    
    // insert node
    while(!inserted){
      for(size_t dim = 0; dim < N; ++dim){ // cycle over dimensions
	if(nodeData[dim] < anchor[dim] + std::pow(0.5, iteration)){
	  // try insert node at position
	  if(tree.insert(nodeData, current->getKey(), LinkDirection::LEFT))  // O(1) operation
	    inserted = true;                     // stop searching for location
	  else
	    current = current->getChildren()[0]; // move to left child
	}else{
	  // try insert node at position
	  if(tree.insert(nodeData, current->getKey(), LinkDirection::RIGHT)) // O(1) operation
	    inserted = true;                     // stop searching for location
	  else{
	    current = current->getChildren()[1]; // move to right child
	    anchor[dim] += std::pow(0.5, iteration);
	  }
	}
      }
      // virtually perform an half split of the hipercube
      iteration++;
    }
  }
  
  // construction ended
  return;
}

#endif // __TREE_SEARCH_H__

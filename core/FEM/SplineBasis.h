#ifndef __SPLINE_BASIS_H__
#define __SPLINE_BASIS_H__

#include "../utils/DataStructures/Tree.h"
#include <queue>
#include <utility>
using fdaPDE::core::node_ptr;
using fdaPDE::core::LinkDirection;
#include <cstddef>
using fdaPDE::core::Tree;
#include <vector>
#include <cmath>
#include <stack>
#include <set>
#include <functional>
#include <map>

/* Let u_0, u_1, ..., u_n n distinct knots. Call U = [u_0, u_1, ..., u_n] knot vector. Define B_ij the i-th spline basis of order j.
   By Cox-de Boor formula splines can be recursively defined as

   N_i0(x) = 1 if x \in [u_i, u_i+1) 0 otherwise
   N_ij(x) = [(x-u_i)/(u_i+j - u_i)]*N_i,j-1(x) + [(u_i+j+1 - x)/(u_i+j+1 - u_i+1)]*N_i+1,j-1(x)

             |---------------------|              |-------------------------------|
                       w1                                         w2
	     
   The above recursive formula can be efficiently stored using a binary tree structure. As such a B-Spline will be stored as a collection
   of well defined binary trees. Each node of the tree will contains the knots required for computing weights w1 and w2. Evaluation will take place
   using a DFS visit with a proper defined functor.
*/

// a node of a spline tree
class SplineNode {
private:
  double k1_, k2_; // knots u_i   and u_i+j
  double d1_ = 0;  // for efficiency reasons precompute 1/(u_i+j - u_i)

public:
  // constructor
  SplineNode() = default;
  SplineNode(double k1, double k2) : k1_(k1), k2_(k2), d1_(1/(k2_ - k1_)){ }
  
  std::pair<double, double> getKnotSpan() const { // for a spline of order j get the knot span [u_i, u_i+j)
    return std::make_pair(k1_, k2_);
  }

  double operator()(double x) const { // evaluate the weight for this node
    return (x - k1_)*d1_;
  }

  bool contains(double x) const {
    if(k1_ < k2_) return (x >= k1_ && x < k2_) ? true : false;
    else          return (x >= k2_ && x < k1_) ? true : false;
  }
};

// just a functor wrapping a Tree<SplineNode>
class Spline {
private:
  Tree<SplineNode> spline_; // the actual spline

  // associate to a pair of non-negative integers a unique integer
  inline int cantorPairingFunction(int x, int y) const {
    return (x+y)*(x+y+1)/2 + y;
  }
  inline int cantorPairingFunction(std::pair<int, int> pair) const {
    return cantorPairingFunction(pair.first, pair.second);
  }
  
public:
  // wraps an already existing spline tree in a spline object
  Spline(const Tree<SplineNode>& spline) : spline_(spline) {}

  // build a spline tree given the knot vector and a pair of indexes i,j,
  // where i denotes the knot the spline refers to and j the spline order
  Spline(const std::vector<double>& knotVector, int i, int j) {
    // build a pair from (i,j)
    typedef std::pair<int, int> pairID;
    pairID rootID = std::make_pair(i,j-1);
    
    // insert root in spline tree
    spline_ = Tree<SplineNode>(SplineNode()/*, cantorPairingFunction(rootID)*/);
    
    // use a queue structure to assist the spline build process
    std::queue<pairID> queue{};
    std::queue<unsigned int> queueID{};
    queue.push(rootID); // push root to queue
    queueID.push(spline_.getRoot()->getKey());
    
    // spline construction
    while(!queue.empty()){
      pairID currentPair = queue.front();
      queue.pop();

      int nodeID = queueID.front(); // identifier of this node
      queueID.pop();

      // child pair identifiers
      pairID leftChild  = std::make_pair(currentPair.first,   currentPair.second - 1);
      pairID rightChild = std::make_pair(currentPair.first+1, currentPair.second - 1);
      
      // build left spline node
      SplineNode leftNode  = SplineNode(knotVector[leftChild.first], knotVector[leftChild.first + currentPair.second]);
      auto left_ptr = spline_.insert(leftNode/*, leftChildID*/, nodeID, LinkDirection::LEFT);
      
      // build right spline node
      SplineNode rightNode = SplineNode(knotVector[rightChild.first + currentPair.second], knotVector[rightChild.first]);
      auto right_ptr = spline_.insert(rightNode/*, rightChildID*/, nodeID, LinkDirection::RIGHT);      
      
      // push child nodes to queue for later processing if children nodes are not leaf
      if(currentPair.second - 1 > 0){
	queue.push(leftChild); queue.push(rightChild);
	queueID.push(left_ptr->getKey()); queueID.push(right_ptr->getKey());
      }
    }
    return;
  }
  
  // evaluates the spline at a given point
  double operator()(double) const;

  // get internal data structure
  Tree<SplineNode> getTree() const { return spline_; }
};

double Spline::operator()(double x) const {
  // perform a slight variation of the DFS search to collect the spline evaluation
  std::set<unsigned int> visitedNodes{};     // set containing the IDs of visited nodes
  node_ptr<SplineNode> currentNode = spline_.getRoot();
  
  double result  = 0;
  std::map<int, double> partials{};
  partials[currentNode->getKey()] = 1;
  
  double partial = 1;
  LinkDirection lastLink;
  
  // start visit
  while(visitedNodes.size() < spline_.getNumberOfNodes()){ // repeat until tree not completely explored
    if(currentNode->isLeaf()){ // leaf node, this is an N_i0 node
      std::pair<double, double> support = currentNode->getData().getKnotSpan();

      if(currentNode->getData().contains(x)){ // add result only if N_i0(x) = 1
	result += partial;
      }
      visitedNodes.insert(currentNode->getKey()); // leaf visited
      
      // back to previous node 
      currentNode = currentNode->getFather();

      partial = partials[currentNode->getKey()];
    }else{
      // search for child to visit
      std::size_t i = 0;
      while(i < 2){
	auto childNode = currentNode->getChildren()[i];
	
	if(childNode != nullptr &&
	   visitedNodes.find(childNode->getKey()) == visitedNodes.end()){ // node still not marked as visited

          lastLink = static_cast<LinkDirection>(i);             // keep track of last move
	  partial *= childNode->getData()(x); // update result
	  partials[childNode->getKey()] = partial;

          currentNode = childNode;                              // move to child
	  visitedNodes.insert(currentNode->getKey());           // insert child as visited node
	  break;
	}
	++i;
      }
      // if all child visited mark this node as visited and go back to father
      if(i == 2){

        visitedNodes.insert(currentNode->getKey());
	// back to previous node
	if(currentNode->getFather() != nullptr){ // only root has nullptr father
	  currentNode = currentNode->getFather();
	  partial = partials[currentNode->getKey()];
	}
      }
    }
  }
  return result;
}

// spline basis implementation
class SplineBasis{
private:
  std::vector<double> knotsVector_{};
  std::vector<Spline> basis_{}; // the spline basis. Each tree represent a single spline

public:
  // constructor
  SplineBasis(const std::vector<double>& knotsVector, int order) : knotsVector_(knotsVector) {
    // build spline basis
    for(std::size_t k = 0; k < knotsVector_.size() - order; ++k){ // create spline at each iteration
      basis_.push_back(Spline(knotsVector_, k, order));
    }
  }
  Spline& operator[](std::size_t i) { return basis_[i]; }
};

#endif // __SPLINE_BASIS_H__

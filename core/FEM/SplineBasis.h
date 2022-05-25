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
    spline_ = Tree<SplineNode>(SplineNode());
    
    // use a queue structure to assist the spline build process
    std::queue<std::pair<unsigned int, pairID>> queue{};
    //    std::queue<unsigned int> queueID{};
    queue.push(std::make_pair(spline_.getRoot()->getKey(),rootID)); // push root to queue
    //    queueID.push(spline_.getRoot()->getKey());
    
    // spline construction
    while(!queue.empty()){
      std::pair<unsigned int, pairID> currentNode = queue.front();
      queue.pop();

      int nodeID = currentNode.first;           // identifier of this node
      pairID currentPair = currentNode.second;  // pair (i,j) of spline node N_ij(x)
      
      // build left spline node
      SplineNode leftNode  = SplineNode(knotVector[currentPair.first],
					knotVector[currentPair.first + currentPair.second]);
      auto left_ptr = spline_.insert(leftNode, nodeID, LinkDirection::LEFT);
      
      // build right spline node
      SplineNode rightNode = SplineNode(knotVector[currentPair.first + 1 + currentPair.second],
					knotVector[currentPair.first + 1]);
      auto right_ptr = spline_.insert(rightNode, nodeID, LinkDirection::RIGHT);      
      
      // push child nodes to queue for later processing if children nodes are not leafs
      if(currentPair.second - 1 > 0){
	queue.push(std::make_pair(left_ptr->getKey(),
				  std::make_pair(currentPair.first,     currentPair.second - 1)));
	queue.push(std::make_pair(right_ptr->getKey(),
				  std::make_pair(currentPair.first + 1, currentPair.second - 1)));
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

  struct ResultCollector{
    double result_ = 0;
    std::map<int, double> partialMap_{};
    double partialProduct_ = 1;
    double inputPoint_;
    
    ResultCollector(const node_ptr<SplineNode>& root_ptr, double inputPoint) : inputPoint_(inputPoint) {
      partialMap_[root_ptr->getKey()] = 1;
    };

    // add result only if N_i0(x) = 1, that is if x is contained in the knot span [u_i, u_i+1)
    void leafAction(const node_ptr<SplineNode>& currentNode) {
      if(currentNode->getData().contains(inputPoint_))
	result_ += partialProduct_;

      // restore partialProduct to the one of the father
      partialProduct_ = partialMap_[currentNode->getFather()->getKey()];
      return;
    }
    
    void firstVisitAction(const node_ptr<SplineNode>& currentNode) {
      // update partial product
      partialProduct_ *= currentNode->getData()(inputPoint_);
      partialMap_[currentNode->getKey()] = partialProduct_;
      return;
    }

    void lastVisitAction(const node_ptr<SplineNode>& currentNode) {
      // restore to partialProduct of the father
      partialProduct_ = partialMap_[currentNode->getFather()->getKey()];
      return;
    }
  };
  
  // perform a DFS search to collect the spline evaluation
  ResultCollector resultObject = ResultCollector(spline_.getRoot(), x);
  spline_.DFS(resultObject);
  return resultObject.result_;;
}

// spline basis implementation
class SplineBasis{
private:
  std::vector<double> knotsVector_{};  // vector of knots
  std::vector<Spline> basis_{};        // the spline basis.

public:
  // constructor
  SplineBasis(const std::vector<double>& knotsVector, int order) : knotsVector_(knotsVector) {
    // build spline basis
    for(std::size_t k = 0; k < knotsVector_.size() - order; ++k){ // create spline at each iteration
      basis_.push_back(Spline(knotsVector_, k, order));
    }
  }

  // return i-th element of the basis
  Spline& operator[](std::size_t i) { return basis_[i]; }
};

#endif // __SPLINE_BASIS_H__

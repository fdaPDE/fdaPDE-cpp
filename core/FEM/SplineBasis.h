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
#include <tuple>

/* Let u_0, u_1, ..., u_n n distinct knots. Call U = [u_0, u_1, ..., u_n] knot vector. Define B_ij the i-th spline basis of order j.
   By Cox-de Boor formula splines can be recursively defined as

   N_i0(x) = 1 if x \in [u_i, u_i+1) 0 otherwise
   N_ij(x) = [(x-u_i)/(u_i+j - u_i)]*N_i,j-1(x) + [(u_i+j+1 - x)/(u_i+j+1 - u_i+1)]*N_i+1,j-1(x)
             |---------------------|              |-------------------------------|
                       wL                                         wR
	     
   The above recursive formula can be efficiently stored using a binary tree structure. For example a cubic spline will be represented by the following
   recursion tree:

                          N_02
	   wL*N_01                    wR*N_11
   wL*N_00         wR*N_10    wL*N_10         wR*N_20
   
   As such a B-Spline will be stored as a collection of well defined binary trees. Each node of the tree will contains the
   knots required for computing weights w1 and w2. Evaluation will take place using a DFS visit with a proper defined functor.
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
  // the actual spline. This is a vector since spline derivatives are stored as a a vector of lower degree spline trees
  std::vector<std::tuple<double, Tree<SplineNode>, int, int> > spline_{}; 
  int i_, j_;
  std::vector<double> knotVector_;

  // evaluate a single tree in the spline_ vector
  double eval(const Tree<SplineNode>& tree, double x) const;
  // extract i-th subtree
  Tree<SplineNode>& operator[](std::size_t i) { return std::get<1>(spline_[i]); } 
  
public:
  // wraps an already existing spline tree in a spline object
  Spline(const Tree<SplineNode>& spline, int i, int j) { spline_.push_back(std::make_tuple(1, spline, i, j)); }

  // build a spline tree given the knot vector and a pair of indexes i,j,
  // where i denotes the knot the spline refers to and j the spline order
  Spline(const std::vector<double>& knotVector, int i, int j);

  // initialize a spline from a vector of spline trees and associated weights (used by gradient() method)
  Spline(const std::vector<std::tuple<double, Tree<SplineNode>, int, int> >& spline) : spline_(spline) {}
  
  double operator()(double) const;    // evaluates the spline vector at a given point
  Spline gradient() const;            // compute derivative of spline. This is another spline object
};

class SplineBasis{
private:
  std::vector<double> knotsVector_{};  // vector of knots
  std::vector<Spline> basis_{};        // the spline basis.

public:
  // constructor
  SplineBasis(const std::vector<double>& knotsVector, int order) {
    
    // pad the knot vector to obtain a full basis for the whole knot span [u_0, u_n]
    for(std::size_t i = 0; i < order; ++i) 
      knotsVector_.push_back(knotsVector[0]);
    knotsVector_.insert(knotsVector_.end(), knotsVector.begin(), knotsVector.end());
    for(std::size_t i = 0; i < order; ++i)
      knotsVector_.push_back(knotsVector[knotsVector.size()-1]);

    // build spline basis
    for(std::size_t k = 0; k < knotsVector_.size() - order; ++k){ // create spline at each iteration
      basis_.push_back(Spline(knotsVector_, k, order));
    }
  }

  Spline& operator[](std::size_t i) { return basis_[i]; }   // return i-th element of the basis
  int getNumberOfBasis() const { return basis_.size(); }    // return the number of basis elements
};

#include "SplineBasis.tpp"

#endif // __SPLINE_BASIS_H__

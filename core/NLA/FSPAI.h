#ifndef __FSPAI_H__
#define __FSPAI_H__

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Cholesky>

#include <Eigen/src/Cholesky/LLT.h>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <ostream>
#include <unordered_set>
#include <utility>
#include <iostream>

#include "../utils/Symbols.h"

// an implementation of the Factorized Sparse Approximate Inverse algorithm with sparsity pattern update.
// FSPAI assumes that the square, n x n, sparse matrix A of which we want to compute the inverse is SPD, in this sense there exists a lower
// triangular matrix L_A such that A = L_A.transpose()*L_A. FSPAI finds an approximate inverse for the Cholesky factor L_A of matrix A
// while keeping L_A sparse (in general indeed the inverse of a sparse matrix is not generally sparse, i.e. it could be dense)

// This FSPAI implementation is based on the minimization of the K-condition number of matrix A

class FSPAI{
private:
  const Eigen::SparseMatrix<double>& A_;    // const reference to target matrix
  Eigen::Index n_;                          // dimension of matrix A_

  std::vector<std::unordered_set<std::size_t>> J_; // sparsity pattern of approximate inverse

  void updateSparsityPattern();
  
public:
  // constructor
  FSPAI(const Eigen::SparseMatrix<double>& A) : A_(A), n_(A.rows()) {

    // initialize the sparsity pattern to the identity matrix
    J_.resize(n_);
    for(std::size_t i = 0; i < n_; ++i){
      J_[i].insert(i);
    }
  }

  // a sparsity pattern is a set of indexes corresponding to non-zero entries of a matrix
  typedef std::unordered_set<std::size_t> column_sparsity_pattern;

  // build system matrix A(p1, p2) given sparsity patterns p1 and p2. The result is a |p1| x |p2| dense matrix
  DMatrix extractMatrix(const column_sparsity_pattern& p1,
			const column_sparsity_pattern& p2){
    DMatrix result;
    result.resize(p1.size(), p2.size());

    // (*it1.first, *it2.first) is each time a pair of coordinates in the cross product  X tildeJk
    for(auto it1 = std::make_pair(p1.begin(), 0);
	it1.first != p1.end();
	it1.first++, it1.second++){
      for(auto it2 = std::make_pair(p2.begin(), 0);
	  it2.first != p2.end();
	  it2.first++, it2.second++){

	//std::cout << *it1.first << " - " << *it2.first << " : " << A_.coeff(*it1.first, *it2.first) << std::endl;
	
	result(it1.second, it2.second) = A_.coeff(*it1.first, *it2.first);
      }
    }
    return result;
  }

  // build vector A(p, k) given column k (fixed) and sparsity pattern p. The result is a |p| x 1 dense column
  DVector extractVector(unsigned col, const column_sparsity_pattern& p){
    DVector result;
    result.resize(p.size(), 1);

    for(auto it1 = std::make_pair(p.begin(), 0);
	it1.first != p.end();
	it1.first++, it1.second++){

      result[it1.second] = A_.coeff(*it1.first, col);
    }
    return result;
  }
  
  // compute the Factorize Sparse Approximate Inverse of A using its K-condition number minimization
  Eigen::SparseMatrix<double> compute(unsigned alpha, unsigned beta, double epsilon){

    Eigen::SparseMatrix<double> result;
    result.resize(n_, n_);

    std::list<Eigen::Triplet<double>> tripetList;
    
    // cycle over each dimension of matrix A
    for(std::size_t k = 0; k < n_; ++k){
      Eigen::Matrix<double, Eigen::Dynamic, 1> Lk;
      Lk.resize(n_, 1);
      Lk.fill(0);
      
      // perform alpha steps of sparsity pattern update along column k
      for(std::size_t s = 0; s < alpha; ++s){
	column_sparsity_pattern tildeJk = J_[k];	
	tildeJk.erase(k); 
	
	if(tildeJk.empty()){
	  // skip computation if tildeJk is empty. No linear system to solve here, just compute the diagonal
	  // element of the current approximate inverse
	  double l_kk = 1/(std::sqrt(A_.coeff(k, k)));
	  tripetList.push_back(Eigen::Triplet<double>(k, k, l_kk));
	  Lk[k] = l_kk;
	}else{

	  // solve optimization problem: find best vector in the set of all vectors having current sparsity pattern
	  // minimizing the K-condition number of L^T*A*L. It can be proven that this problem is equivalent to the solution
	  // of a small dense SPD linear system.
	  
	  // thanks to sparsity of matrix A, given the current sparsity pattern tildeJk we extract just
	  // the dense matrix Ak = A(tildeJk, tildeJk) as well as the dense vector bk = A_k(tildeJk)
	  DMatrix Ak = extractMatrix(tildeJk, tildeJk);  // define system matrix
	  DVector bk = extractVector(k, tildeJk);        // define rhs of linear system

          Eigen::Matrix<double, Eigen::Dynamic, 1> yk;
	  yk.resize(Ak.rows(), 1);
	  
	  // solve linear system AK*yk = bk using Cholesky factorization
	  Eigen::LLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> choleskySolver;
	  choleskySolver.compute(Ak);
	  yk = choleskySolver.solve(bk);
	
	  // compute diagonal entry l_kk
	  double l_kk = A_.coeff(k,k) - bk.transpose().dot(yk);
	  l_kk = 1/(std::sqrt(l_kk));
	  Lk[k] = l_kk;
	  
	  for(auto it = std::make_pair(tildeJk.begin(), 0);
	      it.first != tildeJk.end();
	      ++it.first, ++it.second){

	    Lk[*it.first] = -l_kk*yk[it.second];
	  }	  	  
	  
	  // update current solution
	  tripetList.push_back(Eigen::Triplet<double>(k,k,l_kk));
	  for(auto it = std::make_pair(tildeJk.begin(), 0); it.first != tildeJk.end(); ++it.first, ++it.second){
	    tripetList.push_back(Eigen::Triplet<double>(*it.first, k, Lk[*it.first]));
	  }	  
	}

	// search for best update of sparsity pattern
	if(s != alpha){ // maximum number of sparsity updates not reached

          // computation of candidate indexes for sparsity pattern update
	  std::unordered_map<std::size_t, double> hatJk{};
	  for(std::size_t j = k+1; j < n_; ++j){
	    double v = 0;
	    for(auto it = std::make_pair(J_[k].begin(), 0); it.first != J_[k].end(); ++it.first, ++it.second){
              // compute A(j, jK)*Lk(jK)
	      v += A_.coeff(j, *it.first)*Lk[*it.first];
	    }
	    // index j is a possible candidate for updating the sparsity pattern of the solution along column k
	    if(v != 0 && J_[k].find(j) == J_[k].end()){
	      hatJk.insert(std::make_pair(j, v));
	    }
	  }

	  // compute average score
	  double tau_k = 0;
	  std::map<std::size_t, double> tauMap;
	  double max_tau = 0; // tau > 0, ok to use 0 as starting treshold

	  for(auto it = hatJk.begin(); it != hatJk.end(); ++it){
	    double tau_jk = std::pow(it->second, 2)/A_.coeff(it->first, it->first);
	    tauMap.insert(std::make_pair(it->first, tau_jk));
	    tau_k += tau_jk;
	    if(tau_jk > max_tau)
	      max_tau = tau_jk;
	  }
	  // the end iterator points to the maximum tau_jk computed (O(1) operation)
	  if(max_tau >= epsilon){ // best improvement higher than treshold
	    tau_k /= hatJk.size();

            // select most promising first beta entries
	    for(std::size_t idx = 0; idx < beta && !tauMap.empty(); ++idx){
	      // find maximum element
	      auto max = tauMap.begin();
	      for(auto it = tauMap.begin(); it != tauMap.end(); ++it){
		if(it->second > max->second) max = it; 
	      }
	      
              if(max->second > tau_k) // update sparsity pattern
		J_[k].insert(max->first);

              tauMap.erase(max);
	    }
	  }
	}
      }
    }

    result.setFromTriplets(tripetList.begin(), tripetList.end());
    return result;
  }

  
};

#endif // __FSPAI_H__

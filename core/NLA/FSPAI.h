#ifndef __FSPAI_H__
#define __FSPAI_H__

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Cholesky>

#include <map>
#include <unordered_set>
#include <vector>

#include "../utils/Symbols.h"

// an implementation of the Factorized Sparse Approximate Inverse algorithm with sparsity pattern update.
// FSPAI assumes that the square, n x n, sparse matrix A of which we want to compute the inverse is SPD, in this sense there exists a lower
// triangular matrix L_A such that A = L_A.transpose()*L_A. FSPAI finds an approximate inverse for the Cholesky factor L_A of matrix A
// while keeping L_A sparse (in general indeed the inverse of a sparse matrix is not generally sparse, i.e. it could be dense)

// This FSPAI implementation is based on the minimization of the K-condition number of matrix A. A must be SPD
class FSPAI{
private:
  // a sparsity pattern is a set of indexes corresponding to non-zero entries of the matrix
  typedef std::unordered_set<std::size_t> column_sparsity_pattern;

  const Eigen::SparseMatrix<double>& A_;    // const reference to target matrix
  Eigen::SparseMatrix<double> L_;           // the sparse approximate inverse of the Cholesky factor of A_

  // internal status data members
  Eigen::Index n_;                              // dimension of square sparse matrix A_
  std::vector<column_sparsity_pattern> J_;      // sparsity pattern of approximate inverse
  Eigen::Matrix<double, Eigen::Dynamic, 1> Lk_; // the k-th column of the approximate inverse of the cholesky factor L_
      
  // build system matrix A(p1, p2) given sparsity patterns p1 and p2. The result is a |p1| x |p2| dense matrix
  DMatrix extractMatrix(const column_sparsity_pattern& p1,
			const column_sparsity_pattern& p2) const {
    DMatrix result;
    result.resize(p1.size(), p2.size());

    // (*it1.first, *it2.first) is each time a pair of coordinates in the cross product  X tildeJk
    for(auto it1 = std::make_pair(p1.begin(), 0);
	it1.first != p1.end();
	it1.first++, it1.second++){
      for(auto it2 = std::make_pair(p2.begin(), 0);
	  it2.first != p2.end();
	  it2.first++, it2.second++){

	result(it1.second, it2.second) = A_.coeff(*it1.first, *it2.first);
      }
    }
    return result;
  }

  // build vector A(p, k) given column k (fixed) and sparsity pattern p. The result is a |p| x 1 dense column
  DVector extractVector(unsigned col, const column_sparsity_pattern& p) const {
    DVector result;
    result.resize(p.size(), 1);

    for(auto it1 = std::make_pair(p.begin(), 0);
	it1.first != p.end();
	it1.first++, it1.second++){

      result[it1.second] = A_.coeff(*it1.first, col);
    }
    return result;
  }

  // update approximate inverse of column k
  void updateApproximateInverse(const Eigen::Index& k, const DVector& bk, const DVector& yk,
				const column_sparsity_pattern& tildeJk){
    // compute diagonal entry l_kk
    double l_kk = A_.coeff(k,k) - bk.transpose().dot(yk);
    l_kk = 1/(std::sqrt(l_kk));
    Lk_[k] = l_kk;

    // update other entries according to current sparsity pattern tildeJk
    for(auto it = std::make_pair(tildeJk.begin(), 0);
	it.first != tildeJk.end();
	++it.first, ++it.second){

      Lk_[*it.first] = -l_kk*yk[it.second];
    }
    return;
  }

  // select candidate indexes for sparsity pattern update for column k
  std::unordered_map<std::size_t, double> selectCandidates(const Eigen::Index& k) const {
    std::unordered_map<std::size_t, double> result{};
    for(std::size_t j = k+1; j < n_; ++j){
      double v = 0;
      // compute A(j, jK)*Lk(jK)
      for(auto it = J_[k].begin(); it != J_[k].end(); ++it){
	v += A_.coeff(j, *it)*Lk_[*it]; 
      }
      // index j is a possible candidate for updating the sparsity pattern of the solution along column k
      if(v != 0 && J_[k].find(j) == J_[k].end()){
	result.insert(std::make_pair(j, v));
      }
    }
    return result;
  }
  
public:
  // constructor
  FSPAI(const Eigen::SparseMatrix<double>& A) : A_(A), n_(A.rows()) {

    // initialize the sparsity pattern to the identity matrix
    J_.resize(n_);
    for(std::size_t i = 0; i < n_; ++i){
      J_[i].insert(i);
    }

    // pre-allocate memory
    L_.resize(n_, n_);
    Lk_.resize(n_, 1); 
  }

  // getters

  // returns the Cholesky factor of matrix A_
  const Eigen::SparseMatrix<double>& getL() const { return L_; }
  // returns the approximate inverse of A_
  Eigen::SparseMatrix<double> getInverse()  const { return L_*L_.transpose(); }
  
  // compute the Factorize Sparse Approximate Inverse of A using its K-condition number minimization
  // alpha:   number of sparsity pattern updates to compute for each column k of A_
  // beta:    number of indexes to augment the sparsity pattern of Lk_ per update step
  // epsilon: do not consider an entry of A_ as valid if it causes a reduction to its K-condition number lower than epsilon
  void compute(unsigned alpha, unsigned beta, double epsilon){
    
    // eigen requires a list of triplet to construct a sparse matrix in an efficient way
    std::list<Eigen::Triplet<double>> tripetList; 
    
    // cycle over each column of the sparse matrix A_
    for(std::size_t k = 0; k < n_; ++k){
      Lk_.fill(0); // reset column vector Lk_
      
      // perform alpha steps of inverse update along column k
      for(std::size_t s = 0; s < alpha; ++s){
	column_sparsity_pattern tildeJk = J_[k]; // extract sparsity pattern
	tildeJk.erase(k); 
	
	if(tildeJk.empty()){
	  // skip computation if tildeJk is empty. No linear system to solve here, just compute the diagonal
	  // element of the current approximate inverse
	  double l_kk = 1/(std::sqrt(A_.coeff(k, k)));
	  tripetList.push_back(Eigen::Triplet<double>(k, k, l_kk));
	  Lk_[k] = l_kk;
	}else{
	  // we must find the best vector fixed its sparsity pattern minimizing the K-condition number of L^T*A*L.
	  // It can be proven that this problem is equivalent to the solution of a small dense SPD linear system
	  //         A(tildeJk, tildeJk)*yk = Ak(tildeJk)
	  
	  DMatrix Ak = extractMatrix(tildeJk, tildeJk);  // define system matrix:        A(tildeJk, tildeJk)
	  DVector bk = extractVector(k, tildeJk);        // define rhs of linear system: Ak(tildeJk)

          Eigen::Matrix<double, Eigen::Dynamic, 1> yk;
	  yk.resize(Ak.rows(), 1);
	  
	  // solve linear system A(tildeJk, tildeJk)*yk = Ak(tildeJk) using Cholesky factorization (system is SPD)
	  Eigen::LLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> choleskySolver;
	  choleskySolver.compute(Ak);
	  yk = choleskySolver.solve(bk); // compute yk = A(tildeJk, tildeJk)^{-1}*Ak(tildeJk)

	  // update approximate inverse
	  updateApproximateInverse(k, bk, yk, tildeJk);
	  
	  // update current solution
	  tripetList.push_back(Eigen::Triplet<double>(k,k,Lk_[k]));
	  for(auto it = std::make_pair(tildeJk.begin(), 0);
	      it.first != tildeJk.end();
	      ++it.first, ++it.second){
	    
	    tripetList.push_back(Eigen::Triplet<double>(*it.first, k, Lk_[*it.first]));
	  }	  
	}

	// search for best update of the sparsity pattern
	if(s != alpha){
          // computation of candidate indexes for sparsity pattern update
	  std::unordered_map<std::size_t, double> hatJk = selectCandidates(k);

	  // use average value heuristic for selection of best candidates
	  double tau_k = 0;   // the average improvement to the K-condition number
	  double max_tau = 0; // tau > 0, ok to use 0 as starting treshold for maximum value search
	  
	  for(auto it = hatJk.begin(); it != hatJk.end(); ++it){
	    double tau_jk = std::pow(it->second, 2)/A_.coeff(it->first, it->first);
	    hatJk[it->first] = tau_jk;
	    tau_k += tau_jk;

	    // update maximum value
	    if(tau_jk > max_tau) max_tau = tau_jk;
	  }

	  // if the best improvement is higher than accetable treshold
          if(max_tau >= epsilon){
	    tau_k /= hatJk.size();

            // select most promising first beta entries according to average heuristic
	    for(std::size_t idx = 0; idx < beta && !hatJk.empty(); ++idx){
	      // find maximum
	      auto max = hatJk.begin();
	      for(auto it = hatJk.begin(); it != hatJk.end(); ++it){
		if(it->second > max->second) max = it; 
	      }
	      
              if(max->second > tau_k) // sparsity pattern update
		J_[k].insert(max->first);

	      // erase current maximum to start a new maximum search
              hatJk.erase(max);
	    }
	  }
	}
      }
    }

    // build final result
    L_.setFromTriplets(tripetList.begin(), tripetList.end());
    return;
  }
  
};

#endif // __FSPAI_H__

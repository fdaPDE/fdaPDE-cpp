#ifndef __K_FOLD_CV__
#define __K_FOLD_CV__

#include <algorithm>
#include <vector>
#include <thread> // multithreading support
#include "../core/utils/Symbols.h"
#include "../core/utils/DataStructures/BlockFrame.h"
#include "../models/ModelTraits.h"

namespace fdaPDE{
namespace calibration{

  // general implementation of KFold Cross Validation
  class KFoldCV {
  private:  
    std::size_t K_; // number of folds

    std::vector<double> avg_scores_; // mean CV score for each \lambda
    std::vector<double> std_scores_; // CV score standard deviation for each \lambda
    DMatrix<double> scores_;         // matrix of CV scores
    DVector<double> optimum_;        // optimal smoothing parameter
    
    // split data in K folds
    std::pair<BlockView<Sparse, double, int>, BlockView<Sparse, double, int>>
    split(const BlockFrame<double, int>& data, std::size_t i) {
      std::size_t n = data.rows();      // number of data points
      std::size_t m = std::floor(n/K_); // number of data per fold
    
      // create test-train index sets, as function of test fold i
      std::vector<std::size_t> testIdx(m);
      std::vector<std::size_t> trainIdx(n-m);
      for(std::size_t j = 0; j < n; ++j){
	if(j >= m*i && j < m*(i+1))
	  testIdx[j-m*i] = j;
	else
	  trainIdx[j >= m*(i+1) ? j-m : j] = j;
      }
      // create views
      BlockView<Sparse, double, int> train = data(trainIdx);
      BlockView<Sparse, double, int> test  = data(testIdx);
      return std::make_pair(train, test);
    }
    
  public:
    // constructor
    KFoldCV() = default;
    KFoldCV(std::size_t K) : K_(K) {};

    // selects best smoothing parameter according to a K-fold cross validation strategy using the output of functor F as CV score
    // F receives, in this order: current smoothing parameter, train set and test set
    void compute
    (const std::vector<DVector<double>>& lambdas, const BlockFrame<double, int>& data, 
     const std::function<double(DVector<double>, BlockFrame<double, int>, BlockFrame<double, int>)>& F,
     bool randomize = true){ // if true a randomization of the data is performed before split
      // reserve space for CV scores
      scores_.resize(K_, lambdas.size());
      BlockFrame<double, int> data_;
      if(randomize) // perform a first shuffling of the data if required	
	data_ = data.shuffle();
      else data_ = data;
      
      // cycle over all folds
      for(std::size_t fold = 0; fold < K_; ++fold){
	// create train test partition
	std::pair<BlockView<Sparse, double, int>, BlockView<Sparse, double, int>> train_test = split(data_, fold);
	// decouple training from testing
	BlockFrame<double, int> train = train_test.first.extract();
	BlockFrame<double, int> test = train_test.second.extract();
	
	// fixed a data split, cycle over all lambda values
	for(std::size_t j = 0; j < lambdas.size(); ++j)
	  scores_.coeffRef(fold, j) = F(lambdas[j], train, test); // compute CV score
      }
      // reserve space for storing results
      avg_scores_.clear(); std_scores_.clear(); // clear possible previous execution
      avg_scores_.reserve(lambdas.size());
      std_scores_.reserve(lambdas.size());
      for(std::size_t j = 0; j < lambdas.size(); ++j){
	// record the average score and its standard deviation computed among the K_ runs
	double avg_score = 0;
	for(std::size_t i = 0; i < K_; ++i) avg_score += scores_(i,j);
	avg_score /= K_;
	double std_score = 0;
	for(std::size_t i = 0; i < K_; ++i) std_score += std::pow(scores_(i,j) - avg_score, 2);
	std_score = std::sqrt(std_score/(K_ - 1)); // use unbiased sample estimator for standard deviation
	// store results
	avg_scores_.push_back(avg_score);
	std_scores_.push_back(std_score);
      }

      // store optimal lambda according to given metric F
      std::vector<double>::iterator opt_score = std::min_element(avg_scores_.begin(), avg_scores_.end());
      optimum_ = lambdas[std::distance(avg_scores_.begin(), opt_score)];
    }
    
    // getters
    std::vector<double> avg_scores() const { return avg_scores_; } // mean CV score vector
    std::vector<double> std_scores() const { return std_scores_; } // CV score standard deviation vector
    DMatrix<double> scores() const { return scores_; } // CV scores
    DVector<double> optimum() const { return optimum_; } // optimal smoothing level according to provided CV index
    // setters
    void set_K(std::size_t K) { K_ = K; }
  };
  
}}

#endif // __K_FOLD_CV__

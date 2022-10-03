#ifndef __K_FOLD_CV__
#define __K_FOLD_CV__

#include "../core/utils/Symbols.h"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <vector>
#include <memory>

// just a pair of values
template <typename T>
struct DataSet {
  std::shared_ptr<T> data;
  std::vector<std::size_t> indices;

  // constructor
  DataSet() = default;
  DataSet(const std::shared_ptr<T>& data_, const std::vector<std::size_t>& indices_)
    : data(data_), indices(indices_) {};
};

template <typename M>
class KFoldCV {
private:  
  M& model_;      // model to validate
  std::size_t K_; // number of folds

  std::vector<double> avg_scores_{}; // vector recording the mean score of the model for each value of lambda
  std::vector<double> std_scores_{}; // vector recording the standard deviation of the recordered scores
  DMatrix<double> scores_{};         // matrix of scores (one column for each explored lambda value, one row for each fold)
  
  template <typename T> using TrainTestSet = std::pair<DataSet<T>, DataSet<T>>;
  
  // split data in K folds
  template <typename T>
  std::pair<DataSet<T>, DataSet<T>> split(const T& data, std::size_t i) {
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

    // create train and test sets
    DataSet<T> train
      (std::make_shared<T>(data(trainIdx, Eigen::all)),
       trainIdx);
    DataSet<T> test
      (std::make_shared<T>(data(testIdx,  Eigen::all)),
       testIdx);

    return std::make_pair(train, test);
  }
    
public:
  // constructor
  KFoldCV(M& model, std::size_t K) : model_(model), K_(K) {};

  // select best smoothing parameter according to a K-fold cross validation strategy using the output of functor F as model score
  template <typename F>
  SVector<1> compute(const std::vector<SVector<1>>& lambdas, const F& scoreFunctor){    
    // reserve space for storing scores
    scores_.resize(K_, lambdas.size());
    // cycle over all folds (execute just K_ splits of the data, very expensive operation)
    for(std::size_t fold = 0; fold < K_; ++fold){
      // create train test partition
      TrainTestSet<DMatrix<double>> W_ = split(*model_.W(), fold);
      TrainTestSet<DVector<double>> z_ = split(*model_.z(), fold);      
      // fixed a data split, cycle over all lambda values
      for(std::size_t j = 0; j < lambdas.size(); ++j){
	M m(model_);  // create a fresh copy of current model (O(1) operation)
	m.setLambda(lambdas[j][0]);  // set current lambda

	// fit the model on training set
	m.setObservations(*z_.first.data, z_.first.indices);
	m.setCovariates(*W_.first.data);
	m.smooth(); // fit the model

	// evaluate model score on the fold left out (test set)
	scores_.coeffRef(fold, j) = scoreFunctor(m, z_.second.data, W_.second.data, W_.second.indices);
      }
    }
    // reserve space for storing results
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

    // return optimal lambda according to given metric F
    std::vector<double>::iterator opt_score = std::min_element(avg_scores_.begin(), avg_scores_.end());
    return lambdas[std::distance(avg_scores_.begin(), opt_score)];
  }

  // getters
  std::vector<double> avg_scores() const { return avg_scores_; }
  std::vector<double> std_scores() const { return std_scores_; }
  DMatrix<double> scores() const { return scores_; } // table of all scores (K_ x |lambdas| matrix)
  
};

#endif // __K_FOLD_CV__

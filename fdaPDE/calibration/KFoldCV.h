#ifndef __K_FOLD_CV__
#define __K_FOLD_CV__

#include "../core/utils/Symbols.h"
#include <cmath>
#include <cstddef>
#include <Eigen/Core>
#include <iterator>
#include <utility>
#include <vector>
#include <memory>

template <typename M>
class KFoldCV {
private:

  // a dataset is a pair made by: the actual data and a vector of indices corresponding to data locations
  template <typename T> using dataset = std::pair<std::shared_ptr<T>, std::vector<std::size_t>>;
  
  // split data in K folds returning a train test composed by all the folds but the i-th one, which is used as test set
  template <typename T>
  std::pair<dataset<T>, dataset<T>> trainTestSplit(const T& data, std::size_t i) {
    std::size_t n = data.rows();      // number of data points
    std::size_t m = std::floor(n/K_); // number of data per fold

    // create test-train idx sets, depending on fold of interest i
    std::vector<std::size_t> testIdx(m);
    std::vector<std::size_t> trainIdx(n-m);
    for(std::size_t j = 0; j < n; ++j){
      if(j >= m*i && j < m*(i+1))
	testIdx[j-m*i] = j;
      else
	trainIdx[j >= m*(i+1) ? j-m : j] = j;
    }

    // create train set as the concatenation of K-1 folds
    std::shared_ptr<T> train = std::make_shared<T>(data(trainIdx, Eigen::all));
    // create test set using the K-th fold left out
    std::shared_ptr<T> test  = std::make_shared<T>(data(testIdx , Eigen::all));
    
    return std::make_pair(std::make_pair(train, trainIdx), std::make_pair(test, testIdx));
  }
  
  // compute RMSE: \sqrt{\frac{norm(z - \hat z)^2/}{n}}
  double RMSE(const DVector<double>& z, const DVector<double>& z_hat) const {
    return std::sqrt((z - z_hat).squaredNorm()/z.rows());
  }

  M& model_;
  std::size_t K_;
  
public:
  // constructor
  KFoldCV(M& model, std::size_t K) : model_(model), K_(K) {};

  // select best smoothing parameter according to a K-fold cross validation strategy
  SVector<1> compute(const std::vector<SVector<1>>& lambdas){
    // vector of RMSEs
    std::vector<double> RMSEs{};

    // cycle over all lambdas
    for(SVector<1> lambda : lambdas){
      M m(model_);             // create a copy of current model
      double rmse = 0;         // keep track of average RMSE

      // cycle over all folds
      for(std::size_t f = 0; f < K_; ++f){
	// create train test partition
	std::pair<dataset<DMatrix<double>>, dataset<DMatrix<double>>> train_test_W = trainTestSplit(*model_.W(), f);
	std::pair<dataset<DVector<double>>, dataset<DVector<double>>> train_test_z = trainTestSplit(*model_.z(), f);

	// set train datasets
	m.setObservations(*train_test_z.first.first, train_test_z.first.second);
	m.setCovariates(*train_test_W.first.first);
	m.setLambda(lambda[0]);  // set current lambda
	m.smooth(); // fit the model

	// compute predicted values
	DVector<double> predictions;
	predictions.resize(train_test_z.second.first->rows());
	for(std::size_t i = 0; i < predictions.rows(); ++i){
	  predictions[i] = m.predict(train_test_W.second.first->row(i), train_test_W.second.second[i]);
	}
	// evaluate RMSE
	rmse += RMSE(*train_test_z.second.first, predictions);
      }
      // record the average RMSE computed among the K_ runs
      RMSEs.push_back(rmse/K_);
    }

    // return optimal lambda according to KFoldCV
    std::vector<double>::iterator minRMSE = std::min_element(RMSEs.begin(), RMSEs.end());
    return lambdas[std::distance(RMSEs.begin(), minRMSE)];
  }

  
};

#endif // __K_FOLD_CV__

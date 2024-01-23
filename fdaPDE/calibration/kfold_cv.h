// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __K_FOLD_CV__
#define __K_FOLD_CV__

#include <fdaPDE/utils.h>
#include <fdaPDE/linear_algebra.h>

#include <algorithm>
#include <vector>
using fdapde::core::BlockFrame;
using fdapde::core::BinaryVector;
using fdapde::Dynamic;

#include "calibration_base.h"

namespace fdapde {
namespace calibration {
  
// general implementation of KFold Cross Validation
class KCV : public CalibratorBase<KCV> {
   private:
    // algorithm's parameters
    std::size_t K_;      // number of folds
    std::size_t seed_;   // seed used for BlockFrame shuffling
    bool shuffle_;       // whether to shuffle data before splitting into folds

    DVector<double> avg_scores_;   // mean CV score for each \lambda
    DVector<double> std_scores_;   // CV score standard deviation for each \lambda
    DMatrix<double> scores_;       // matrix of CV scores
    DVector<double> optimum_;      // optimal smoothing parameter

    // produces a bit_mask identifying the i-th fold to be supplied to the model via mask_obs
    using TrainTestPartition = std::pair<BinaryVector<Dynamic>, BinaryVector<Dynamic>>;
    TrainTestPartition split(const BlockFrame<double, int>& data, std::size_t i) {
        std::size_t n = data.rows();          // number of data points
        std::size_t m = std::floor(n / K_);   // number of data per fold

        // create test-train index sets, as function of test fold i
        BinaryVector<Dynamic> test_mask(n);
	BinaryVector<Dynamic> train_mask(n);
        for (std::size_t j = 0; j < n; ++j) {
            if (j >= m * i && j < m * (i + 1)) {
                test_mask.set(j);
            } else {
                train_mask.set(j);
            }
        }
        return std::make_pair(train_mask, test_mask);
    }
   public:
    // constructor
    KCV() = default;
    KCV(std::size_t K, bool shuffle = true) : K_(K), seed_(std::random_device()()), shuffle_(shuffle) {};
    KCV(std::size_t K, int seed, bool shuffle = true) :
        K_(K), seed_((seed == fdapde::random_seed) ? std::random_device()() : seed), shuffle_(shuffle) {};

    // selects best smoothing parameter of model according to a K-fold cross validation strategy
    template <typename ModelType, typename ScoreType>
    DVector<double> fit(ModelType& model, const std::vector<DVector<double>>& lambdas, ScoreType cv_score) {
        // reserve space for CV scores
        scores_.resize(K_, lambdas.size());
        if (shuffle_) {   // perform a first shuffling of the data if required
            model.set_data(model.data().shuffle(seed_));
	}
        // cycle over all tuning parameters
        for (std::size_t j = 0; j < lambdas.size(); ++j) {
            for (std::size_t fold = 0; fold < K_; ++fold) {   // fixed a tuning parameter, cycle over all data splits
                // compute train-test partition and evaluate CV score
                TrainTestPartition partition_mask = split(model.data(), fold);
                scores_.coeffRef(fold, j) = cv_score(lambdas[j], partition_mask.first, partition_mask.second);
            }
        }
        // reserve space for storing results
        avg_scores_ = DVector<double>::Zero(lambdas.size());
        std_scores_ = DVector<double>::Zero(lambdas.size());
        for (std::size_t j = 0; j < lambdas.size(); ++j) {
            // record the average score and its standard deviation computed among the K_ runs
            double avg_score = 0;
            for (std::size_t i = 0; i < K_; ++i) avg_score += scores_(i, j);
            avg_score /= K_;
            double std_score = 0;
            for (std::size_t i = 0; i < K_; ++i) std_score += std::pow(scores_(i, j) - avg_score, 2);
            std_score = std::sqrt(std_score / (K_ - 1));   // use unbiased sample estimator for standard deviation
            // store results
            avg_scores_[j] = avg_score;
            std_scores_[j] = std_score;
        }

        // store optimal lambda according to given metric F
        Eigen::Index opt_score;
        avg_scores_.minCoeff(&opt_score);
        optimum_ = lambdas[opt_score];
	return optimum_;
    }

    // getters
    const DVector<double>& avg_scores() const { return avg_scores_; }   // mean CV score vector
    const DVector<double>& std_scores() const { return std_scores_; }   // CV score standard deviation vector
    const DMatrix<double>& scores() const { return scores_; }           // CV scores
    const DVector<double>& optimum() const { return optimum_; }         // optimal tuning parameter
    // setters
    void set_n_folds(std::size_t K) { K_ = K; }
};

}   // namespace calibration
}   // namespace fdapde

#endif   // __K_FOLD_CV__

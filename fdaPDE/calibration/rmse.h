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

#ifndef __RMSE_H__
#define __RMSE_H__

#include <fdaPDE/utils.h>
#include <fdaPDE/linear_algebra.h>

#include "../models/regression/regression_type_erasure.h"
using fdapde::core::BlockFrame;
using fdapde::core::BinaryVector;
using fdapde::models::RegressionView;

namespace fdapde {
namespace calibration {

// functor implementing the Root Mean Squared Error (RMSE)
class RMSE {
   private:
    RegressionView<void> model_;
    double alpha_; 
    double eps_; 
   public:
    RMSE() = default;
    template <typename ModelType> RMSE(const ModelType& model, double alpha, double eps) : model_(model), alpha_(alpha), eps_(eps) {};
    // template <typename ModelType> RMSE(const ModelType& model) : model_(model) {};
    template <typename ModelType> void set_model(const ModelType& model) { model_ = model; }

    double operator()(
      const DVector<double>& lambda, const BinaryVector<fdapde::Dynamic>& train_mask,
      const BinaryVector<fdapde::Dynamic>& test_mask) {
        model_.set_lambda(lambda);
        // fit model on train set
        model_.set_mask(test_mask);   // discard test set from training phase
        model_.init();
        model_.solve();

        // compute RMSE over all testing observation which are not missing
        BinaryVector<fdapde::Dynamic> rmse_mask = test_mask & ~model_.nan_mask();
        // RMSE evaluation
        double rmse = 0;
        std::size_t n = 0;   // cardinality of test set

        for (std::size_t i = 0, sz = rmse_mask.size(); i < sz; ++i) {
            if (rmse_mask[i]) {                                    // not a missing value
                double hat_y = model_.Psi().row(i) * model_.f();   // non-parametric field evaluation at i-th location
                if (model_.has_covariates()) { hat_y += model_.X().row(i) * model_.beta(); }

                    // rmse += std::pow(model_.y()(i, 0) - hat_y, 2);

                    // M 
                    rmse += pinball(model_.y()(i, 0) - hat_y, eps_);  

                n++;
            }
        }

        // return std::sqrt(rmse / n);   // \sqrt{\frac{norm(y - \hat y)^2/}{n}}

        return rmse / n;   // M ATT: è la pinball, non rmse
   
    }


    // smoothed pinball version
    double pinball(double x, double eps) const {   // quantile check function
        return (alpha_ - 1) * x + eps * fdapde::log1pexp(x / eps);
    };
};



}   // namespace calibration
}   // namespace fdapde

#endif   // __RMSE_H__

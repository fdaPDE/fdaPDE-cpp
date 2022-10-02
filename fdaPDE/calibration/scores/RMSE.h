#ifndef __RMSE_H__
#define __RMSE_H__

#include <vector>
#include <memory>
#include "../core/utils/Symbols.h"

// functor implementing Root Mean Squared Error performance measure
struct RMSE {

  // compute RMSE: \sqrt{\frac{norm(z - \hat z)^2/}{n}} 
  template <typename M>
  double operator()(const M& model, std::shared_ptr<DVector<double>> z_test,
		    std::shared_ptr<DMatrix<double>> W_test, const std::vector<std::size_t>& obs_indices) const {
    // compute predicted values \hat z
    DVector<double> z_hat(z_test->rows());
    for(std::size_t i = 0; i < z_hat.rows(); ++i){
      z_hat[i] = model.predict(W_test->row(i), obs_indices[i]);
    }

    return std::sqrt((*z_test - z_hat).squaredNorm()/z_test->rows());
  }
  
};

#endif // __RMSE_H__

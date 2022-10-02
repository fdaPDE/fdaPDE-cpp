#ifndef __R2_H__
#define __R2_H__

#include <vector>
#include <memory>
#include "../core/utils/Symbols.h"

// functor implementing the R^2 statistics (coefficient of determination)
struct R2 {

  // compute R^2: 1 - \frac{SSres}{SStot}
  // SSres = \sum_{i=1}^n {z_i - \hat z_i}^2
  // SStot = \sum_{i=1}^n {z_i - \frac{1}{n}*\sum_{j=1}^n z_j}
  template <typename M>
  double operator()(const M& model, std::shared_ptr<DVector<double>> z_test,
		    std::shared_ptr<DMatrix<double>> W_test, const std::vector<std::size_t>& obs_indices) const {
    // compute predicted values \hat z
    DVector<double> z_hat(z_test->rows());
    for(std::size_t i = 0; i < z_hat.rows(); ++i){
      z_hat[i] = model.predict(W_test->row(i), obs_indices[i]);
    }
    // compute average of observed data
    double z_avg = 0;
    for(std::size_t i = 0; i < z_hat.rows(); ++i) z_avg += (*z_test)[i];
    z_avg /= z_test->rows();

    // compute SSres
    double SSres = (*z_test - z_hat).squaredNorm();
    // compute SStot
    double SStot = 0;
    for(std::size_t i = 0; i < z_hat.rows(); ++i) SStot += std::pow((*z_test)[i] - z_avg, 2);

    // need to return a value such that the lower the better
    return SSres/SStot - 1;
  }
  
};

#endif // __R2_H__

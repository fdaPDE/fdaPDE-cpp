#ifndef __MAPE_H__
#define __MAPE_H__

#include <vector>
#include <memory>
#include "../core/utils/Symbols.h"

// functor implementing the Mean Absolute Percentage Error evaluation metric
class MAPE {
private:
  double epsilon = 1e-20;
public:

  // compute MAPE: \frac{1}{n} * \sum_{i=1}^n \frac{|z_i - \hat z_i|}{\max{\epsilon, |z_i|}}
  // with \epsilon a very small yet strictly positive value to avoid undefined results when z_i is zero
  template <typename M>
  double operator()(const M& model, std::shared_ptr<DVector<double>> z_test,
		    std::shared_ptr<DMatrix<double>> W_test, const std::vector<std::size_t>& obs_indices) const {
    // compute predicted values \hat z
    DVector<double> z_hat(z_test->rows());
    for(std::size_t i = 0; i < z_hat.rows(); ++i){
      z_hat[i] = model.predict(W_test->row(i), obs_indices[i]);
    }
    
    // compute sum of percentage absolute errors
    double abs_err = 0;
    for(std::size_t i = 0; i < z_hat.rows(); ++i)
      // \frac{|z_i - \hat z_i|}{\max{\epsilon, |z_i|}}
      abs_err += std::abs((*z_test)[i] - z_hat[i])/std::max(epsilon, std::abs((*z_test)[i]));
    
    // need to return a value such that the lower the better
    return abs_err/z_test->rows();
  }
  
};

#endif // __MAPE_H__

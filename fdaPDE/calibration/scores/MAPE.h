#ifndef __MAPE_H__
#define __MAPE_H__

#include <vector>
#include <memory>
#include "../core/utils/Symbols.h"
#include "../core/utils/DataStructures/BlockFrame.h"

namespace fdaPDE{
namespace calibration{

  // functor implementing the Mean Absolute Percentage Error evaluation metric
  class MAPE {
  private:
    double epsilon = 1e-20;
  public:

    // compute MAPE: \frac{1}{n} * \sum_{i=1}^n \frac{|z_i - \hat z_i|}{\max{\epsilon, |z_i|}}
    // with \epsilon a very small yet strictly positive value to avoid undefined results when z_i is zero
    template <typename M>
    double operator()(const M& model, const BlockFrame<double, int>& test) const {
      // compute predicted values \hat z
      const DMatrix<double>& z_test = test.get<double>(STAT_MODEL_Z_BLK);
      std::size_t n = z_test.rows();

      DVector<double> z_hat(n);
      for(std::size_t i = 0; i < n; ++i){
	z_hat[i] = model.predict(test.get<double>(STAT_MODEL_W_BLK).row(i), test.get<int>(STAT_MODEL_I_BLK)(i,0));
      }
    
      // compute sum of percentage absolute errors
      double abs_err = 0;
      for(std::size_t i = 0; i < n; ++i){
	// \frac{|z_i - \hat z_i|}{\max{\epsilon, |z_i|}}
	abs_err += std::abs(z_test(i,0) - z_hat[i])/std::max(epsilon, std::abs(z_test(i,0)));
      }
      return abs_err/n;
    }  
  };
}}

#endif // __MAPE_H__

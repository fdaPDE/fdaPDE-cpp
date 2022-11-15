#ifndef __RMSE_H__
#define __RMSE_H__

#include <vector>
#include <memory>
#include "../core/utils/Symbols.h"
#include "../core/utils/DataStructures/BlockFrame.h"

namespace fdaPDE{
namespace calibration{

  // functor implementing Root Mean Squared Error performance measure
  struct RMSE {

    // compute RMSE: \sqrt{\frac{norm(y - \hat y)^2/}{n}} 
    template <typename M>
    double operator()(const M& model, const BlockFrame<double, int>& test) const {
      // compute predicted values \hat y
      const DMatrix<double>& y_test = test.get<double>(STAT_MODEL_Y_BLK);
      std::size_t n = y_test.rows();

      DVector<double> y_hat(n);
      for(std::size_t i = 0; i < n; ++i){
	y_hat[i] = model.predict(test.get<double>(STAT_MODEL_W_BLK).row(i), test.get<int>(STAT_MODEL_I_BLK)(i,0));
      }
      // \sqrt{\frac{norm(y - \hat y)^2/}{n}} 
      return std::sqrt((y_test - y_hat).squaredNorm()/n);
    }
  };
}}
  
#endif // __RMSE_H__

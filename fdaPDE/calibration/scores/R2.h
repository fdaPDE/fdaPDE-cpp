#ifndef __R2_H__
#define __R2_H__

#include <vector>
#include <memory>
#include "../core/utils/Symbols.h"
#include "../core/utils/DataStructures/BlockFrame.h"

namespace fdaPDE{
namespace calibration{

  // functor implementing the R^2 statistics (coefficient of determination)
  struct R2 {

    // compute R^2: 1 - \frac{SSres}{SStot}
    // SSres = \sum_{i=1}^n {z_i - \hat z_i}^2
    // SStot = \sum_{i=1}^n {z_i - \frac{1}{n}*\sum_{j=1}^n z_j}^2
    template <typename M>
    double operator()(const M& model, const BlockFrame<double, int>& test) const {
      // compute predicted values \hat z
      const DMatrix<double>& z_test = test.get<double>(STAT_MODEL_Z_BLK);
      std::size_t n = z_test.rows();
    
      DVector<double> z_hat(n);
      for(std::size_t i = 0; i < n; ++i){
	z_hat[i] = model.predict(test.get<double>(STAT_MODEL_W_BLK).row(i), test.get<int>(STAT_MODEL_I_BLK)(i,0));
      }
      // compute average of observed data
      double z_avg = 0;
      for(std::size_t i = 0; i < n; ++i) z_avg += z_test(i,0);
      z_avg /= n;

      // compute SSres
      double SSres = (z_test - z_hat).squaredNorm();
      // compute SStot
      double SStot = 0;
      for(std::size_t i = 0; i < n; ++i) SStot += std::pow(z_test(i,0) - z_avg, 2);

      // need to return a value such that the lower the better
      return SSres/SStot - 1;
    }  
  };
}}

#endif // __R2_H__

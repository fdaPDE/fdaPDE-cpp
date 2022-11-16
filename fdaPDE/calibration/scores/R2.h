#ifndef __R2_H__
#define __R2_H__

#include <vector>
#include <memory>
#include "../../core/utils/Symbols.h"
#include "../../core/utils/DataStructures/BlockFrame.h"

namespace fdaPDE{
namespace calibration{

  // functor implementing the R^2 statistics (coefficient of determination)
  struct R2 {

    // compute R^2: 1 - \frac{SSres}{SStot}
    // SSres = \sum_{i=1}^n {y_i - \hat y_i}^2
    // SStot = \sum_{i=1}^n {y_i - \frac{1}{n}*\sum_{j=1}^n y_j}^2
    template <typename M>
    double operator()(const M& model, const BlockFrame<double, int>& test) const {
      // compute predicted values \hat y
      const DMatrix<double>& y_test = test.get<double>(OBSERVATIONS_BLK);
      std::size_t n = y_test.rows();
    
      DVector<double> y_hat(n);
      for(std::size_t i = 0; i < n; ++i){
	y_hat[i] = model.predict(test.get<double>(DESIGN_MATRIX_BLK).row(i), test.get<int>(INDEXES_BLK)(i,0));
      }
      // compute average of observed data
      double y_avg = 0;
      for(std::size_t i = 0; i < n; ++i) y_avg += y_test(i,0);
      y_avg /= n;

      // compute SSres
      double SSres = (y_test - y_hat).squaredNorm();
      // compute SStot
      double SStot = 0;
      for(std::size_t i = 0; i < n; ++i) SStot += std::pow(y_test(i,0) - y_avg, 2);

      // need to return a value such that the lower the better
      return SSres/SStot - 1;
    }  
  };
}}

#endif // __R2_H__

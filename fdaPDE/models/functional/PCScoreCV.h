#ifndef __R2_H__
#define __R2_H__

#include <vector>
#include <memory>
#include "../../core/utils/Symbols.h"
#include "../../core/utils/DataStructures/BlockFrame.h"

namespace fdaPDE{
namespace models{

  // functor implementing the CV statistics for selection of \lambda parameters in PCA
  struct PCScoreCV {

    // compute \norm_F(Y - s*f^T) where \norm_F{} is the Frobenius norm
    // s are the predicted scores while f the computed loadings
    template <typename M>
    double operator()(const M& m, const BlockFrame<double, int>& test) const {
      // compute scores from estimated loadings
      DVector<double> scores = test.get<double>(OBSERVATIONS_BLK)*m.f(); // s = Y*f
      double norm = m.f().squaredNorm() + (m.f().transpose()*m.R0()*m.f()).coeff(0,0); // \norm{f} + P(f)
      // normalize scores
      scores = scores/norm; // s/(\norm{f} + P(f))

      // return CV statistic \norm_F(Y - s*f^T)
      return (test.get<double>(OBSERVATIONS_BLK) - scores*m.f().transpose()).squaredNorm();
    }  
  };
  
}}

#endif // __R2_H__

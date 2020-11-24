#ifndef __K_FOLD_CV_L2_ERROR_H__
#define __K_FOLD_CV_L2_ERROR_H__

// This file implements the cross validation error based on the L2 norm useful for the Density Estimation problem

//! @brief A class to compute the L2 error during cross-validation.
template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class KfoldCV_L2_error{
  private:
    // A member to acess  data problem methods
    const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dataProblem_;

  public:
    //! A constructor.
    KfoldCV_L2_error(const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dp): dataProblem_(dp) {};
    //! A call operator to compute the L2 error.
    Real operator()(const SpMat& Psi, const VectorXr& g);

};

#include "K_Fold_CV_L2_Error_imp.h"

#endif

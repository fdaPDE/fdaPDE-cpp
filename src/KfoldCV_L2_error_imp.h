#ifndef _KFOLDCV_L2_ERROR_IMP_H_
#define _KFOLDCV_L2_ERROR_IMP_H_

template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
Real KfoldCV_L2_error<Integrator_noPoly, ORDER,mydim,ndim>::operator() (const SpMat& Psi, const VectorXr& g)
{
  Real integral = dataProblem_.FEintegrate_exponential(2.*g);
  Real test = (Psi*g).array().exp().sum();

  return (integral - 2./Psi.rows() *test);
}

#endif

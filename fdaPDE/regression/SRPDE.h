#ifndef __SRPDE_H__
#define __SRPDE_H__

#include <memory>
// CORE imports
#include "../core/utils/Symbols.h"
#include "../core/FEM/PDE.h"
using fdaPDE::core::FEM::PDE;
#include "../core/NLA/SparseBlockMatrix.h"
using fdaPDE::core::NLA::SparseBlockMatrix;
#include "../core/NLA/SMW.h"
using fdaPDE::core::NLA::SMW;
// calibration module imports
#include "../calibration/iGCV.h"
// regression module imports
#include "iStatModel.h"
#include "Internals.h"

template <unsigned int M, unsigned int N, unsigned int K, typename E, typename B>
class SRPDE : public iStatModel<M,N,K,E,B>, public iGCV {
public:
  // constructor
  SRPDE() = default;
  SRPDE(const PDE<M,N,K,E,B>& pde, double lambda)
    : iStatModel<M,N,K,E,B>(pde, lambda) {};
  IMPORT_STAT_MODEL_SYMBOLS(M,N,K,E,B);

  // iStatModel interface implementation
  virtual void smooth();                  // finds a solution to the smoothing problem
  virtual DVector<double> fitted() const; // computes the fitted values \hat z
  virtual double predict(const DVector<double>& covs, const std::size_t loc) const;

  // iGCV interface implementation
  virtual std::shared_ptr<DMatrix<double>> T(); // T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
  
  virtual ~SRPDE() = default;
};

// template argument deduction guide
template <unsigned int M, unsigned int N, unsigned int K, typename E, typename B>
SRPDE(const PDE<M,N,K,E,B>& pde_, double lambda_) -> SRPDE<M,N,K,E,B>;

#include "SRPDE.tpp"

#endif // __SRPDE_H__

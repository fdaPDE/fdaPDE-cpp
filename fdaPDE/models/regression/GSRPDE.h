#ifndef __GSRPDE_H__
#define __GSRPDE_H__

namespace fdaPDE {
namespace models {
namespace regression {

  template <unsigned int M, unsigned int N, unsigned int K, typename E, typename F, typename B, typename I, typename S>
  class GSRPDE : public iRegressionModel<GSRPDE<M,N,K,E,F,B,I,S>>, public iGCV {
  public:
    // constructor
    GSRPDE() = default;
    GSRPDE(const PDE<M,N,K,E,F,B,I,S>& pde, double lambda)
      : iStatModel<M,N,K,E,F,B,I,S>(pde, lambda) {};
    IMPORT_REGRESSION_MODEL_SYMBOLS(GSRPDE<M,N,K,E,F,B,I,S>);

    // iStatModel interface implementation
    virtual void solve();                   // finds a solution to the smoothing problem

    // iRegressionModel interface implementation
    virtual DMatrix<double> fitted() const; // computes the fitted values \hat z
    virtual double predict(const DVector<double>& covs, const std::size_t loc) const;

    // iGCV interface implementation
    virtual std::shared_ptr<DMatrix<double>> T(); // T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
  
    virtual ~SRPDE() = default;
  };

  
}}}

#endif // __GSRPDE_H__

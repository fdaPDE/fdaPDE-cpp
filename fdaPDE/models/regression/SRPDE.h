#ifndef __SRPDE_H__
#define __SRPDE_H__

#include <memory>
#include <type_traits>
#include "../../core/utils/Symbols.h"
#include "../../core/FEM/PDE.h"
using fdaPDE::core::FEM::PDEBase;
#include "../ModelBase.h"
#include "../ModelMacros.h"
#include "../../core/NLA/SparseBlockMatrix.h"
using fdaPDE::core::NLA::SparseBlockMatrix;
#include "../../core/NLA/SMW.h"
using fdaPDE::core::NLA::SMW;
#include "../../calibration/iGCV.h"
using fdaPDE::calibration::iGCV;
#include "../SamplingDesign.h"
#include "RegressionBase.h"
using fdaPDE::models::RegressionBase;

namespace fdaPDE{
namespace models{
  
  template <typename PDE, typename SamplingDesign>
  class SRPDE : public RegressionBase<SRPDE<PDE, SamplingDesign>>, public iGCV {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
  private:
    typedef RegressionBase<SRPDE<PDE, SamplingDesign>> Base;   
    SparseBlockMatrix<double,2,2> A_{}; // system matrix of non-parametric problem (2N x 2N matrix)
    fdaPDE::SparseLU<SpMatrix<double>> invA_; // factorization of matrix A
    DVector<double> b_{}; // right hand side of problem's linear system (1 x 2N vector)
  public:
    IMPORT_REGRESSION_SYMBOLS;
    using Base::lambdaS; // smoothing parameter in space
    // constructor
    SRPDE() = default;
    SRPDE(const PDE& pde) : Base(pde) {};
    
    // ModelBase implementation
    void init_model();    // update model object in case of **structural** changes in its definition
    virtual void solve(); // finds a solution to the smoothing problem
    
    // iGCV interface implementation
    virtual const DMatrix<double>& T(); // T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
    virtual const DMatrix<double>& Q(); // Q = W(I - H) = W - W*X*(X^T*W*X)^{-1}X^T*W
    // returns the euclidian norm of y - \hat y
    virtual double norm(const DMatrix<double>& obs, const DMatrix<double>& fitted) const;

    // getters
    const SparseBlockMatrix<double,2,2>& A() const { return A_; }
    const fdaPDE::SparseLU<SpMatrix<double>>& invA() const { return invA_; }
    
    virtual ~SRPDE() = default;
  };
  template <typename PDE_, typename SamplingDesign_>
  struct model_traits<SRPDE<PDE_, SamplingDesign_>> {
    typedef PDE_ PDE;
    typedef SpaceOnly regularization;
    typedef SamplingDesign_ sampling;
    typedef MonolithicSolver solver;
    static constexpr int n_lambda = 1;
  };

  // srpde trait
  template <typename Model>
  struct is_srpde { static constexpr bool value = is_instance_of<Model, SRPDE>::value; };
  
#include "SRPDE.tpp"
}}
    
#endif // __SRPDE_H__

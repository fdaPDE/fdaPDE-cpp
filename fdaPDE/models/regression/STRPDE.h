#ifndef __STRPDE_H__
#define __STRPDE_H__

#include <memory>
#include <type_traits>
// CORE imports
#include "../../core/utils/Symbols.h"
#include "../../core/FEM/PDE.h"
using fdaPDE::core::FEM::PDEBase;
#include "../space_time/SplineBasis.h"
#include "../../core/NLA/SparseBlockMatrix.h"
using fdaPDE::core::NLA::SparseBlockMatrix;
#include "../../core/NLA/SMW.h"
using fdaPDE::core::NLA::SMW;
#include "../../core/NLA/KroneckerProduct.h"
using fdaPDE::core::NLA::SparseKroneckerProduct;
using fdaPDE::core::NLA::Kronecker;
#include "../../core/utils/DataStructures/BlockVector.h"
using fdaPDE::BlockVector;
// calibration module imports
#include "../../calibration/iGCV.h"
using fdaPDE::calibration::iGCV;
// regression module imports
#include "RegressionBase.h"
using fdaPDE::models::RegressionBase;
#include "../ModelTraits.h"
using fdaPDE::models::Gaussian;

namespace fdaPDE{
namespace models{

  // base class for STRPDE model
  template <typename PDE, typename TimeRegularization, Sampling SamplingDesign, SolverType Solver> class STRPDE;

  // implementation of STRPDE for separable space-time regularization
  template <typename PDE, Sampling SamplingDesign>
  class STRPDE<PDE, SpaceTimeSeparableTag, SamplingDesign, SolverType::Monolithic>
    : public RegressionBase<STRPDE<PDE, SpaceTimeSeparableTag, SamplingDesign, SolverType::Monolithic>>, public iGCV {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
  private:
    typedef SpaceTimeSeparableTag TimeRegularization;
    typedef RegressionBase<STRPDE<PDE, TimeRegularization, SamplingDesign, SolverType::Monolithic>> Base;
    SpMatrix<double> A_{}; // system matrix of non-parametric problem (2N x 2N matrix)
    fdaPDE::SparseLU<SpMatrix<double>> invA_; // factorization of matrix A
    DVector<double> b_{};  // right hand side of problem's linear system (1 x 2N vector)

    // matrices related to woodbury decomposition
    DMatrix<double> U_{};
    DMatrix<double> V_{};
  public:
    // import commonly defined symbols from base
    IMPORT_REGRESSION_SYMBOLS;
    using Base::lambdaS; // smoothing parameter in space
    using Base::lambdaT; // smoothing parameter in time
    using Base::Pt;      // time penalization matrix: [Pt_]_{ij} = \int_{[0,T]} (\phi_i)_tt*(\phi_j)_tt
    using Base::Rt;      // time mass matrix: [Rt_]_{ij} = \int_{[0,T]} \phi_i*\phi_j
    // constructor
    STRPDE() = default;
    STRPDE(const PDE& pde, const DMatrix<double>& time) : Base(pde, time) {};
    
    // ModelBase interface implementation
    void init_model();    // update model object in case of **structural** changes in its definition
    virtual void solve(); // finds a solution to the smoothing problem

    // iGCV interface implementation
    virtual const DMatrix<double>& T() { // T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
      // compute value of R = R1^T*R0^{-1}*R1, cache for possible reuse
      if(R_.size() == 0){
	invR0_.compute(R0());
	R_ = R1().transpose()*invR0_.solve(R1());
      }
      // compute and store matrix T for possible reuse
      if(!hasCovariates()) // case without covariates, Q is the identity matrix
	T_ = PsiTD()*W()*Psi()   + lambdaS()*R_;
      else // general case with covariates
	T_ = PsiTD()*lmbQ(Psi()) + lambdaS()*R_;
      return T_;
    }; 
    virtual const DMatrix<double>& Q() { // Q = W(I - H) = W - W*X*(X^T*W*X)^{-1}X^T*W
      if(Q_.size() == 0){ // Q is computed on request since not needed in general
	// compute Q = W(I - H) = W - W*X*(X*W*X^T)^{-1}*X^T*W
	Q_ = W()*(DMatrix<double>::Identity(n_obs(), n_obs()) - X()*invXtWX().solve(X().transpose()*W()));
      }
      return Q_;
    };
    // returns the euclidian norm of y - \hat y
    virtual double norm(const DMatrix<double>& op1, const DMatrix<double>& op2) const {
      return (op1 - op2).squaredNorm();
    }

    // getters
    const SpMatrix<double>& A() const { return A_; }
    const fdaPDE::SparseLU<SpMatrix<double>>& invA() const { return invA_; }
    const DMatrix<double>& U() const { return U_; }
    const DMatrix<double>& V() const { return V_; }
    
    virtual ~STRPDE() = default;
  };
  template <typename PDE_, Sampling SamplingDesign>
  struct model_traits<STRPDE<PDE_, SpaceTimeSeparableTag, SamplingDesign, SolverType::Monolithic>> {
    typedef PDE_ PDE;
    typedef SpaceTimeSeparableTag RegularizationType;
    typedef SplineBasis<3> TimeBasis; // use cubic B-splines
    static constexpr Sampling sampling = SamplingDesign;
    static constexpr SolverType solver = SolverType::Monolithic;
    static constexpr int n_lambda = 2;
    typedef Gaussian DistributionType;
  };

  // implementation of STRPDE for parabolic space-time regularization, monolithic solver
  template <typename PDE, Sampling SamplingDesign>
  class STRPDE<PDE, SpaceTimeParabolicTag, SamplingDesign, SolverType::Monolithic>
    : public RegressionBase<STRPDE<PDE, SpaceTimeParabolicTag, SamplingDesign, SolverType::Monolithic>>/*, public iGCV*/ {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
  private:
    typedef SpaceTimeParabolicTag TimeRegularization;
    typedef RegressionBase<STRPDE<PDE, SpaceTimeParabolicTag, SamplingDesign, SolverType::Monolithic>> Base;
    SpMatrix<double> A_{}; // system matrix of non-parametric problem (2N x 2N matrix)
    fdaPDE::SparseLU<SpMatrix<double>> invA_; // factorization of matrix A
    DVector<double> b_{};  // right hand side of problem's linear system (1 x 2N vector)

    // matrices related to woodbury decomposition
    DMatrix<double> U_{};
    DMatrix<double> V_{};    
  public:
    // import commonly defined symbols from base
    IMPORT_REGRESSION_SYMBOLS;
    using Base::lambdaS; // smoothing parameter in space
    using Base::lambdaT; // smoothing parameter in time
    using Base::L;       // [L]_{ii} = 1/DeltaT for i \in {1 ... m} and [L]_{i,i-1} = -1/DeltaT for i \in {1 ... m-1}
    using Base::n_time;  // number of time instants m defined over [0,T]
    using Base::s;       // initial condition
    // constructor
    STRPDE() = default;
    STRPDE(const PDE& pde, const DMatrix<double>& time) : Base(pde, time) {};
    
    // ModelBase interface implementation
    void init_model();    // update model object in case of **structural** changes in its definition
    virtual void solve(); // finds a solution to the smoothing problem
    
    // iGCV interface implementation
    // virtual const DMatrix<double>& T(); // T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
    // virtual const DMatrix<double>& Q(); // Q = W(I - H) = W - W*X*(X^T*W*X)^{-1}X^T*W
    // returns the euclidian norm of y - \hat y
    // virtual double norm(const DMatrix<double>& obs, const DMatrix<double>& fitted) const;
    
    virtual ~STRPDE() = default;
  };
  template <typename PDE_, Sampling SamplingDesign>
  struct model_traits<STRPDE<PDE_, SpaceTimeParabolicTag, SamplingDesign, SolverType::Monolithic>> {
    typedef PDE_ PDE;
    typedef SpaceTimeParabolicTag RegularizationType;
    static constexpr Sampling sampling = SamplingDesign;
    static constexpr SolverType solver = SolverType::Monolithic;
    static constexpr int n_lambda = 2;
    typedef Gaussian DistributionType;
  };
  
  // implementation of STRPDE for parabolic space-time regularization, monolithic solver
  template <typename PDE, Sampling SamplingDesign>
  class STRPDE<PDE, SpaceTimeParabolicTag, SamplingDesign, SolverType::Iterative>
    : public RegressionBase<STRPDE<PDE, SpaceTimeParabolicTag, SamplingDesign, SolverType::Iterative>>/*, public iGCV*/ {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
  private:
    typedef SpaceTimeParabolicTag TimeRegularization;
    typedef RegressionBase<STRPDE<PDE, SpaceTimeParabolicTag, SamplingDesign, SolverType::Iterative>> Base;
    SpMatrix<double> A_{}; // system matrix of non-parametric problem (2N x 2N matrix)
    fdaPDE::SparseLU<SpMatrix<double>> invA_; // factorization of matrix A
    DVector<double> b_{};  // right hand side of problem's linear system (1 x 2N vector)
    
    // the functional minimized by the iterative scheme
    // J(f,g) = \sum_{k=1}^m (z^k - \Psi*f^k)^T*(z^k - \Psi*f^k) + \lambda_S*(g^k)^T*(g^k)
    double J(const DMatrix<double>& f, const DMatrix<double>& g) const;
    // internal solve routine used by the iterative method
    void solve(std::size_t t, BlockVector<double>& f_new, BlockVector<double>& g_new) const;
  public:
    // import commonly defined symbols from base
    IMPORT_REGRESSION_SYMBOLS;
    using Base::lambdaS;   // smoothing parameter in space
    using Base::lambdaT;   // smoothing parameter in time
    using Base::n_time;    // number of time instants m defined over [0,T]
    using Base::tol_;      // tolerance on std::abs((Jnew - Jold)/Jnew)
    using Base::max_iter_; // maximum number of allowed iterations before forced stop
    using Base::DeltaT;    // distance between two time instants
    // constructor
    STRPDE() = default;
    STRPDE(const PDE& pde, const DMatrix<double>& time) : Base(pde, time) {};
    
    // ModelBase interface implementation
    void init_model() { return; }
    virtual void solve(); // finds a solution to the smoothing problem

    // iGCV interface implementation
    // virtual const DMatrix<double>& T(); // T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
    // virtual const DMatrix<double>& Q(); // Q = W(I - H) = W - W*X*(X^T*W*X)^{-1}X^T*W
    // returns the euclidian norm of y - \hat y
    // virtual double norm(const DMatrix<double>& obs, const DMatrix<double>& fitted) const;
    
    virtual ~STRPDE() = default;
  };
  template <typename PDE_, Sampling SamplingDesign>
  struct model_traits<STRPDE<PDE_, SpaceTimeParabolicTag, SamplingDesign, SolverType::Iterative>> {
    typedef PDE_ PDE;
    typedef SpaceTimeParabolicTag RegularizationType;
    static constexpr Sampling sampling = SamplingDesign;
    static constexpr SolverType solver = SolverType::Iterative;
    static constexpr int n_lambda = 2;
    typedef Gaussian DistributionType;
  };
  
#include "STRPDE.tpp"
}}
#endif // __STRPDE_H__

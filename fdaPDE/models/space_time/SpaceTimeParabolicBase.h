#ifndef __SPACE_TIME_PARABOLIC_BASE_H__
#define __SPACE_TIME_PARABOLIC_BASE_H__

#include "../../core/utils/Symbols.h"
#include "../../core/NLA/KroneckerProduct.h"
#include "../../core/FEM/operators/BilinearFormTraits.h"
using fdaPDE::core::FEM::is_parabolic;
#include "../../core/FEM/PDE.h"
using fdaPDE::core::NLA::SparseKroneckerProduct;
using fdaPDE::core::NLA::Kronecker;
#include "../SpaceTimeBase.h"
using fdaPDE::models::SpaceTimeBase;

namespace fdaPDE{
namespace models{

  // base class for parabolic regularization solved using either a moholitic or iterative solution strategy
  template <typename Model, typename Solver> class SpaceTimeParabolicBase;

  // base class for parabolic regularization, monholitic solver
  template <typename Model>
  class SpaceTimeParabolicBase<Model, MonolithicSolver> : public SpaceTimeBase<Model> {
    static_assert(is_parabolic<typename model_traits<Model>::PDE::BilinearFormType>::value,
		  "you have asked for parabolic regularization but using a non-parabolic differential operator");
  private:
    // let m the number of time points
    DMatrix<double> s_;    // N x 1 initial condition vector
    DMatrix<double> u_;    // discretized forcing [1/DeltaT * (u_1 + R_0*s) \ldots u_n]
    SpMatrix<double> Im_;  // m x m sparse identity matrix (assembled once and cached for reuse)
    SpMatrix<double> L_;   // m x m matrix associated with the derivation in time
    double DeltaT_;        // time step (assumes equidistant points in time)
    SpMatrix<double> R0_;  // Im \kron R0 (R0: discretization of the identity operator)
    SpMatrix<double> R1_;  // Im \kron R1 (R1: discretization of the differential operator L in the regularizing PDE)

    // quantites related to discretization of the penalty term
    SpMatrix<double> penT_; // discretization of the time derivative: L \kron R0
    fdaPDE::SparseLU<SpMatrix<double>> invR0_; // factorization of Im \kron R0
    SpMatrix<double> pen_;  // discretized regularizing term: (Im \kron R1 + L \kron R0)^T*(I_m \kron R0)^{-1}*(Im \kron R1 + L \kron R0)
  public:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    typedef typename model_traits<Model>::regularization TimeRegularization; // regularization in time
    typedef SpaceTimeBase<Model> Base;
    using Base::pde_;  // regularizing term in space
    using Base::model; // underlying model object
    using Base::time_; // time interval [0,T]
    using Base::lambdaS; // smoothing parameter in space
    using Base::lambdaT; // smoothing parameter in time
    
    // constructor
    SpaceTimeParabolicBase() = default;
    SpaceTimeParabolicBase(const PDE& pde, const DVector<double>& time)
      : SpaceTimeBase<Model>(pde, time) {}
    // init data structure related to parabolic regularization
    void init_regularization() {
      std::size_t m_ = time_.rows(); // number of time points
      DeltaT_ = time_[1] - time_[0]; // time step (assuming equidistant points)
      
      // assemble once the m x m identity matrix and cache for fast access
      Im_.resize(m_,m_);
      Im_.setIdentity();

      // assemble matrix associated with derivation in time L_
      // [L_]_{ii} = 1/DeltaT for i \in {1 ... m} and [L_]_{i,i-1} = -1/DeltaT for i \in {1 ... m-1} 
      std::vector<fdaPDE::Triplet<double>> tripletList;
      tripletList.reserve(2*m_);
      // start assembly loop
      double invDeltaT = 1.0/DeltaT_;
      tripletList.emplace_back(0,0, invDeltaT);
      for(std::size_t i = 1; i < m_; ++i){
	tripletList.emplace_back(i,   i,  invDeltaT);
	tripletList.emplace_back(i, i-1, -invDeltaT);
      }
      // finalize construction
      L_.resize(m_,m_);
      L_.setFromTriplets(tripletList.begin(), tripletList.end());
      L_.makeCompressed();
      // compute tensorized matrices
      R0_ = Kronecker(Im_, pde_->R0());
      R1_ = Kronecker(Im_, pde_->R1());
    } 
    
    // getters
    const SpMatrix<double>& R0()  const { return R0_; }
    const SpMatrix<double>& R1()  const { return R1_; }
    // matrices proper of separable regularization
    const SpMatrix<double>& L() const { return L_; }
    // return discretized force corrected by initial conditions
    const DMatrix<double>& u() {
      if(is_empty(u_)){ // compute once and cache
	u_ = pde_->force();
	// correct first n rows of discretized force as (u_1 + R0*s/DeltaT)
	u_.block(0,0, model().n_basis(),1) += (1.0/DeltaT_)*(pde_->R0()*s_);
      }
      return u_;
    }
    const DMatrix<double>& s() { return s_; } // initial condition
    double DeltaT() const { return DeltaT_; }

    // computes and cache matrices (Im \kron R0)^{-1} and L \kron R0
    // returns an expression encoding \lambdaS*((Im \kron R1 + \lambda_T*(L \kron R0))^T*(I_m \kron R0)^{-1}*(Im \kron R1 + \lambda_T*(L \kron R0)))
    auto pen() {
      if(is_empty(pen_)){ // compute once and cache result
	invR0_.compute(R0());
	penT_ = Kronecker(L_, pde_->R0());	
      }
      return lambdaS()*(R1() + lambdaT()*penT_).transpose()*invR0_.solve(R1() + lambdaT()*penT_);
    }
    
    // setters
    void setInitialCondition(const DMatrix<double>& s) { s_ = s; }
    
    // destructor
    virtual ~SpaceTimeParabolicBase() = default;  
  };

  // base class for parabolic regularization, iterative solver
  template <typename Model>
  class SpaceTimeParabolicBase<Model, IterativeSolver> : public SpaceTimeBase<Model> {
    static_assert(is_parabolic<typename model_traits<Model>::PDE::BilinearFormType>::value,
		  "you have asked for parabolic regularization but using a non-parabolic differential operator");
  protected:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    typedef typename model_traits<Model>::regularization TimeRegularization; // regularization in time
    typedef SpaceTimeBase<Model> Base;
    using Base::pde_;  // regularizing term in space
    using Base::model; // underlying model object
    using Base::time_; // time interval [0,T]

    DMatrix<double> s_; // N x 1 initial condition vector
    DMatrix<double> u_; // discretized forcing [1/DeltaT * (u_1 + R_0*s) \ldots u_n]
    double DeltaT_;     // time step (assumes equidistant points in time)

    // quantities related to iterative scheme
    double tol_ = 1e-4;         // tolerance used as stopping criterion
    std::size_t max_iter_ = 50; // maximum number of allowed iterations
  public:
    // constructor
    SpaceTimeParabolicBase() = default;
    SpaceTimeParabolicBase(const PDE& pde, const DVector<double>& time)
      : SpaceTimeBase<Model>(pde, time) {}
    // init data required for iterative solution of parabolic regularization
    void init_regularization() {
      // compute time step (assuming equidistant points)
      DeltaT_ = time_[1] - time_[0];
      u_ = pde_->force(); // compute forcing term
      // correct first n rows of discretized force as (u_1 + R0*s/DeltaT)
      u_.block(0,0, model().n_basis(),1) += (1.0/DeltaT_)*(pde_->R0()*s_);
    }
    
    // getters
    const SpMatrix<double>& R0()  const { return pde_->R0(); }    // mass matrix in space
    const SpMatrix<double>& R1()  const { return pde_->R1(); }    // discretization of differential operator L
    DMatrix<double> u(std::size_t k) const { // discretization of forcing term u at time k
      return u_.block(model().n_basis()*k,0, model().n_basis(),1); }
    double DeltaT() const { return DeltaT_; } 
    const DMatrix<double>& s() const { return s_; }
    
    // setters
    void setInitialCondition(const DMatrix<double>& s) { s_ = s; }
    void setTolerance(double tol) { tol_ = tol; }
    void setMaxIterations(std::size_t max_iter) { max_iter_ = max_iter; }
    
    // destructor
    virtual ~SpaceTimeParabolicBase() = default;  
  };
  
}}

#endif // __SPACE_TIME_PARABOLIC_BASE_H__

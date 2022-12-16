#ifndef __SPACE_TIME_PARABOLIC_BASE_H__
#define __SPACE_TIME_PARABOLIC_BASE_H__

#include "../../core/utils/Symbols.h"
#include "../../core/NLA/KroneckerProduct.h"
#include "../../core/FEM/operators/BilinearFormTraits.h"
using fdaPDE::core::FEM::is_parabolic;
#include <cstddef>
using fdaPDE::core::NLA::SparseKroneckerProduct;
using fdaPDE::core::NLA::Kronecker;
#include "../SpaceTimeBase.h"
using fdaPDE::models::SpaceTimeBase;
#include "TimeAssembler.h"

namespace fdaPDE{
namespace models{

  // base class for parabolic regularization
  template <typename Model>
    class SpaceTimeParabolicBase : public SpaceTimeBase<Model> {
    static_assert(is_parabolic<typename model_traits<Model>::PDE>::value,
		  "you have asked for parabolic smoothing but using a non-parabolic differential operator");
  private:
    // let m the number of time points
    DVector<double> s_;   // N x 1 initial condition vector
    DMatrix<double> u_;   // discretized forcing [1/DeltaT * (u_1 + R_0*s) \ldots u_n]
    SpMatrix<double> Im_; // m x m sparse identity matrix (assembled once and cached for reuse)
    SpMatrix<double> L_;  // m x m matrix associated with the derivation in time
    double DeltaT_;       // time step (assumes equidistant points in time)
  public:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    typedef typename model_traits<Model>::RegularizationType TimeRegularization; // regularization in time
    typedef SpaceTimeBase<Model> Base;
    using Base::pde_;  // regularizing term in space
    using Base::model; // underlying model object
    
    // constructor
    SpaceTimeParabolicBase(const PDE& pde, const DVector<double>& time)
      : SpaceTimeBase<Model>(pde, time) {
      std::size_t m_ = time.rows(); // number of time points
      DeltaT_ = time[1] - time[0]; // time step (assuming equidistant points)
      
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
      for(std::size_t i = 1; i < m_-1; ++i){
	tripletList.emplace_back(i,   i,  invDeltaT);
	tripletList.emplace_back(i, i-1, -invDeltaT);
      }
      tripletList.emplace_back(m_,m_, invDeltaT);
      // finalize construction
      L_.setFromTriplets(tripletList.begin(), tripletList.end());
      L_.makeCompressed();
    }
    // getters
    SparseKroneckerProduct<> R0()  const { return Kronecker(Im_, pde_->R0()); }
    SparseKroneckerProduct<> R1()  const { return Kronecker(Im_, pde_->R1()); }
    SparseKroneckerProduct<> Psi() const { return Kronecker(Im_, model().Psi_); }
    // matrices proper of separable regularization
    const SpMatrix<double>& L() const { return L_; }
    
    // return discretized force corrected by initial conditions
    const DMatrix<double>& u() {
      u_ = pde_->force();
      u_.block(0,0, model().n_basis(),1) /= DeltaT_;
      // in case of initial condition add R0*s/DeltaT to the first block
      if(hasInitialCondition())
	u_.block(0,0, model().n_basis(),1) += (1.0/DeltaT_)*(pde_->R0()*s_);
      
      return u_;
    }

    void setInitialCondition(const DVector<double>& s) { s_ = s; }
    bool hasInitialCondition() const { return s_.size() != 0; }
    
    // destructor
    virtual ~SpaceTimeParabolicBase() = default;  
  };
  
}}

#endif // __SPACE_TIME_PARABOLIC_BASE_H__

#ifndef __SPACE_TIME_SEPARABLE_BASE_H__
#define __SPACE_TIME_SEPARABLE_BASE_H__

#include "../../core/utils/Symbols.h"
#include "../../core/NLA/KroneckerProduct.h"
using fdaPDE::core::NLA::SparseKroneckerProduct;
using fdaPDE::core::NLA::Kronecker;
#include "../SpaceTimeBase.h"
using fdaPDE::models::SpaceTimeBase;
#include "TimeAssembler.h"

namespace fdaPDE{
namespace models{

  // base class for separable regularization
  template <typename Model>
    class SpaceTimeSeparableBase : public SpaceTimeBase<Model> {
  private:
    // let \phi_i the i-th basis function in time
    SpMatrix<double> Rt_;  // mass matrix in time: [Rt_]_{ij} = \int_{[0,T]} \phi_i*\phi_j
    SpMatrix<double> Pt_;  // penalty matrix in time: [Pt_]_{ij} = \int_{[0,T]} (\phi_i)_tt*(\phi_j)_tt
    SpMatrix<double> Phi_; // [Phi_]_{ij} = \phi_i(t_j)
    DVector<double> u_;    // stacked discretized forcing [u_1 \ldots u_n, \ldots, u_1 \ldots u_n]
    
    typedef typename model_traits<Model>::TimeBasis TimeBasis; // basis used for discretization in time
    TimeBasis basis_;
    DMatrix<double> time_locations_; // time instants t_1, ..., t_m. used only if provided from BlockFrame
    SpMatrix<double> pen_; // discretized regularizing term: \lambda_S*((R1^T*R0^{-1}*R1) \kron Rt) + \lambda_T*(R0 \kron Pt)
  public:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    typedef typename model_traits<Model>::RegularizationType TimeRegularization; // regularization in time
    typedef SpaceTimeBase<Model> Base;
    using Base::pde_;    // regularizing term in space Lf = u
    using Base::model;   // underlying model object
    using Base::time_;   // time interval [0,T]
    using Base::lambdaS; // smoothing parameter in space \lambda_S
    using Base::lambdaT; // smoothing parameter in time \lambda_T
    
    // constructor
    SpaceTimeSeparableBase() = default;
    SpaceTimeSeparableBase(const PDE& pde, const DVector<double>& time)
      : SpaceTimeBase<Model>(pde, time), basis_(time) {}
    // init data structure related to separable regularization
    void init_regularization() {
      // compute \Phi = [\Phi]_{ij} = \phi_i(t_j) using the provided basis function
      if(!model().data().hasBlock(TIME_LOCATIONS_BLK)){
	// assume time instants t_1, ..., t_m equal to time nodes
	Phi_ = basis_.eval(time_);
      }else{
	// extract time locations t_1, ..., t_m from BlockFrame and evaluate basis on them
	time_locations_ = model().data().template extract_unique<double>(TIME_LOCATIONS_BLK);
	Phi_ = basis_.eval(time_locations_);
      }
      // discretize operators in time
      TimeAssembler<TimeBasis> assembler(time_);
      Rt_ = assembler.assemble(TimeMass<TimeBasis>());    // mass matrix in time
      Pt_ = assembler.assemble(TimePenalty<TimeBasis>()); // time penalty matrix
    } 
    // getters
    SparseKroneckerProduct<> R0()  const { return Kronecker(Rt_, pde_->R0()); }
    SparseKroneckerProduct<> R1()  const { return Kronecker(Rt_, pde_->R1()); }
    SparseKroneckerProduct<> Psi() const { return Kronecker(Phi_, model().Psi_); }
    // matrices proper of separable regularization
    const SpMatrix<double>& Rt()  const { return Rt_; }
    const SpMatrix<double>& Pt()  const { return Pt_; }
    const SpMatrix<double>& Phi() const { return Phi_; }
    inline std::size_t n_time() const { // overloads Base::n_time()
      if(!model().data().hasBlock(TIME_LOCATIONS_BLK)) return time_.rows();
      return time_locations_.rows();
    } 

    // return stacked version of discretized forcing field
    const DVector<double>& u() {
      // compute result dimensions.
      std::size_t M = basis_.size();
      std::size_t N = Base::n_basis();
      if(u_.size() == 0){ // first time this is computed
	u_.resize(N*M);
	// in separable regularization PDE doesn't depend on time. stack forcing term m times
	for(std::size_t i = 0; i < M; ++i) u_.segment(i*N, N) = pde_->force();
      }
      return u_;
    }

    // computes and returns \lambda_S*((R1^T*R0^{-1}*R1) \kron Rt) + \lambda_T*(R0 \kron Pt)
    const SpMatrix<double>& pen() {
      if(pen_.size() == 0){ // compute once and cache result
	fdaPDE::SparseLU<SpMatrix<double>> invR0_;
	invR0_.compute(pde_->R0());
	SpMatrix<double> Ps_ = pde_->R1().transpose()*invR0_.solve(pde_->R1()); // R1^T*R0^{-1}*R1
	// \lambda_S*((R1^T*R0^{-1}*R1) \kron Rt) + \lambda_T*(R0 \kron Pt)
	pen_ = lambdaS()*Kronecker(Ps_, Rt_) + lambdaT()*Kronecker(pde_->R0(), Pt_);
      }
      return pen_;
    }
    
    // destructor
    virtual ~SpaceTimeSeparableBase() = default;  
  };
  
}}

#endif // __SPACE_TIME_SEPARABLE_BASE_H__

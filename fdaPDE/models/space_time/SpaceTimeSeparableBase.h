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

  // base class for separable regularization solved using either a moholitic or iterative solution strategy
  template <typename Model, typename Solver> class SpaceTimeSeparableBase;
  
  // base class for separable regularization
  template <typename Model>
    class SpaceTimeSeparableBase<Model, MonolithicSolver> : public SpaceTimeBase<Model> {
  private:
    // let \phi_i the i-th basis function in time
    SpMatrix<double> Rt_;  // mass matrix in time: [Rt_]_{ij} = \int_{[0,T]} \phi_i*\phi_j
    SpMatrix<double> Pt_;  // penalty matrix in time: [Pt_]_{ij} = \int_{[0,T]} (\phi_i)_tt*(\phi_j)_tt
    SpMatrix<double> Phi_; // [Phi_]_{ij} = \phi_i(t_j)
    DVector<double> u_;    // stacked discretized forcing [u_1 \ldots u_n, \ldots, u_1 \ldots u_n]
    SpMatrix<double> R0_;  // Rt_ \kron R0 (R0: discretization of the identity operator)
    SpMatrix<double> R1_;  // Rt_ \kron R1 (R1: discretization of the differential operator L in the regularizing PDE)
    
    typedef typename model_traits<Model>::TimeBasis TimeBasis; // basis used for discretization in time
    TimeBasis basis_;
    DVector<double> time_locations_; // time instants t_1, ..., t_m. used only if supplied via set_temporal_locations()
    SpMatrix<double> penS_; // discretization of space regularization: (R1^T*R0^{-1}*R1) \kron Rt
    SpMatrix<double> penT_; // discretization of time regularization:  (R0 \kron Pt)
  public:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    typedef typename model_traits<Model>::regularization TimeRegularization; // regularization in time
    typedef SpaceTimeBase<Model> Base;
    using Base::pde_;    // regularizing term in space Lf = u
    using Base::model;   // underlying model object
    using Base::time_;   // time interval [0,T]
    using Base::lambdaS; // smoothing parameter in space
    using Base::lambdaT; // smoothing parameter in time
    
    // constructor
    SpaceTimeSeparableBase() = default;
    SpaceTimeSeparableBase(const PDE& pde, const DVector<double>& time)
      : SpaceTimeBase<Model>(pde, time) {}
    // init data structure related to separable regularization
    void init_regularization() {
      basis_ = TimeBasis(time_);
      // compute \Phi = [\Phi]_{ij} = \phi_i(t_j) using the provided basis function
      if(is_empty(time_locations_))
	Phi_ = basis_.eval(time_); // assume time instants t_1, ..., t_m equal to time nodes
      else Phi_ = basis_.eval(time_locations_); // evaluate on given time locations

      // discretize operators in time
      TimeAssembler<TimeBasis> assembler(time_);
      Rt_ = assembler.assemble(TimeMass<TimeBasis>());    // mass matrix in time
      Pt_ = assembler.assemble(TimePenalty<TimeBasis>()); // time penalty matrix
      // compute tensorized matrices
      R0_ = Kronecker(Rt_, pde_->R0());
      R1_ = Kronecker(Rt_, pde_->R1());
    }
    // setters
    void set_temporal_locations(const DVector<double> time_locations) { time_locations_ = time_locations; }
    // getters
    const SpMatrix<double>& R0()  const { return R0_; }
    const SpMatrix<double>& R1()  const { return R1_; }
    // matrices proper of separable regularization
    const SpMatrix<double>& Rt()  const { return Rt_; }
    const SpMatrix<double>& Pt()  const { return Pt_; }
    const SpMatrix<double>& Phi() const { return Phi_; }
    inline std::size_t n_temporal_locs() const { return is_empty(time_locations_) ? time_.rows() : time_locations_.rows(); }
    const DVector<double>& time_locs() const { return is_empty(time_locations_) ? time_ : time_locations_; }

    // return stacked version of discretized forcing field
    const DVector<double>& u() {
      if(is_empty(u_)){ // compute once and cache
	std::size_t M = basis_.size();
	std::size_t N = Base::n_basis();
	u_.resize(N*M);
	// in separable regularization PDE doesn't depend on time. stack forcing term m times
	for(std::size_t i = 0; i < M; ++i) u_.segment(i*N, N) = pde_->force();
      }
      return u_;
    }

    // computes and cache matrices (R1^T*R0^{-1}*R1) \kron Rt and R0 \kron Pt.
    // returns an expression encoding \lambda_S*((R1^T*R0^{-1}*R1) \kron Rt) + \lambda_T*(R0 \kron Pt)
    auto pen() {
      if(is_empty(penS_)){ // compute once and cache result
	fdaPDE::SparseLU<SpMatrix<double>> invR0_;
	invR0_.compute(pde_->R0());
	penS_ = Kronecker(pde_->R1().transpose()*invR0_.solve(pde_->R1()), Rt_); // (R1^T*R0^{-1}*R1) \kron Rt
	penT_ = Kronecker(pde_->R0(), Pt_); // (R0 \kron Pt)
      }
      return lambdaS()*penS_ + lambdaT()*penT_;
    }
    
    // destructor
    virtual ~SpaceTimeSeparableBase() = default;  
  };
  
}}

#endif // __SPACE_TIME_SEPARABLE_BASE_H__

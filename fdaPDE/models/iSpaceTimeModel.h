#ifndef __I_SPACE_TIME_MODEL__
#define __I_SPACE_TIME_MODEL__

#include <cstddef>
#include <memory>
#include <Eigen/LU>
#include <type_traits>

#include "../core/utils/Symbols.h"
#include "../core/utils/Traits.h"
#include "iStatModel.h"
using fdaPDE::models::iStatModel;
#include "space_time/TimeAssembler.h"
using fdaPDE::models::TimeAssembler;
using fdaPDE::models::TimeMass;
using fdaPDE::models::TimePenalty;
#include "../core/NLA/KroneckerProduct.h"
using fdaPDE::core::NLA::SparseKroneckerProduct;
using fdaPDE::core::NLA::Kronecker;

namespace fdaPDE {
namespace models {
  
  // abstract base interface for any *space-time* fdaPDE statistical model. Uses CRTP pattern. Here is handled the different
  // type of time penalization, either separable or parabolic
  template <typename Model>
  class iSpaceTimeModel : public iStatModel<Model> {
    // check Model refers to the space-time case
    static_assert(!std::is_same<typename model_traits<Model>::RegularizationType, SpaceOnly>::value);
  protected:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    typedef typename model_traits<Model>::RegularizationType TimeRegularization; // type of time regularization
    typedef iStatModel<Model> Base;

    double lambdaS_; // smoothing parameter in space
    double lambdaT_; // smoothing parameter in time
    DVector<double> time_; // time domain [0, T]
  public:
    // import symbols from general stat model base
    using Base::pde_;
    // constructor
    iSpaceTimeModel() = default;
    iSpaceTimeModel(const PDE& pde, const DVector<double>& time)
      : iStatModel<Model>(pde), time_(time) {};
    // copy constructor, copy only pde object (as a consequence also the space domain) and time domain
    iSpaceTimeModel(const iSpaceTimeModel& rhs) { pde_ = rhs.pde_; time_ = rhs.time_;}

    // smoothing parameters
    void setLambdaS(double lambdaS) { lambdaS_ = lambdaS; }
    inline double lambdaS() const { return lambdaS_; }
    void setLambdaT(double lambdaT) { lambdaT_ = lambdaT; }    
    inline double lambdaT() const { return lambdaT_; }
    
    // destructor
    virtual ~iSpaceTimeModel() = default;  
  };

  // base class for separable regularization
  template <typename Model>
  class iSpaceTimeSeparableModel : public iSpaceTimeModel<Model> {
  private:
    SpMatrix<double> Rt_;  // mass matrix in time: [Rt_]_{ij} = \int_{[0,T]} \phi_i*\phi_j
    SpMatrix<double> Pt_;  // penalty matrix in time: [Pt_]_{ij} = \int_{[0,T]} (\phi_i)_tt*(\phi_j)_tt
    SpMatrix<double> Phi_; // [Phi_]_{ij} = \phi_i(t_j)

    DVector<double> y_; // vector of stacked observations [y_{11}, \ldots, y_{n1}, \ldots, y_{n1}, \ldots, y_{nm}]
    DVector<double> u_; // vector of stacked discretized force [u_1 \ldots u_n, \ldots, u_1 \ldots u_n] m times
    SpMatrix<double> PsiTD_; // matrix \Psi corrected for areal observations (stores \Psi^T*D if D \neq I)
    
    typedef typename model_traits<Model>::TimeBasis TimeBasis; // basis used for discretization in time
    TimeBasis basis_;
    DVector<double> time_; // interval [0,T]
  public:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    typedef typename model_traits<Model>::RegularizationType TimeRegularization; // type of time regularization
    typedef iSpaceTimeModel<Model> Base;
    using Base::Psi_; // matrix of spatial basis evaluations
    using Base::pde_; // regularizing term in space
    using Base::df_;  // block frame of data
    
    // constructor
    iSpaceTimeSeparableModel(const PDE& pde, const DVector<double>& time)
      : iSpaceTimeModel<Model>(pde, time), time_(time), basis_(time) {
      // compute \Phi = [\Phi]_{ij} = \phi_i(t_j) using the provided basis function
      Phi_ = basis_.eval(time);
      if(this->sampling() == SamplingStrategy::Areal){
	// in case of areal data store \Psi^T corrected for areal observations as (\Phi \kron \Psi)^T * (I_m \kron D)
	// being D diagonal I_m \kron D = diag(d_11 d_22 \ldots d_nn d_11 \ldots d_nn)
	std::size_t m = time.rows();
	std::size_t n = this->locs();
	DVector<double> D(n*m);
	for(std::size_t i = 0; i < m; ++i) D.segment(i*n, n) = this->D_.diagonal();
	// compute and store result
	PsiTD_ = Psi().transpose()*D.asDiagonal();
      }      
      // discretize operators in time
      TimeAssembler<TimeBasis> assembler(time);
      Rt_ = assembler.assemble(TimeMass<TimeBasis>()); // mass matrix in time
      Pt_ = assembler.assemble(TimePenalty<TimeBasis>()); // time penalty matrix
    }

    SparseKroneckerProduct<> R0()  const { return Kronecker(Rt_, pde_->R0()); }
    SparseKroneckerProduct<> R1()  const { return Kronecker(Rt_, pde_->R1()); }
    SparseKroneckerProduct<> Psi() { return Kronecker(Phi_, this->__Psi()); }
    // getters to raw matrices relative to separable time penalization
    const SpMatrix<double>& Rt()  const { return Rt_; }
    const SpMatrix<double>& Pt()  const { return Pt_; }
    const SpMatrix<double>& Phi() const { return Phi_; }

    // overload of iStatModel::u to directly return a stacked version of the discretized forcing field
    const DVector<double>& u() {
      // compute result dimensions.
      std::size_t M = basis_.size();
      std::size_t N = this->nbasis();
      if(u_.size() == 0){ // first time this is computed
	u_.resize(N*M);
	// in separable case the PDE used for space penalization is not function of time.
	// repeat discretized forcing term m times 
	for(std::size_t i = 0; i < M; ++i)
	  u_.segment(i*N, N) = pde_->force();
      }
      return u_;
    }
    // returns the block (\Phi \kron \Psi)^T*(I_m \kron D) as eigen expression, if D = I returns (\Phi \kron \Psi)^T
    auto PsiTD() { return this->sampling() == SamplingStrategy::Areal ? PsiTD_ : Psi().transpose(); }; 
    const DVector<double>& time_domain() const { return time_; } // time interval [0,T]
    
    // destructor
    virtual ~iSpaceTimeSeparableModel() = default;  
  };
  
}}

#endif // __I_SPACE_ONLY_MODEL__

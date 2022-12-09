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
#include "space_time/SpaceTime.h"
using fdaPDE::models::Phi;
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
    double lambdaS_; // smoothing parameter in space
    double lambdaT_; // smoothing parameter in time
    DVector<double> time_; // time domain [0, T]
  public:
    // import symbols from general stat model base
    IMPORT_STAT_MODEL_SYMBOLS(Model);
    // constructor
    iSpaceTimeModel() = default;
    iSpaceTimeModel(const PDE& pde, const DVector<double>& time) : iStatModel<Model>(pde), time_(time) {};
    // copy constructor, copy only pde object (as a consequence also the space domain) and time domain
    iSpaceTimeModel(const iSpaceTimeModel& rhs) { pde_ = rhs.pde_; time_ = rhs.time_;}

    // smoothing parameters
    void setLambdaS(double lambdaS) { lambdaS_ = lambdaS; }
    inline double lambdaS() const { return lambdaS_; }
    void setLambdaT(double lambdaT) { lambdaT_ = lambdaT; }    
    inline double lambdaT() const { return lambdaT_; }

    // an efficient implementation of left multiplication by \Psi
    // DMatrix<double> lmbPsi(const DMatrix<double>& x) const;
    // auto PsiTD() const { // returns the block \Psi^T*D as eigen expression, if D = I returns \Psi^T
    //   return sampling() == SamplingStrategy::Areal ? PsiTD_ : Psi_.transpose(); }; 
    
    // destructor
    virtual ~iSpaceTimeModel() = default;  
  };

  // import all symbols from iStatModel interface in derived classes
#define IMPORT_SPACE_TIME_MODEL_SYMBOLS( ... )			 \
  /* first import general stat model symbols */			 \
  IMPORT_STAT_MODEL_SYMBOLS(__VA_ARGS__)			 \
  using iSpaceTimeModel<__VA_ARGS__>::lambdaS;			 \
  using iSpaceTimeModel<__VA_ARGS__>::lambdaT;			 \
  using iSpaceTimeModel<__VA_ARGS__>::time_;			 \
  
  template <typename Model>
  class iSpaceTimeSeparableModel : public iSpaceTimeModel<Model> {
  private:
    SpMatrix<double> Rt_;  // mass matrix in time: [Rt_]_{ij} = \int_{[0,T]} \phi_i*\phi_j
    SpMatrix<double> Pt_;  // penalty matrix in time: [Pt_]_{ij} = \int_{[0,T]} (\phi_i)_tt*(\phi_j)_tt
    SpMatrix<double> Phi_; // [Phi_]_{ij} = \phi_i(t_j)
  public:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    typedef typename model_traits<Model>::RegularizationType TimeRegularization; // type of time regularization
    typedef typename model_traits<Model>::TimeBasis TimeBasis; // basis used for discretization in time
    IMPORT_SPACE_TIME_MODEL_SYMBOLS(Model);
    using iStatModel<Model>::Psi_;
    // constructor
    iSpaceTimeSeparableModel(const PDE& pde, const DVector<double>& time) : iSpaceTimeModel<Model>(pde, time) {
	// compute space-time matrices for separable case
	Phi_ = PhiAssembler<TimeBasis>().compute(time_);
	TimeAssembler<TimeBasis> assembler(time_);
	Rt_ = assembler.assemble(TimeMass<TimeBasis>());
	Pt_ = assembler.assemble(TimePenalty<TimeBasis>());
    }

    SparseKroneckerProduct<SpMatrix<double>, SpMatrix<double>> R0()  const { // \tilde R0 = Rt \kron R0
      return Kronecker(Rt_, pde_->R0()); }
    SparseKroneckerProduct<SpMatrix<double>, SpMatrix<double>> R1()  const { // \tilde R1 = Rt \kron R1
      return Kronecker(Rt_, pde_->R1()); }
    SparseKroneckerProduct<SpMatrix<double>, SpMatrix<double>> Psi() const { // \tilde \Psi = \Psi \kron \Phi 
      return Kronecker(Psi_, Phi_); }
    // getters to raw matrices relative to separable time penalization
    const SpMatrix<double>& Rt()  const { return Rt_; }
    const SpMatrix<double>& Pt()  const { return Pt_; }
    const SpMatrix<double>& Phi() const { return Phi_; }
    const DMatrix<double>&  u()  const { return pde_->force(); }

    // destructor
    virtual ~iSpaceTimeSeparableModel() = default;  
  };

#define IMPORT_SPACE_TIME_SEPARABLE_MODEL_SYMBOLS( ... )	 \
  /* first import general stat model symbols */			 \
  IMPORT_SPACE_TIME_MODEL_SYMBOLS(__VA_ARGS__)			 \
  using iSpaceTimeSeparableModel<__VA_ARGS__>::R0;		 \
  using iSpaceTimeSeparableModel<__VA_ARGS__>::R1;		 \
  using iSpaceTimeSeparableModel<__VA_ARGS__>::Psi;		 \
  using iSpaceTimeSeparableModel<__VA_ARGS__>::Rt;		 \
  using iSpaceTimeSeparableModel<__VA_ARGS__>::Pt;		 \
  using iSpaceTimeSeparableModel<__VA_ARGS__>::Phi;		 \
  using iSpaceTimeSeparableModel<__VA_ARGS__>::u;		 \
  
}}

#endif // __I_SPACE_ONLY_MODEL__

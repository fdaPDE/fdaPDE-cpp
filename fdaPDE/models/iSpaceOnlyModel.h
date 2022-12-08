#ifndef __I_SPACE_ONLY_MODEL__
#define __I_SPACE_ONLY_MODEL__

#include <cstddef>
#include <memory>
#include <Eigen/LU>
#include <type_traits>

#include "../core/utils/Symbols.h"
#include "../core/utils/Traits.h"
#include "iStatModel.h"
using fdaPDE::models::iStatModel;

namespace fdaPDE {
namespace models {

  // abstract base interface for any *space-only* fdaPDE statistical model. Uses CRTP pattern
  template <typename Model>
  class iSpaceOnlyModel : public iStatModel<Model> {
  protected:
    typedef typename model_traits<Model>::PDE PDE; // PDE used for regularization in space
    double lambda_; // smoothing parameter
  public:
    // import symbols from general stat model base
    IMPORT_STAT_MODEL_SYMBOLS(Model);
    using iStatModel<Model>::Psi_;
    using iStatModel<Model>::PsiTD_;
    
    // constructor
    iSpaceOnlyModel() = default;
    iSpaceOnlyModel(const PDE& pde)
      : iStatModel<Model>(pde) {};
    // copy constructor, copy only pde object (as a consequence also the problem domain)
    iSpaceOnlyModel(const iSpaceOnlyModel& rhs) { pde_ = rhs.pde_; }

    void setLambda(double lambda) { lambda_ = lambda; } 
    double lambda() const { return lambda_; } // smoothing parameter \lambda
    // pointers to FEM related quantites
    const SpMatrix<double>& R0() const { return pde_->R0(); }
    const SpMatrix<double>& R1() const { return pde_->R1(); }
    const DMatrix<double>&  u()  const { return pde_->force(); }
    const SpMatrix<double>& Psi() { return this->__Psi(); }; // pointer to n x N sparse matrix \Psi. This computes \Psi if not available
    // an efficient implementation of left multiplication by \Psi
    DMatrix<double> lmbPsi(const DMatrix<double>& x) const;
    auto PsiTD() const { // returns the block \Psi^T*D as eigen expression, if D = I returns \Psi^T
      return sampling() == SamplingStrategy::Areal ? PsiTD_ : Psi_.transpose(); }; 
    
    // destructor
    virtual ~iSpaceOnlyModel() = default;  
  };

template <typename E>
DMatrix<double> iSpaceOnlyModel<E>::lmbPsi(const DMatrix<double>& x) const {
  // compute dimensions of resulting matrix
  std::size_t n = Psi_.rows();
  std::size_t m = x.cols();
  // preallocate space for n x m result matrix
  DMatrix<double> result(n,m);
  // if data are sampled at mesh nodes (or a subset of them) then \Psi is a permutation matrix
  if(dataAtNodes()){
    // just permute input matrix columns
    for(std::size_t k = 0; k < Psi_.outerSize(); ++k){
      for (SpMatrix<double>::InnerIterator it(Psi_, k); it; ++it){
	result.row(it.row()) = x.row(it.col());
      }
    }
  }else{
    // in the general case no optimization can be put in place
    result = Psi_*x;
  }
  return result;
}
  
  // import all symbols from iStatModel interface in derived classes
#define IMPORT_SPACE_ONLY_MODEL_SYMBOLS(Model)			 \
  /* first import general stat model symbols */			 \
  IMPORT_STAT_MODEL_SYMBOLS(Model)				 \
  using iSpaceOnlyModel<Model>::lambda;				 \
  using iSpaceOnlyModel<Model>::Psi;				 \
  using iSpaceOnlyModel<Model>::R0;				 \
  using iSpaceOnlyModel<Model>::R1;				 \
  using iSpaceOnlyModel<Model>::u;				 \
  using iSpaceOnlyModel<Model>::PsiTD;				 \
  
}}

#endif // __I_SPACE_ONLY_MODEL__

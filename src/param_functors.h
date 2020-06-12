/*
 * param_functors.hpp
 *
 *  Created on: Jun 2, 2015
 *      Author: eardi
 */
#ifndef PARAM_FUNCTORS_H_
#define PARAM_FUNCTORS_H_

template <UInt ORDER, UInt mydim, UInt ndim>
class FiniteElement;

#include "pde_expression_templates.h"

template <UInt ndim, bool is_SV = false>
struct Diffusion{
  using diffusion_matr = Eigen::Matrix<Real,ndim,ndim>;

	Diffusion(const diffusion_matr& K) :
			K_(K) {}

  #ifdef R_VERSION_
  Diffusion(SEXP RGlobalVector){
    K_ = Eigen::Map<const diffusion_matr>(&REAL(RGlobalVector)[0]);
  }
  #endif

  template<UInt ORDER, UInt mydim, UInt WRONG>
  auto operator() (const FiniteElement<ORDER,mydim,WRONG>& fe_, UInt iq, UInt i, UInt j) const -> decltype(fe_.stiff_impl(iq, i, j, diffusion_matr())) {
    static_assert(ndim==WRONG, "ERROR! INCOMPATIBLE DIMENSIONS OF DIFFUSION PARAMETER AND FINITE ELEMENT! See param_functors.h");
    return fe_.stiff_impl(iq, i, j, K_);
  }

private:
  diffusion_matr K_;
};


template <UInt ndim>
struct Diffusion<ndim, true>{
  using diffusion_matr = Eigen::Matrix<Real,ndim,ndim>;
	using diff_matr_container = std::vector<diffusion_matr, Eigen::aligned_allocator<diffusion_matr> >;

	Diffusion(const diff_matr_container& K):
			K_(K){}

  #ifdef R_VERSION_
	Diffusion(SEXP RGlobalVector){
		UInt num_int_nodes = Rf_length(RGlobalVector)/(ndim*ndim);
		K_.reserve(num_int_nodes);
		for(UInt i=0; i<num_int_nodes; ++i)
				K_.push_back(Eigen::Map<const diffusion_matr>(&REAL(RGlobalVector)[ndim*ndim*i]));
	}
	#endif

  template<UInt ORDER, UInt mydim, UInt WRONG>
  auto operator() (const FiniteElement<ORDER,mydim,WRONG>& fe_, UInt iq, UInt i, UInt j) const -> decltype(fe_.stiff_impl(iq, i, j, diffusion_matr())) {
    static_assert(ndim==WRONG, "ERROR! INCOMPATIBLE DIMENSIONS OF DIFFUSION PARAMETER AND FINITE ELEMENT! See param_functors.h");
    UInt globalIndex=fe_.getGlobalIndex(iq);
    return fe_.stiff_impl(iq, i, j, K_[globalIndex]);
  }

private:
  diff_matr_container K_;
};


template <UInt ndim, bool is_SV = false>
struct Advection{
	using advection_vec = Eigen::Matrix<Real,ndim,1>;

	Advection(const advection_vec& beta) :
		beta_(beta){}

  #ifdef R_VERSION_
  Advection(SEXP RGlobalVector){
    beta_ = Eigen::Map<const advection_vec>(&REAL(RGlobalVector)[0]);
  }
  #endif

  template<UInt ORDER, UInt mydim, UInt WRONG>
  auto operator() (const FiniteElement<ORDER,mydim,WRONG>& fe_, UInt iq, UInt i, UInt j) const -> decltype(fe_.grad_impl(iq, i, j, advection_vec())) {
    static_assert(ndim==WRONG, "ERROR! INCOMPATIBLE DIMENSIONS OF ADVECTION PARAMETER AND FINITE ELEMENT! See param_functors.h");
    return fe_.grad_impl(iq, i, j, beta_);
  }

  EOExpr<const Advection&> dot(const EOExpr<Grad>& grad) const {
    typedef EOExpr<const Advection&> ExprT;
    return ExprT(*this);
  }

private:
  advection_vec beta_;
};

template<UInt ndim>
struct Advection<ndim, true>{
	using advection_vec = Eigen::Matrix<Real,ndim,1>;
	using adv_vec_container = std::vector<advection_vec, Eigen::aligned_allocator<advection_vec> >;

	Advection(const adv_vec_container& beta):
		beta_(beta){}

	#ifdef R_VERSION_
	Advection(SEXP RGlobalVector){
		UInt num_int_nodes = Rf_length(RGlobalVector)/ndim;
		beta_.reserve(num_int_nodes);
		for(UInt i=0; i<num_int_nodes; ++i)
			beta_.push_back(Eigen::Map<const advection_vec>(&REAL(RGlobalVector)[ndim*i]));
	}
	#endif

  template<UInt ORDER, UInt mydim, UInt WRONG>
  auto operator() (const FiniteElement<ORDER,mydim,WRONG>& fe_, UInt iq, UInt i, UInt j) const -> decltype(fe_.grad_impl(iq, i, j, advection_vec())) {
    static_assert(ndim==WRONG, "ERROR! INCOMPATIBLE DIMENSIONS OF ADVECTION PARAMETER AND FINITE ELEMENT! See param_functors.h");
    UInt globalIndex=fe_.getGlobalIndex(iq);
    return fe_.grad_impl(iq, i, j, beta_[globalIndex]);
  }

  EOExpr<const Advection&> dot(const EOExpr<Grad>& grad) const {
    typedef EOExpr<const Advection&> ExprT;
    return ExprT(*this);
  }

private:
  adv_vec_container beta_;
};

struct Reaction{
	Reaction(const std::vector<Real>& c):
		c_(c) {}

#ifdef R_VERSION_
	Reaction(SEXP RGlobalVector){
		UInt num_int_nodes = Rf_length(RGlobalVector);
		c_.reserve(num_int_nodes);
		for(UInt i=0; i<num_int_nodes; ++i)
			c_.push_back(REAL(RGlobalVector)[i]);
	}
	#endif

  template<UInt ORDER, UInt mydim, UInt ndim>
  auto operator() (const FiniteElement<ORDER, mydim, ndim>& fe_, UInt iq, UInt i, UInt j) const -> decltype(0.*fe_.mass_impl(iq, i, j)){
    UInt globalIndex=fe_.getGlobalIndex(iq);
    return c_[globalIndex]*fe_.mass_impl(iq, i, j);
  }

  EOExpr<const Reaction&> operator* (const EOExpr<Mass>&  mass) const {
      typedef EOExpr<const Reaction&> ExprT;
      return ExprT(*this);
  }


private:
  std::vector<Real> c_;
};




struct ForcingTerm{

  ForcingTerm() = default;
	ForcingTerm(const std::vector<Real>& u) :
		u_(u) {}

	#ifdef R_VERSION_
	ForcingTerm(SEXP RGlobalVector){
		UInt num_int_nodes = Rf_length(RGlobalVector);
		u_.reserve(num_int_nodes);
		for(UInt i=0; i<num_int_nodes; ++i)
			u_.push_back(REAL(RGlobalVector)[i]);
	}
	#endif

	Real operator[](UInt globalNodeIndex) const {return u_[globalNodeIndex];}
private:
  std::vector<Real> u_;
};


#endif /* PARAM_FUNCTORS_H_ */

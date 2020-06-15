/*
 * param_functors.hpp
 *
 *  Created on: Jun 2, 2015
 *      Author: eardi
 */
#ifndef PARAM_FUNCTORS_H_
#define PARAM_FUNCTORS_H_

#include "pde_expression_templates.h"

// Forward declaration!
template <UInt ORDER, UInt mydim, UInt ndim>
class FiniteElement;

// Convenience enum for options
enum class PDEParameterOptions{Constant, SpaceVarying};

template<PDEParameterOptions OPTION>
struct Diffusion;

template<>
struct Diffusion<PDEParameterOptions::Constant>{

  Diffusion(const Real* const K_ptr) :
    K_ptr_(K_ptr) {}

  #ifdef R_VERSION_
  Diffusion(SEXP RGlobalVector) :
    K_ptr_(REAL(RGlobalVector)) {}
  #endif

  template<UInt ORDER, UInt mydim, UInt ndim>
  Real operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq, UInt i, UInt j) const {
    using EigenMap2Diff_matr = Eigen::Map<const Eigen::Matrix<Real,ndim,ndim> >;
    return fe_.stiff_impl(iq, i, j, EigenMap2Diff_matr(K_ptr_));
  }
private:
  const Real* const K_ptr_;
};

template<>
struct Diffusion<PDEParameterOptions::SpaceVarying>{

  Diffusion(const Real* const K_ptr) :
    K_ptr_(K_ptr) {}

  #ifdef R_VERSION_
  Diffusion(SEXP RGlobalVector) :
    K_ptr_(REAL(RGlobalVector)) {}
  #endif

  template<UInt ORDER, UInt mydim, UInt ndim>
  Real operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq, UInt i, UInt j) const {
    using EigenMap2Diff_matr = Eigen::Map<const Eigen::Matrix<Real,ndim,ndim> >;
    const UInt index = fe_.getGlobalIndex(iq)*ndim*ndim;
    return fe_.stiff_impl(iq, i, j, EigenMap2Diff_matr(&K_ptr_[index]));
  }
private:
  const Real* const K_ptr_;
};


template<PDEParameterOptions OPTION>
struct Advection;

template<>
struct Advection<PDEParameterOptions::Constant>{

  Advection(const Real* const b_ptr) :
    b_ptr_(b_ptr) {}

  #ifdef R_VERSION_
  Advection(SEXP RGlobalVector) :
    b_ptr_(REAL(RGlobalVector)) {}
  #endif

  template<UInt ORDER, UInt mydim, UInt ndim>
  Real operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq, UInt i, UInt j) const {
    using EigenMap2Adv_vec = Eigen::Map<const Eigen::Matrix<Real,ndim,1> >;
    return fe_.grad_impl(iq, i, j, EigenMap2Adv_vec(b_ptr_));
  }

  EOExpr<const Advection&> dot(const EOExpr<Grad>& grad) const {
    typedef EOExpr<const Advection&> ExprT;
    return ExprT(*this);
  }

private:
  const Real* const b_ptr_;
};

template<>
struct Advection<PDEParameterOptions::SpaceVarying>{

  Advection(const Real* const b_ptr) :
    b_ptr_(b_ptr) {}

  #ifdef R_VERSION_
  Advection(SEXP RGlobalVector) :
    b_ptr_(REAL(RGlobalVector)) {}
  #endif

  template<UInt ORDER, UInt mydim, UInt ndim>
  Real operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq, UInt i, UInt j) const {
    using EigenMap2Adv_vec = Eigen::Map<const Eigen::Matrix<Real,ndim,1> >;
    const UInt index = fe_.getGlobalIndex(iq)*ndim;
    return fe_.grad_impl(iq, i, j, EigenMap2Adv_vec(&b_ptr_[index]));
  }

  EOExpr<const Advection&> dot(const EOExpr<Grad>& grad) const {
    typedef EOExpr<const Advection&> ExprT;
    return ExprT(*this);
  }

private:
  const Real* const b_ptr_;
};


struct Reaction{

	Reaction(const Real* const  c_ptr) :
		c_ptr_(c_ptr) {}

  #ifdef R_VERSION_
	Reaction(SEXP RGlobalVector) :
    c_ptr_(REAL(RGlobalVector)) {}
	#endif

  template<UInt ORDER, UInt mydim, UInt ndim>
  Real operator() (const FiniteElement<ORDER, mydim, ndim>& fe_, UInt iq, UInt i, UInt j) const {
    const UInt index = fe_.getGlobalIndex(iq);
    return c_ptr_[index]*fe_.mass_impl(iq, i, j);
  }

  EOExpr<const Reaction&> operator* (const EOExpr<Mass>&  mass) const {
      typedef EOExpr<const Reaction&> ExprT;
      return ExprT(*this);
  }

private:
  const Real* const c_ptr_;
};


struct ForcingTerm{

  ForcingTerm() :
    u_ptr_(nullptr) {}
	ForcingTerm(const Real* const u_ptr) :
		u_ptr_(u_ptr) {}

	#ifdef R_VERSION_
	ForcingTerm(SEXP RGlobalVector) :
    u_ptr_(REAL(RGlobalVector)) {}
	#endif

  template<UInt ORDER, UInt mydim, UInt ndim>
  Real integrate (const FiniteElement<ORDER, mydim, ndim>& fe_, UInt i) const {
    const UInt index = fe_.getGlobalIndex(0);
    return fe_.forcing_integrate(i, &u_ptr_[index]);
  }

private:
  const Real* const u_ptr_;
};


#endif /* PARAM_FUNCTORS_H_ */

#ifndef __PARAM_FUNCTORS_H__
#define __PARAM_FUNCTORS_H__

#include "Pde_Expression_Templates.h"

// Forward declaration!
template <UInt ORDER, UInt mydim, UInt ndim>
class FiniteElement;

// Convenience enum for options
enum class PDEParameterOptions{Constant, SpaceVarying};

template<PDEParameterOptions OPTION>
class Diffusion{
public:

  Diffusion(const Real* const K_ptr) :
    K_ptr_(K_ptr) {}

  Diffusion(SEXP RGlobalVector) :
    K_ptr_(REAL(RGlobalVector)) {}

  template<UInt ORDER, UInt mydim, UInt ndim>
  Real operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq, UInt i, UInt j) const;

private:
  const Real* const K_ptr_;
};

template<>
template<UInt ORDER, UInt mydim, UInt ndim>
Real Diffusion<PDEParameterOptions::Constant>::operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq, UInt i, UInt j) const {
  using EigenMap2Diff_matr = Eigen::Map<const Eigen::Matrix<Real,ndim,ndim> >;

  return fe_.stiff_impl(iq, i, j, EigenMap2Diff_matr(K_ptr_));
}

template<>
template<UInt ORDER, UInt mydim, UInt ndim>
Real Diffusion<PDEParameterOptions::SpaceVarying>::operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq, UInt i, UInt j) const {
  using EigenMap2Diff_matr = Eigen::Map<const Eigen::Matrix<Real,ndim,ndim> >;

  const UInt index = fe_.getGlobalIndex(iq) * EigenMap2Diff_matr::SizeAtCompileTime;
  return fe_.stiff_impl(iq, i, j, EigenMap2Diff_matr(&K_ptr_[index]));
}

template<PDEParameterOptions OPTION>
class Advection{
public:

  Advection(const Real* const b_ptr) :
    b_ptr_(b_ptr) {}

  Advection(SEXP RGlobalVector) :
    b_ptr_(REAL(RGlobalVector)) {}

  template<UInt ORDER, UInt mydim, UInt ndim>
  Real operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq, UInt i, UInt j) const;

  EOExpr<const Advection&> dot(const EOExpr<Grad>& grad) const {
    typedef EOExpr<const Advection&> ExprT;
    return ExprT(*this);
  }

private:
  const Real* const b_ptr_;
};

template<>
template<UInt ORDER, UInt mydim, UInt ndim>
Real Advection<PDEParameterOptions::Constant>::operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq, UInt i, UInt j) const {
  using EigenMap2Adv_vec = Eigen::Map<const Eigen::Matrix<Real,ndim,1> >;

  return fe_.grad_impl(iq, i, j, EigenMap2Adv_vec(b_ptr_));
}

template<>
template<UInt ORDER, UInt mydim, UInt ndim>
Real Advection<PDEParameterOptions::SpaceVarying>::operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq, UInt i, UInt j) const {
  using EigenMap2Adv_vec = Eigen::Map<const Eigen::Matrix<Real,ndim,1> >;

  const UInt index = fe_.getGlobalIndex(iq) * EigenMap2Adv_vec::SizeAtCompileTime;
  return fe_.grad_impl(iq, i, j, EigenMap2Adv_vec(&b_ptr_[index]));
}


class Reaction{
public:

	Reaction(const Real* const  c_ptr) :
		c_ptr_(c_ptr) {}

	Reaction(SEXP RGlobalVector) :
    c_ptr_(REAL(RGlobalVector)) {}

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


class ForcingTerm{
public:

  ForcingTerm() :
    u_ptr_(nullptr) {}
	ForcingTerm(const Real* const u_ptr) :
		u_ptr_(u_ptr) {}

	ForcingTerm(SEXP RGlobalVector) :
    u_ptr_(REAL(RGlobalVector)) {}

  template<UInt ORDER, UInt mydim, UInt ndim>
  Real integrate (const FiniteElement<ORDER, mydim, ndim>& fe_, UInt i) const {
    const UInt index = fe_.getGlobalIndex(0);
    return fe_.forcing_integrate(i, &u_ptr_[index]);
  }

private:
  const Real* const u_ptr_;
};


#endif /* PARAM_FUNCTORS_H_ */

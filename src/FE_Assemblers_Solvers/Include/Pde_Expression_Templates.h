#ifndef __EXPRESSION_H__
#define __EXPRESSION_H__

//Forward declarations!
template <UInt ORDER, UInt mydim, UInt ndim>
class FiniteElement;

enum class PDEParameterOptions;

template<PDEParameterOptions OPTION>
struct Diffusion;


struct Stiff {
  template<UInt ORDER, UInt mydim, UInt ndim>
  Real operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq, UInt i, UInt j) const {
    return fe_.stiff_impl(iq, i, j);
  }
};

struct Grad {
  template<UInt ORDER, UInt mydim, UInt ndim>
  Real operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq, UInt i, UInt j) const {
    return fe_.grad_impl(iq, i, j);
  }
};

struct Mass {
  template<UInt ORDER, UInt mydim, UInt ndim>
  Real operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq, UInt i, UInt j) const {
    return fe_.mass_impl(iq, i, j);
  }
};


template<typename A>
class EOExpr{
	  //! "A" is a generic type
	  A a_;
public:
	//! A constructor.
	/*!
	 * \param object is a constant reference to a generic operator.
	 */
	 EOExpr(const A& a) : a_(a) {}

	 template<UInt ORDER, UInt mydim, UInt ndim>
 	 Real operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq, UInt i, UInt j) const {
		 return a_(fe_, iq, i, j);
	 }
};

template<>
class EOExpr<Stiff>{
	  //! "A" is a generic type
	  Stiff a_;
public:
	//! A constructor.
	/*!
	 * \param object is a constant reference to a generic operator.
	 */
	 EOExpr(const Stiff& a) : a_(a) {}

   template<PDEParameterOptions OPTION>
   EOExpr<const Diffusion<OPTION>&> operator[] (const Diffusion<OPTION>& K){
     typedef EOExpr<const Diffusion<OPTION>&> ExprT;
     return ExprT(K);
   }

	 template<UInt ORDER, UInt mydim, UInt ndim>
 	 Real operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq, UInt i, UInt j) const {
		 return a_(fe_, iq, i, j);
	 }
};


//composition of two wrappers (operator)
//composition of two wrappers (operator)
template<typename A, typename B, typename Op>
class EOBinOp{
		//! "A" is a generic type.
		/*!
		 * Stores the first operand.
		 */
		A a_;
		//! "B" is a generic type.
		/*!
		 * Stores the second operand.
		 */
		B b_;
public:
	//! A constructor.
	/*!
	 * \param a is a constant reference to a generic type.
	 * \param b is a constant reference to a generic type.
	 */
		EOBinOp(const A& a ,const B& b) : a_(a), b_(b) {}
	 //! A definition of operator () taking two arguments.
	 /*!
     * \param i is an unsigned int
     * \param j is an unsigned int
     * applies the generic operation defined by the type Op to the two generic objects a_, b_;
     * returns a type P variable
	 */
		template<UInt ORDER,UInt mydim,UInt ndim>
 		Real operator () (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq, UInt i, UInt j) const {
		  return Op::apply(a_(fe_, iq, i, j), b_(fe_, iq, i, j));
	  }
};

template<class B, class Op>
class EOBinOp<Real, B, Op>{
	Real M_a;
	B M_b;
public:
	EOBinOp(Real a, const B& b) : M_a(a), M_b(b) {};

	template<UInt ORDER, UInt mydim, UInt ndim>
  Real operator () (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq, UInt i, UInt j) const {
		return Op::apply(M_a, M_b(fe_, iq, i, j));
	}
};

//wrappers addition
//! A ETWAdd class: Expression Template Wrapper Addition
/*!
 * Class that defines Addition operation, following:
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 */

struct EOAdd{

	//! A static inline method taking two arguments.
	/*!
	 *The actual addition operation
	 * \param a is of P type, first addend
	 * \param b is of P type, second addend
	 */
	template<class T, class V>
	static auto apply(const T& a, const V& b) -> decltype(a+b) {
		return a+b;
	}
};

//multiplication by real scalar
//! A ETWMult class: Expression Template Wrapper Multiplication.
/*!
 * Class that defines Multiplication operation, following:
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 */

struct EOMult{
	 //! A static inline method taking two arguments.
	/*!
	 * The actual multiplication operation.
	 * \param a is of P type, first operand
	 * \param b is a Real, second operand
	 */
	template<class T>
 	static auto apply(Real a, const T &b) -> decltype(a*b) {
 		return a*b;
 	}
};

//operator +
//! Overloading of operator +.
/*!
 * Following:
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 * Takes two arguments:
 * \param a is const reference ETWrapper<P, A>
 * \param b is const reference ETWrapper<P, A>
 * \return a ETWrapper<P,ETWBinOp<P, ETWrapper<P,A>, ETWrapper<P, B>, ETWAdd<P> > which is resolved at compile time.
 */
template<typename A, typename B>
EOExpr<EOBinOp<EOExpr<A>, EOExpr<B>, EOAdd> >
operator + (const EOExpr<A>&  a, const EOExpr<B>&  b){

	  typedef EOBinOp<EOExpr<A>, EOExpr<B>, EOAdd > ExprT;
	  return EOExpr<ExprT> (ExprT(a,b));
}

template<typename B>
EOExpr<EOBinOp<Real, EOExpr<B>, EOMult> >
operator * (Real a, const EOExpr<B>& b){
	  typedef EOBinOp<Real, EOExpr<B>, EOMult> ExprT;
	  return EOExpr<ExprT> (ExprT(a,b));
}




#endif

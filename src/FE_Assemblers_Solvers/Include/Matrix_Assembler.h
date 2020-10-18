#ifndef __MATRIX_ASSEMBLER_H__
#define __MATRIX_ASSEMBLER_H__


#include "../../FdaPDE.h"
#include "Finite_Element.h"
#include "../../Mesh/Include/Mesh_Objects.h"
#include "Param_Functors.h"
#include "Spline.h"
#include "../../Mesh/Include/Mesh.h"
//! A Stiff class: a class for the stiffness operator.

class Stiff{
  private:
  public:
	//! A definition of operator () taking three arguments.
    /*!
     * Evaluates the stiffness operator (i,j) of the current planar finite element.
     * \param currentfe_ is an object of class FiniteElement<Integrator, ORDER,2,2>, current planar finite element
     * \param i is an unsigned int, current finite element local index
     * \param j is an unsigned int, current finite element local index
     * returns a double.
     */
	template<class Integrator, UInt ORDER>
	inline Real operator() (FiniteElement<Integrator, ORDER,2,2>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
	{
	   	Real s = 0;
	   	for (UInt icoor=0; icoor < 2; icoor ++)
	   	{
	   		s += currentfe_.invTrJPhiDerMaster(i, icoor, iq)*currentfe_.invTrJPhiDerMaster(j, icoor, iq);

	   	}
	   	return s;
	}

	//! A definition of operator () taking three arguments.
    /*!
     * Evaluates the stiffness operator (i,j) of the current superficial finite element.
     * \param currentfe_ is an object of class FiniteElement<Integrator, ORDER,2,3>, current finite element
     * \param i is an unsigned int, current finite element local index
     * \param j is an unsigned int, current finite element local index
     * returns a double.
     */

	template<class Integrator, UInt ORDER>
	inline Real operator() (FiniteElement<Integrator, ORDER,2,3>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
	{
		Real s = 0;

		Eigen::Matrix<Real,2,1> grad_phi_i;
		Eigen::Matrix<Real,2,1> grad_phi_j;

		grad_phi_i(0) = currentfe_.phiDerMaster(i, 0, iq);
		grad_phi_i(1) = currentfe_.phiDerMaster(i, 1, iq);
		grad_phi_j(0) = currentfe_.phiDerMaster(j, 0, iq);
		grad_phi_j(1) = currentfe_.phiDerMaster(j, 1, iq);

		s = grad_phi_i.dot(currentfe_.metric()*grad_phi_j);

	   	return s;
	}


	//! A definition of operator () taking three arguments.
    /*!
     * Evaluates the stiffness operator (i,j) of the current superficial finite element.
     * \param currentfe_ is an object of class FiniteElement<Integrator, ORDER,3,3>, current finite element
     * \param i is an unsigned int, current finite element local index
     * \param j is an unsigned int, current finite element local index
     * returns a double.
     */

	template<class Integrator, UInt ORDER>
	inline Real operator() (FiniteElement<Integrator, ORDER,3,3>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
	{
	   	Real s = 0;
	   	for (UInt icoor=0; icoor < 3; icoor ++)
	   	{
	   		s += currentfe_.invTrJPhiDerMaster(i, icoor, iq)*currentfe_.invTrJPhiDerMaster(j, icoor, iq);

	   	}
	   	return s;
	}


};

template <class Type>
class StiffAnys{
};

template <>
class StiffAnys<Eigen::Matrix<Real,2,2>>{
  private:

	const Eigen::Matrix<Real,2,2>& K_;
  public:
	//! A constructor.

    StiffAnys(const Eigen::Matrix<Real,2,2>& K): K_(K){};

     //! A definition of operator () taking four arguments.
    /*!
     * Evaluates the product of: the derivative of basis(i) with respect to coordinate ic1 and the derivative of basis(j) with respect
     * to coordinate ic2 ,on current finite elemente.
     * \param i is an unsigned int, current finite element local index
     * \param j is an unsigned int, current finite element local index
     * \param ic1 is an unsigned int, the variable respect whom the derivative is take: ic1=0 abscissa, ic1=1 ordinata
     * \param ic1 is an unsigned int, the variable respect whom the derivative is take: ic1=0 abscissa, ic1=1 ordinata
     * returns a double.
     */

    template<class Integrator, UInt ORDER>
	inline Real operator() (FiniteElement<Integrator, ORDER,2,2>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
	{
	   	Real s = 0;
	   	for (UInt icoor=0; icoor < 2; icoor ++)
	   	{
	   		s += currentfe_.invTrJPhiDerMaster(i, 0, iq)*K_(0,icoor)*currentfe_.invTrJPhiDerMaster(j, icoor, iq) +
			currentfe_.invTrJPhiDerMaster(i, 1, iq)*K_(1,icoor)*currentfe_.invTrJPhiDerMaster(j, icoor, iq);

	   		//s += currentfe_.invTrJPhiDerMaster(i, icoor, iq)*currentfe_.invTrJPhiDerMaster(j, icoor, iq);

	   	}
	   	return s;
	}

    template<class Integrator, UInt ORDER>
	inline Real operator() (FiniteElement<Integrator, ORDER,2,3>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0){return 0;}

template<class Integrator, UInt ORDER>
	inline Real operator() (FiniteElement<Integrator, ORDER,3,3>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0){return 0;}

};

template <>
class StiffAnys<Diffusivity>{
  private:

	const Diffusivity& K_;
  public:
	//! A constructor.

	StiffAnys(const Diffusivity& K): K_(K){};


     //! A definition of operator () taking four arguments.
    /*!
     * Evaluates the product of: the derivative of basis(i) with respect to coordinate ic1 and the derivative of basis(j) with respect
     * to coordinate ic2 ,on current finite elemente.
     * \param i is an unsigned int, current finite element local index
     * \param j is an unsigned int, current finite element local index
     * \param ic1 is an unsigned int, the variable respect whom the derivative is take: ic1=0 abscissa, ic1=1 ordinata
     * \param ic1 is an unsigned int, the variable respect whom the derivative is take: ic1=0 abscissa, ic1=1 ordinata
     * returns a double.
     */

    template <class Integrator, UInt ORDER>
	inline Real operator() (FiniteElement<Integrator, ORDER,2,2>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
	{
	   	Real s = 0;
	   	for (UInt icoor=0; icoor < 2; icoor ++)
	   	{
	   		UInt globalIndex = currentfe_.getGlobalIndex(iq);
	   		s += currentfe_.invTrJPhiDerMaster(i, 0, iq)*K_(globalIndex)(0,icoor)*currentfe_.invTrJPhiDerMaster(j, icoor, iq) +
			currentfe_.invTrJPhiDerMaster(i, 1, iq)*K_(globalIndex)(1,icoor)*currentfe_.invTrJPhiDerMaster(j, icoor, iq);

	   		//s += currentfe_.invTrJPhiDerMaster(i, icoor, iq)*currentfe_.invTrJPhiDerMaster(j, icoor, iq);
	   	}
	   	return s;
	}

    template <class Integrator, UInt ORDER>
	inline Real operator() (FiniteElement<Integrator, ORDER,2,3>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0){return 0;}

    template <class Integrator, UInt ORDER>
	inline Real operator() (FiniteElement<Integrator, ORDER,3,3>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0){return 0;}

};

//! A Mass class: a class for the mass operator.
class Mass{
	private:

	public:

    //! A definition of operator () taking three arguments.
    /*!
     * Evaluates the mass operator (i,j) of the current finite element.
     * \param currentfe_ is an object of class FiniteElement<Integrator, ORDER,2,2>, current planar finite element
     * \param i is an unsigned int, current finite element local index
     * \param j is an unsigned int, current finite element local index
     * returns a double.
     */
	template <class Integrator ,UInt ORDER>
    inline Real operator() (FiniteElement<Integrator, ORDER,2,2>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
    {
    	return currentfe_.phiMaster(i,iq)*  currentfe_.phiMaster(j,iq);
    };

    //! A definition of operator () taking three arguments.
    /*!
     * Evaluates the mass operator (i,j) of the current finite element.
     * \param currentfe_ is an object of class FiniteElement<Integrator, ORDER,2,3>, current planar finite element
     * \param i is an unsigned int, current finite element local index
     * \param j is an unsigned int, current finite element local index
     * returns a double.
     */

	template <class Integrator ,UInt ORDER>
    inline Real operator() (FiniteElement<Integrator, ORDER,2,3>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
    {
    	return currentfe_.phiMaster(i,iq)*  currentfe_.phiMaster(j,iq);
    };


    //! A definition of operator () taking three arguments.
    /*!
     * Evaluates the mass operator (i,j) of the current finite element.
     * \param currentfe_ is an object of class FiniteElement<Integrator, ORDER,3,3>, current planar finite element
     * \param i is an unsigned int, current finite element local index
     * \param j is an unsigned int, current finite element local index
     * returns a double.
     */
    template <class Integrator ,UInt ORDER>
    inline Real operator() (FiniteElement<Integrator, ORDER,3,3>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
    {
    	return currentfe_.phiMaster(i,iq)*  currentfe_.phiMaster(j,iq);
    };

};

//! A vGrad class: a class for the the vectorial Gradient operator.

class Grad{
	private:

	public:
    //! A definition of operator () taking three arguments.
    /*!
     * Evaluates the component ic of the vGrad operator (i,j) on the current finite elemente.
     * \param i is an unsigned int, current finite element local index
     * \param j is an unsigned int, current finite element local index
     * \param ic is an unsigned int, vGrad component to be evaluated
     * returns a double.
     */
	 template<class Integrator, UInt ORDER>
     inline Real operator() (FiniteElement<Integrator,ORDER,2,2>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
     {
    	 return currentfe_.phiMaster(i,iq)*currentfe_.invTrJPhiDerMaster(j,ic,iq);
     }


     // AGGIUNGERE NUOVI METODI PER MESH SUPERFICIALI:
     //E' UNA SORTA DI DUMMY. DA IMPLEMENTARE SERIAMENTE PER ndim=3

     	 template<class Integrator, UInt ORDER>
     inline Real operator() (FiniteElement<Integrator, ORDER,2,3>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0){return 0;};

     template<class Integrator, UInt ORDER>
     inline Real operator() (FiniteElement<Integrator, ORDER,3,3>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0){return 0;};


};


class TimeMass{
	private:
    //! A reference to FiniteElement<Integrator>
    /*!
     * Stores a reference to the finite element where the mass operator is evaluated.
     */
	//FiniteElement<Integrator, ORDER>& currentfe_;
	public:
	//! A constructor.
	/*!
	 \param fe is a reference to FiniteElement<Integrator>
	 */
	//Mass(FiniteElement<Integrator, ORDER> & fe ):currentfe_(fe){};
    //! A definition of operator () taking two arguments.
    /*!
     * Evaluates the mass operator (i,j) of the current finite elemente.
     * \param i is an unsigned int, current finite element local index
     * \param j is an unsigned int, current finite element local index
     * returns a double.
     */

    template <class Integrator, UInt DEGREE, UInt ORDER_DERIVATIVE>
    inline Real operator() (Spline<Integrator, DEGREE, ORDER_DERIVATIVE>& spline_, UInt i, UInt j, Real u)
    {
        return spline_.BasisFunctionDerivative(DEGREE, ORDER_DERIVATIVE, i, u) * spline_.BasisFunctionDerivative(DEGREE, ORDER_DERIVATIVE, j, u);
    }

};

//generic template class wrapper
//! A ETWrapper class: Expression Template Wrapper.
/*!
 * Class that mimic the behaviour of a generic operator defined above: following
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 */



template<typename A>
class EOExpr{
	private:
	  //! "A" is a generic type
	  A a_;
	public:
	//! A constructor.
	/*!
	 * \param object is a constant reference to a generic operator.
	 */
	  EOExpr(const A& a):a_(a){};
	 //! A definition of operator () which takes two arguments.
     /*!
     * Masks the behaviour of the correspondent operator in the above classes.
     * \param i is an unsigned int
     * \param j is an unsigned int
     * returns a P variable.
     */
	  //P operator() (UInt i, UInt j) {return a_(i,j);}
	 //! A definition of operator () which takes three arguments.
     /*!
     * Masks the behaviour of the correspondent operator in the above classes.
     * \param i is an unsigned int
     * \param j is an unsigned int
     * \param ic is an unsigned int
     * returns a P variable.
     */
	  EOExpr<StiffAnys<Eigen::Matrix<Real,2,2>> >  operator[] (const Eigen::Matrix<Real,2,2>& K)
      {
		  typedef EOExpr<StiffAnys<Eigen::Matrix<Real,2,2> > > ExprT;
		  StiffAnys<Eigen::Matrix<Real,2,2> > anys(K);
    	  return ExprT(anys);
    	  //StiffAnys<Eigen::Matrix<Real,2,2> > a(K);
      }

	  EOExpr<StiffAnys<Diffusivity> > operator[] (const Diffusivity& K)
	  {
		  typedef EOExpr<StiffAnys<Diffusivity> > ExprT;
		  StiffAnys<Diffusivity> anys(K);
		  return ExprT(anys);
		  //return EOExpr<P,A>(A(K));
	  }

	  template<typename Integrator, UInt ORDER,UInt mydim, UInt ndim>
      Real operator() (FiniteElement<Integrator, ORDER,mydim,ndim>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
      {
          return a_(currentfe_, i,j,iq,ic);
      }

    template<typename Integrator, UInt DEGREE, UInt ORDER_DERIVATIVE>
      Real operator() (Spline<Integrator, DEGREE, ORDER_DERIVATIVE>& spline_, UInt i, UInt j, Real u)
      {
          return a_(spline_, i, j, u);
      }

	 //! A definition of operator () which takes four arguments.
     /*!
     * Masks the behaviour of the correspondent operator in the above classes.
     * \param i is an unsigned int
     * \param j is an unsigned int
     * \param ic is an unsigned int
     * returns a P variable.
     */
	  //P operator() (UInt i, UInt j, UInt ic1, UInt ic2) { return a_(i,j,ic1,ic2);}
	};

//composition of two wrappers (operator)
//! A ETWBinOp class: Expression Template Wrapper Binary Operation
/*!
 * Class that implements an abstract binary operation defined by Op between two ETWrappers, following:
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 */
template<typename A, typename B, typename Op>
class EOBinOp{
	private:
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
		EOBinOp(const A& a ,const B& b): a_(a),b_(b){};
	 //! A definition of operator () taking two arguments.
	 /*!
     * \param i is an unsigned int
     * \param j is an unsigned int
     * applies the generic operation defined by the type Op to the two generic objects a_, b_;
     * returns a type P variable
	 */
		template<typename Integrator, UInt ORDER,UInt mydim,UInt ndim>
	  Real operator () (FiniteElement<Integrator, ORDER,mydim,ndim>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
	  {
		  return Op::apply(a_(currentfe_,i,j, iq, ic),b_(currentfe_, i,j, iq, ic));
	  }
	};

template<class B, class Op>
class EOBinOp<Real, B, Op>
{
	Real M_a;
	B M_b;
public:
	EOBinOp(Real a, const B& b):M_a(a),M_b(b) {};

	template<typename Integrator, UInt ORDER,UInt mydim, UInt ndim>
	inline Real operator()(FiniteElement<Integrator, ORDER,mydim,ndim>& currentfe_, int i, int j, int iq, int ic = 0)
	{
		return Op::apply(M_a,M_b(currentfe_, i, j, iq, ic));
	}
};

template<class B, class Op>
class EOBinOp<Function, B, Op>
{
	const Function& M_a;
	B M_b;
public:
	EOBinOp(const Function& a, const B& b):M_a(a),M_b(b) {};

	template<typename Integrator, UInt ORDER,UInt mydim, UInt ndim>
	inline Real operator()(FiniteElement<Integrator, ORDER,mydim,ndim>& currentfe_, int i, int j, int iq, int ic = 0)
	{
		UInt globalIndex = currentfe_.getGlobalIndex(iq);
		return Op::apply(M_a(globalIndex),M_b(currentfe_, i, j, iq, ic));
	}
};

//wrappers addition
//! A ETWAdd class: Expression Template Wrapper Addition
/*!
 * Class that defines Addition operation, following:
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 */

class EOAdd{
	public:
	//! A constructor.
	EOAdd(){}
	//! A stastic inline method taking two arguments.
	/*!
	 *The actual addition operation
	 * \param a is of P type, first addend
	 * \param b is of P type, second addend
	 */
	static inline Real apply(Real a, Real b){ return (a+b); }
};

//multiplication by real scalar
//! A ETWMult class: Expression Template Wrapper Multiplication.
/*!
 * Class that defines Multiplication operation, following:
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 */

class EOMult{
	public:
	//! A constructor
	EOMult(){}
	 //! A stastic inline method taking two arguments.
	/*!
	 * The actual multiplication operation.
	 * \param a is of P type, first operand
	 * \param b is a Real, second operand
	 */
	  static inline Real apply( Real a, Real b){ return (a*b);}
	  //not needed since in ETRBinOp I did "Op::apply(a_(i,j),b_)"
	  //static inline P apply(Real b, const P&a){ return (a*b);}
	};

//Dot
template<class A, class B>
class EODotProd
{
	A a_;
	B b_;
	public:
	EODotProd(const A& a, const B& b):a_(a), b_(b) {};
	template<typename Integrator, UInt ORDER, UInt ndim>
	inline Real operator()(FiniteElement<Integrator, ORDER,2,ndim>& currentfe_, int i, int j, int iq, int ic = 0)
	{
		Real s = 0.;
		for(ic = 0; ic < 2; ++ic)
			s += a_(ic) * b_(currentfe_, i, j, iq, ic);
	return s;
	}

	template<typename Integrator, UInt ORDER>
	inline Real operator()(FiniteElement<Integrator, ORDER,3,3>& currentfe_, int i, int j, int iq, int ic = 0)
	{
		Real s = 0.;
		for(ic = 0; ic < 3; ++ic)
			s += a_(ic) * b_(currentfe_, i, j, iq, ic);
	return s;
	}
};

//Dot
template<class B>
class EODotProd<Function, B>
{
	const Function& M_a;
	B M_b;
	public:
	EODotProd(const Function& a, const B& b):M_a(a), M_b(b) {};
	template<typename Integrator, UInt ORDER, UInt ndim>
	inline Real operator()(FiniteElement<Integrator, ORDER,2,ndim>& currentfe_, int i, int j, int iq, int ic = 0)
	{
		Real s = 0.;
		UInt globalIndex = currentfe_.getGlobalIndex(iq);
		for(ic = 0; ic < 2; ic++)
			s += M_a(globalIndex, ic) * M_b(currentfe_, i, j, iq, ic);
	return s;
	}

	template<typename Integrator, UInt ORDER>
	inline Real operator()(FiniteElement<Integrator, ORDER,3,3>& currentfe_, int i, int j, int iq, int ic = 0)
	{
		Real s = 0.;
		UInt globalIndex = currentfe_.getGlobalIndex(iq);
		for(ic = 0; ic < 3; ic++)
			s += M_a(globalIndex, ic) * M_b(currentfe_, i, j, iq, ic);
	return s;
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
EOExpr<EOBinOp<EOExpr<A>, EOExpr<B>, EOAdd > >
operator + (const EOExpr<A>&  a, const EOExpr<B>&  b){

	  typedef EOBinOp<EOExpr<A>, EOExpr<B>, EOAdd > ExprT;
	  return EOExpr<ExprT> (ExprT(a,b));
}

template<typename B>
EOExpr<EOBinOp<Function, EOExpr<B>, EOMult > >
operator * (const Function&  a, const EOExpr<B>&  b){

	  typedef EOBinOp<Function, EOExpr<B>, EOMult> ExprT;
	  return EOExpr<ExprT> (ExprT(a,b));
}

template<typename B>
EOExpr<EOBinOp<Real, EOExpr<B>, EOMult > >
operator * (Real a, const EOExpr<B>&  b){

	  typedef EOBinOp<Real, EOExpr<B>, EOMult > ExprT;
	  return EOExpr<ExprT> (ExprT(a,b));
}

template<typename B>
EOExpr<EODotProd<Eigen::Matrix<Real,2,1>, EOExpr<B> > >
dot(const Eigen::Matrix<Real,2,1>& a, const EOExpr<B>&  b){

	  typedef EODotProd<Eigen::Matrix<Real,2,1>, EOExpr<B> > ExprT;
	  return EOExpr<ExprT> (ExprT(a,b));
}

template<typename B>
EOExpr<EODotProd<Eigen::Matrix<Real,3,1>, EOExpr<B> > >
dot(const Eigen::Matrix<Real,3,1>& a, const EOExpr<B>&  b){

	  typedef EODotProd<Eigen::Matrix<Real,3,1>, EOExpr<B> > ExprT;
	  return EOExpr<ExprT> (ExprT(a,b));
}

template<typename B>
EOExpr<EODotProd<Function, EOExpr<B> > >
dot(const Function& a, const EOExpr<B>&  b){

	  typedef EODotProd<Function, EOExpr<B> > ExprT;
	  return EOExpr<ExprT> (ExprT(a,b));
}



//!A Assmbler class: discretize a generic differential operator in a sparse matrix
//template<UInt mydim, UInt ndim>
class Assembler{
	private:
	public:
	  //! A constructor
	  //Assembler (){};
	  //! A template member taking three arguments: discretize differential operator
	  /*!
	   * \param oper is a template expression : the differential operator to be discretized.
	   * \param mesh is const reference to a MeshHandler<ORDER,2,2>: the mesh where we want to discretize the operator.
	   * \param fe is a const reference to a FiniteElement
	   * stores the discretization in SPoper_mat_
	   */

	  //Return triplets vector
	  template<UInt ORDER, typename Integrator, typename A>
	  static void operKernel(EOExpr<A> oper,const MeshHandler<ORDER,2,2>& mesh,
	  	                     FiniteElement<Integrator, ORDER,2,2>& fe, SpMat& OpMat);

    template<UInt DEGREE, UInt ORDER_DERIVATIVE, typename Integrator, typename A>
    static void operKernel(EOExpr<A> oper, Spline<Integrator, DEGREE, ORDER_DERIVATIVE>& spline, SpMat& OpMat);

	  template<UInt ORDER, typename Integrator>
	  static void forcingTerm(const MeshHandler<ORDER,2,2>& mesh, FiniteElement<Integrator, ORDER,2,2>& fe, const ForcingTerm& u, VectorXr& forcingTerm);

	  //! A template member taking three arguments: discretize differential operator
	  /*!
	   * \param oper is a template expression : the differential operator to be discretized.
	   * \param mesh is const reference to a MeshHandler<ORDER,2,3>: the mesh where we want to discretize the operator.
	   * \param fe is a const reference to a FiniteElement
	   * stores the discretization in SPoper_mat_
	   */

	  template<UInt ORDER, typename Integrator, typename A>
	  static void operKernel(EOExpr<A> oper,const MeshHandler<ORDER,2,3>& mesh,
	  	                     FiniteElement<Integrator, ORDER,2,3>& fe, SpMat& OpMat);

	  template<UInt ORDER, typename Integrator>
	  static void forcingTerm(const MeshHandler<ORDER,2,3>& mesh, FiniteElement<Integrator, ORDER,2,3>& fe, const ForcingTerm& u, VectorXr& forcingTerm);


	  template<UInt ORDER, typename Integrator, typename A>
	  static void operKernel(EOExpr<A> oper,const MeshHandler<ORDER,3,3>& mesh,
	  	                     FiniteElement<Integrator, ORDER,3,3>& fe, SpMat& OpMat);

	  template<UInt ORDER, typename Integrator>
	  static void forcingTerm(const MeshHandler<ORDER,3,3>& mesh, FiniteElement<Integrator, ORDER,3,3>& fe, const ForcingTerm& u, VectorXr& forcingTerm);


};


#include "Matrix_Assembler_imp.h"

#endif

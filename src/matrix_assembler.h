#ifndef MATRIX_ASSEMBLER_H_
#define MATRIX_ASSEMBLER_H_

#include "param_functors.h"

#include "mesh_objects.h"

template <UInt ORDER, UInt mydim, UInt ndim>
class MeshHandler;

template <UInt ORDER, UInt mydim, UInt ndim>
class FiniteElement;

template<typename A>
class EOExpr;


//!A Assembler class: discretize a generic differential operator in a sparse matrix
//template<UInt mydim, UInt ndim>
struct Assembler{
  //! A template member taking three arguments: discretize differential operator
  /*!
   * \param oper is a template expression : the differential operator to be discretized.
   * \param mesh is const reference to a MeshHandler<ORDER,2,2>: the mesh where we want to discretize the operator.
   * \param fe is a const reference to a FiniteElement
   * stores the discretization in SPoper_mat_
   */

  //Return triplets vector
  template<UInt ORDER, UInt mydim, UInt ndim, typename A>
  static void operKernel(EOExpr<A> oper,const MeshHandler<ORDER,mydim,ndim>& mesh,
  	                     FiniteElement<ORDER,mydim,ndim>& fe, SpMat& OpMat);

  template<UInt ORDER, UInt mydim, UInt ndim>
  static void forcingTerm(const MeshHandler<ORDER,mydim,ndim>& mesh, FiniteElement<ORDER,mydim,ndim>& fe, const ForcingTerm& u, VectorXr& forcingTerm);

  template<UInt DEGREE, UInt ORDER_DERIVATIVE, typename Integrator, typename A>
  static void operKernel(EOExpr<A> oper, Spline<Integrator, DEGREE, ORDER_DERIVATIVE>& spline, SpMat& OpMat);

};

#include "matrix_assembler_imp.h"

#endif

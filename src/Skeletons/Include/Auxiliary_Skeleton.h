#ifndef __AUXILIARY_SKELETON_H__
#define __AUXILIARY_SKELETON_H__

#include "../../FdaPDE.h"
#include "../../FE_Assemblers_Solvers/Include/Finite_Element.h"
#include "../../FE_Assemblers_Solvers/Include/Matrix_Assembler.h"
#include "../../Mesh/Include/Mesh.h"

template<UInt ORDER, UInt mydim, UInt ndim>
SEXP get_integration_points_skeleton(SEXP Rmesh)
{
	using Integrator = typename FiniteElement<ORDER, mydim, ndim>::Integrator;
	using meshElement = typename MeshHandler<ORDER, mydim, ndim>::meshElement;

	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);

	SEXP result;
	PROTECT(result=Rf_allocVector(REALSXP, ndim * Integrator::NNODES * mesh.num_elements()));
	for(UInt i=0; i<mesh.num_elements(); ++i){
		meshElement el = mesh.getElement(i);
		for(UInt l = 0; l < Integrator::NNODES; ++l){
			Point<ndim> p{el.getM_J() * Integrator::NODES[l].eigenView()};
			p += el[0];
			for(UInt j=0; j < ndim; ++j)
				REAL(result)[j * mesh.num_elements() * Integrator::NNODES + i * Integrator::NNODES + l] = p[j];
		}
	}

	UNPROTECT(1);
	return(result);
}

template<UInt ORDER, UInt mydim, UInt ndim, typename A>
SEXP get_FEM_Matrix_skeleton(SEXP Rmesh, EOExpr<A> oper)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);

	FiniteElement<ORDER, mydim, ndim> fe;

	SpMat AMat;
	Assembler::operKernel(oper, mesh, fe, AMat);

	//Copy result in R memory
	SEXP result;
	result = PROTECT(Rf_allocVector(VECSXP, 2));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(INTSXP, AMat.nonZeros() , 2));
	SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, AMat.nonZeros()));

	int *rans = INTEGER(VECTOR_ELT(result, 0));
	Real  *rans2 = REAL(VECTOR_ELT(result, 1));
	UInt i = 0;
	for (UInt k=0; k < AMat.outerSize(); ++k)
		{
			for (SpMat::InnerIterator it(AMat,k); it; ++it)
			{
				//std::cout << "(" << it.row() <<","<< it.col() <<","<< it.value() <<")\n";
				rans[i] = 1+it.row();
				rans[i + AMat.nonZeros()] = 1+it.col();
				rans2[i] = it.value();
				i++;
			}
		}
	UNPROTECT(1);
	return(result);
}

#endif

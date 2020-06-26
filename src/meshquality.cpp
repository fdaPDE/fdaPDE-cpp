#define R_VERSION_

#include "fdaPDE.h"
#include "mesh_objects.h"
#include "mesh.h"

#include <cmath>


template<UInt mydim>
Eigen::Matrix<Real, mydim, mydim> constructW();

template<>
Eigen::Matrix<Real, 2, 2> constructW(){
	Eigen::Matrix<Real, 2, 2> W;
	W<< 0, 			.5,
		1, 	 std::sqrt(3)/2;
	return W;
}

template<UInt mydim>
Eigen::Matrix<Real, 3, 3> constructW(){
	Eigen::Matrix<Real, 3, 3> W;
	W<< 1, 			   .5,  			.5,
		0, 			   .5, 			   -.5,
		0, 1/std::sqrt(2), 1./std::sqrt(2);
	return W;

}


extern "C" {

template<UInt ORDER, UInt mydim, UInt ndim>
void meshQuality(SEXP Rmesh, SEXP Rtol){

	const Real tol = *REAL(Rtol);

	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);

	Eigen::Matrix<Real, mydim, mydim> W=constructW();
	Eigen::Matrix<Real, mydim, mydim> invW=W.inverse(); 

	for(UInt i=0; i<mesh.num_elements(); ++i){
		auto el=mesh.getElement(i);
		if(mydim/(el.getM_J().lazyProduct(invW).squaredNorm() * W.lazyProduct(el.getM_invJ()).squaredNorm()) < tol*tol)
			Rprintf("WARNING: Element %d is badly shaped!\n", i+1);
	}

}


}

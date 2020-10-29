#ifndef __FINITE_ELEMENT_IMP_H__
#define __FINITE_ELEMENT_IMP_H__

// Template auxiliary function forward declaration
// This function evaluate ^phi_i's at quadrature nodes
template<UInt NBASES, UInt mydim>
Eigen::Matrix<Real, NBASES, 1> reference_eval_point(const Point<mydim> &node);

// This function evaluate nabla ^phi_i's at quadrature nodes
template<UInt NBASES, UInt mydim>
Eigen::Matrix<Real, NBASES, mydim> reference_eval_der_point(const Point<mydim> &node);

//! FEData implementation
template <UInt ORDER, UInt mydim, UInt ndim>
FiniteElementData<ORDER, mydim, ndim>::FiniteElementData()
{
	setPhi();
	setPhiDer();
}

template <UInt ORDER, UInt mydim, UInt ndim>
FiniteElementData<ORDER, mydim, ndim>::~FiniteElementData()
{
}

template <UInt ORDER, UInt mydim, UInt ndim>
void FiniteElementData<ORDER, mydim, ndim>::updateElement(const Element<NBASES,mydim,ndim> &t)
{
	t_ = t;
	setElementPhiDer();
}

template <UInt ORDER, UInt mydim, UInt ndim>
void FiniteElementData<ORDER, mydim, ndim>::setPhi()
{
	for(UInt iq=0; iq<Integrator::NNODES; ++iq)
		referencePhi.row(iq) = reference_eval_point<NBASES,mydim>(Integrator::NODES[iq]).transpose();
}

template <UInt ORDER, UInt mydim, UInt ndim>
void FiniteElementData<ORDER, mydim, ndim>::setPhiDer()
{
	// .template is needed! See Eigen documentation regarding
	// the template and typename keywords in C++
	for(UInt iq=0; iq<Integrator::NNODES; ++iq)
		referencePhiDer.template block<mydim, NBASES>(0, NBASES*iq) = reference_eval_der_point<NBASES,mydim>(Integrator::NODES[iq]).transpose();
}

template <UInt ORDER, UInt mydim, UInt ndim>
void FiniteElementData<ORDER, mydim, ndim>::setElementPhiDer()
{
	// we need J^(-T) nabla( phi)
	for (UInt iq=0; iq < Integrator::NNODES; ++iq)
		elementPhiDer.template block<ndim, NBASES>(0, NBASES*iq).noalias() = t_.getM_invJ().transpose() * referencePhiDer.template block<mydim, NBASES>(0, NBASES*iq);
}



// Templates for auxiliary functions
// This function evaluates the basis function on the reference element
// at the quadrature nodes

template<>
inline Eigen::Matrix<Real, 3, 1> reference_eval_point(const Point<2> &node){
	return makeBaryCoord(node.eigenConstView());
}

template<>
inline Eigen::Matrix<Real, 4, 1> reference_eval_point(const Point<3> &node){
	return makeBaryCoord(node.eigenConstView());
}

template<>
inline Eigen::Matrix<Real, 6, 1> reference_eval_point<6,2>(const Point<2> &node){
	Eigen::Matrix<Real, 6, 1> phi;
	phi << (1-node[0]-node[1]) * (1-2*node[0]-2*node[1]),
										 node[0] * (2*node[0]-1),
										 node[1] * (2*node[1]-1),
									 4*node[0] * node[1],
								   4*node[1] * (1-node[0]-node[1]),
									 4*node[0] * (1-node[0]-node[1]);
	return phi;
}

template<>
inline Eigen::Matrix<Real, 10, 1> reference_eval_point<10,3>(const Point<3> &node){
	Eigen::Matrix<Real, 10, 1> phi;
	phi << (1-node[0]-node[1]-node[2]) * (1-2*node[0]-2*node[1]-2*node[2]),
														 node[0] * (2*node[0]-1),
														 node[1] * (2*node[1]-1),
														 node[2] * (2*node[2]-1),
													 4*node[0] * (1-node[0]-node[1]-node[2]),
													 4*node[1] * (1-node[0]-node[1]-node[2]),
													 4*node[2] * (1-node[0]-node[1]-node[2]),
													 4*node[0] * node[1],
													 4*node[1] * node[2],
													 4*node[2] * node[0];
	return phi;
}

// This function evaluates the ndim-gradient of basis function on the reference element
// at the quadrature nodes
template<>
inline Eigen::Matrix<Real, 3,2> reference_eval_der_point(const Point<2> &node){
	Eigen::Matrix<Real,3,2> B1;
	B1.row(0).setConstant(-1);
	B1.bottomRows(2).setIdentity();
	return B1;
}

template<>
inline Eigen::Matrix<Real, 4,3> reference_eval_der_point(const Point<3> &node){
	Eigen::Matrix<Real,4,3> B1;
	B1.row(0).setConstant(-1);
	B1.bottomRows(3).setIdentity();
	return B1;
}

template<>
inline Eigen::Matrix<Real, 6,2> reference_eval_der_point<6,2>(const Point<2> &node){
	Eigen::Matrix<Real,6,2> B2;
	B2 << 1-4*(1-node[0]-node[1]), 1-4*(1-node[0]-node[1]),
										4*node[0]-1, 											 0,
															0, 						 4*node[1]-1,
											4*node[1],  						 4*node[0],
										 -4*node[1], 4*(1-node[0]-2*node[1]),
				4*(1-2*node[0]-node[1]), 						  -4*node[0];
	return B2;
}

template<>
inline Eigen::Matrix<Real, 10,3> reference_eval_der_point<10,3>(const Point<3> &node){
	Eigen::Matrix<Real,10,3> B2;


	B2 << 1-4*(1-node[0]-node[1]-node[2]), 	1-4*(1-node[0]-node[1]-node[2]), 	1-4*(1-node[0]-node[1]-node[2]),
							  4*node[0]-1, 		 			 			  0,			                      0,
										0, 						4*node[1]-1, 								  0,
										0,								  0,						4*node[2]-1,
		  4*(1-2*node[0]-node[1]-node[2]),  				     -4*node[0], 						 -4*node[0],
							   -4*node[1],  4*(1-node[0]-2*node[1]-node[2]), 						 -4*node[1],
							   -4*node[2], 			  			 -4*node[2],    4*(1-node[0]-node[1]-2*node[2]),
								4*node[1],						  4*node[0], 					  			  0,
										0, 						  4*node[2], 						  4*node[1],
								4*node[2],								  0,						  4*node[0];
	return B2;
}

#endif

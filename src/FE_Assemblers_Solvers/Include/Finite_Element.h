#ifndef __FINITE_ELEMENT_H__
#define __FINITE_ELEMENT_H__

#include "../../FdaPDE.h"
#include "Integration.h"
#include "../../Mesh/Include/Mesh_Objects.h"

// This is an abstract base class that wraps Element objects
// It stores the data needed by a triangular or tetrahedral finite element
template <UInt ORDER, UInt mydim, UInt ndim>
class FiniteElementData{
	static_assert((ORDER==1 || ORDER==2) &&
								(mydim==2 || mydim==3) &&
								 mydim <= ndim,
								 "ERROR! TRYING TO INSTANTIATE FINITE ELEMENT WITH WRONG NUMBER OF NODES AND/OR DIMENSIONS! See finite_element.h");
public:
	// This type encodes an appropriate quadrature rule depending on order and dimension
	using Integrator = typename SpaceIntegratorHelper::Integrator<ORDER,mydim>;

	// Number of basis function on the element
	static constexpr UInt NBASES = how_many_nodes(ORDER,mydim);

	// A default constructor
	FiniteElementData();

	// No move/copy constructors/operations
	// Since the class acts as a wrapper they are not needed!
	// Moreover this guarantees that the class can't be passed by value
	// (it is less error prone)
	FiniteElementData(const FiniteElementData&) = delete;
	FiniteElementData(FiniteElementData&&) = delete;
	FiniteElementData& operator=(const FiniteElementData&) = delete;
	FiniteElementData& operator=(FiniteElementData&&) = delete;

	// Pure virtual destructor making the class an ABC
	// Note: this has negligible runtime cost for this class
	virtual ~FiniteElementData()=0;

	// A member that accepts a new element to wrap
	void updateElement(const Element<NBASES,mydim,ndim>& t);

	// Overloaded subscript operator
	// It returns the i-th point of the underlying element
	// Note: read-only access because changing a point in an element
	// currently invalidates its state (see: "mesh_objects.h")
	const Point<ndim>& operator[] (UInt i) const {return t_[i];}

	// Members returning the area/volume of the underlying element
	Real getMeasure() const {return t_.getMeasure();}
	Real getArea() const {return t_.getMeasure();}
	Real getVolume() const {return t_.getMeasure();}

	// A member returning the ID of the underlying element
	Real getId() const {return t_.getId();}

	// A member returning the global index of a quadrature node
	UInt getGlobalIndex(UInt iq) const {return Integrator::NNODES * t_.getId() + iq;}

protected:
	// The underlying element
	Element<NBASES,mydim,ndim> t_;
	// A matrix Phi
	// Phi(i,iq) is ^phi_i(node_iq), i.e. the i-th basis function evaluated at the
	// iq-th quadrature node on the reference element
	Eigen::Matrix<Real, Integrator::NNODES, NBASES> referencePhi;
	// A block matrix PhiDer
	// Each block iq is made of nabla ^phi(node_iq), i.e. it stores the gradients
	// of the basis functions evaluated at the quadrature nodes on the reference element
	Eigen::Matrix<Real, mydim, NBASES*Integrator::NNODES> referencePhiDer;
	// A block matrix elementPhiDer
	// Each block iq is made of nabla phi(node_iq), i.e. it stores the gradients
	// of the basis functions evaluated at the quadrature nodes on the underlying element
	Eigen::Matrix<Real, ndim, NBASES*Integrator::NNODES> elementPhiDer;

	// Members initializing the corresponding matrices at construction
	void setPhi();
	void setPhiDer();
	// A member updating elementPhiDer after each element update
	void setElementPhiDer();

};


// This class implements all the needed methods to assemble the FE matrix

template <UInt ORDER, UInt mydim, UInt ndim>
class FiniteElement : public FiniteElementData<ORDER, mydim, ndim> {
	
public:
	using Integrator = typename FiniteElementData<ORDER, mydim, ndim>::Integrator;
	static constexpr UInt NBASES = FiniteElementData<ORDER, mydim, ndim>::NBASES;

	FiniteElement()=default;

	Real stiff_impl(UInt iq, UInt i, UInt j) const {
		return this->elementPhiDer.col(iq*NBASES+i).dot(this->elementPhiDer.col(iq*NBASES+j));
	}

	template <class ArgType>
	Real stiff_impl(UInt iq, UInt i, UInt j, const Eigen::MatrixBase<ArgType>& K) const {
		static_assert(ArgType::RowsAtCompileTime==ndim && ArgType::ColsAtCompileTime==ndim,
			"ERROR! WRONG DIMENSIONS OF THE INPUT! See finite_element.h");
		return this->elementPhiDer.col(iq*NBASES+i).dot(K.lazyProduct(this->elementPhiDer.col(iq*NBASES+j)));
	}

	Real grad_impl(UInt iq, UInt i, UInt j) const {
		return this->referencePhi(iq,i) * this->elementPhiDer(0,iq*NBASES+j);
	}
	
	template <class ArgType>
	Real grad_impl(UInt iq, UInt i, UInt j, const Eigen::MatrixBase<ArgType>& b) const {
		static_assert(ArgType::RowsAtCompileTime==ndim && ArgType::ColsAtCompileTime==1,
			"ERROR! WRONG DIMENSIONS OF THE INPUT! See finite_element.h");
		return this->referencePhi(i,iq) * b.dot(this->elementPhiDer.col(iq*NBASES+j));
	}

	Real mass_impl(UInt iq, UInt i, UInt j) const {
		return this->referencePhi(iq,i) * this->referencePhi(iq,j);
	}

	Real forcing_integrate(UInt i, const Real* const local_u) const {
		using EigenMap2Forcing_vec = Eigen::Map<const Eigen::Matrix<Real, Integrator::NNODES, 1> >;
		return this->referencePhi.col(i).dot(EigenMap2Forcing_vec(local_u).cwiseProduct(EigenMap2Forcing_vec(&Integrator::WEIGHTS[0])));
	}

};

#include "Finite_Element_imp.h"

#endif

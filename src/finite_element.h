#ifndef __FINITE_ELEMENT_HPP__
#define __FINITE_ELEMENT_HPP__

#include "fdaPDE.h"
#include "integration.h"
#include "mesh_objects.h"

//!  This class implements all properties of a Triangular or Tetrahedral Finite Element
/*!
 * This class is the most important one of the entire code
 * and implements everything needed by a triangular or tetrahedral finite elemnt
 * 
 * It takes as a template parameter a class that implements the method used
 * for determining the mass, stiff and grad matrices
*/

template <class Integrator ,UInt ORDER, UInt mydim, UInt ndim>
class FiniteElement{
};


template <class Integrator ,UInt ORDER>
class FiniteElement<Integrator, ORDER, 2,2>{
private:
	Element<3*ORDER,2,2> reference_;
	Element<3*ORDER,2,2> t_;
	Eigen::Matrix<Real,3*ORDER, Integrator::NNODES> phiMapMaster_;
	//Numero basi locali x Num coordinate x numero nodi integrazione
	Eigen::Matrix<Real,3*ORDER, Integrator::NNODES*2> phiDerMapMaster_;
	Eigen::Matrix<Real,3*ORDER, Integrator::NNODES*2> invTrJPhiDerMapMaster_;
	
	void setPhiMaster();
	void setPhiDerMaster();
	void setInvTrJPhiDerMaster();

public:

	//! This is an empty constructor
    /*!
        For efficiency and Expression Templates organization of the
        code, the use of this class is based on the updateElement class
    */
	FiniteElement();
	
	//! A member updating the Finite Element properties
    /*!
      \param t an element from which to update the finite element properties
    */
	void updateElement(const Element<3*ORDER,2,2> &t);
	
	Real getAreaReference()
	{
		return reference_.getArea();
	}

	Real getDet()
	{
		return t_.getDetJ();
	}

	Point coorQuadPt(UInt iq)
	{
		return Point(t_.getM_J()(0,0)*Integrator::NODES[iq][0] + t_.getM_J()(0,1)*Integrator::NODES[iq][1] + t_[0][0],
				t_.getM_J()(1,0)*Integrator::NODES[iq][0] + t_.getM_J()(1,1)*Integrator::NODES[iq][1] + t_[0][1]);
	}
	
	UInt getGlobalIndex(UInt iq)
	{
		return Integrator::NNODES * t_.getId() + iq;
	}

	//Returns \hat{phi}
	Real phiMaster(UInt i, UInt iq) const;

	//Returns \nabla \hat{phi}
	Real phiDerMaster(UInt i, UInt ic, UInt iq) const;

	//Returns J^{-1} \nabla \hat{phi}
	Real invTrJPhiDerMaster(UInt i, UInt ic, UInt iq) const;
};


template <class Integrator ,UInt ORDER>
class FiniteElement<Integrator, ORDER, 2,3>{
private:
	Element<3*ORDER,2,2> reference_;
	Element<3*ORDER,2,3> t_;
	Eigen::Matrix<Real,3*ORDER, Integrator::NNODES> phiMapMaster_;
	//Numero basi locali x Num coordinate x numero nodi integrazione
	Eigen::Matrix<Real,3*ORDER, Integrator::NNODES*2> phiDerMapMaster_;
	//Eigen::Matrix<Real,3*ORDER, Integrator::NNODES*2> invTrJPhiDerMapMaster_;
	Eigen::Matrix<Real,2,2> metric_;
	
	void setPhiMaster();
	void setPhiDerMaster();
	//void setInvTrJPhiDerMaster();

public:

	//! This is an empty constructor
    /*!
        For efficiency and Expression Templates organization of the
        code, the use of this class is based on the updateElement class
    */
	FiniteElement();
	
	//! A member updating the Finite Element properties
    /*!
      \param t an element from which to update the finite element properties
    */
	void updateElement(const Element<3*ORDER,2,3> &t);
	
	Real getAreaReference()
	{
		return reference_.getArea();
	}

	Real getDet()
	{
		return t_.getDetJ();
	}

	Point coorQuadPt(UInt iq)
	{
		return Point(t_.getM_J()(0,0)*Integrator::NODES[iq][0] + t_.getM_J()(0,1)*Integrator::NODES[iq][1] + t_[0][0],
				t_.getM_J()(1,0)*Integrator::NODES[iq][0] + t_.getM_J()(1,1)*Integrator::NODES[iq][1] + t_[0][1],
				t_.getM_J()(2,0)*Integrator::NODES[iq][0] + t_.getM_J()(2,1)*Integrator::NODES[iq][1] + t_[0][2]);
	}
	
	UInt getGlobalIndex(UInt iq) 
	{
		return Integrator::NNODES * t_.getId() + iq;
	} 

	//Returns \hat{phi}
	Real phiMaster(UInt i, UInt iq) const;

	//Returns \nabla \hat{phi}
	Real phiDerMaster(UInt i, UInt ic, UInt iq) const;

	//Returns J^{-1} \nabla \hat{phi}
	//Real invTrJPhiDerMaster(UInt i, UInt ic, UInt iq) const;
	
	Eigen::Matrix<Real,2,2> metric()const {return metric_;};

};



//Implementazione FiniteElement con mydim=3 e ndim=3
template <class Integrator ,UInt ORDER>
class FiniteElement<Integrator, ORDER, 3,3>{
private:
	Element<6*ORDER-2,3,3> reference_;
	Element<6*ORDER-2,3,3> t_;
	Eigen::Matrix<Real,6*ORDER-2, Integrator::NNODES> phiMapMaster_;
	//Numero basi locali x Num coordinate x numero nodi integrazione
	Eigen::Matrix<Real,6*ORDER-2, Integrator::NNODES*3> phiDerMapMaster_;
	Eigen::Matrix<Real,6*ORDER-2, Integrator::NNODES*3> invTrJPhiDerMapMaster_;
	Eigen::Matrix<Real,3,3> metric_;
	
	void setPhiMaster();
	void setPhiDerMaster();
	void setInvTrJPhiDerMaster();
	

public:

	//! This is an empty constructor
    /*!
        For efficiency and Expression Templates organization of the
        code, the use of this class is based on the updateElement class
    */
	FiniteElement();
	
	//! A member updating the Finite Element properties
    /*!
      \param t an element from which to update the finite element properties
    */
	void updateElement(const Element<6*ORDER-2,3,3> &t);
	
	Real getVolumeReference()
	{
		return reference_.getVolume();
	}

	Real getDet()
	{
		return t_.getDetJ();
	}

	Point coorQuadPt(UInt iq)
	{
		return Point(t_.getM_J()(0,0)*Integrator::NODES[iq][0] + t_.getM_J()(0,1)*Integrator::NODES[iq][1]+t_.getM_J()(0,2)*Integrator::NODES[iq][2] + t_[0][0],
				t_.getM_J()(1,0)*Integrator::NODES[iq][0] + t_.getM_J()(1,1)*Integrator::NODES[iq][1]+t_.getM_J()(1,2)*Integrator::NODES[iq][2] + t_[0][1],
				t_.getM_J()(2,0)*Integrator::NODES[iq][0] + t_.getM_J()(2,1)*Integrator::NODES[iq][1]+t_.getM_J()(2,2)*Integrator::NODES[iq][2] + t_[0][2]);
	}
	
	UInt getGlobalIndex(UInt iq) 
	{
		return Integrator::NNODES * t_.getId() + iq;
	} 

	//Returns \hat{phi}
	Real phiMaster(UInt i, UInt iq) const;

	//Returns \nabla \hat{phi}
	Real phiDerMaster(UInt i, UInt ic, UInt iq) const;

	//Returns J^{-1} \nabla \hat{phi}
	Real invTrJPhiDerMaster(UInt i, UInt ic, UInt iq) const;
	
	Eigen::Matrix<Real,3,3> metric()const {return metric_;};

};

#include "finite_element_imp.h"

#endif

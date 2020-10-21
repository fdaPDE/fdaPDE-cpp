#ifndef __EVALUATOR_H__
#define __EVALUATOR_H__

#include <iostream>
#include <algorithm>

#include "../../FdaPDE.h"
#include "Finite_Element.h"
#include "Matrix_Assembler.h"
#include "../../Mesh/Include/Mesh.h"

//!  A class for the evaluation of the solution, given the coefficients of the global basis
/*!
 * This class, given a vector of coordinates evaluates the solution, fist locating the point
 * and secondly evaluating it respect to the coefficients of the basis
 * It depends on a template parameter that specifies the Order of the initializing mesh
*/

template <UInt ORDER, UInt mydim, UInt ndim>
class Evaluator
{
};


template <UInt ORDER>
class Evaluator<ORDER,2,2>
{
	public:
		//! A constructor. It initializes the constructor given a mesh object.
		Evaluator(const MeshHandler<ORDER,2,2>& mesh): mesh_(mesh){};

		//! A member that computes the evaluation of a Point in a mesh, given the bases' coefficients.
		/*!
		\param X a pointer to the x coordinates to evaluate.
		\param Y a pointer to the y coordinates to evaluate.
		\param length a unsigned integer containing the number of points to evaluate.
		\param coef a pointer to the vector of coefficients of the solution, the value in position i
		is associated to the basis \phi(i)
		\param fast a boolean that specifies if the algorithm is completely based on the walking
				algorithm (can miss locations in case of non convex structures)
		\param result a double pointer to an already allocated memory space, where the evaluations
		will be stored
		*/
		void eval(Real* X, Real *Y, UInt length, const Real *coef, bool redundancy, Real* result, std::vector<bool>& isinside);
		void evalWithInfo(Real* X, Real *Y, UInt length, const Real *coef, bool redundancy, Real* result, std::vector<bool>& isinside, const std::vector<UInt> & element_id, Real **barycenters);

		//! A member that computes the integral over regions divided by the measure of the region in a mesh,
		//  given the bases' coefficients.
		/*!
		\param incidenceMatrix a nRegions*nElements array telling which triangle compose each region.
		\param nRegions an unsigned integer containing the number of regions.
		\param nElements an unsigned integer containing the number of triangles.
		\param coef a pointer to the vector of coefficients of the solution, the value in position i
		is associated to the basis \phi(i).
		\param result a double pointer to an already allocated memory space, where the evaluations
		will be stored.
		*/
		void integrate(UInt** incidenceMatrix, UInt nRegions, UInt nElements, const Real *coef, Real* result);

	private:
		const MeshHandler<ORDER,2,2> &mesh_;

};

template <UInt ORDER>
class Evaluator<ORDER,2,3>
{
	public:
		//! A constructor. It initializes the constructor given a mesh object.
		Evaluator(const MeshHandler<ORDER,2,3>& mesh): mesh_(mesh){};

		//! A member that computes the evaluation of a Point in a mesh, given the bases' coefficients.
		/*!
		\param X a pointer to the x coordinates to evaluate.
		\param Y a pointer to the y coordinates to evaluate.
		\param Z a pointer to the x coordinates to evaluate.
		\param length a unsigned integer containing the number of points to evaluate.
		\param coef a pointer to the vector of coefficients of the solution, the value in position i
		is associated to the basis \phi(i)
		\param fast a boolean that specifies if the algorithm is completely based on the walking
				algorithm (can miss locations in case of non convex structures)
		\param result a double pointer to an already allocated memory space, where the evaluations
		will be stored
		*/
		void eval(Real* X, Real *Y, Real *Z, UInt length, const Real *coef, bool redundancy, Real* result, std::vector<bool>& isinside);
		void evalWithInfo(Real* X, Real *Y, Real *Z, UInt length, const Real *coef, bool redundancy, Real* result, std::vector<bool>& isinside, const std::vector<UInt> & element_id, Real **barycenters);

		//! A member that computes the integral over regions divided by the measure of the region in a mesh,
		//  given the bases' coefficients.
		/*!
		\param incidenceMatrix a nRegions*nElements array telling which triangle compose each region.
		\param nRegions an unsigned integer containing the number of regions.
		\param nElements an unsigned integer containing the number of triangles.
		\param coef a pointer to the vector of coefficients of the solution, the value in position i
		is associated to the basis \phi(i).
		\param result a double pointer to an already allocated memory space, where the evaluations
		will be stored.
		*/
		void integrate(UInt** incidenceMatrix, UInt nRegions, UInt nElements, const Real *coef, Real* result);

	private:
		const MeshHandler<ORDER,2,3> &mesh_;

};


//Implementazione evaluator mydim=3 ndim=3

template <UInt ORDER>
class Evaluator<ORDER,3,3>
{
	public:
		//! A constructor. It initializes the constructor given a mesh object.
		Evaluator(const MeshHandler<ORDER,3,3>& mesh): mesh_(mesh){};

		//! A member that computes the evaluation of a Point in a mesh, given the bases' coefficients.
		/*!
		\param X a pointer to the x coordinates to evaluate.
		\param Y a pointer to the y coordinates to evaluate.
		\param Z a pointer to the x coordinates to evaluate.
		\param length a unsigned integer containing the number of points to evaluate.
		\param coef a pointer to the vector of coefficients of the solution, the value in position i
		is associated to the basis \phi(i)
		\param fast a boolean that specifies if the algorithm is completely based on the walking
				algorithm (can miss locations in case of non convex structures)
		\param result a double pointer to an already allocated memory space, where the evaluations
		will be stored
		*/
		void eval(Real* X, Real *Y, Real *Z, UInt length, const Real *coef, bool redundancy, Real* result, std::vector<bool>& isinside);
		void evalWithInfo(Real* X, Real *Y, Real *Z, UInt length, const Real *coef, bool redundancy, Real* result, std::vector<bool>& isinside, const std::vector<UInt> & element_id, Real **barycenters);
		//! A member that computes the integral over regions divided by the measure of the region in a mesh,
		//  given the bases' coefficients.
		/*!
		\param incidenceMatrix a nRegions*nElements array telling which tetrahedron compose each region.
		\param nRegions an unsigned integer containing the number of regions.
		\param nElements an unsigned integer containing the number of tetrahedrons.
		\param coef a pointer to the vector of coefficients of the solution, the value in position i
		is associated to the basis \phi(i).
		\param result a double pointer to an already allocated memory space, where the evaluations
		will be stored.
		*/
		void integrate(UInt** incidenceMatrix, UInt nRegions, UInt nElements, const Real *coef, Real* result);

	private:
		const MeshHandler<ORDER,3,3> &mesh_;

};


#include "Evaluator_imp.h"

#endif

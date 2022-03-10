#ifndef __EVALUATOR_H__
#define __EVALUATOR_H__

#include <iostream>
#include <algorithm>

#include "../../FdaPDE.h"
#include "Finite_Element.h"
#include "Matrix_Assembler.h"
#include "../../Mesh/Include/Mesh.h"

//!Template class for function evaluation on given locations
/*!
 \tparam ORDER UInt representing order of the element of the mesh
 \tparam mydim UInt representing the mesh space size, 1 for 1.5D, 2 for 2D and 2.5D, 3 for 3D
 \tparam ndim UInt representing the space size, 2 for 1.5D and 2D, 3 for 2.5D and 3D
 \sa Mesh
*/
template<UInt ORDER, UInt mydim, UInt ndim>
class Evaluator
{
public:
    //! A constructor. It initializes the constructor given a mesh object.
    Evaluator(const MeshHandler<ORDER,mydim,ndim>& mesh): mesh_(mesh){};

    //! A member computing the evaluation of a Point in a mesh, given the bases' coefficients.
    /*!
    \param locations a RNumericMatrix containing the coordinates of points to evaluate.
    \param coef a RNumericMatrix containing the coefficients of the solution, the value in position i
    is associated to the basis \phi(i).
    \param redundancy a boolean to enforce verbose search for Walking algorithm. If it is true and the walking algorithm fails a new search is performed based on the naive algorithm
    \param result a RNumericMatrix to an already protected memory space in which the evaluations
    will be stored
    */
    template<bool isManifold=(ndim!=mydim)>
    typename std::enable_if<isManifold,void>::type
    eval(const RNumericMatrix& locations,const RNumericMatrix& coef, bool redundancy,RNumericMatrix& result, std::vector<bool>& isinside);

    template<bool isManifold=(ndim!=mydim)>
    typename std::enable_if<!isManifold,void>::type
    eval(const RNumericMatrix& locations,const RNumericMatrix& coef, bool redundancy,RNumericMatrix& result, std::vector<bool>& isinside);

    void evalWithInfo(const RNumericMatrix& locations,const RNumericMatrix& coef, bool redundancy,RNumericMatrix& result, std::vector<bool>& isinside, const RIntegerMatrix& element_id,const RNumericMatrix& barycenters);

    //! A member that computes the integral over regions divided by the measure of the region in a mesh,
    //!  given the bases' coefficients.
    /*!
    \param incidenceMatrix a nRegions*nElements RIntegerMatrix telling which triangle compose each region.
    \param coef a pointer to the vector of coefficients of the solution, the value in position i
    is associated to the basis \phi(i).
    \param result a RNumericMatrix to an already protected memory space, where the evaluations
    will be stored.
    */
    void integrate(const RIntegerMatrix& incidenceMatrix, const RNumericMatrix& coef, RNumericMatrix& result);

private:
    const MeshHandler<ORDER,mydim,ndim> & mesh_;
};



#include "Evaluator_imp.h"

#endif

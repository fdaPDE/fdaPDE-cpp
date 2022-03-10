#ifndef __MESH_OBJECTS_HPP__
#define __MESH_OBJECTS_HPP__


#include "../../FdaPDE.h"

#include "Point.h"

using Id=int;

/*! \brief A function to evaluate factorial at compile time
 *
 *  This function evaluates the factorial of a given UInt n
 *  It is constexpr so that it can be evaluated at compile time
 *  \param n a const UInt
 */
constexpr UInt factorial(const UInt n) {
    return n ? (n * factorial(n - 1)) : 1;
}

/*! \brief A function to compute the number of nodes in an element
 *
 *  This function compute the number of nodes in an element
 *  It is constexpr so that it can be evaluated at compile time
 *  \param ORDER a const UInt representing the order of the element
 *  \param mydim a const UInt representing the dimension of the element
 */

constexpr UInt how_many_nodes(const UInt ORDER, const UInt mydim) {
    return factorial(ORDER+mydim)/(factorial(ORDER)*factorial(mydim));
}


//! This class implements an Element (i.e. a triangle or tetrahedron)
//!  This class implements a Triangle as an objects composed by three or six nodes.
/*!
 *  The first three nodes represent the vertices, the others the internal nodes,
 *  following this enumeration: !IMPORTANT! different from Sangalli code!
 *
 * 		        3
 * 			    *
 * 		     /    \
 * 		  5 *	   * 4
 * 		  /	        \
 * 		 *_____*_____*
 * 		1	     6      2
*/


//!  This class implements a Tetrahedron as an objects composed by four or ten nodes, embedded in a 3-dimensional space.
// For NNODES=10 the edges are ordered like so: (1,2), (1,3), (1,4), (2,3), (3,4), (2,4)
// The midpoints are also expected to follow this convention!

// Note: only order 1 and 2 are implemented at the moment

template <UInt NNODES, UInt mydim, UInt ndim>
class Element : public Identifier {

  static_assert((mydim==2 || mydim==3) &&
                 mydim <= ndim &&
                 (NNODES==how_many_nodes(1,mydim) || NNODES==how_many_nodes(2,mydim)),
                 "ERROR! TRYING TO INSTANTIATE ELEMENT WITH WRONG NUMBER OF NODES AND/OR DIMENSIONS! See mesh_objects.h");

public:
  // Note: this macro is needed to avoid alignment issues!
  // See Eigen documentation "Structures Having Eigen Members"
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  using elementPoints = std::array<Point<ndim>,NNODES>;
  using iterator = typename elementPoints::iterator;
  using const_iterator = typename elementPoints::const_iterator;

  static constexpr UInt numVertices=mydim+1;
  static constexpr UInt numSides=mydim+1;
  static constexpr UInt myDim=mydim; 
  static constexpr UInt nDim=ndim; 

  // Note: these don't really mean anything, they're just here for compatibility
  // with the adtree implementation
  static constexpr UInt dp() {return ndim;}
  static constexpr UInt dt() {return 2*ndim;}
  static constexpr UInt coordsize() {return (ndim+1)*ndim;}

  //! This constructor creates an "empty" Element, with an Id Not Valid
  Element()=default;

	//! This constructor creates an Element, given its Id and an array containing the Points
	Element(UInt id, const elementPoints& points) :
					Identifier(id), points_(points) {computeProperties();}

	//! Overloaded subscript operator
  // Note: only the const version is available because any change in points_
  // would require a call to computeProperties() to keep the element in a valid state
	const Point<ndim>& operator[](UInt i) const {return points_[i];}

  //! Define begin and end iterators (this also gives "ranged for" for free)
  // Note: only the const version is available because any change in points_
  // would require a call to computeProperties() to keep the element in a valid state
  const_iterator begin() const {return points_.begin();}
  const_iterator end() const {return points_.end();}

	const Eigen::Matrix<Real,ndim,mydim>& getM_J() const {return M_J_;}
	const Eigen::Matrix<Real,mydim,ndim>& getM_invJ() const {return M_invJ_;}

  //! A member function that computes the barycentric coordinates of a given Point
	Eigen::Matrix<Real,mydim+1,1> getBaryCoordinates(const Point<ndim>& point) const;

  //! A member function that tests if a given Point is located inside the Element
  bool isPointInside(const Point<ndim>& point) const;

  //! A member function that returns beyond which edge/face a given Point lies
  // This function returns -1 if the point is inside the element
  int getPointDirection(const Point<ndim>&) const;

  //! Some members returning the area/volume of the element
  // Note: all three are available for convenience and compatibility reasons with existing code
  Real getMeasure() const {return element_measure;}
  Real getArea() const {return getMeasure();}
  Real getVolume() const {return getMeasure();}

  //! A member to evaluate a function at a point given the function's coefficients
  // on the element's basis functions
  // Note: this function assumes that the point is inside the element!
  Real evaluate_point(const Point<ndim>&, const Eigen::Matrix<Real,NNODES,1>&) const;
  //! A member to evaluate integrals on the element
  Real integrate(const Eigen::Matrix<Real,NNODES,1>&) const;

	//! Overload the << operator to easily print element info (note: define it in class
  // to avoid a forward declaration)
  template <UInt nnodes, UInt MYDIM, UInt NDIM>
  friend std::ostream& operator<<(std::ostream& os, const Element<nnodes, MYDIM, NDIM>& el);

private:
  // Data members
	elementPoints points_;
  // A matrix encoding a linear transformation from barycentric coordinates to
  // standard coordinates (modulo a translation)
	Eigen::Matrix<Real,ndim,mydim> M_J_;
  // A matrix which is the inverse of M_J_
	Eigen::Matrix<Real,mydim,ndim> M_invJ_;
  // A Real storing the area/volume of the element
  Real element_measure;

  // A member initializing M_J_, M_invJ_ and element_measure at construction
  void computeProperties();

  Real evaluate_point(const Eigen::Matrix<Real, mydim+1, 1>&, const Eigen::Matrix<Real,NNODES,1>&) const;

};


//! Partial template specialization for surface elements
//!  This class implements a Triangle as an objects composed by three or six nodes, embedded in a 3-dimensional space
/*!
 *  The first three nodes represent the vertices, the others the internal nodes,
 *  following this enumeration: !IMPORTANT! different from Sangalli code!
 *
 * 			    3
 * 			    *
 * 		     /    \
 * 		  5 *	   * 4
 * 		  /	        \
 * 		 *_____*_____*
 * 		1	   6	  2
*/

template <UInt NNODES>
class Element<NNODES, 2, 3> : public Identifier {

  static_assert(NNODES==3 || NNODES==6,
                 "ERROR! TRYING TO INSTANTIATE SURFACE ELEMENT WITH WRONG NUMBER OF NODES! See mesh_objects.h");

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  using elementPoints = std::array<Point<3>, NNODES>;
  using iterator = typename elementPoints::iterator;
  using const_iterator = typename elementPoints::const_iterator;

  static constexpr UInt numVertices=3;
  static constexpr UInt numSides=3;
  static constexpr UInt myDim=2; 
  static constexpr UInt nDim=3; 


  // Note: these don't really mean anything, they're just here for compatibility
  // with the adtree implementation
  static constexpr UInt dp() {return 3;}
  static constexpr UInt dt() {return 6;}
  static constexpr UInt coordsize() {return 9;}

  //! This constructor creates an "empty" Element, with an Id Not Valid
  Element()=default;

	//! This constructor creates an Element, given its Id and an std array with the Points
	Element(UInt id, const elementPoints& points) :
					Identifier(id), points_(points) {computeProperties();}

  //! Overloaded subscript operator
  // Note: only the const version is available because any change in points_
  // would require a call to computeProperties() to keep the element in a valid state
	const Point<3>& operator[](UInt i) const {return points_[i];}

  //! Define begin and end iterators (this also gives "ranged for" for free)
  // Note: only the const version is available because any change in points_
  // would require a call to computeProperties() to keep the element in a valid state
  const_iterator begin() const {return points_.begin();}
  const_iterator end() const {return points_.end();}

	const Eigen::Matrix<Real,3,2>& getM_J() const {return M_J_;}
	const Eigen::Matrix<Real,2,3>& getM_invJ() const {return M_invJ_;} //this is actually the pseudoinverse!

  //! A member function that computes the barycentric coordinates of a given point
	Eigen::Matrix<Real,3,1> getBaryCoordinates(const Point<3>& point) const;

  //! A member function that tests if a Point is located inside the Element
  bool isPointInside(const Point<3>& point) const;

  //! Some members returning the area/volume of the element
  // Note: all three are available for convenience and compatibility reasons with existing code
  Real getMeasure() const {return element_measure;}
  Real getArea() const {return getMeasure();}
  Real getVolume() const {return getMeasure();}

  //! A member to evaluate functions at a point inside the element
  Real evaluate_point(const Point<3>&, const Eigen::Matrix<Real,NNODES,1>&) const;
  //! A member to evaluate integrals on the element
  Real integrate(const Eigen::Matrix<Real,NNODES,1>&) const;

  // A member function that projects a 3D point (XYZ coordinates!) onto the element
  // Note: if the projection lies outside the element the function returns
  // the closest point on the boundary of the element instead
  Point<3> computeProjection(const Point<3>&) const;

	//! Overload the << operator to easily print element info (note: define it in class
  // to avoid a forward declaration)

  template <UInt nnodes, UInt MYDIM, UInt NDIM>
  friend std::ostream& operator<<(std::ostream& os, const Element<nnodes, MYDIM, NDIM>& el);

private:
  // Data members
	elementPoints points_;
  // A matrix encoding a linear transformation from barycentric coordinates to
  // standard coordinates (modulo a translation)
	Eigen::Matrix<Real,3,2> M_J_;
  // A matrix which is the pseudoinverse of M_J_
	Eigen::Matrix<Real,2,3> M_invJ_;
  // A Real storing the area/volume of the element
  Real element_measure;

  // A member initializing M_J_, M_invJ_ and element_measure at construction
  void computeProperties();

  Real evaluate_point(const Eigen::Matrix<Real, 3, 1>&, const Eigen::Matrix<Real,NNODES,1>&) const;
};

//! Partial template specialization for line elements
//!  This class implements an Edge as an objects composed by two or three nodes, embedded in a 2-dimensional space
/*!
 *  The first two nodes represent the vertices, the other one the internal node,
 *  following this enumeration:
 *
 *
 * 		 *_____*_____*
 * 		1	   3	  2
*/

template <UInt NNODES>
class Element<NNODES, 1, 2> : public Identifier {

    static_assert(NNODES==2 || NNODES==3,
    "ERROR! TRYING TO INSTANTIATE LINE ELEMENT WITH WRONG NUMBER OF NODES! See mesh_objects.h");

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using elementPoints = std::array<Point<2>, NNODES>;
    using iterator = typename elementPoints::iterator;
    using const_iterator = typename elementPoints::const_iterator;

    static constexpr UInt numVertices=2;
    static constexpr UInt numSides=2;
    static constexpr UInt myDim=1;
    static constexpr UInt nDim=2;


    // Note: these don't really mean anything, they're just here for compatibility
    // with the adtree implementation
    static constexpr UInt dp() {return 2;} //  physical dim         (?!)
    static constexpr UInt dt() {return 4;} //  2*physical dim       (?!)
    static constexpr UInt coordsize() {return 4;} // nvertex*ndim   (?!)

    //! This constructor creates an "empty" Element, with an Id Not Valid
    Element()=default;

    //! This constructor creates an Element, given its Id and an std array with the Points
    Element(UInt id, const elementPoints& points) :
            Identifier(id), points_(points) {computeProperties();}

    //! Overloaded subscript operator
    // Note: only the const version is available because any change in points_
    // would require a call to computeProperties() to keep the element in a valid state
    const Point<2>& operator[](UInt i) const {return points_[i];}

    //! Define begin and end iterators (this also gives "ranged for" for free)
    // Note: only the const version is available because any change in points_
    // would require a call to computeProperties() to keep the element in a valid state
    const_iterator begin() const {return points_.begin();}
    const_iterator end() const {return points_.end();}

    const Eigen::Matrix<Real,2,1>& getM_J() const {return M_J_;}
    const Eigen::Matrix<Real,1,2>& getM_invJ() const {return M_invJ_;} //this is actually the pseudoinverse!

    //! A member function that computes the barycentric coordinates of a given point
    Eigen::Matrix<Real,2,1> getBaryCoordinates(const Point<2>& point) const;

    //! A member function that tests if a Point is located inside the Element
    bool isPointInside(const Point<2>& point) const;

    //! Some members returning the area/volume of the element
    // Note: all three are available for convenience and compatibility reasons with existing code
    // Add ->Real getLength()const{return element_measure;}  ?
    Real getMeasure() const {return element_measure;}
    Real getArea() const {return getMeasure();}
    Real getVolume() const {return getMeasure();}

    //! A member to evaluate functions at a point inside the element
    Real evaluate_point(const Point<2>&, const Eigen::Matrix<Real,NNODES,1>&) const;
    //! A member to evaluate integrals on the element
    Real integrate(const Eigen::Matrix<Real,NNODES,1>&) const;

    // A member function that projects a 2D point (XY coordinates!) onto the element
    // Note: if the projection lies outside the element the function returns
    // the closest point on the boundary of the element instead
    Point<2> computeProjection(const Point<2>&) const;

    //! Overload the << operator to easily print element info (note: define it in class
    // to avoid a forward declaration)

    template <UInt nnodes, UInt MYDIM, UInt NDIM>
    friend std::ostream& operator<<(std::ostream& os, const Element<nnodes, MYDIM, NDIM>& el);

private:
    // Data members
    elementPoints points_;
    // A matrix encoding a linear transformation from barycentric coordinates to
    // standard coordinates (modulo a translation)
    Eigen::Matrix<Real,2,1> M_J_;
    // A matrix which is the pseudoinverse of M_J_
    Eigen::Matrix<Real,1,2> M_invJ_;
    // A Real storing the area/volume of the element
    Real element_measure;

    // A member initializing M_J_, M_invJ_ and element_measure at construction
    void computeProperties();

    Real evaluate_point(const Eigen::Matrix<Real, 2, 1>&, const Eigen::Matrix<Real,NNODES,1>&) const;
};


#include "Mesh_Objects_imp.h"
#endif

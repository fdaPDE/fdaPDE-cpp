#ifndef MESH_H_
#define MESH_H_

// Needed for std::enable_if
#include <type_traits>

// Note: how_many_nodes constexpr function is defined in mesh_objects.h
// Also Point and Element
#include "mesh_objects.h"


template <UInt ORDER, UInt mydim, UInt ndim>
class MeshHandler{
  static_assert((ORDER==1 || ORDER==2) &&
								(mydim==2 || mydim==3) &&
								 mydim <= ndim,
								 "ERROR! TRYING TO INSTANTIATE MESH_HANDLER WITH WRONG NUMBER OF NODES AND/OR DIMENSIONS! See mesh.h");
public:
	using meshElement = Element<how_many_nodes(ORDER,mydim),mydim,ndim>;

  // A constructor
  // Note: this constructor is never actually used in practice outside of debug
	MeshHandler(Real* points, UInt* sides, UInt* elements, UInt* neighbors, UInt num_nodes, UInt num_sides, UInt num_elements):
			points_(points), sides_(sides), elements_(elements), neighbors_(neighbors), num_nodes_(num_nodes), num_sides_(num_sides), num_elements_(num_elements) {};

  //! A constructor.
    /*!
      * This constructor permits the initialization of the mesh from an R object
      * constructed with the TriLibrary (our R wrapper for the Triangle library)
      * in 2D (in 2.5D and 3D R functions can produce a compatible object if the
      * triangulation is already available)
    */

	#ifdef R_VERSION_
	MeshHandler(SEXP Rmesh);
	#endif

	//! A member returning the number of nodes in the mesh
	UInt num_nodes() const {return num_nodes_;}

	//! A member returning the number of elements in the mesh
	UInt num_elements() const {return num_elements_;}

	//! A member returning the number of distinct sides
  // (edges for mydim=2,faces for mydim=3) in the mesh
  // All three implemented for convenience
	UInt num_sides() const {return num_sides_;}
  UInt num_edges() const {return num_sides();}
  UInt num_faces() const {return num_sides();}


  const Real& nodes(const UInt i, const UInt j) const {return points_[i+j*num_nodes_];}

  const UInt& elements(const UInt i, const UInt j) const {return elements_[i+j*num_elements_];}

  const UInt& sides(const UInt i, const UInt j) const {return sides_[i+j*num_sides_];}
  const UInt& edges(const UInt i, const UInt j) const {return sides(i,j);}
  const UInt& faces(const UInt i, const UInt j) const {return sides(i,j);}

  const UInt& neighbors(const UInt i, const UInt j) const {return neighbors_[i+j*num_elements_];}



	//! A member returning the ndim-dimensional Point with the specified ID
	Point<ndim> getPoint(UInt id) const;

	//! A member returning an Element with the specified ID
	meshElement getElement(UInt id) const;

  //! A member returning the area/volume of a given element of the mesh
  Real elementMeasure(UInt id) const {return getElement(id).getMeasure();}

  //! A member returning the "number"-th neighbor of element id_element,
  // i.e. the neighbor opposite the "number"-th vertex of the element id_element
  // Note: this function returns an empty element if the neighbor lies outside the boundary
	meshElement getNeighbors(UInt id_element, UInt number) const;

	void printPoints(std::ostream &);
	void printElements(std::ostream &);
	void printNeighbors(std::ostream &);

	//! A normal member returning the element on which a point is located
		/*!
		 * This method implements a simply research between all the elements of the mesh
		*/
	meshElement findLocationNaive(const Point<ndim>&) const;

  //! A normal member returning the element on which a point is located
    /*!
    * This method implements a Visibility Walk Algorithm (further details in: Walking in a triangulation, Devillers et al)
    */
  // Note: this method is guaranteed to work only on convex domains
  // It does not work for manifold data!
  // We make sure that it is not available for manifold data at compile time using enable_if
  template <UInt m=mydim, UInt n=ndim>
  typename std::enable_if<n==m && n==ndim && m==mydim, meshElement>::type
  findLocationWalking(const Point<ndim>&, const meshElement&) const;


private:
	#ifdef R_VERSION_
	SEXP mesh_;
	#endif
	Real *points_;
	UInt *sides_;
	UInt *elements_;
	UInt *neighbors_;

	UInt num_nodes_, num_sides_, num_elements_;

};

#include "mesh_imp.h"

#endif

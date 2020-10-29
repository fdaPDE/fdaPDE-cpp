#ifndef __MESH_H__
#define __MESH_H__

#include <set>
#include <memory>

#include "../../FdaPDE.h"
// Note: how_many_nodes constexpr function is defined in mesh_objects.h
// Also Point and Element
#include "Mesh_Objects.h"
#include "AD_Tree.h"

template <UInt ORDER, UInt mydim, UInt ndim>
class MeshHandler{
  static_assert((ORDER==1 || ORDER==2) &&
                (mydim==2 || mydim==3) &&
                 mydim <= ndim,
                 "ERROR! TRYING TO INSTANTIATE MESH_HANDLER WITH WRONG NUMBER OF NODES AND/OR DIMENSIONS! See mesh.h");
public:
  using meshElement = Element<how_many_nodes(ORDER,mydim),mydim,ndim>;


  //! A constructor.
    /*!
      * This constructor permits the initialization of the mesh from an R object
      * constructed with the TriLibrary (our R wrapper for the Triangle library)
      * in 2D (in 2.5D and 3D R functions can produce a compatible object if the
      * triangulation is already available)
    */

  MeshHandler(SEXP Rmesh, UInt search=1);

  MeshHandler(const MeshHandler&) = delete;
  MeshHandler(MeshHandler&&) = delete;
  MeshHandler& operator=(const MeshHandler&) = delete;
  MeshHandler& operator=(MeshHandler&&) = delete;

  //! A member returning the number of nodes in the mesh
  UInt num_nodes() const {return points_.nrows();}

  //! A member returning the number of elements in the mesh
  UInt num_elements() const {return elements_.nrows();}

  //! A member returning the number of distinct sides
  // (edges for mydim=2,faces for mydim=3) in the mesh
  // All three implemented for convenience
  UInt num_edges() const {return sides_.nrows();}
  UInt num_faces() const {return sides_.nrows();}


  const Real& nodes(const UInt i, const UInt j) const {return points_(i,j);}

  const UInt& elements(const UInt i, const UInt j) const {return elements_(i,j);}

  const UInt& edges(const UInt i, const UInt j) const {return sides_(i,j);}
  const UInt& faces(const UInt i, const UInt j) const {return sides_(i,j);}

  const UInt& neighbors(const UInt i, const UInt j) const {return neighbors_(i,j);}

  //! A member returning the ndim-dimensional Point with the specified ID
  Point<ndim> getPoint(const UInt id) const;

  //! A member returning an Element with the specified ID
  meshElement getElement(const UInt id) const;

  //! A member returning the area/volume of a given element of the mesh
  Real elementMeasure(const UInt id) const {return getElement(id).getMeasure();}

  UInt getSearch() const {return search_;}

  bool hasTree() const {return tree_ptr_.get()!=nullptr;}
  const ADTree<meshElement>& getTree() const {return *tree_ptr_;}
  //! A member returning the "number"-th neighbor of element id_element,
  // i.e. the neighbor opposite the "number"-th vertex of the element id_element
  // Note: this function returns an empty element if the neighbor lies outside the boundary
  meshElement getNeighbors(const UInt id_element, const UInt number) const;

  void printPoints(std::ostream &) const;
  void printElements(std::ostream &) const;
  void printNeighbors(std::ostream &) const;
  void printTree(std::ostream & out) const;


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
  // We make sure that it is not available for manifold data at compile time using static_assert

  meshElement findLocationWalking(const Point<ndim>&, const meshElement&) const;

  meshElement findLocationTree(const Point<ndim>& point) const;

private:

  const RNumericMatrix points_;
  const RIntegerMatrix sides_;
  const RIntegerMatrix elements_;
  const RIntegerMatrix neighbors_;

  const UInt search_;

  std::unique_ptr<const ADTree<meshElement> > tree_ptr_;

};

#include "Mesh_imp.h"

#endif

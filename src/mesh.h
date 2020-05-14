#ifndef MESH_H_
#define MESH_H_

#include "fdaPDE.h"
#include "mesh_objects.h"
#include "bounding_box.h"
#include "tree_header.h"
#include "domain.h"
#include "treenode.h"
#include "exception_handling.h"
#include "adtree.h"
#include <math.h>

using std::vector;

template <UInt ORDER,UInt mydim, UInt ndim>
class MeshHandler{
};

//!  2D MESH:
//!  This class gives an object-oriented reading interface to the output of the library Triangle (Jonathan Richard Shewchuk).
/*!
 * The template parameters specify the order of its elements.
 * The aim of this class is to do not introduce any initialization overhead,
 * beacuse it will be called many time during the execution of a R script
*/
template <UInt ORDER>
class MeshHandler<ORDER,2,2> {
public:
  typedef int UInt;
  //! A constructor.
    /*!
      * The constructor permits the initialization of the mesh from an R object
      * constructed with the TriLibrary (our R wrapper for the Triangle library)
    */

  MeshHandler(Real* points, UInt* edges, UInt* triangles, UInt* neighbors, UInt num_nodes, UInt num_edges, UInt num_triangles):
    points_(points), edges_(edges), elements_(triangles), neighbors_(neighbors), num_nodes_(num_nodes), num_edges_(num_edges), num_elements_(num_triangles) 
    {
      search_=2;
      ADTree<Element<3*ORDER,2,2>> tmp(points_, elements_, num_nodes_, num_elements_);
      tree_ = tmp;  
    };

  #ifdef R_VERSION_
  MeshHandler(SEXP Rmesh, UInt search_=2); //default search_=tree
  #endif

  ~MeshHandler(){};

  //! A normal member returning an unsigned integer value.
    /*!
      \return The number of nodes in the mesh
    */
    UInt num_nodes() const {return num_nodes_;}

  //! A normal member returning an unsigned integer value.
    /*!
      \return The number of elements in the mesh
    */
    UInt num_elements() const {return num_elements_;}

    //! A normal member returning an unsigned integer value.
    /*!
      \return The number of edges in the mesh
    */
    UInt num_edges() const {return num_edges_;}

    //! A normal member returning a Point
    /*!
     * \param id an Id argument
      \return The point with the specified id
    */
    Point getPoint(Id id);

    //! A normal member returning an Edge
    /*!
     * \param id an Id argument
      \return The edge with the specified id
    */
    Edge getEdge(Id id);

   //! A normal member setting an Element
        /*!
         * \param id an Id argument
          \return The element with order coerent to that of the mesh with the specified id
        */
    //void setElement(Element<3*ORDER,2,2>& tri, Id id) const;

    //! A normal member returning an Element
    /*!
     * \param id an Id argument 
      \return The element with order coerent to that of the mesh with the specified id
    */ 
    Element<3*ORDER,2,2>  getElement(Id id) const;

    //The "number" neighbor of element i is opposite the "number" corner of element i
    //! A normal member returning the Neighbors of a element
    /*!
     * \param id the id of the element
     * \param number the number of the vertex
      \return The element that has as an edge the one opposite to the specified
      vertex
    */
    Element<3*ORDER,2,2> getNeighbors(Id id_element, UInt number) const;

    //! A normal member returning the ADTree
    /*!
     *  \return The ADTree, the nodes contain the index of the triangle in the mesh
    */ 
    const ADTree<Element<3*ORDER,2,2>> &  getTree() const {return tree_;};

    void printPoints(std::ostream & out);
    void printEdges(std::ostream & out);
    void printElements(std::ostream & out);
    void printNeighbors(std::ostream & out);
    void printTree(std::ostream & out);

     //! A normal member returning the element on which a point is located
    /*!
     * This method implements a simply research between all the elements of the mesh
     * \param point the point we want to locate
      \return The element that contains the point
    */
    Element<3*ORDER,2,2> findLocationNaive(Point point) const;

     //! A normal member returning the element on which a point is located
    /*!
     * This method implements a Visibility Walk Algorithm (further details in: Walking in a triangulation, Devillers et al)
     * \param point the point we want to locate
     * \param starting_element Element that specifies the poposed starting points for the walking algorithm
      \return The element that contains the point
    */
    Element<3*ORDER,2,2> findLocationWalking(const Point& point, const Element<3*ORDER,2,2>& starting_element) const;


     //! A normal member returning the triangle on which a point is located
    /*!
     * This method implements a ADTree algorithm
     * \param point the point we want to locate
      \return The triangle that contains the point
    */ 
    Element<3*ORDER,2,2> findLocationTree(const Point& point) const;

    //! A normal member returning the area of an Element
    /*!
     * \param id an Id argument
      \return The volume of the element with the given id
    */
    Real elementMeasure(Id id) const;
  UInt getSearch() const {return search_;};                    

private:
  #ifdef R_VERSION_
  SEXP mesh_;
  #endif
  Real *points_;
  UInt *edges_;
  UInt *elements_;
  UInt *neighbors_;

  UInt *border_edges; //contains the list of edge_id at the border
  UInt num_nodes_, num_edges_, num_elements_;
  UInt search_;
  ADTree<Element<3*ORDER,2,2>> tree_; // adtree associated to the mesh

};


//!  SURFACE MESH:
//!  This class gives an object-oriented reading interface to the mesh object passed from R
/*!
 * The template parameters specify the order of its elements.
*/
template <UInt ORDER>
class MeshHandler<ORDER,2,3> {
public:
  typedef int UInt;
  //! A constructor.
    
    MeshHandler(Real* points, UInt* triangles, UInt num_nodes, UInt num_triangles):
      points_(points), elements_(triangles), num_nodes_(num_nodes), num_elements_(num_triangles) {
        search_=2;
        ADTree<Element<3*ORDER,2,3>> tmp(points_, elements_, num_nodes_, num_elements_);
        tree_ = tmp;  
      };
  
    //! A constructor.
    /*!
      * The constructor permits the initialization of the mesh from a .csv file, useful for
      * debugging purposes
    */
  
    MeshHandler(std::string &filename){
       if(filename.find(".csv") != std::string::npos){
          importfromCSV(filename);
       }
    }
  

    void importfromCSV(std::string &filename);
  
    //! A constructor.
    /*!
      * The constructor permits the initialization of the mesh from an R object
    */
    #ifdef R_VERSION_
  MeshHandler(SEXP Rmesh, UInt search_=2); //default search_=tree
  #endif

  ~MeshHandler(){};

    //! A normal member returning an unsigned integer value.
    /*!
      \return The number of nodes in the mesh
    */
    UInt num_nodes() const {return num_nodes_;}

    //! A normal member returning an unsigned integer value.
    /*!
      \return The number of elements in the mesh
    */
    UInt num_elements() const {return num_elements_;}

    //! A normal member returning a Point
    /*!
     * \param id an Id argument
      \return The point with the specified id
    */
    Point getPoint(Id id);

    //! A normal member returning an Element
    /*!
     * \param id an Id argument
      \return The element with order coerent to that of the mesh with the specified id
    */
    Element<3*ORDER,2,3>  getElement(Id id) const;

    //! A normal member returning the ADTree
    /*!
     *  \return The ADTree, the nodes contain the index of the triangle in the mesh
    */ 
    const ADTree<Element<3*ORDER,2,3>> &  getTree() const {return tree_;};

    void printPoints(std::ostream & out);
    void printElements(std::ostream & out);
   

     //! A normal member returning the element on which a point is located
    /*!
     * This method implements a simply research between all the elements of the mesh
     * \param point the point we want to locate
      \return The element that contains the point
    */
    Element<3*ORDER,2,3> findLocationNaive(Point point) const;

     //! A normal member returning the triangle on which a point is located
    /*!
     * This method implements a ADTree algorithm
     * \param point the point we want to locate
      \return The triangle that contains the point
    */ 
    Element<3*ORDER,2,3> findLocationTree(const Point& point) const;

    //! A normal member returning the area of an Element
    /*!
     * \param id an Id argument
      \return The volume of the element with the given id
    */
    Real elementMeasure(Id id) const;
  UInt getSearch() const {return search_;};                    


private:
  #ifdef R_VERSION_
  SEXP mesh_;
  #endif

    Real *points_;
    UInt *elements_;
    

  UInt num_nodes_, num_elements_;
  UInt search_;
  ADTree<Element<3*ORDER,2,3>> tree_; //adtree associated to the mesh

};


//!  VOLUME MESH:
//!  This class gives an object-oriented reading interface to the mesh object passed from R
/*!
 * The template parameters specify the order of its elements.
*/
template <UInt ORDER>
class MeshHandler<ORDER,3,3> {
public:
  typedef int UInt;
  //! A constructor.
    
    MeshHandler(Real* points, UInt* tetrahedrons, UInt num_nodes, UInt num_tetrahedrons):
      points_(points), elements_(tetrahedrons), num_nodes_(num_nodes), num_elements_(num_tetrahedrons) {
        search_=2;
        ADTree<Element<6*ORDER-2,3,3>> tmp(points_, elements_, num_nodes_, num_elements_);
        tree_ = tmp;  
       };
  
  //! A constructor.
    /*!
      * The constructor permits the initialization of the mesh from an R object
    */
    #ifdef R_VERSION_
  MeshHandler(SEXP Rmesh, UInt search_=2); //default search_=tree
  #endif

  ~MeshHandler(){};

  //! A normal member returning an unsigned integer value.
    /*!
      \return The number of nodes in the mesh
    */
    UInt num_nodes() const {return num_nodes_;}

  //! A normal member returning an unsigned integer value.
    /*!
      \return The number of elements in the mesh
    */
    UInt num_elements() const {return num_elements_;}

    //! A normal member returning a Point
    /*!
     * \param id an Id argument
      \return The point with the specified id
    */
    Point getPoint(Id id);

    //! A normal member returning an Element
    /*!
     * \param id an Id argument
      \return The element with order coerent to that of the mesh with the specified id
    */
    Element<6*ORDER-2,3,3>  getElement(Id id) const;

    //! A normal member returning the ADTree
    /*!
     *  \return The ADTree, the nodes contain the index of the triangle in the mesh
    */ 
    const ADTree<Element<6*ORDER-2,3,3>> &  getTree() const {return tree_;};

    void printPoints(std::ostream & out);
    void printElements(std::ostream & out);
   

     //! A normal member returning the element on which a point is located
    /*!
     * This method implements a simply research between all the elements of the mesh
     * \param point the point we want to locate
      \return The element that contains the point
    */
    Element<6*ORDER-2,3,3> findLocationNaive(Point point) const;

  //! A normal member returning the triangle on which a point is located
    /*!
     * This method implements a ADTree algorithm
     * \param point the point we want to locate
      \return The triangle that contains the point
    */ 
    Element<6*ORDER-2,3,3> findLocationTree(const Point& point) const;

    //! A normal member returning the volume of an Element
    /*!
     * \param id an Id argument
      \return The volume of the element with the given id
    */
    Real elementMeasure(Id id) const;
  UInt getSearch() const {return search_;};                    


private:
  #ifdef R_VERSION_
  SEXP mesh_;
  #endif

  Real *points_;
  UInt *elements_;

  UInt num_nodes_, num_elements_;
  UInt search_;      
  ADTree<Element<6*ORDER-2,3,3>> tree_; //adtree associated to the mesh
};


#include "mesh_imp.h"

#endif
/**
 *	\file treenode.hpp
 *	\author Alessandra Cardani
 *	\author Pigoli Davide
 *  \author Prada Daniele
 */

#ifndef TREENODE_H_
#define TREENODE_H_

#include "fdaPDE.h"
//#include "mesh_objects.h"
//#include "bounding_box.h"

/**	\class TreeNode
 * 	\brief Class defining a tree node.
* \param T: template parameter, is the original shape, in the treenode are stored the bbox of the original object and an index that identify that object
 */
template<UInt ndim>
class TreeNode{
protected:
  // Position of the father node. (It's used by the algorithm for deleting a tree node.)
  //`int father_;

  ///Bounding Box of the object to be stored
  Box<ndim> box_;

  /// Positions of left and right children.
  std::array<int, 2> children_;

  /// The id of the Element that create the Box
  //Be careful! This is not the id of the Treenode but id of Element
  Id id_ = std::numeric_limits<UInt>::max();

public:
  /**	Default constructor.
   *
   *	It's fundamental in creating a vector of TreeNode objects.
   */
  TreeNode() :
    box_(), children_() {}
  /**	Another constructor.
   *
   * 	T is the shape of the id. It is needed for Box constructor. It works with Element or Box.
   */

  TreeNode(Id const id, const Box<ndim>& box) :
    box_(box), id_(id), children_() {}


  // constructor in case there is already tree information
  TreeNode(Id const id, Box<ndim> const & box, std::array<int, 2> const & children) :
    box_(box), id_(id), children_(children) {}


  /**	\brief Sets a child.
   *
   * 	\param[in] flag Index of the child to be set. \n
   * 					If:
   * 					<ul>
   * 					<li> flag = 0, set the left child
   * 					<li> flag = 1, set the right child
   * 					</ul>
   * 	\param[in] child The child node.
   */
  inline void setchild(short int const & flag, int const & child) { children_[flag] = child; }
  /// Returns the father.
  //inline int getfather() const { return father_; }
  /**	\brief Returns a child.
   *
   * 	\param[in] flag Index of the child to be returned. \n
   * 					If:
   * 					<ul>
   * 					<li> flag = 0, return the left child
   * 					<li> flag = 1, return the right child
   * 					</ul>
   */
  inline int getchild(short int const & flag) const { return children_[flag]; }
  ///	Sets the id stored in the node.
  inline void setid(Id id) { id_=id; }
  /// Gets id stored in the node.
  inline Id getid() const { return id_; }
  /// Returns a reference to box_.
  inline Box<ndim> & box() { return box_; }
  inline Box<ndim> const & box() const { return box_; }

  ///print information about the treenode
  void print(std::ostream & out) const;
};

#include "treenode_imp.h"

#endif /* TREENODE_H_ */

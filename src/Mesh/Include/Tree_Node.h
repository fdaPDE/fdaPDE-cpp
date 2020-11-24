/**
 *	\file treenode.hpp
 *	\author Alessandra Cardani
 *	\author Pigoli Davide
 *  \author Prada Daniele
 */

#ifndef __TREE_NODE_H__
#define __TREE_NODE_H__

#include "../../FdaPDE.h"
//#include "mesh_objects.h"
//#include "bounding_box.h"

/**	\class TreeNode
 * 	\brief Class defining a tree node.
* \param T: template parameter, is the original shape, in the treenode are stored the bbox of the original object and an index that identify that object
 */
template<class T>
class TreeNode{
protected:
  // Position of the father node. (It's used by the algorithm for deleting a tree node.)
  //`int father_;

  ///Bounding Box of the object to be stored
  Box<T::dp()> box_;

  /// Positions of left and right children.
  int children_[2];

  /// The id of the Element that create the Box
  //Be careful! This is not the id of the Treenode but id of Element
  Id id_;

  //dò per scontato che ci sia un oggetto mesh che contenga tutte le forme (salvate in qualche modo) e che le identifichi attraverso un id di tipo Uint!
  //generalizzando si potrebbe usare un puntatore a Shape (parametro template della forma generica, può essere un triangolo, o l'id stesso se si vuole tornare al caso precedente), in teoria poi non devo distruggere la memoria perchè la forma è salvata in una struttra mesh che deve rimanere inalterata
  //Shape * id_;
  //un'altra possibilità è salvare la forma anche nella struttura mesh con un shared_ptr e mettere anche qui uno shared_ptr<Shape>


public:
  /**	Default constructor.
   *
   *	It's fundamental in creating a vector of TreeNode objects.
   */
  TreeNode(): box_() { //father_(0),
    children_[0] = 0;
    children_[1] = 0;
    id_ = std::numeric_limits<UInt>::max();
  }
  /**	Another constructor.
   *
   * 	T is the shape of the id. It is needed for Box constructor. It works with Element or Box.
   */
  TreeNode(Id const id, T shape): box_(shape) { //father_(0),
    children_[0] = 0;
    children_[1] = 0;
    id_ = id;
  }

  TreeNode(Id const id, const Box<T::dp()>& shape): box_(shape) { //father_(0),
    children_[0] = 0;
    children_[1] = 0;
    id_ = id;
  }


  // constructor in case there is already tree information
  TreeNode(Box<T::dp()> const & box, Id const & id, int const & left_child, int const & right_child):
    box_(box), id_(id) {
    children_[0] = left_child;
    children_[1] = right_child;
    //father_(0)

  }
  /// Sets the father.
  //inline void setfather(int const & ifth) { father_ = ifth; }

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
  ///	Sets the coordinates of the bbox stored in the node.
  inline void setcoords(std::vector<Real> const & data) { box_.set(data); }
  /// Gets the i-th coordinate of the bounding box of the object stored in the node.
  inline Real getcoord(int const & i) const { return box_[i]; }
  ///	Sets the id stored in the node.
  inline void setid(Id id) { id_=id; }
  /// Gets id stored in the node.
  inline Id getid() const { return id_; }
  /// Returns a reference to box_.
  inline Box<T::dp()> & getbox() { return box_; }
  ///print information about the treenode
  void print(std::ostream & out) const;
};

#include "Tree_Node_imp.h"

#endif /* TREENODE_H_ */

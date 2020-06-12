/**
 *	\file adtree.h
 *	\author Cardani Alessandra
 *	\author Pigoli Davide
 *	\author Prada Daniele
 */

#ifndef ADTREE_H_
#define ADTREE_H_



#include "fdaPDE.h"
#include "tree_header.h"
#include "mesh_objects.h"
#include "bounding_box.h"
#include "domain.h"
#include "treenode.h"
#include "exception_handling.h"


/**	\class ADTree
 * 	\brief Alternating binary range searching tree.
 *	\param Shape: template parameter, the original shape
 */
template<UInt ndim>
class ADTree {
protected:
  /** The header.
   *
   *  It contains general information about the tree.
   */
  TreeHeader<ndim> header_;
  /// Vector of tree nodes.
  std::vector<TreeNode<ndim> > data_;
  /** \brief Adds a point to the tree.
   * 	It throws:
   * 	<ul>
   * 	<li> a TreeDomainError exception if the point is out of domain;
   * 	<li> a TreeAlloc exception if there is no more space in the tree to add the node;
   * 	<li> a LevRuntimeError if you exceed the limit set for the tree levels due to the inclusion of the node.
   * 	</ul>
   */
  int adtrb(Id shapeid, std::vector<Real> const & coords);
  /// Handles a TreeDomainError exception.
  int handledomerr(Id shapeid, std::vector<Real> const & coords);
  /// Handles a TreeAlloc exception.
  int handletreealloc(Id shapeid, std::vector<Real> const & coords);
  /// Handles a LevRuntimeError exception.
  int handleleverr(Id shapeid, std::vector<Real> const & coords);
  /** Searches dimension associated to a given level.
   *
   * 	\param[in] lev The given level.
   * 	\param[in] dim The number of dimensions used for the search.
   */
  inline int searchdim(int const & lev, int const & dim) const {
    return (lev % dim);
  }
  /** Finds delta associated to division at a given level.
   *
   * 	\param[in] lev The given level.
   * 	\param[in] dim The number of dimensions used for the search.
   */
  inline double delta(int const & lev, int const & dim) const {
    return std::pow(0.5, int(lev/dim)+1);
  }
public:
  /**	A default constructor.
   *
   *	It initializes the tree header and reserve a suitable number of memory
   *	locations to store tree nodes. It doesn't fill the tree.
   */
  ADTree() : header_() {
    /*
     * The first element in the tree nodes vector is the head.
     * It stores the address of the tree root (i.e. the first node in the tree).
     * If it stores 0 it means that the tree is empty.
     */
   	data_.reserve(header_.gettreeloc()+1);

    // Id, obj are arbitrary parameters. Remember that data_[0] is the head, not a tree node.
  	data_.emplace_back();
	};

  /**	A base constructor.
   *
   *	It initializes the tree header and reserve a suitable number of memory
   *	locations to store tree nodes. It doesn't fill the tree.
   */
  ADTree(TreeHeader<ndim> const & header) : header_(header) {
    data_.reserve(header_.gettreeloc()+1);
    data_.emplace_back();
  }

  // constructor in case there is already tree information
  ADTree(TreeHeader<ndim> const & header, std::vector<TreeNode<ndim> > const & data) : header_(header), data_(data) {}

  /** It fills all the locations of the tree. Object's coordinates are stored to perform searching operations.
   * 	See mesh_handler to verify what points and triangle must contain.
   */
  template<UInt ORDER, UInt mydim>
  ADTree(const MeshHandler<ORDER,mydim,ndim>&);

  /// Returns a reference to the tree header.
  inline TreeHeader<ndim>& gettreeheader() const { return header_; }
  /** Adds a node to the tree.
   * 	It calls the handlers of the exceptions that can be thrown by adtrb().
   *
   * 	\param[in] coords Coordinates of the point.
   *
   *	The location of the current node in the tree is returned.
   */
  int addtreenode(Id shapeid, std::vector<Real> const & coords);
  /** Gets out the informations stored at a given node.
   *
   * 	\param[in] loc Location of the searched node.
   * 	\param[out] coord Bounding box coordinates of the object stored in the loc-th location.
   *  \param[out] id Id of searched node.
   */
  inline void gettri(int const & loc, std::vector<Real> & coord, Id & id);
  /** Gets out the node stored at a given location.
   *
   * 	\param[in] loc Location of the searched node.
   */
  inline TreeNode<ndim> gettreenode(int const & loc) const {return data_[loc];}

  /** Finds all points or (bounding) boxes that intersect a given box.
   *
   * 	\param[in] region Box where searching described by the representative point obtained through a corner transformation.
   * 	\param[out] found Indices of the found elements.
   *
   * 	This function returns true if it has completed successfully, false otherwise.
   */
  bool search(std::vector<Real> const & region, std::set<int> & found) const;

  /// Gets the the Id of the original object of the p-th treenode.
  inline Id pointId (int const & p) const { return data_[p].getid(); }
  /// Outputs informations contained in the tree header.
  template<UInt NDIM>
  friend std::ostream & operator<<(std::ostream & ostr, ADTree<NDIM> const & myadt);

};

#include "adtree_imp.h"

#endif /* ADTREE_H_ */

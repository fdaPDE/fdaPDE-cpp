/**
 *	\file tree_header.h
 *	\author Cardani Alessandra
 *	\author Pigoli Davide
 *	\author Prada Daniele
 */

#ifndef TREE_HEADER_H_
#define TREE_HEADER_H_

#include "../../FdaPDE.h"
#include "Mesh_Objects.h"
#include "Bounding_Box.h"
#include "Domain.h"
#include "Tree_Node.h"
#include "Exception_Handling.h"



/**	\class Header
 * 	\brief It contains general information about the tree.
* 	\param	T The template parameter expresses the original objects for the construction of the tree, it is used to store Domain<T>,
				T::dp() gives the physical dimension of the object and T::dt() gives the search dimension \n
 * 				The tree always store Box<NDIMP> and the index to find the original object (Triangle...).
 				If we want to store a point, we have to create a box from point.
 				(for example repeating the same coordinates for Min_point and Max_point)
 				In the simplest case, or to change the structure adding the poin case everywhere (treenode...)
 *
 * 	Default copy constructor and default destructor work fine.
 */
template<class T>
class TreeHeader {
protected:
	// Tree memory locations.
	int tree_loc_;

	// Tree levels.
	int tree_lev_;

	// Number of physical space dimensions (typically 2 or 3).
	int ndimp_;

	// Number of dimensions used for the search (2*ndimp because we use box).
	int ndimt_;

	//	Number of logical locations currently used in the tree. Initialized to 0.
	int nele_;

	/**	@name Tree indices
	 *	The use of iava and iend avoids the necessity of initializing the stack of available locations.
	 * In fact, the stack contains only locations that have been previously deleted from the tree.
	 */
	//@{
	/**	Next available location in the stack. Initialized to 1.
	 *
	 *	The stack of available nodes is a linked list with all nodes that have been previously deleted. If iava=0 the stack is empty.
	 */
	int iava_;

	/**	Next available location in the yet not allocated part of the tree ("tree free store"). Initialized to 1.
	 *
	 *	The "tree free store" is the yet unassigned portion of the vector storing the tree.
	 *	It effectively acts as a free storage area. If iend_=1, the free storage is at his maximum size. iend never decreases.
	 *	It may increase during insertions and will always remain unaltered after deletions.
	 *	It isn't equal to nele.
	 */
	int iend_;
	//@}

	// Tree's domain.
	Domain<T> tree_domain_;

	/**	A protected constructor.
	 *
	 *	\param[in] ntree Tree dimension needed.
	 *	\param[in] d Tree's domain.
	 *
	 *	Public function createtreeheader must be used to create an object of Header class.
	 *	This avoids the creation of a tree with more memory locations than a fixed limit.
	 */
	TreeHeader(int const & ntree, Domain<T> const & d);
	/// Tries to set the number of tree memory locations (throws a LocLengthError exception if nt is out of range).
	void stml(int const & nt);
public:
	/** Default constructor.
	 *
	 * 	It's fundamental in creating an ADTree object from a MeshFile::ff2dmesh or a MeshFile::ff3dmesh object.
	 */
	TreeHeader():tree_loc_(0), tree_lev_(0), ndimp_(T::dp()), ndimt_(T::dt()), nele_(0), iava_(1), iend_(1), tree_domain_(){}

	// constructor in case there is already tree information
	TreeHeader(int const & tree_loc, int const & tree_lev, int const & ndimp, int const & ndimt,
		int const & nele, int const & iava, int const & iend, Domain<T> const & tree_domain):
			tree_loc_(tree_loc), tree_lev_(tree_lev), ndimp_(ndimp), ndimt_(ndimt),
			nele_(nele), iava_(iava), iend_(iend), tree_domain_(tree_domain){}

	/// Gets the number of tree memory locations.
	inline int gettreeloc() const { return tree_loc_; }
	/// Sets the number of tree memory locations (handles a LocLengthError exception).
	void settreeloc(int const & nt);
	/// Gets the number of tree levels.
	inline int gettreelev() const { return tree_lev_; }
	/// Sets the number of tree levels.
	inline void settreelev(int const & nl) { tree_lev_ = nl; }
	/// Gets the number of physical space dimension.
	inline int getndimp() const { return ndimp_; }
	/// Gets the number of dimensions used for the search.
	inline int getndimt() const { return ndimt_; }
	///	Gets the number of logical locations currently used in the tree.
	inline int getnele() const { return nele_; }
	/// Sets the number of logical locations currently used in the tree.
	inline void setnele(int const & ne) { nele_ = ne; }
	/// Gets the next available location in the stack.
	inline int getiava() const { return iava_; }
	/// Sets the next available location in the stack.
	inline void setiava(int const & ia) { iava_ = ia; }
	/// Gets the next available location in the tree free store.
	inline int getiend() const { return iend_; }
	/// Sets the next available location in the tree free store.
	inline void setiend(int const & ie) { iend_ = ie; }
	/// Gets the i-th coordinate of the origin of the tree's (domain) bounding box. --> queste falle friend
	inline double domainorig(int const & i) const { return tree_domain_.orig(i); }
	/// Gets the i-th scaling factor of the tree's bounding box.
	inline double domainscal(int const & i) const { return tree_domain_.scal(i); }
	/**	Output operator.
	 *
	 *	It prints out:
	 *	- number of tree memory locations;
	 *	- number of tree levels;
	 *	- number of physical space dimensions;
	 *	- number of pieces of information carried by the tree;
	 *	- number of dimensions used for the search;
	 *	- number of logical locations currently used in the tree;
	 *	- tree domain.
	 */
	template<class S>
	friend std::ostream & operator<<(std::ostream &, TreeHeader<S> const &);

	/**	\fn Header<S> createtreeheader(int const & nt, int const & nk, Domain<S> const & d)
	 *	\brief Creates a tree header handling the exception condition in which tree maximum dimension is exceeded.
	 *	\param[in] nt Tree dimension needed.
	 *	\param[in] nk Number of pieces of information carried by the tree.
	 *	\param[in] d Tree's domain.
	 */
	template<class S>
	friend TreeHeader<S> createtreeheader(int const & nt, Domain<S> const & d);
};

#include "Tree_Header_imp.h"

#endif /* TREE_HEADER_H_ */

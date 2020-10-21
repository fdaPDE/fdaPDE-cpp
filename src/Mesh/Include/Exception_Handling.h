/**
 *	\file exception_handling.h
 *	\author Cardani Alessandra
 *	\author Pigoli Davide
 *	\author Prada Daniele
 */

#ifndef __EXCEPTION_HANDLING_H__
#define __EXCEPTION_HANDLING_H__

#include "../../FdaPDE.h"
#include "Mesh_Objects.h"
#include "Bounding_Box.h"

/**	\class TreeException
 *
 * 	All exceptions thrown from the adtree library are derived from the base class TreeException.
 */
template<class Shape>
class TreeException {};

/**	\class TreeLogicError
 *
 * 	This class is used to report a logic error concerning the tree, i.e. an error that, at least
 * 	in theory, could be avoided by the program; for example, by performing additional tests of function arguments.
 */
template<class Shape>
class TreeLogicError: public TreeException<Shape> {};

/**	\class TreeAlloc
 *
 * 	This class is used to report the exception condition in which there is no more space in the tree to store a new node.
 */
template<class Shape>
class TreeAlloc: public TreeLogicError<Shape> {};

/**	\class TreeDomainError
 *
 * 	This class is used to report the exception condition in which an object is out of domain.
 */
template<class Shape>
class TreeDomainError: public TreeLogicError<Shape> {
protected:
	/// Number of logical locations currently used in the tree plus 1.
	int nelep1;
	/// Size of the array storing the coordinates of the object which is out of domain.
	int csize;
	/// Coordinates of the object which is out of domain.
	std::vector<Real> outobj;
public:
	/**	A constructor.
	 *
	 */
	TreeDomainError(int const & np1, int const & size, std::vector<Real> const coord);
	/// Returns number of logical locations currently used in the tree plus 1.
	inline int getnelep1() const {  return nelep1; }
	/**	Output operator.
	 *
	 * 	It outputs coordinates of the object which is out of domain.
	 */
	template<class T>
	friend std::ostream & operator<<(std::ostream & ostr, TreeDomainError<T> const & de);
};

template<class Shape>
TreeDomainError<Shape>::TreeDomainError(int const & np1, int const & size, std::vector<Real> const coord):
	nelep1(np1), csize(size) {
	outobj = coord;
}

template<class Shape>
std::ostream & operator<<(std::ostream & ostr, TreeDomainError<Shape> const & de) {
	for(int i = 0; i < de.csize; ++i) {
		ostr << "Coordinate " << i+1 << ": " << de.outobj[i] << std::endl;
	}

	return ostr;
}


/**	\class TreeLengthError
 *
 *  This class is used to report an attempt to do something that exceeds a maximum allowable size
 *  concerning the tree.
 */
template<class Shape>
class TreeLengthError: public TreeLogicError<Shape> {};

/**	\class LocLengthError
 * 	\brief This class is used to report an attempt to build a tree with more memory locations than a fixed limit.
 */
template<class Shape>
class LocLengthError: public TreeLengthError<Shape> {
protected:
	/// Maximum tree memory locations available.
	int max_tree_loc;
	/// Tree memory locations needed.
	int tree_loc;
public:
	/**	A constructor.
	 *
	 * 	Initialize the number of the required tree locations.
	 */
	LocLengthError(int const & mtree, int const & ntree): max_tree_loc(mtree), tree_loc(ntree) {}
	/// Gets the maximum number of tree memory locations available.
	inline int const getmaxtreeloc() const { return max_tree_loc; }
	/// Gets the number of tree memory locations needed.
	inline int const gettreeloc() const { return tree_loc; }
};

/**	\class TreeRuntimeError
 *

 *	Exceptions derived from runtime_error are provided to report events that are beyond the
 *	scope of a program and are not easily avoidable.
 */
template<class Shape>
class TreeRuntimeError: public TreeException<Shape> {};

/**	\class LevRuntimeError
 * 	\brief This class is used to report an attempt to build a tree with more levels than a fixed limit.
 */
template<class Shape>
class LevRuntimeError: public TreeRuntimeError<Shape> {
protected:
	/// Maximum tree levels.
	static int max_tree_lev;
public:
	/// Sets the maximum number of tree levels.
	inline static void setmaxtreelev(int const & mlev) { max_tree_lev = mlev; }
	/// Gets the maximum number of tree levels.
	inline static int getmaxtreelev() { return max_tree_lev; }
};

template<class Shape>
int LevRuntimeError<Shape>::max_tree_lev = 1.e3;


#endif /* EXCEPTION_HANDLING_H_ */

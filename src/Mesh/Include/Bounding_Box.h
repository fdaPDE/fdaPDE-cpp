
#ifndef __BOUNDING_BOX_H__
#define __BOUNDING_BOX_H__


#include "../../FdaPDE.h"
#include "Mesh_Objects.h"

/**	\class Box
 * 	\brief Policy used when boxes have to be stored in the tree.
 */
/**
 *	\template parameter NDIMP is the physical dimension of Box: 2 -> 2D, 3 -> 3D
*/
template<int NDIMP>
class Box {
protected:
	/** A vector of rectangle corner coordinates.
	 * 	First NDIMP values are the coordinates of the rectangle corner with minimum coordinates,
	 *  followed by the coordinates of the opposite corner. (2D: xmin, ymin, xmax, ymax)
	 */
	std::vector<Real> x_;

public:
	/**	Default constructor.
	 *
	 *	It's fundamental in creating a vector of Box objects.
	 *  Sets coordinate values to zero.
	 */
	Box();

	/**	Another constructor.
	 *
	 *	\param[in] coord A vector of input coordinates.
	 *					 Make sure that its size is equal to the size of protected member x_.
	 */
	// constructor in case there is already tree information
	Box(std::vector<Real> const & coord);


	template <UInt NNODES,int NDIME,int NDIMPP>
	Box(Element<NNODES,NDIME,NDIMPP> const & element);


	/// Returns the i-th coordinate value.
	inline Real operator[](int const & i) { return x_[i]; }
	inline Real operator[](int const & i) const { return x_[i]; }
	/// Returns the number of physical space dimension.
	inline static constexpr int dp() { return NDIMP; }
	/// Returns the number of dimensions used for the search (2*NDIMP).
	inline static constexpr int dt() { return 2*NDIMP; }
	/// Returns the size of coordinate array. (composed of min and max)
	inline static constexpr int coordsize() { return 2*NDIMP; }
	/** Sets coordinate values.
	 *
	 *	\param[in] data A vector of new coordinates.
	 *					Make sure that its size is equal to the size of protected member x_.
	 */
	void set(std::vector<Real> const & data);

	/** Gets coordinate values.
	 */
	std::vector<Real> get() const {return x_; };

	/** print minimum box point and maximum box point
	*/
	void print(std::ostream & out) const;
};

#include "Bounding_Box_imp.h"

#endif /* BOUNDING_BOX_H_ */


#ifndef BOUNDING_BOX_H_
#define BOUNDING_BOX_H_


#include "fdaPDE.h"
#include "mesh_objects.h"

/**	\class Box
 * 	\brief Policy used when boxes have to be stored in the tree.
 */
/**
 *	\template parameter NDIMP is the physical dimension of Box: 2 -> 2D, 3 -> 3D
*/
template<UInt ndim>
class Box {
protected:
	/** A vector of rectangle corner coordinates.
	 * 	First NDIMP values are the coordinates of the rectangle corner with minimum coordinates,
	 *  followed by the coordinates of the opposite corner. (2D: xmin, ymin, xmax, ymax)
	 */
	Point<ndim> minPoint_;
	Point<ndim> maxPoint_;

public:
	/**	Default constructor.
	 *
	 *	It's fundamental in creating a vector of Box objects.
	 *  Sets coordinate values to zero.
	 */
	Box()=default;

	/**	Another constructor.
	 *
	 *	\param[in] coord A vector of input coordinates.
	 *					 Make sure that its size is equal to the size of protected member x_.
	 */
	// constructor in case there is already tree information
	Box(const Point<ndim>& min, const Point<ndim>& max) :
			minPoint_(min), maxPoint_(max) {}

	Box(std::vector<Points<ndim> > const & points);


	Point<ndim>& minPoint() {return minPoint_;}
	const Point<ndim>& minPoint() const {return minPoint_;}

	Point<ndim>& maxPoint() {return maxPoint_;}
	const Point<ndim>& maxPoint() const {return maxPoint_;}
	
 	/** Sets coordinate values.
	 *
	 *	\param[in] data A vector of new coordinates.
	 *					Make sure that its size is equal to the size of protected member x_.

	/** print minimum box point and maximum box point
	*/
	void print(std::ostream & out) const;
};

#include "bounding_box_imp.h"

#endif /* BOUNDING_BOX_H_ */

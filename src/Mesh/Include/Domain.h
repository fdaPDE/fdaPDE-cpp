/**
 *	\file domain.h
 *	\author Alessandra Cardani
 *	\author Pigoli Davide
 *  \author Prada Daniele
 */

#ifndef __DOMAIN_H__
#define __DOMAIN_H__

#include "../../FdaPDE.h"


/**	\class Domain
 * 	\brief It defines geometric limits of a set of points.
 * template parameter Shape is the original geometric figure (ex: Point, Triangle, Tetrahedron, Box),
 * it's useful to extract the physical dimension
 */
template<class Shape>
class Domain {
protected:
	/// Origin of the object's bounding box = min(coord(*,1:number of points)).
	// ex) 2D: xmin, ymin, xmin, ymin
	std::vector<Real> origin_;
	/// Scaling factors = 1./(max(coord(*,1:number of points)) - min(coord(*,1:number of points))).
	// ex) 2D: xscale, yscale, xscale, yscale
	std::vector<Real> scalingfactors_;
	/// Tolerance being applied to the object's bounding box.
	static Real tolerance_;
	/// Minimum difference between coordinates allowed.
	static Real mindiff_;
public:
	/** Default constructor.
	 */
	Domain();
	/**	Another constructor.
	 *
	 * 	\param[in] coord Point coordinates organized in a matrix (coord[i, j] is the i-th coordinate of the j-th point). \n
	 * 	This matrix is created through a vector of vectors in order to use standard algorithms.
	 *
	 * 	It finds geometric limits of a set of points first. Then adds a tolerance.
	 *  Repeats the limits if the tree dimension is 2 * physical space dimension.
	 *  This is an useful trick. For example, when you have to scale dimensions.
	 */

	// constructor in case there is already tree information
	Domain(std::vector<Real> const & origin, std::vector<Real> const & scalingfactors):origin_(origin), scalingfactors_(scalingfactors) {};
	Domain(std::vector<std::vector<Real> > const & coord);

	/// Sets the tolerance being applied to the object's bounding box.
	inline static void settolerance(Real const & tol) { tolerance_ = tol; }
	/// Gets the tolerance being applied to the object's bounding box.
	inline static Real gettolerance() { return tolerance_; }
	/// Sets the minimum difference between coordinates allowed.
	inline static void setmindiff(Real const & md) { mindiff_ = md; }
	/// Gets the minimum difference between coordinates allowed.
	inline static Real getmindiff() { return mindiff_; }
	/// Gets the i-th coordinate of the origin of the object's bounding box.
	inline Real orig(int const & i) const { return origin_[i]; }
	/// Gets the i-th scaling factor of the object's bounding box.
	inline double scal(int const & i) const { return scalingfactors_[i]; }
	/// Gets the size of origin_ vector
	inline int getoriginsize() { return int(origin_.size()); }
	/// Gets the size of scalingfactors_  vector
	inline int getscalingsize() { return int(scalingfactors_.size()); }


	/**	Output operator.
	 *
	 * 	It outputs the bounding box.
	 */
	template<class S>
	friend std::ostream & operator<<(std::ostream &, Domain<S> const &);
};

#include "Domain_imp.h"

#endif /* DOMAIN_H_ */

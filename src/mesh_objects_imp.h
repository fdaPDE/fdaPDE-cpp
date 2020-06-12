//#include "mesh_objects.hpp"
#ifndef __MESH_OBJECTS_IMP_HPP__
#define __MESH_OBJECTS_IMP_HPP__

// Member functions for class Element

// This function is called to construct elements in 2D and 3D
template <UInt NNODES, UInt mydim, UInt ndim>
void Element<NNODES,mydim,ndim>::computeProperties()
{
	// Note: with Eigen Map one can use non-eigen types as if they were eigen types
	// It is rather useful and no copies are involved
	// A map to const ensures that no change is performed on the underlying object
	using EigenMap2Const_t = Eigen::Map<const Eigen::Matrix<Real,ndim,1> >;

	{
		// Note: the first point in the element is taken as a reference point
		// Note: while this is an arbitrary choice other parts of the code rely on it
		// so any change needs to be considered carefully
		EigenMap2Const_t basePointCoord(&points_[0][0]);
		// The columns of M_J_ are P_i - P_0, i=1,...,mydim
		for (int i=0; i<mydim; ++i)
			M_J_.col(i) = EigenMap2Const_t(&points_[i+1][0])-basePointCoord;
	}
	// NOTE: for small (not bigger than 4x4) matrices eigen directly calculates
	// determinants and inverses, it is very efficient (and much less error prone)!
	M_invJ_ = M_J_.inverse();

	// Area/Volume of the element is the absolute value of det(M_J_)/ndim!
	element_measure = std::abs(M_J_.determinant())/factorial(ndim);

}


template <UInt NNODES, UInt mydim, UInt ndim>
Eigen::Matrix<Real,mydim+1,1> Element<NNODES,mydim,ndim>::getBaryCoordinates(const Point<ndim> &point) const
{
	// Note: with Eigen Map one can use non-eigen types as if they were eigen types
	// It is rather useful and no copies are involved
	// A map to const ensures that no change is performed on the underlying object
	using EigenMap2Const_t = Eigen::Map<const Eigen::Matrix<Real,ndim,1> >;

	Eigen::Matrix<Real,mydim+1,1> lambda;

	// .template is needed! See Eigen documentation regarding
	// the template and typename keywords in C++
	// lambda = M_invJ_ * (P - P_0) where P is the given point and P_0 is the first point of the element
	lambda.template tail<mydim>().noalias() = M_invJ_ * (EigenMap2Const_t(&point[0])-EigenMap2Const_t(&points_[0][0]));

	// The barycentric coordinate corresponding to P_0 can be computed from the others
  lambda(0) = 1 - lambda.template tail<mydim>().sum();

	return lambda;

}

template <UInt NNODES, UInt mydim, UInt ndim>
bool Element<NNODES,mydim,ndim>::isPointInside(const Point<ndim>& point) const
{
	static constexpr Real eps = std::numeric_limits<Real>::epsilon(),
		 tolerance = 10 * eps;

	Eigen::Matrix<Real,mydim+1,1> lambda = getBaryCoordinates(point);

	// Point is inside if all barycentric coords are positive!
	return (-tolerance<=lambda.array()).all();

}

template <UInt NNODES, UInt mydim, UInt ndim>
int Element<NNODES,mydim,ndim>::getPointDirection(const Point<ndim>& point) const
{
	static constexpr Real eps = std::numeric_limits<Real>::epsilon(),
		 tolerance = 10 * eps;

	Eigen::Matrix<Real,mydim+1,1> lambda = getBaryCoordinates(point);

	//Find the minimum barycentric coord
	UInt min_index;
	lambda.minCoeff(&min_index);
	// If the minimum barycentric coord is negative then the point lies in the direction
	// of the opposing edge/face, else the point is inside the element
	return (lambda[min_index] < -tolerance)	? min_index : -1;
}

// Implementation of function evaluation at a point inside the element
template <UInt NNODES, UInt mydim, UInt ndim>
inline Real Element<NNODES,mydim,ndim>::evaluate_point(const Point<ndim>& point, const Eigen::Matrix<Real,NNODES,1>& coefficients) const
{
  return coefficients.dot(getBaryCoordinates(point));
}

// Full specialization for order 2 in 2D
// Note: needs to be declared inline because it is defined in a header file!
// These formulas come from the book "The Finite Element Method: its Basis and Fundamentals" by Zienkiewicz, Taylor and Zhu
template <>
inline Real Element<6,2,2>::evaluate_point(const Point<2>& point, const Eigen::Matrix<Real,6,1>& coefficients) const
{
	Eigen::Matrix<Real,3,1> lambda = getBaryCoordinates(point);
  return coefficients[0] * lambda[0] * (2*lambda[0]-1) +
         coefficients[1] * lambda[1] * (2*lambda[1]-1) +
         coefficients[2] * lambda[2] * (2*lambda[2]-1) +
         coefficients[3] * 4 * lambda[1] * lambda[2] +
         coefficients[4] * 4 * lambda[2] * lambda[0] +
         coefficients[5] * 4 * lambda[0] * lambda[1];
}

// Full specialization for order 2 in 3D
// MEMO: this works assuming edges are ordered like so: (1,2), (1,3), (1,4), (2,3), (3,4), (2,4)
template <>
inline Real Element<10,3,3>::evaluate_point(const Point<3>& point, const Eigen::Matrix<Real,10,1>& coefficients) const
{
 Eigen::Matrix<Real,4,1> lambda = getBaryCoordinates(point);
 return coefficients[0] * lambda[0] * (2*lambda[0]-1) +
        coefficients[1] * lambda[1] * (2*lambda[1]-1) +
        coefficients[2] * lambda[2] * (2*lambda[2]-1) +
        coefficients[3] * lambda[3] * (2*lambda[3]-1) +
        coefficients[4] * 4 * lambda[1] * lambda[0] +
        coefficients[5] * 4 * lambda[2] * lambda[0] +
        coefficients[6] * 4 * lambda[3] * lambda[0] +
        coefficients[7] * 4 * lambda[1] * lambda[2] +
        coefficients[8] * 4 * lambda[2] * lambda[3] +
        coefficients[9] * 4 * lambda[3] * lambda[1];
}

// Implementation of integration on the element
template <UInt NNODES, UInt mydim, UInt ndim>
inline Real Element<NNODES,mydim,ndim>::integrate(const Eigen::Matrix<Real,NNODES,1>& coefficients) const
{
	// A 1-node quadrature formula is enough in the linear case
	return getMeasure() * coefficients.mean();
}

// Full specialization for order 2 in 2D
template <>
inline Real Element<6,2,2>::integrate(const Eigen::Matrix<Real,6,1>& coefficients) const
{
	// A 3-node quadrature formula is needed for the quadratic case on the triangle
	// Nodes: (2/3, 1/6, 1/6), (1/6, 2/3, 1/6), (1/6, 2/3, 1/6)
  // Weights: 1/3, 1/3, 1/3
	// Note: this is actually 9*basis_fun (hence we divide by 3*9=27 in the formula below)
	static constexpr Real basis_fun[]={2, -1, -1, 1, 4, 4, //quadratic basis functions evaluated at node 0
																		 -1, 2, -1, 4, 1, 4, //quadratic basis functions evaluated at node 1
																		 -1, -1, 2, 4, 4, 1  //quadratic basis functions evaluated at node 2
																		};

	return getMeasure()/27 * coefficients.replicate<3,1>().dot(Eigen::Map<const Eigen::Matrix<Real,18,1> >(basis_fun));
}


// Full specialization for order 2 in 3D
// MEMO: this works assuming edges are ordered like so: (1,2), (1,3), (1,4), (2,3), (3,4), (2,4)
template <>
inline Real Element<10,3,3>::integrate(const Eigen::Matrix<Real,10,1>& coefficients) const
{
	// A 4-node quadrature formula is needed for the quadratic case on the tetrahedron
	// See Zienkiewicz, Taylor and Zhu p. 182 for the node coordinates
	// Weights: 1/4, 1/4, 1/4, 1/4
	static constexpr Real basis_fun[]={0.1, -0.1, -0.1, -0.1, 0.323606797749979, 0.323606797749979, 0.323606797749979, 0.076393202250021, 0.076393202250021, 0.076393202250021,
																		 -0.1, 0.1, -0.1, -0.1, 0.323606797749979, 0.076393202250021, 0.076393202250021, 0.323606797749979, 0.076393202250021, 0.323606797749979,
																		 -0.1, -0.1, 0.1, -0.1, 0.076393202250021, 0.323606797749979, 0.076393202250021, 0.323606797749979, 0.323606797749979, 0.076393202250021,
																		 -0.1, -0.1, -0.1, 0.1, 0.076393202250021, 0.076393202250021, 0.323606797749979, 0.076393202250021, 0.323606797749979, 0.323606797749979
																		};

	return getMeasure() * 0.25*coefficients.replicate<4,1>().dot(Eigen::Map<const Eigen::Matrix<Real,40,1> >(basis_fun));
}


// Member functions for class Element (Surface element specialization)

template <UInt NNODES>
void Element<NNODES,2,3>::computeProperties()
{
	// Note: with Eigen Map one can use non-eigen types as if they were eigen types
	// It is rather useful and no copies are involved
	// A map to const ensures that no change is performed on the underlying object
	using EigenMap2Const_t = Eigen::Map<const Eigen::Matrix<Real,3,1> >;

	{
		// Note: the first point in the element is taken as a reference point
		// Note: while this is an arbitrary choice other parts of the code rely on it
		// so any change needs to be considered carefully
		EigenMap2Const_t basePointCoord(&points_[0][0]);
		// The columns of M_J_ are P_i - P_0, i=1,2
		for (int i=0; i<2; ++i)
			M_J_.col(i) = EigenMap2Const_t(&points_[i+1][0])-basePointCoord;
	}

	// NOTE: for small (not bigger than 4x4) matrices eigen directly calculates
	// determinants and inverses, it is very efficient!
	M_invJ_.noalias() = (M_J_.transpose()*M_J_).inverse() * M_J_.transpose();
	// Area of 3D triangle is half the norm of cross product of two sides!
	element_measure = .5 * M_J_.col(0).cross(M_J_.col(1)).norm();
}


template <UInt NNODES>
Eigen::Matrix<Real,3,1> Element<NNODES,2,3>::getBaryCoordinates(const Point<3> &point) const
{
	// Note: with Eigen Map one can use non-eigen types as if they were eigen types
	// It is rather useful and no copies are involved
	// A map to const ensures that no change is performed on the underlying object
	using EigenMap2Const_t = Eigen::Map<const Eigen::Matrix<Real,3,1> >;

	Eigen::Matrix<Real,3,1> lambda;

	// .template is needed! See Eigen documentation regarding
	// the template and typename keywords in C++
	// lambda = M_invJ_ * (P - P_0) where P is the given point and P_0 is the first point of the element
	lambda.template tail<2>().noalias() = M_invJ_ * (EigenMap2Const_t(&point[0])-EigenMap2Const_t(&points_[0][0]));

	// The barycentric coordinate corresponding to P_0 can be computed from the others
  lambda(0) = 1 - lambda.template tail<2>().sum();

	return lambda;

}

// Note: this function is more expensive for manifold data because one must check
// that the point actually lies on the same plane as the 3D triangle
template <UInt NNODES>
bool Element<NNODES,2,3>::isPointInside(const Point<3>& point) const
{
	static constexpr Real eps = std::numeric_limits<Real>::epsilon(),
		 tolerance = 10 * eps;

		 // Note: with Eigen Map one can use non-eigen types as if they were eigen types
	 	// It is rather useful and no copies are involved
	 	// A map to const ensures that no change is performed on the underlying object
	 	using EigenMap2Const_t = Eigen::Map<const Eigen::Matrix<Real,3,1> >;


	 	Eigen::Matrix<Real,3,1> lambda = getBaryCoordinates(point);

		// If the point's projection onto the 3D triangle lies outside the triangle
	 	// (i.e. there is at least one negative barycentric coordinate)
		// then the point does not belong to the triangle
	 	if ((-tolerance>lambda.array()).any())
	 		return false;


	 	Eigen::Matrix<Real,3,3> A;
	 	A << M_J_, (EigenMap2Const_t(&point[0])-EigenMap2Const_t(&points_[0][0]));

	 	Eigen::FullPivHouseholderQR<Eigen::Matrix<Real,3,3> > qr(A);
	 	qr.setThreshold(tolerance);
	 	// Note: the point is inside the element if it lies on the same plane as the 3D triangle
	 	// (i.e. A is rank deficient) AND its projection onto the 3D triangle lies inside the triangle
	 	// (i.e. all barycentric coordinates are positive)
	  	// Note: fullPivHouseholderQr is as fast as ColPivHouseholderQR for such small matrices
	  	// but this is optimized for rank computations (see eigen documentation)
	 	return !qr.isInvertible();

}

template <UInt NNODES>
Point<3> Element<NNODES,2,3>::computeProjection(const Point<3>& point) const
{
	// Note: no need for tolerances here because the projection is continuous
	Eigen::Matrix<Real,3,1> lambda = getBaryCoordinates(point);
	// Convention: (+,+,+) means that all lambda are positive and so on
  // For visual reference: (remember that edges are numbered wrt the opposing node)
  //
	//\ (-,-,+)|
  //  \      |
  //    \    |
  //      \  |
  //        \|
  //         |3
  //         | \
  //         |   \
  // (+,-,+) |     \       (-,+,+)
  //         |       \
  //         |         \
  //         |           \
  //         |  (+,+,+)    \
  //  _____1 |_______________\ 2 _____________
  //         |                 \
  // (+,-,-) |    (+,+,-)        \  (-,+,-)
  //         |                     \

	// If (+,-,-) the projection lies beyond node 1
	// Simply return node 1 (same for the others)
	if(lambda(0)>0 && lambda(1)<0 && lambda(2)<0)
		return points_[0];
	else if (lambda(0)<0 && lambda(1)>0 && lambda(2)<0)
		return points_[1];
	else if (lambda(0)<0 && lambda(1)<0 && lambda(2)>0)
		return points_[2];

	Eigen::Matrix<Real,3,1> coords3D;
	// If (+,+,+) the projection lies inside the element
	// So just convert back to 3D coords
	if(lambda(0)>0 && lambda(1)>0 && lambda(2)>0)
		coords3D = M_J_ * lambda.tail<2>();
	// If (+,+,-) the projection lies beyond edge 3
	// Simply scale it back on the edge and convert
	else if(lambda(0)>0 && lambda(1)>0)
		coords3D = M_J_.col(0)/lambda.head<2>().sum();
	// If (+,-,+) the projection lies beyond edge 2
	else if (lambda(0)>0 && lambda(2)>0)
		coords3D = M_J_.col(1)/(lambda(0)+lambda(2));
	// If (-,+,+) the projection lies beyond edge 1
	else
		coords3D = M_J_ * lambda.tail<2>()/lambda.tail<2>().sum();

	// Translate back by the first point and return
	return Point<3>(coords3D)+=points_[0];
}


// Implementation of function evaluation at a point inside the surface element
template <>
inline Real Element<3,2,3>::evaluate_point(const Point<3>& point, const Eigen::Matrix<Real,3,1>& coefficients) const
{
  return coefficients.dot(getBaryCoordinates(point));
}

 // Full specialization for order 2 in 2.5D
template <>
inline Real Element<6,2,3>::evaluate_point(const Point<3>& point, const Eigen::Matrix<Real,6,1>& coefficients) const
{
	Eigen::Matrix<Real,3,1> lambda = getBaryCoordinates(point);
	return coefficients[0] * lambda[0] * (2*lambda[0]-1) +
         coefficients[1] * lambda[1] * (2*lambda[1]-1) +
         coefficients[2] * lambda[2] * (2*lambda[2]-1) +
         coefficients[3] * 4 * lambda[1]*lambda[2] +
         coefficients[4] * 4 * lambda[2]*lambda[0] +
         coefficients[5] * 4 * lambda[0]*lambda[1];
}


// Implementation of integration on the surface element
template<>
inline Real Element<3,2,3>::integrate(const Eigen::Matrix<Real,3,1>& coefficients) const
{
	// A 1-node quadrature formula is enough in the linear case
	return getMeasure() * coefficients.mean();
}

// Full specialization for order 2 in 2.5D
template <>
inline Real Element<6,2,3>::integrate(const Eigen::Matrix<Real,6,1>& coefficients) const
{
	// A 3-node quadrature formula is needed for the quadratic case on the triangle
	// Nodes: (2/3, 1/6, 1/6), (1/6, 2/3, 1/6), (1/6, 2/3, 1/6)
	// Weights: 1/3, 1/3, 1/3
	// Note: this is actually 9*basis_fun (hence we divide by 3*9=27 in the formula below)

	static constexpr Real basis_fun[]={2, -1, -1, 1, 4, 4,
																		 -1, 2, -1, 4, 1, 4,
																		 -1, -1, 2, 4, 4, 1
																		};

	return getMeasure()/27 * coefficients.replicate<3,1>().dot(Eigen::Map<const Eigen::Matrix<Real,18,1> >(basis_fun));
}

template <UInt nnodes, UInt MYDIM, UInt NDIM>
std::ostream& operator<<(std::ostream& os, const Element<nnodes, MYDIM, NDIM>& el){
	os<< el.getId() <<":";
	for (const auto &p : el)
		os<<" "<<p.getId();
	return os<<std::endl;
}


#endif

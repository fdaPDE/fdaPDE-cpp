// this file contains a set of all required specialization to extend the Mesh code to the manifold case
// (surface elements and linear network elements)

// specialization for surface elements (2.5D manifolds)
template <>
bool Element<2,3>::contains(const SVector<3> &x) const {
  // for surface we should check if the point x belongs to the plane spanned by the 3 vertices of
  // the triangle, if so check if its barycentric coordinates are positive.
  // we can equivalently proceed by computing the orthogonal projection of the point over the plane containing
  // the element and then check if the projection belongs to the element (if the point belongs to the plane
  // the orthogonal projection coincides with the starting point itself) using the barycentric coords of the projection

  SVector<3> a = (coords[1] - coords[0]);
  SVector<3> b = (coords[2] - coords[0]);
  SVector<3> normalVector = a.cross(b).normalized();

  // compute projection of x over the plane
  SVector<3> projection = x - ((x- coords[0]).dot(normalVector))*normalVector;
  
  // get barycentric coordinates of input point
  SVector<4> baryCoord = computeBarycentricCoordinates(projection);

  // use Eigen visitor to check for positiveness of elements
  return (baryCoord.array() >= 0).all();
}

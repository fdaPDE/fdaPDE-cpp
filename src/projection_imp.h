#ifndef _PROJECTION_IMP_H
#define _PROJECTION_IMP_H

#include <cmath>
#include <utility>
#include <limits>


template<UInt ORDER>
UInt projection<ORDER,2,3>::getMaxCoor(const Point& p) const
{
  std::vector<Real> v(3);

  v[0] = std::abs(p[0]);
  v[1] = std::abs(p[1]);
  v[2] = std::abs(p[2]);

  return std::distance(v.begin(), std::max_element(v.begin(), v.end()));
}

template<UInt ORDER>
Real projection<ORDER,2,3>::computeDistance(const Point& p, const Point& q) const
{
  return std::sqrt((p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]) + (p[2]-q[2])*(p[2]-q[2]));
}

template<UInt ORDER>
Real projection<ORDER,2,3>::getAreaTriangle2d(const Point& A, const Point& B, const Point& C) const
{
  Real det = (B[0]-A[0])*(C[1]-A[1])-(C[0]-A[0])*(B[1]-A[1]);
  return (0.5*det);
}

template<UInt ORDER>
std::vector<UInt> projection<ORDER,2,3>::computeNodePatch(UInt id_node) const
{
  std::vector<UInt> patch_node;
  constexpr UInt Nodes = 3*ORDER;

  for(UInt t=0; t<mesh_.num_elements(); t++){
    Element<Nodes, 2, 3> current_element = mesh_.getElement(t);
    bool exit = false;
    for(UInt i=0; i<Nodes && exit==false ; i++){
      if(current_element[i].id() == id_node){
        patch_node.push_back(t);
        exit = true;
      }
    }
  }
  return patch_node;
}

template<UInt ORDER>
std::pair<Point, Real> projection<ORDER,2,3>::project(const Element<3*ORDER,2,3>& triangle ,const Point& P) const 
{
  // only for triangles
  Point A = triangle[0],
        B = triangle[1],
        C = triangle[2];

  // Compute normal to the plane defined by the triangle
  Real nx = (B[1]-A[1])*(C[2]-A[2]) - (B[2]-A[2])*(C[1]-A[1]);
  Real ny = (B[2]-A[2])*(C[0]-A[0]) - (B[0]-A[0])*(C[2]-A[2]);
  Real nz = (B[0]-A[0])*(C[1]-A[1]) - (B[1]-A[1])*(C[0]-A[0]);
  VectorXr n(3); // ndim
  n << nx, ny, nz;
  n.normalize(); 
  nx = n[0];
  ny = n[1];
  nz = n[2];

  // Compute the "d" coefficient of the plane equation ax+by+cz=d, where (a,b,c) is the normale vector n
  Real d = A[0]*nx + A[1]*ny + A[2]*nz;
  // Compute the projection
  Real t = d - (P[0]*nx + P[1]*ny + P[2]*nz);
  // pt = p + t*n;
  Point pt(P[0]+t*nx, P[1]+t*ny, P[2]+t*nz);

  if(triangle.isPointInside(pt)){
    return (std::pair<Point, Real>(pt, computeDistance(P, pt)));
  } 

  // x-y plane
  UInt z = getMaxCoor(Point(nx,ny,nz));
  UInt x = (z+1) % 3;
  UInt y = (z+2) % 3;

  Point q(pt[x],pt[y]);
  Point a(A[x],A[y]);
  Point b(B[x],B[y]);
  Point c(C[x],C[y]);

  // If the projection does not belong to the triangle:

  // vertex

  Eigen::Matrix<Real,3,1> lambda = triangle.getBaryCoordinates(pt);
  if(lambda[0]>0 && lambda[1]<0 && lambda[2]<0){
    return (std::pair<Point, Real>(A, computeDistance(P, A)));
  }
  if(lambda[0]<0 && lambda[1]>0 && lambda[2]<0){
    return (std::pair<Point, Real>(B, computeDistance(P, B)));
  }
  if(lambda[0]<0 && lambda[1]<0 && lambda[2]>0){
    return (std::pair<Point, Real>(C, computeDistance(P, C)));
  }

  // edge

  if(lambda[0]>0 && lambda[1]>0 && lambda[2]<0){
    VectorXr pt_A(3), B_A(3), proj_AB(3), vecA(3);
    pt_A << pt[0]-A[0], pt[1]-A[1], pt[2]-A[2];
    B_A << B[0]-A[0], B[1]-A[1], B[2]-A[2];
    vecA << A[0], A[1], A[2];
    proj_AB = (pt_A.dot(B_A) / B_A.dot(B_A))*B_A + vecA;
    Point p03d(proj_AB[0], proj_AB[1], proj_AB[2]);

    Point pt2d(pt[x], pt[y]);
    Point p02d(proj_AB[x], proj_AB[y]);

    if(triangle.isPointInside(p03d)){
      return(std::pair<Point,Real>(p03d, computeDistance(P, p03d)));
    }
    else{
      Real distA = computeDistance(p02d, Point(A[x], A[y]));
      Real distB = computeDistance(p02d, Point(B[x], B[y]));
      std::pair<Point, Real> ret;
      distA < distB ? ret = std::make_pair(A, computeDistance(P, A)) : ret = std::make_pair(B, computeDistance(P, B));
      return ret;
    }
  }

  if(lambda[0]>0 && lambda[1]<0 && lambda[2]>0){
    VectorXr pt_A(3), C_A(3), proj_AC(3), vecA(3);
    pt_A << pt[0]-A[0], pt[1]-A[1], pt[2]-A[2];
    C_A << C[0]-A[0], C[1]-A[1], C[2]-A[2];
    vecA << A[0], A[1], A[2];
    proj_AC = (pt_A.dot(C_A) / C_A.dot(C_A))*C_A + vecA;
    Point p03d(proj_AC[0], proj_AC[1], proj_AC[2]);

    Point pt2d(pt[x], pt[y]);
    Point p02d(proj_AC[x], proj_AC[y]);

    if(triangle.isPointInside(p03d)){
      return(std::pair<Point,Real>(p03d, computeDistance(P, p03d)));
    }
    else{
      Real distA = computeDistance(p02d, Point(A[x], A[y]));
      Real distC = computeDistance(p02d, Point(C[x], C[y]));
      std::pair<Point, Real> ret;
      distA < distC ? ret=std::make_pair(A, computeDistance(P, A)) : ret=std::make_pair(C, computeDistance(P, C));
      return ret;
    }
  }

  if(lambda[0]<0 && lambda[1]>0 && lambda[2]>0){
    VectorXr pt_B(3), C_B(3), proj_BC(3), vecB(3);
    pt_B << pt[0]-B[0], pt[1]-B[1], pt[2]-B[2];
    C_B << C[0]-B[0], C[1]-B[1], C[2]-B[2];
    vecB << B[0], B[1], B[2];
    proj_BC = (pt_B.dot(C_B) / C_B.dot(C_B))*C_B + vecB;
    Point p03d(proj_BC[0], proj_BC[1], proj_BC[2]);

    Point pt2d(pt[x], pt[y]);
    Point p02d(proj_BC[x], proj_BC[y]);

    if(triangle.isPointInside(p03d)){
      return(std::pair<Point,Real>(p03d, computeDistance(P, p03d)));
    }
    else{
      Real distB = computeDistance(p02d, Point(B[x], B[y]));
      Real distC = computeDistance(p02d, Point(C[x], C[y]));
      std::pair<Point, Real> ret;
      distB < distC ? ret=std::make_pair(B, computeDistance(P, B)) : ret=std::make_pair(C, computeDistance(P, C));
      return ret;
    }
  }

  // centroid
  Point G = Point(0.333*(A[0] + B[0] + C[0]), 0.333*(A[1] + B[1] + C[1]), 0.333*(A[2] + B[2] + C[2]));
  return(std::pair<Point,Real>(G, computeDistance(P, G)));
}

template<UInt ORDER>
std::vector<Point> projection<ORDER,2,3>::computeProjection()
{
  std::vector<Point> res;
  res.resize(getNumPoints());

  // (1): nearest node to each point
  std::vector<UInt> nearest_node;
  nearest_node.reserve(getNumPoints());

  for(UInt d=0; d<getNumPoints(); d++){
    std::pair<Real, UInt> min_dist (computeDistance(deData_[d], mesh_.getPoint(0)), 0);
    Real dist;
    for(UInt n=1; n<mesh_.num_nodes(); n++){
      dist = computeDistance(deData_[d], mesh_.getPoint(n));
      if(dist < min_dist.first) min_dist = std::make_pair(dist, n);
    }
    nearest_node.push_back(min_dist.second);
  }

  // (2): patch of the nearest node to each point
  std::vector<std::vector<UInt>> patch_element; // for each node -> id of elements that contains the node
  patch_element.reserve(getNumPoints());
  for(UInt i=0; i<getNumPoints(); i++){
    patch_element.push_back(computeNodePatch(nearest_node[i]));
  }

  // (3): I project each data in the patch of the nearest node
  for(UInt i=0; i<getNumPoints(); i++){
    Real dist = std::numeric_limits<Real>::max();
    Point new_datum;
    for(auto elem: patch_element[i]){
      constexpr UInt Nodes = 3*ORDER;
      Element<Nodes, 2, 3> element = mesh_.getElement(elem);
      std::pair<Point, Real> proj = project(element, deData_[i]);

      if(proj.second < dist){
        res[i] = proj.first;

        dist = proj.second;
      }
    }

  }

  return res;
}



#endif

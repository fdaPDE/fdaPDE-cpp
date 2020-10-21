#ifndef __PROJECTION_H__
#define __PROJECTION_H__

#include "../../Mesh/Include/Mesh.h"

template <UInt ORDER,UInt mydim, UInt ndim>
class projection{
};

template <UInt ORDER>
class projection<ORDER,2,2>{
private:
  MeshHandler<ORDER,2,2> mesh_;
  const std::vector<Point> & deData_; // the points to be projected
  UInt num_points;

public:
  projection(MeshHandler<ORDER,2,2> m, const std::vector<Point> & d): mesh_(m), deData_(d), num_points(d.size()) {};

  std::vector<Point> computeProjection() {return deData_;}
};


template<UInt ORDER>
class projection<ORDER,2,3>{
private:
  MeshHandler<ORDER,2,3> mesh_;
  const std::vector<Point> & deData_; // the points to be projected
  UInt num_points;

  UInt getNumPoints() const {return num_points;}
  UInt getMaxCoor(const Point& ) const;
  Real computeDistance(const Point&, const Point&) const;
  Real getAreaTriangle2d(const Point&, const Point&, const Point&) const;
  std::vector<UInt> computeNodePatch(UInt ) const;

  std::pair<Point, Real> project(const Element<3*ORDER,2,3>& , const Point& ) const;

public:
  projection(MeshHandler<ORDER,2,3> m, const std::vector<Point> & d): mesh_(m), deData_(d), num_points(d.size()) {};

  std::vector<Point> computeProjection();

};


template <UInt ORDER>
class projection<ORDER,3,3>{
private:
  MeshHandler<ORDER,3,3> mesh_;
  const std::vector<Point> & deData_; // the points to be projected
  UInt num_points;

public:
  projection(MeshHandler<ORDER,3,3> m, const std::vector<Point> & d): mesh_(m), deData_(d), num_points(d.size()) {};

  std::vector<Point> computeProjection() {return deData_;}
};


#include "Projection_imp.h"

#endif

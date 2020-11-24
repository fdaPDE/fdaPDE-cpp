#ifndef __PROJECTION_H__
#define __PROJECTION_H__

#include "../../Mesh/Include/Mesh.h"

template <UInt ORDER,UInt mydim, UInt ndim>
class projection{
};

template <UInt ORDER>
class projection<ORDER,2,2>{
private:
  const MeshHandler<ORDER,2,2>& mesh_;
  const std::vector<Point<2> >& deData_; // the points to be projected
  UInt num_points;

public:
  projection(const MeshHandler<ORDER,2,2>& m, const std::vector<Point<2> > & d): mesh_(m), deData_(d), num_points(d.size()) {};

  std::vector<Point<2> > computeProjection() {return deData_;}
};


template<UInt ORDER>
class projection<ORDER,2,3>{
private:
  const MeshHandler<ORDER,2,3>& mesh_;
  const std::vector<Point<3> > & deData_; // the points to be projected
  const UInt num_points;
  
  std::vector<UInt> computeNodePatch(UInt ) const;

public:
  projection(const MeshHandler<ORDER,2,3>& m, const std::vector<Point<3> > & d): mesh_(m), deData_(d), num_points(d.size()) {};

  std::vector<Point<3> > computeProjection();

};


template <UInt ORDER>
class projection<ORDER,3,3>{
private:
  const MeshHandler<ORDER,3,3>& mesh_;
  const std::vector<Point<3> > & deData_; // the points to be projected
  UInt num_points;

public:
  projection(const MeshHandler<ORDER,3,3>& m, const std::vector<Point<3> > & d): mesh_(m), deData_(d), num_points(d.size()) {};

  std::vector<Point<3> > computeProjection() {return deData_;}
};


#include "Projection_imp.h"

#endif

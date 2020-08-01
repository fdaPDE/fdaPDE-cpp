#ifndef _PROJECTION_IMP_H
#define _PROJECTION_IMP_H

#include <cmath>
#include <utility>
#include <limits>
#include "mesh_objects.h"


template<UInt ORDER>
std::vector<UInt> projection<ORDER,2,3>::computeNodePatch(UInt id_node) const
{
  constexpr UInt NNODES = how_many_nodes(ORDER,2);
  std::vector<UInt> patch_node;

  for(UInt t=0; t<mesh_.num_elements(); ++t){
    for(UInt i=0; i<NNODES; ++i){
      if(mesh_.elements(t,i) == id_node){
        patch_node.push_back(t);
        break;
      }
    }
  }
  return patch_node;
}


template<UInt ORDER>
std::vector<Point<3> > projection<ORDER,2,3>::computeProjection(){

  constexpr UInt NNODES = how_many_nodes(ORDER,2);

  std::vector<Point<3> > res;
  res.reserve(num_points);

  for(UInt i=0; i<num_points; ++i){
    Real dist2 = std::numeric_limits<Real>::max();
    UInt pos_min;
    for(UInt n=0; n<mesh_.num_nodes(); ++n){
      Point<3> current_point = mesh_.getPoint(n);
      if(deData_[i].dist2(current_point) < dist2){
        pos_min = n;
        dist2 = deData_[i].dist2(current_point);
      }
    }
  
    res.push_back(mesh_.getPoint(pos_min));

    for(auto elem : computeNodePatch(pos_min)){
      Element<NNODES, 2, 3> element = mesh_.getElement(elem);
      Point<3> proj = element.computeProjection(deData_[i]);

      if(proj.dist2(deData_[i]) < dist2){
        res[i] = proj;
        dist2 = proj.dist2(deData_[i]);
      }
    }
  }
  return res;
}



#endif

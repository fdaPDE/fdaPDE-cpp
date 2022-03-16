#ifndef __PROJECTION_IMP_H__
#define __PROJECTION_IMP_H__

#include <cmath>
#include <utility>
#include <limits>
#include "../../Mesh/Include/Mesh_Objects.h"

template<UInt ORDER, UInt mydim ,UInt ndim>
template<bool isManifold>
typename std::enable_if<isManifold,std::vector<UInt> >::type
projection<ORDER,mydim,ndim>::computeNodePatch(UInt id_node) const{
    constexpr UInt NNODES = how_many_nodes(ORDER,mydim);
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

template<UInt ORDER, UInt mydim ,UInt ndim>
template<bool isManifold>
typename std::enable_if<isManifold,std::vector<Point<ndim>> >::type
projection<ORDER,mydim,ndim>::computeProjection(){

    constexpr UInt NNODES = how_many_nodes(ORDER,mydim);

    std::vector<Point<ndim> > res;
    res.reserve(num_points);

    for(UInt i=0; i<num_points; ++i){
        Real dist2 = std::numeric_limits<Real>::max();
        UInt pos_min;

        for(UInt n=0; n<mesh_.num_nodes(); ++n){
            Point<ndim> current_point = mesh_.getPoint(n);
            if(deData_[i].dist2(current_point) < dist2){
                pos_min = n;
                dist2 = deData_[i].dist2(current_point);
            }
        }

        res.push_back(mesh_.getPoint(pos_min));

        for(auto elem : this->computeNodePatch(pos_min)){

            Point<ndim> proj = mesh_.getElement(elem).computeProjection(deData_[i]);

            if(proj.dist2(deData_[i]) < dist2){
                res[i] = proj;
                dist2 = proj.dist2(deData_[i]);
            }
        }
    }
return res;
}

#endif

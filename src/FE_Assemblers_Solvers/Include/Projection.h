#ifndef __PROJECTION_H__
#define __PROJECTION_H__

#include "../../Mesh/Include/Mesh.h"

//! Template class for points projection on a given mesh.
/*!
\tparam ORDER UInt representing order of the element of the mesh
\tparam mydim UInt representing the mesh space size, 1 for 1.5D, 2 for 2D and 2.5D, 3 for 3D
\tparam ndim UInt representing the space size, 2 for 1.5D and 2D, 3 for 2.5D and 3D
\sa Mesh
*/
template<UInt ORDER, UInt mydim, UInt ndim>
class projection{
private:
    const MeshHandler<ORDER,mydim,ndim>& mesh_;
    const std::vector<Point<ndim> > & deData_; // the points to be projected
    const UInt num_points;

    //! SFINAE based method: find all the elements to which \p id_node belongs.
    /*!
    It is available only for not Manifold meshes.
    \param id_node global node number
    */
    template<bool isManifold=(ndim!=mydim)>
    typename std::enable_if<isManifold, std::vector<UInt> >::type
    computeNodePatch(UInt id_node) const;

public:
    //! A constructor. It initializes the constructor given a mesh object and a RMatrix storing the coordinates of points to project.
    projection(const MeshHandler<ORDER,mydim,ndim>& m, const std::vector<Point<ndim> > & d): mesh_(m), deData_(d), num_points(d.size()) {};

    //! SFINAE based method: Computing projections of the stored points
    /*!
    It is only defined for Manifold meshes
    */
    template<bool isManifold=(ndim!=mydim)>
    typename std::enable_if<isManifold, std::vector<Point<ndim>> >::type
    computeProjection();

    //! SFINAE based method: declaration-only for not-manifold meshes
    template<bool isManifold=(ndim!=mydim)>
    typename std::enable_if<!isManifold, std::vector<Point<ndim>> >::type
    computeProjection();
};

#include "Projection_imp.h"

#endif

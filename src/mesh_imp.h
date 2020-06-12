#ifndef MESH_IMP_H_
#define MESH_IMP_H_

#ifdef R_VERSION_
template <UInt ORDER, UInt mydim, UInt ndim>
MeshHandler<ORDER,mydim,ndim>::MeshHandler(SEXP mesh)
{
	mesh_ 		= mesh;
	points_ 	= REAL(VECTOR_ELT(mesh_, 0));
	sides_ 		= INTEGER(VECTOR_ELT(mesh_, 6));
	elements_  = INTEGER(VECTOR_ELT(mesh_, 3));
	neighbors_  = INTEGER(VECTOR_ELT(mesh_, 8));

	num_nodes_ = INTEGER(Rf_getAttrib(VECTOR_ELT(mesh_, 0), R_DimSymbol))[0];
	num_sides_ = INTEGER(Rf_getAttrib(VECTOR_ELT(mesh_, 6), R_DimSymbol))[0];
	num_elements_ = INTEGER(Rf_getAttrib(VECTOR_ELT(mesh_, 3), R_DimSymbol))[0];

}
#endif

template <UInt ORDER, UInt mydim, UInt ndim>
inline Point<ndim> MeshHandler<ORDER,mydim,ndim>::getPoint(UInt id) const
{
	return Point<ndim>(id, points_, num_nodes_);
}

template <UInt ORDER, UInt mydim, UInt ndim>
typename MeshHandler<ORDER,mydim,ndim>::meshElement MeshHandler<ORDER,mydim,ndim>::getElement(UInt id) const
{
	typename meshElement::elementPoints elPoints;
	for (int j=0; j<how_many_nodes(ORDER,mydim); ++j)
		elPoints[j] = getPoint(elements(id,j));
	return meshElement(id,elPoints);
}

template <UInt ORDER, UInt mydim, UInt ndim>
typename MeshHandler<ORDER,mydim,ndim>::meshElement MeshHandler<ORDER,mydim,ndim>::getNeighbors(UInt id_element, UInt number) const
{
	UInt id_neighbor{neighbors(id_element, number)};
	//return empty element if "neighbor" not present (out of boundary!)
	return (id_neighbor==-1) ? meshElement() : getElement(id_neighbor);
}

template <UInt ORDER, UInt mydim, UInt ndim>
typename MeshHandler<ORDER,mydim,ndim>::meshElement MeshHandler<ORDER,mydim,ndim>::findLocationNaive(const Point<ndim>& point) const
{
	for(UInt id=0; id < num_elements_; ++id){
		meshElement current_element{getElement(id)};
		if(current_element.isPointInside(point))
			return current_element;
	}
	return meshElement(); //default element with NVAL ID
}

// Visibility walk algorithm which uses barycentric coordinate [Sundareswara et al]
//Starting triangles usually n^(1/3) points
template <UInt ORDER, UInt mydim, UInt ndim>
template <UInt m, UInt n>																																		//vvvvvvvvv actual return type if enabled
typename std::enable_if<n==m && n==ndim && m==mydim, typename MeshHandler<ORDER,mydim,ndim>::meshElement>::type
MeshHandler<ORDER,mydim,ndim>::findLocationWalking(const Point<ndim>& point, const Element<how_many_nodes(ORDER,mydim),mydim,ndim>& starting_element) const
{
	meshElement current_element{starting_element};
	//Test for found Element, or out of border
	while(current_element.hasValidId() && !current_element.isPointInside(point))
		current_element = getNeighbors(current_element.getId(), current_element.getPointDirection(point));

	return current_element;
}


template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandler<ORDER,mydim,ndim>::printPoints(std::ostream& os)
{
	os<<"# Nodes: "<<num_nodes_<<std::endl;
	for(UInt i=0; i<num_nodes_; ++i)
		os<<getPoint(i);
}


template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandler<ORDER,mydim,ndim>::printElements(std::ostream& os)
{
	os << "# Triangles: "<< num_elements_ <<std::endl;
	for (UInt i = 0; i < num_elements_; ++i )
		os<<getElement(i);
}

template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandler<ORDER,mydim,ndim>::printNeighbors(std::ostream& os)
{
	os << "# Neighbors list: "<< num_elements_ <<std::endl;
	for (UInt i = 0; i < num_elements_; ++i ){
		for( UInt j = 0; j < mydim+1; ++j)
			os<<neighbors(i,j)<<" ";
		os<<std::endl;
	}
}

#endif

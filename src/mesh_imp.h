#ifndef MESH_IMP_H_
#define MESH_IMP_H_


template <UInt ORDER, UInt mydim, UInt ndim>
MeshHandler<ORDER,mydim,ndim>::MeshHandler(Real* points, UInt* sides, UInt* elements, UInt* neighbors, UInt num_nodes, UInt num_sides, UInt num_elements, UInt search) :
		points_(points), sides_(sides), elements_(elements), neighbors_(neighbors),
			num_nodes_(num_nodes), num_sides_(num_sides), num_elements_(num_elements),
				search_(search) {
					tree_ptr_.reset(new ADTree<meshElement>(points_, elements_, num_nodes_, num_elements_));
				}


#ifdef R_VERSION_
template <UInt ORDER, UInt mydim, UInt ndim>
MeshHandler<ORDER,mydim,ndim>::MeshHandler(SEXP Rmesh, UInt search) :
	points_(REAL(VECTOR_ELT(Rmesh, 0))), sides_(INTEGER(VECTOR_ELT(Rmesh, 6))),
		elements_(INTEGER(VECTOR_ELT(Rmesh, 3))), neighbors_(INTEGER(VECTOR_ELT(Rmesh, 8))),
			num_nodes_(INTEGER(Rf_getAttrib(VECTOR_ELT(Rmesh, 0), R_DimSymbol))[0]),
				num_sides_(INTEGER(Rf_getAttrib(VECTOR_ELT(Rmesh, 6), R_DimSymbol))[0]),
					num_elements_(INTEGER(Rf_getAttrib(VECTOR_ELT(Rmesh, 3), R_DimSymbol))[0]),
					 	search_(search) {
							if(XLENGTH(Rmesh)==11 && search==2)
								tree_ptr_.reset(new ADTree<meshElement>(points_, elements_, num_nodes_, num_elements_));
							else if (search==2)
								tree_ptr_.reset(new ADTree<meshElement>(Rmesh));
						}
#endif

template <UInt ORDER, UInt mydim, UInt ndim>
inline Point<ndim> MeshHandler<ORDER,mydim,ndim>::getPoint(const UInt id) const
{
	return Point<ndim>(id, points_, num_nodes_);
}

template <UInt ORDER, UInt mydim, UInt ndim>
typename MeshHandler<ORDER,mydim,ndim>::meshElement MeshHandler<ORDER,mydim,ndim>::getElement(const UInt id) const
{
	typename meshElement::elementPoints elPoints;
	for (int j=0; j<how_many_nodes(ORDER,mydim); ++j)
		elPoints[j] = getPoint(elements(id,j));
	return meshElement(id,elPoints);
}

template <UInt ORDER, UInt mydim, UInt ndim>
typename MeshHandler<ORDER,mydim,ndim>::meshElement MeshHandler<ORDER,mydim,ndim>::getNeighbors(const UInt id_element, const UInt number) const
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
typename MeshHandler<ORDER,mydim,ndim>::meshElement MeshHandler<ORDER,mydim,ndim>::findLocationWalking(const Point<ndim>& point, const meshElement& starting_element) const
{
	static_assert(mydim==ndim, "ERROR! WALKING SEARCH CANNOT BE USED ON MANIFOLD MESHES! See mesh_imp.h");
	meshElement current_element{starting_element};
	//Test for found Element, or out of border
	while(current_element.hasValidId() && !current_element.isPointInside(point))
		current_element = getNeighbors(current_element.getId(), current_element.getPointDirection(point));

	return current_element;
}

template <UInt ORDER, UInt mydim, UInt ndim>
typename MeshHandler<ORDER,mydim,ndim>::meshElement MeshHandler<ORDER,mydim,ndim>::findLocationTree(const Point<ndim>& point) const {
	std::set<int> found;
	std::vector<Real> region;
	region.reserve(2*ndim);

	for (UInt i=0; i<ndim; ++i){
		region.push_back(point[i]);
	}
	for (UInt i=0; i<ndim; ++i){
		region.push_back(point[i]);
	}

	if(!tree_ptr_->search(region, found)) {
		return meshElement();
	}

	for (const auto &i : found) {
		const UInt index = tree_ptr_->pointId(i);
	  meshElement tmp = getElement(index);
		if(tmp.isPointInside(point)) {
			return tmp;
		}
	}
	return meshElement();
}


template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandler<ORDER,mydim,ndim>::printPoints(std::ostream& os) const
{
	os<<"# Nodes: "<<num_nodes_<<std::endl;
	for(UInt i=0; i<num_nodes_; ++i)
		os<<getPoint(i);
}


template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandler<ORDER,mydim,ndim>::printElements(std::ostream& os) const
{
	os << "# Triangles: "<< num_elements_ <<std::endl;
	for (UInt i = 0; i < num_elements_; ++i )
		os<<getElement(i);
}

template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandler<ORDER,mydim,ndim>::printNeighbors(std::ostream& os) const
{
	os << "# Neighbors list: "<< num_elements_ <<std::endl;
	for (UInt i = 0; i < num_elements_; ++i ){
		for( UInt j = 0; j < mydim+1; ++j)
			os<<neighbors(i,j)<<" ";
		os<<std::endl;
	}
}

template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandler<ORDER,mydim,ndim>::printTree(std::ostream & os) const
{
	os << "# Tree characteristic: " <<std::endl;
	if (tree_ptr_)
		os << *tree_ptr_ << std::endl;
	else
		os << "No tree!" <<std::endl;
}

#endif

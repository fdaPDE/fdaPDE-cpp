#ifndef __POINT_IMP_HPP__
#define __POINT_IMP_HPP__

template <>
constexpr Point<2>::Point(UInt id, UInt bcId, const Real(&coord)[2]) :
    Identifier(id, bcId), coord_({coord[0], coord[1]}) {}

template <>
constexpr Point<3>::Point(UInt id, UInt bcId, const Real(&coord)[3]) :
    Identifier(id, bcId), coord_({coord[0], coord[1], coord[2]}) {}

template <>
inline Point<2>::Point(UInt id, UInt bcId, const EigenCoords& coord) :
  Identifier(id, bcId), coord_({coord[0], coord[1]}) {}

template <>
inline Point<3>::Point(UInt id, UInt bcId, const EigenCoords& coord) :
  Identifier(id, bcId), coord_({coord[0], coord[1], coord[2]}) {}

template <>
inline Point<2>::Point(UInt id, const Real* const points, const UInt num_points) :
		Point(id, {points[id], points[id+num_points]}) {}

template <>
inline Point<3>::Point(UInt id, const Real* const points, const UInt num_points) :
		Point(id, {points[id], points[id+num_points], points[id+2*num_points]}) {}


// This function returns the squared euclidean distance between "this" point and other
template <UInt ndim>
inline Real Point<ndim>::dist2(const Point<ndim> &other) const {
	Real dist2{0.};
	for (UInt i=0; i<ndim; ++i)
		dist2+=(coord_[i]-other[i])*(coord_[i]-other[i]);
	return dist2;
}

// This function returns the euclidean distance between "this" point and other
template <UInt ndim>
inline Real Point<ndim>::dist(const Point<ndim> &other) const {
	return std::sqrt(dist2(other));
}

// Overloaded +=/-= operators
template <UInt ndim>
inline Point<ndim>& Point<ndim>::operator+=(const Point &other){
	for (UInt i=0; i<ndim; ++i)
		coord_[i]+=other[i];
	return *this;
}

template <UInt ndim>
inline Point<ndim>& Point<ndim>::operator-=(const Point &other){
	for (UInt i=0; i<ndim; ++i)
		coord_[i]-=other[i];
	return *this;
}

template <UInt NDIM>
std::ostream& operator<<(std::ostream& os, const Point<NDIM>& p){
  if(p.hasValidId())
    os<<p.getId()<<":";
  for (const auto &c : p.coord_)
    os<<" "<<c;
  return os<<std::endl;
}



#endif

#ifndef __POINT_IMP_H__
#define __POINT_IMP_H__

// Note: has to be done like this because of constexpr!
template<>
constexpr Point<1>::Point(UInt id, UInt bcId, const Real(&coord)[1]) :
        Identifier(id, bcId), coord_({coord[0]}) {}

template<>
constexpr Point<2>::Point(UInt id, UInt bcId, const Real(&coord)[2]) : 
  Identifier(id, bcId), coord_({coord[0],coord[1]}) {}

template<>
constexpr Point<3>::Point(UInt id, UInt bcId, const Real(&coord)[3]) : 
  Identifier(id, bcId), coord_({coord[0],coord[1],coord[2]}) {}


template<UInt ndim>
Point<ndim>::Point(UInt id, UInt bcId, const EigenCoords& coord) : 
  Identifier(id, bcId) {
    for(UInt i=0; i<ndim; ++i)
      coord_[i]=coord[i];
  }

template<UInt ndim>
Point<ndim>::Point(UInt id, const RNumericMatrix& points) : 
  Identifier(id) {
    for(UInt i=0; i<ndim; ++i)
      coord_[i]=points(id, i);
  }


// This function returns the squared euclidean distance between "this" point and other
template <UInt ndim>
inline Real Point<ndim>::dist2(const Point<ndim> &other) const {
	return (this->eigenView()-other.eigenView()).squaredNorm();
}

// This function returns the euclidean distance between "this" point and other
template <UInt ndim>
inline Real Point<ndim>::dist(const Point<ndim> &other) const {
	return std::sqrt(dist2(other));
}

// Overloaded +=/-= operators
template <UInt ndim>
inline Point<ndim>& Point<ndim>::operator+=(const Point &other){
	this->eigenView()+=other.eigenView();
	return *this;
}

template <UInt ndim>
inline Point<ndim>& Point<ndim>::operator-=(const Point &other){
	this->eigenView()-=other.eigenView();
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

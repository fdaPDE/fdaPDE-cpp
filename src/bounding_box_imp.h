#ifndef __BOUNDING_BOX_IMP_H__
#define __BOUNDING_BOX_IMP_H__


template<UInt ndim>
Box<ndim>::Box(std::vector<Points<ndim> > const & points) {

	minPoint_=points[0];
	maxPoint_=points[0];

	for(int i=1; i<points.size(); ++i){
		for(int j=0; j<ndim; ++j){

			if (points[i][j]<minPoint_[j])
				minPoint_[j]=points[i][j];

			if (points[i][j]>maxPoint_[j])
				maxPoint_[j]=points[i][j];
		}
	}
}

template<int ndim>
void Box<ndim>::print(std::ostream & out) const {
	out << "Min Point: "<<minPoint_<<std::endl;
	out << "Max Point: "<<maxPoint_<<std::endl;
}

#endif //__BOUNDING_BOX_IMP_H__

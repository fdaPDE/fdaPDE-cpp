#ifndef DOMAIN__IMP_H_
#define DOMAIN__IMP_H_


template<>
Domain<2>::Domain() :
	origin_(), scalingfactors_({1,1}) {}

template<>
Domain<3>::Domain() :
	origin_(), scalingfactors_({1,1,1}) {}


template<UInt ndim>
Domain<ndim>::Domain(const std::vector<Point<ndim> >& points) {

	origin_=points.front();
	scalingfactors_=points.front().coord();

	for(int i=1; i<points.size(); ++i){
		for(int j=0; j<ndim; ++j){
			if (points[i][j]<origin_[j])
				origin_[j]=points[i][j];
			if (points[i][j]>scalingfactors_[j])
				scalingfactors_[j]=points[i][j];
		}
	}

	for (int j=0; j<ndim; ++j){
		// Add the tolerance.
		Real delta = scalingfactors_[j] - origin_[j];
		origin_[j] -= delta*gettolerance();
		scalingfactors_[j] += delta*gettolerance();

		delta = scalingfactors_[j] - origin_[j];
		scalingfactors_[j] = 1./std::max(delta, getmindiff());
	}
}


template<UInt NDIM>
std::ostream & operator<<(std::ostream & ostr, Domain<NDIM> const & d) {
	ostr << "Domain" << std::endl;
	ostr << "------" << std::endl;

	int dimp = NDIM; //not recommended to use T::dp() because of exception case of Point<2>
	for(int i = 0; i < dimp; ++i)
		ostr << "x" << i+1 << "_min = " << d.origin_[i] << std::endl;
	ostr << std::endl;

	for(int i = 0; i < dimp; ++i)
		ostr << "x" << i+1 << "_max = " << d.origin_[i]+1./d.scalingfactors_[i] << std::endl;
	ostr << "------" << std::endl;

	return ostr;
}




#endif /* DOMAIN_IMP_H_ */

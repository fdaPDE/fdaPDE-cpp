#ifndef __DOMAIN_IMP_H__
#define __DOMAIN_IMP_H__

template<class T>
Real Domain<T>::tolerance_ = 1.e-3;

template<class T>
Real Domain<T>::mindiff_ = std::numeric_limits<Real>::min();

template<class T>
Domain<T>::Domain(){ //default domain: whole domain scale from [0,1]
	origin_.resize(T::dt());
	std::vector<Real> tmp(T::dt(),1);
	scalingfactors_ = tmp;
}

template<class T>
Domain<T>::Domain(std::vector<std::vector<Real> > const & coord) {
	int ndimp = T::dp();
	origin_.resize(T::dt());
	scalingfactors_.resize(T::dt());

	/* Find geometric limits.
	 *
	 * If loops are put outside for loops in order to improve performance.
	 */
	if(ndimp == int(coord.size())) {
		// T is equal to Point<3>, Element<mydim, ndim>, Box<NDIMP>
		if(T::dp() == T::dt()) {
			// T is equal to Point<3>
			for(int i = 0; i < ndimp; ++i) {
				origin_[i] = *(std::min_element(coord[i].begin(), coord[i].end()));
				scalingfactors_[i] = *(std::max_element(coord[i].begin(), coord[i].end()));

				// Add the tolerance.
				Real delta = scalingfactors_[i] - origin_[i];
				origin_[i] -= delta*gettolerance();
				scalingfactors_[i] += delta*gettolerance();

				delta = scalingfactors_[i] - origin_[i];
				scalingfactors_[i] = 1./std::max(delta, getmindiff());
			}
		} else {
			// T is equal to Element<mydim, ndim>, Box<NDIMP>
			for(int i = 0; i < ndimp; ++i) {
				origin_[i] = *(std::min_element(coord[i].begin(), coord[i].end()));
				scalingfactors_[i] = *(std::max_element(coord[i].begin(), coord[i].end()));

				// Add the tolerance.
				double delta = scalingfactors_[i] - origin_[i];
				origin_[i] -= delta*gettolerance();
				scalingfactors_[i] += delta*gettolerance();

				delta = scalingfactors_[i] - origin_[i];
				scalingfactors_[i] = 1./std::max(delta, getmindiff());

				/* Repeat the limits because tree dimension is in fact 2 * physical space dimension
				 * because the tree contains triangle or tetrahedron bounding boxes.
				 */
				origin_[i + ndimp] = origin_[i];
				scalingfactors_[i + ndimp] = scalingfactors_[i];
			}
		}
	} else {
		//Special case of Point<2> because Point has default myDim=3
		if(T::dp() == T::dt()) {
			int realndimp = int(coord.size()); // in case of 2-dimen Point (default dimension is 3 and real dimension is 2)
			//need to resize the vectors of origin_, scalingfactors_
			origin_.resize(realndimp);
			scalingfactors_.resize(realndimp);

			for(int i = 0; i < realndimp; ++i) {
					origin_[i] = *(std::min_element(coord[i].begin(), coord[i].end()));
					scalingfactors_[i] = *(std::max_element(coord[i].begin(), coord[i].end()));

					// Add the tolerance.
					Real delta = scalingfactors_[i] - origin_[i];
					origin_[i] -= delta*gettolerance();
					scalingfactors_[i] += delta*gettolerance();

					delta = scalingfactors_[i] - origin_[i];
					scalingfactors_[i] = 1./std::max(delta, getmindiff());
			}
		}
	}
}




template<class T>
std::ostream & operator<<(std::ostream & ostr, Domain<T> const & d) {
	ostr << "Domain" << std::endl;
	ostr << "------" << std::endl;

	int dimp = d.origin_.size(); //not recommended to use T::dp() because of exception case of Point<2>
	for(int i = 0; i < dimp; ++i)
		ostr << "x" << i+1 << "_min = " << d.origin_[i] << std::endl;
	ostr << std::endl;

	for(int i = 0; i < dimp; ++i)
		ostr << "x" << i+1 << "_max = " << d.origin_[i]+1./d.scalingfactors_[i] << std::endl;
	ostr << "------" << std::endl;

	return ostr;
}




#endif /* DOMAIN_IMP_H_ */

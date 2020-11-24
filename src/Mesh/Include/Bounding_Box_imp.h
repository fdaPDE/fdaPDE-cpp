#ifndef __BOUNDING_BOX_IMP_H__
#define __BOUNDING_BOX_IMP_H__


template<int NDIMP>
Box<NDIMP>::Box() {
	x_.resize(2*NDIMP); //multiply 2 for min, max
	for(int i = 0; i < 2*NDIMP; ++i) {
		x_[i] = 0;
	}
}

template<int NDIMP>
Box<NDIMP>::Box(std::vector<Real> const & coord) {
	if(coord.size()==2*NDIMP){
	x_.resize(2*NDIMP); //multiply 2 for min, max
	for(int i = 0; i < 2*NDIMP; ++i)
		x_[i] = coord[i];
	}
	else if(NDIMP==2 && coord.size()==6) {
			x_.resize(2*NDIMP); //4
			x_[0] = std::min(std::min(coord[0], coord[2]), coord[4]); //min x
			x_[1] = std::min(std::min(coord[1], coord[3]), coord[5]); //min y
			x_[2] = std::max(std::max(coord[0], coord[2]), coord[4]); //max x
			x_[3] = std::max(std::max(coord[1], coord[3]), coord[5]); //max y
			}
	else if(NDIMP==3 && coord.size()==9) {
			x_.resize(3*NDIMP); //6
			x_[0] = std::min(std::min(coord[0], coord[3]), coord[6]); //min x
			x_[1] = std::min(std::min(coord[1], coord[4]), coord[7]); //min y
			x_[2] = std::min(std::min(coord[2], coord[5]), coord[8]); //min z
			x_[3] = std::max(std::max(coord[0], coord[3]), coord[6]); //max x
			x_[4] = std::max(std::max(coord[1], coord[4]), coord[7]); //max y
			x_[5] = std::max(std::max(coord[2], coord[5]), coord[8]); //max z
			}
	else if(NDIMP==3 && coord.size()==12) { //Tetrahedron has 4 nodes
			x_.resize(3*NDIMP); //6
			x_[0] = std::min(std::min(std::min(coord[0], coord[3]), coord[6]), coord[9]); //min x
			x_[1] = std::min(std::min(std::min(coord[1], coord[4]), coord[7]), coord[10]); //min y
			x_[2] = std::min(std::min(std::min(coord[2], coord[5]), coord[8]), coord[11]); //min z
			x_[3] = std::max(std::max(std::max(coord[0], coord[3]), coord[6]), coord[9]); //max x
			x_[4] = std::max(std::max(std::max(coord[1], coord[4]), coord[7]), coord[10]); //max y
			x_[5] = std::max(std::max(std::max(coord[2], coord[5]), coord[8]), coord[11]); //max z
	}
}



template<int NDIMP>
template <UInt NNODES,int NDIME,int NDIMPP>
Box<NDIMP>::Box(Element<NNODES,NDIME,NDIMPP> const & element) {
	// if(typeid(NDIMP) != typeid(2)) {
	 //    std::cout << std::endl << std::endl;
	 //    std::cout << "error! Box<NDIMP>::Box(Element<NNODES,NDIMP> const & element) : bad template parameter" << std::endl << std::endl;
	 //    std::cout << "In order to build the bounding box associated to the element" << std::endl;
	 //    std::cout << "template parameter must be equal to '2' and not to " << typeid(NDIMP).name() << std::endl;
	 //    std::exit(EXIT_FAILURE);
	 //  } else {

	if(NDIME == 2 && NDIMPP == 2) {
		x_.resize(2*NDIMP); //4
		x_[0] = std::min(std::min(element[0][0], element[1][0]), element[2][0]); //min x
		x_[1] = std::min(std::min(element[0][1], element[1][1]), element[2][1]); //min y
		x_[2] = std::max(std::max(element[0][0], element[1][0]), element[2][0]); //max x
		x_[3] = std::max(std::max(element[0][1], element[1][1]), element[2][1]); //max y
  	} else if(NDIME == 2 && NDIMPP == 3) {
  		x_.resize(3*NDIMP); //6
		x_[0] = std::min(std::min(element[0][0], element[1][0]), element[2][0]); //min x
		x_[1] = std::min(std::min(element[0][1], element[1][1]), element[2][1]); //min y
		x_[2] = std::min(std::min(element[0][2], element[1][2]), element[2][2]); //min z
		x_[3] = std::max(std::max(element[0][0], element[1][0]), element[2][0]); //max x
		x_[4] = std::max(std::max(element[0][1], element[1][1]), element[2][1]); //max y
		x_[5] = std::max(std::max(element[0][2], element[1][2]), element[2][2]); //max z
  	} else if(NDIME == 3 && NDIMPP == 3) { //Tetrahedron has 4 nodes
  		x_.resize(3*NDIMP); //6
		x_[0] = std::min(std::min(std::min(element[0][0], element[1][0]), element[2][0]), element[3][0]); //min x
		x_[1] = std::min(std::min(std::min(element[0][1], element[1][1]), element[2][1]), element[3][1]); //min y
		x_[2] = std::min(std::min(std::min(element[0][2], element[1][2]), element[2][2]), element[3][2]); //min z
		x_[3] = std::max(std::max(std::max(element[0][0], element[1][0]), element[2][0]), element[3][0]); //max x
		x_[4] = std::max(std::max(std::max(element[0][1], element[1][1]), element[2][1]), element[3][1]); //max y
		x_[5] = std::max(std::max(std::max(element[0][2], element[1][2]), element[2][2]), element[3][2]); //max z
  }
}



template<int NDIMP>
void Box<NDIMP>::set(std::vector<Real> const & data) {
	for(int i = 0; i < 2*NDIMP; ++i) { //multiply 2 for min, max
		x_[i] = data[i];
	}
}


template<int NDIMP>
void Box<NDIMP>::print(std::ostream & out) const {
	out << "Min_Point:  ( - ";
	for (UInt i=0; i<NDIMP; ++i)
		out<< x_[i] <<" - ";
	out << " ) " << std::endl;
	out << "Max_Point:  ( - ";
	for (UInt i=0; i<NDIMP; ++i)
		out<< x_[i + NDIMP] <<" - ";
	out << " ) " << std::endl;
	out<<std::endl;
}

#endif //__BOUNDING_BOX_IMP_H__

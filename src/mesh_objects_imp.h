//#include "mesh_objects.h"
#ifndef __MESH_OBJECTS_IMP_H__
#define __MESH_OBJECTS_IMP_H__


template <UInt NNODES>
void Element<NNODES,2,2>::computeProperties()
{

	Element<NNODES,2,2> &t = *this;
	Point d1(t[1][0]-t[0][0], t[1][1]-t[0][1]);
	Point d2(t[2][0]-t[0][0], t[2][1]-t[0][1]);   //Point d2 = t[2] - t[0]; reimplementare sottrazione


	M_J_(0,0) = d1[0];			// (x2-x1)
	M_J_(1,0) = d1[1];			// (y2-y1)
	M_J_(0,1) = d2[0];			// (x3-x1)
	M_J_(1,1) = d2[1];			// (y3-y1)

	detJ_ = M_J_(0,0) * M_J_(1,1) - M_J_(1,0) * M_J_(0,1);

	Real idet = 1. / detJ_;

	M_invJ_(0,0) =  idet * M_J_(1,1);	// (y3-y1)	-(x3-x1)
	M_invJ_(1,0) = -idet * M_J_(1,0);	// -(y2-y1) (x2-x1)
	M_invJ_(0,1) = -idet * M_J_(0,1);	//
	M_invJ_(1,1) =  idet * M_J_(0,0);	//	è la trasposta di quella della Sangalli (Ael)

	metric_ = M_invJ_ * M_invJ_.transpose();
}

template <UInt NNODES>
Eigen::Matrix<Real,3,1> Element<NNODES,2,2>::getBaryCoordinates(const Point& point) const{

	Element<NNODES,2,2> t=*this;
	Eigen::Matrix<Real,3,1> lambda;
	Eigen::Matrix<Real,4,1> bary_coef;
	//Real eps = 2.2204e-016,
	//	 tolerance = 10000 * eps;

	//cout<<"primovert  "<<t[0](0)<<endl;
	//cout<<"primovert  "<<t[0](1)<<endl;

	bary_coef[0] = t[0][0]-t[2][0];  //x1-x3
	bary_coef[1] = t[1][0]-t[2][0];  //x2-x3
	bary_coef[2] = t[0][1]-t[2][1];  //y1-y3
	bary_coef[3] = t[1][1]-t[2][1];  //y2-y3
	//cout<<baryccoef<<endl;
	//cout<<baryccoef<<endl;

	Real detT = bary_coef[0]*bary_coef[3]-bary_coef[1]*bary_coef[2];
	bary_coef = bary_coef / detT;
	//cout<<"detT  "<<detT<<endl;

	//Compute barycentric coordinates for the point
	Real x_diff_third = point[0] - t[2][0],
		 y_diff_third = point[1] - t[2][1];

	lambda[0] =  (bary_coef[3] * x_diff_third - bary_coef[1] * y_diff_third),
	lambda[1] = (-bary_coef[2] * x_diff_third + bary_coef[0] * y_diff_third),
	lambda[2] = 1 - lambda[0] - lambda[1];

	return lambda;

}


template <UInt NNODES>
bool Element<NNODES,2,2>::isPointInside(const Point& point) const
{
	Real eps = 2.2204e-016,
		 tolerance = 10 * eps;

	Eigen::Matrix<Real,3,1> lambda = getBaryCoordinates(point);

	return ((-tolerance <= lambda[0]) &&
			(-tolerance <= lambda[1]) &&
			(-tolerance <= lambda[2]) );

}


// TO BE FIXED: if one dir -1, try with others
template <UInt NNODES>
int Element<NNODES,2,2>::getPointDirection(const Point& point) const
{
	Real eps = 2.2204e-016,
		 tolerance = 10 * eps;

	Eigen::Matrix<Real,3,1> lambda = getBaryCoordinates(point);

	//Find the minimum coordinate (if negative stronger straight to the point searched)
	int min_index;
	if(lambda[0] <= lambda[1] && lambda[0] <= lambda[2]) 		min_index = 0;
	else if(lambda[1] <= lambda[0] && lambda[1] <= lambda[2]) min_index = 1;
	else 													min_index = 2;

	if(lambda[min_index] < -tolerance) 	return min_index;
	else 							   	return -1;
}


template <UInt NNODES>
void Element<NNODES,2,2>::print(std::ostream & out) const
{
	out<<"Triangle Id -"<< id_ <<"- " <<"< ";
	for (UInt i=0; i<NNODES; ++i)
		out<<points_[i].getId()<<"  ";
	out << ">" << std::endl;
}





//IMPLEMENTAZIONE myDim=2, nDim=3

template <UInt NNODES>
void Element<NNODES,2,3>::computeProperties()
{

	Element<NNODES,2,3> &t = *this;
	Point d1(t[1][0]-t[0][0], t[1][1]-t[0][1], t[1][2]-t[0][2]);
	Point d2(t[2][0]-t[0][0], t[2][1]-t[0][1], t[2][2]-t[0][2]);


	M_J_(0,0) = d1[0];			// (x2-x1)
	M_J_(1,0) = d1[1];			// (y2-y1)
	M_J_(2,0) = d1[2];			// (z2-z1)
	M_J_(0,1) = d2[0];			// (x3-x1)
	M_J_(1,1) = d2[1];			// (y3-y1)
	M_J_(2,1) = d2[2];			// (z3-z1)


	G_J_=M_J_.transpose()*M_J_;

	detJ_ = G_J_(0,0) * G_J_(1,1) - G_J_(1,0) * G_J_(0,1);

	Real idet = 1. / detJ_;

	//M_invJ_(0,0) =  idet * M_J_(1,1);	// (y3-y1)	-(x3-x1)
	//M_invJ_(1,0) = -idet * M_J_(1,0);	// -(y2-y1) (x2-x1)
	//M_invJ_(0,1) = -idet * M_J_(0,1);	//
	//M_invJ_(1,1) =  idet * M_J_(0,0);	//	è la trasposta di quella della Sangalli (Ael)

	metric_(0,0) =  idet * G_J_(1,1);
	metric_(1,0) = -idet * G_J_(1,0);
	metric_(0,1) = -idet * G_J_(0,1);
	metric_(1,1) =  idet * G_J_(0,0);

}

template <UInt NNODES>
Eigen::Matrix<Real,3,1> Element<NNODES,2,3>::getBaryCoordinates(const Point& point) const{


//	Element<NNODES,2,3> t=*this;            implementation trial using area ratios->does not work!
//	Eigen::Matrix<Real,3,1> lambda;
//	Real detJ_point;
//	Eigen::Matrix<Real,3,2> M_J_point;
//	Eigen::Matrix<Real,2,2> G_J_point;
//
//	for (int k=0; k<2; ++k){
////		Point d1(t[1][0]-point[0], t[1][1]-point[1], t[1][2]-point[2]);
////		Point d2(t[2][0]-point[0], t[2][1]-point[1], t[2][2]-point[2]);
//		
//		Point d1(t[k][0]-point[0], t[k][1]-point[1], t[k][2]-point[2]);
//		Point d2(t[k+1][0]-point[0], t[k+1][1]-point[1], t[k+1][2]-point[2]);
//
//		M_J_point(0,0) = d1[0];			// (x2-x1)
//		M_J_point(1,0) = d1[1];			// (y2-y1)
//		M_J_point(2,0) = d1[2];			// (z2-z1)
//		M_J_point(0,1) = d2[0];			// (x3-x1)
//		M_J_point(1,1) = d2[1];			// (y3-y1)
//		M_J_point(2,1) = d2[2];			// (z3-z1)
//
//
//		G_J_point=M_J_point.transpose()*M_J_point;
//
//		detJ_point = G_J_point(0,0) * G_J_point(1,1) - G_J_point(1,0) * G_J_point(0,1);  
//		lambda[2-k]=std::sqrt(detJ_point)/(2*t.getArea());
//
//	}
//	lambda[0]=1-lambda[2]-lambda[1];
//
//	return lambda;

    Element<NNODES,2,3> t=*this;         // implementation using the linear system-> not perfect but the best implementation so far
    Eigen::Matrix<Real,3,1> lambda;
    
	Eigen::Matrix<Real,3,2> A;
	Eigen::Matrix<Real,3,1> b;
	Eigen::Matrix<Real,2,1> sol;
	Eigen::Matrix<Real,3,1> err;

	A(0,0) = t[1][0]-t[0][0];
	A(0,1) = t[2][0]-t[0][0];
	A(1,0) = t[1][1]-t[0][1];
	A(1,1) = t[2][1]-t[0][1];
	A(2,0) = t[1][2]-t[0][2];
	A(2,1) = t[2][2]-t[0][2];

	b(0) = point[0]-t[0][0];
	b(1) = point[1]-t[0][1];
	b(2) = point[2]-t[0][2];

	sol = A.colPivHouseholderQr().solve(b);
	
	err = A*sol-b;
	
//	Real tolerance = (A(0,0)*A(0,0) + A(1,0)*A(1,0) + A(2,0)*A(2,0) + A(0,1)*A(0,1) + A(1,1)*A(1,1) + A(2,1)*A(2,1))/4;  not clear why this tolerance should be used
//	
//	if((err(0)*err(0) + err(1)*err(1) + err(2)*err(2)) >= tolerance ){
//		
//		#ifdef R_VERSION_
//		Rprintf("Warning: finding barycentric coordinates for this point is ill-conditioned");
//		#else
//		std::cout<<"Warning: finding barycentric coordinates for this point is ill-conditioned";
//		#endif
//	}

	
	lambda(0)=1-sol(0)-sol(1);
	lambda(1)=sol(0);
	lambda(2)=sol(1);
	
	return lambda;

}


// THIS COMMENT IS FROM BERAHA, COSMO: We solve 3 scalar equation in 2 unknowns(u,v)
// u*(P1-P0)+v*(P2-P0)=P-P0
// if the system is solveable, P is in the plane (P1,P2,P0), if in addition
// u,v>=0 and u+v<=1 then P is inside the triangle

template <UInt NNODES>
bool Element<NNODES,2,3>::isPointInside(const Point& point) const
{
	Real eps = 2.2204e-016;
	Real tolerance = 10 * eps;
	//THIS COMMENT IS FROM BERAHA, COSMO First: check consistency trough Rouchè-Capelli theorem 
	Element<NNODES,2,3> t=*this;

	Eigen::Matrix<Real,3,2> A;
	Eigen::Matrix<Real,3,1> b;
	Eigen::Matrix<Real,2,1> sol;
	Eigen::Matrix<Real,3,1> err;

	A(0,0) = t[1][0]-t[0][0];
	A(0,1) = t[2][0]-t[0][0];
	A(1,0) = t[1][1]-t[0][1];
	A(1,1) = t[2][1]-t[0][1];
	A(2,0) = t[1][2]-t[0][2];
	A(2,1) = t[2][2]-t[0][2];

	b(0) = point[0]-t[0][0];
	b(1) = point[1]-t[0][1];
	b(2) = point[2]-t[0][2];

	sol = A.colPivHouseholderQr().solve(b);
	err = A*sol-b;
	
	//Real tolerance = (A(0,0)*A(0,0) + A(1,0)*A(1,0) + A(2,0)*A(2,0) + A(0,1)*A(0,1) + A(1,1)*A(1,1) + A(2,1)*A(2,1))/4;

	if( (err(0)*err(0)<tolerance) && (err(1)*err(1)<tolerance) && (err(2)*err(2)<tolerance) ){
		return((sol(0)+sol(1)<=1+2*eps) && (sol(0)>=0-eps) && (sol(1)>=0-eps));
	} else {
		return 0;
	}
}


/*
template <UInt NNODES>
int Triangle<NNODES,2,3>::getPointDirection(const Point& point) const
{
	//da implementare
	std::cerr<<"ancora da implementare";
}
*/

template <UInt NNODES>
void Element<NNODES,2,3>::print(std::ostream & out) const
{
	out<<"Triangle Id -"<< id_ <<"- " <<"< ";
	for (UInt i=0; i<NNODES; ++i)
		out<<points_[i].getId()<<"  ";
	out << ">" << std::endl;
}



//IMPLEMENTAZIONE myDim=3, nDim=3

template <UInt NNODES>
void Element<NNODES,3,3>::computeProperties()
{

	Element<NNODES,3,3> &t = *this;
	Point d1(t[1][0]-t[0][0], t[1][1]-t[0][1], t[1][2]-t[0][2]);
	Point d2(t[2][0]-t[0][0], t[2][1]-t[0][1], t[2][2]-t[0][2]);
	Point d3(t[3][0]-t[0][0], t[3][1]-t[0][1], t[3][2]-t[0][2]);


	M_J_(0,0) = d1[0];			// (x2-x1)
	M_J_(1,0) = d1[1];			// (y2-y1)
	M_J_(2,0) = d1[2];			// (z2-z1)
	M_J_(0,1) = d2[0];			// (x3-x1)
	M_J_(1,1) = d2[1];			// (y3-y1)
	M_J_(2,1) = d2[2];			// (z3-z1)
	M_J_(0,2) = d3[0];			// (x4-x1)
	M_J_(1,2) = d3[1];			// (y4-y1)
	M_J_(2,2) = d3[2];			// (z4-z1)
	
	Real detMJ_ = M_J_(0,0) * (M_J_(1,1) * M_J_(2,2) - M_J_(1,2) * M_J_(2,1)) -
		M_J_(0,1) * (M_J_(1,0) * M_J_(2,2) - M_J_(1,2) * M_J_(2,0)) +
		M_J_(0,2) * (M_J_(1,0) * M_J_(2,1) - M_J_(1,1) * M_J_(2,0));
	
	
	Real idetMJ = 1. / detMJ_;
	
	
	M_invJ_(0,0) =  idetMJ * (M_J_(1, 1) * M_J_(2, 2) - M_J_(1, 2) * M_J_(2, 1));
	M_invJ_(0,1) =  idetMJ * (M_J_(0, 2) * M_J_(2, 1) - M_J_(0, 1) * M_J_(2, 2));
	M_invJ_(0,2) =  idetMJ * (M_J_(0, 0) * M_J_(2, 2) - M_J_(0, 2) * M_J_(2, 0));
	M_invJ_(1,0) =  idetMJ * (M_J_(1, 2) * M_J_(2, 0) - M_J_(1, 0) * M_J_(2, 2));
	M_invJ_(1,1) =  idetMJ * (M_J_(0, 0) * M_J_(2, 2) - M_J_(0, 2) * M_J_(2, 0));
	M_invJ_(1,2) =  idetMJ * (M_J_(1, 0) * M_J_(0, 2) - M_J_(0, 0) * M_J_(1, 2));
	M_invJ_(2,0) =  idetMJ * (M_J_(1, 0) * M_J_(2, 1) - M_J_(2, 0) * M_J_(1, 1));
	M_invJ_(2,1) =  idetMJ * (M_J_(2, 0) * M_J_(0, 1) - M_J_(0, 0) * M_J_(2, 1));
	M_invJ_(2,2) =  idetMJ * (M_J_(0, 0) * M_J_(1, 1) - M_J_(1, 0) * M_J_(0, 1));
	
	
	G_J_=M_J_.transpose()*M_J_;

	detJ_ = G_J_(0,0) * (G_J_(1,1) * G_J_(2,2) - G_J_(1,2) * G_J_(2,1)) -
		G_J_(0,1) * (G_J_(1,0) * G_J_(2,2) - G_J_(1,2) * G_J_(2,0)) +
		G_J_(0,2) * (G_J_(1,0) * G_J_(2,1) - G_J_(1,1) * G_J_(2,0));

	Real idet = 1. / detJ_;

	metric_(0,0) =  idet * (G_J_(1, 1) * G_J_(2, 2) - G_J_(1, 2) * G_J_(2, 1));
	metric_(0,1) =  idet * (G_J_(0, 2) * G_J_(2, 1) - G_J_(0, 1) * G_J_(2, 2));
	metric_(0,2) =  idet * (G_J_(0, 0) * G_J_(2, 2) - G_J_(0, 2) * G_J_(2, 0));
	metric_(1,0) =  idet * (G_J_(1, 2) * G_J_(2, 0) - G_J_(1, 0) * G_J_(2, 2));
	metric_(1,1) =  idet * (G_J_(0, 0) * G_J_(2, 2) - G_J_(0, 2) * G_J_(2, 0));
	metric_(1,2) =  idet * (G_J_(1, 0) * G_J_(0, 2) - G_J_(0, 0) * G_J_(1, 2));
	metric_(2,0) =  idet * (G_J_(1, 0) * G_J_(2, 1) - G_J_(2, 0) * G_J_(1, 1));
	metric_(2,1) =  idet * (G_J_(2, 0) * G_J_(0, 1) - G_J_(0, 0) * G_J_(2, 1));
	metric_(2,2) =  idet * (G_J_(0, 0) * G_J_(1, 1) - G_J_(1, 0) * G_J_(0, 1));
	
	
	
	Eigen::Matrix<Real,4,4> m;
	m(0,0) = 1;
	m(0,1) = 1;
	m(0,2) = 1;
	m(0,3) = 1;
	m(1,0) = t[0][0];			
	m(1,1) = t[1][0];			
	m(1,2) = t[2][0];			
	m(1,3) = t[3][0];			
	m(2,0) = t[0][1];			
	m(2,1) = t[1][1];			
	m(2,2) = t[2][1];			
	m(2,3) = t[3][1];			
	m(3,0) = t[0][2];			
	m(3,1) = t[1][2];			
	m(3,2) = t[2][2];			
	m(3,3) = t[3][2];			
	
	
	Volume_= 1./6*
	std::abs(m(0,3) * m(1,2) * m(2,1) * m(3,0) - m(0,2) * m(1,3) * m(2,1) * m(3,0) -
         	 m(0,3) * m(1,1) * m(2,2) * m(3,0) + m(0,1) * m(1,3) * m(2,2) * m(3,0) +
         	 m(0,2) * m(1,1) * m(2,3) * m(3,0) - m(0,1) * m(1,2) * m(2,3) * m(3,0) -
         	 m(0,3) * m(1,2) * m(2,0) * m(3,1) + m(0,2) * m(1,3) * m(2,0) * m(3,1) +
         	 m(0,3) * m(1,0) * m(2,2) * m(3,1) - m(0,0) * m(1,3) * m(2,2) * m(3,1) -
         	 m(0,2) * m(1,0) * m(2,3) * m(3,1) + m(0,0) * m(1,2) * m(2,3) * m(3,1) +
         	 m(0,3) * m(1,1) * m(2,0) * m(3,2) - m(0,1) * m(1,3) * m(2,0) * m(3,2) -
         	 m(0,3) * m(1,0) * m(2,1) * m(3,2) + m(0,0) * m(1,3) * m(2,1) * m(3,2) +
         	 m(0,1) * m(1,0) * m(2,3) * m(3,2) - m(0,0) * m(1,1) * m(2,3) * m(3,2) -
         	 m(0,2) * m(1,1) * m(2,0) * m(3,3) + m(0,1) * m(1,2) * m(2,0) * m(3,3) +
         	 m(0,2) * m(1,0) * m(2,1) * m(3,3) - m(0,0) * m(1,2) * m(2,1) * m(3,3) -
         	 m(0,1) * m(1,0) * m(2,2) * m(3,3) + m(0,0) * m(1,1) * m(2,2) * m(3,3));

}


	


template <UInt NNODES>
Eigen::Matrix<Real,4,1> Element<NNODES,3,3>::getBaryCoordinates(const Point& point) const{


	Element<NNODES,3,3> t=*this;
	Eigen::Matrix<Real,4,1> lambda;
	Eigen::Matrix<Real,3,3> M_J_point;
	Eigen::Matrix<Real,3,1> rhs;
	Eigen::Matrix<Real,3,1> sol;
	
	Point d1(t[1][0]-t[0][0], t[1][1]-t[0][1], t[1][2]-t[0][2]);
	Point d2(t[2][0]-t[0][0], t[2][1]-t[0][1], t[2][2]-t[0][2]);
	Point d3(t[3][0]-t[0][0], t[3][1]-t[0][1], t[3][2]-t[0][2]);
	


	M_J_point(0,0) = d1[0];			// (x2-x1)
	M_J_point(1,0) = d1[1];			// (y2-y1)
	M_J_point(2,0) = d1[2];			// (z2-z1)
	M_J_point(0,1) = d2[0];			// (x3-x1)
	M_J_point(1,1) = d2[1];			// (y3-y1)
	M_J_point(2,1) = d2[2];			// (z3-z1)
	M_J_point(0,2) = d3[0];			// (x4-x1)
	M_J_point(1,2) = d3[1];			// (y4-y1)
	M_J_point(2,2) = d3[2];			// (z4-z1)
	
	
	rhs(0)= point[0]-t[0][0];
	rhs(1)= point[1]-t[0][1];
	rhs(2)= point[2]-t[0][2];
	
	sol = M_J_point.colPivHouseholderQr().solve(rhs);
	
	lambda[1]=sol(0);
	lambda[2]=sol(1);
	lambda[3]=sol(2);
		
	lambda[0]=1-lambda[1]-lambda[2]-lambda[3];

	return lambda;

}


template <UInt NNODES>
bool Element<NNODES,3,3>::isPointInside(const Point& point) const
{
	Real eps = 2.2204e-016;
	Real tolerance = 10 * eps;

	Element<NNODES,3,3> t=*this;
	Eigen::Matrix<Real,4,1> bary_coeff = t.getBaryCoordinates(point);
	return -tolerance <= bary_coeff[0] && -tolerance <= bary_coeff[1] && -tolerance <= bary_coeff[2] && -tolerance <= bary_coeff[3];
}


template <UInt NNODES>
void Element<NNODES,3,3>::print(std::ostream & out) const
{
	out<<"Tetrahedron Id -"<< id_ <<"- "<<"< ";
	for (UInt i=0; i<NNODES; ++i)
		out<<points_[i].getId()<<"  ";
	out << ">" << std::endl;
}




#endif
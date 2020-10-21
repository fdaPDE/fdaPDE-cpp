#include "../../FdaPDE.h"
#include "../Include/Mesh_Objects.h"
#include "../Include/Bounding_Box.h"
#include "../Include/Domain.h"
#include "../Include/Tree_Node.h"
#include "../Include/Tree_Header.h"
#include "../Include/AD_Tree.h"
#include "../Include/Mesh.h"

extern "C"{
void Test_Point_c(int *n){

	try{
	//2-dimen
	Point a(10,20);     //Point with coordinates
	Point b(1,2,15,25); //Point with Id, bcId and coordinates

	//3-dimen
	Point c;	//default constructor
	Point d(10,20,30);     //Point with coordinates
	Point e(1,2,15,25,35); //Point with Id, bcId and coordinates

	std::cout << std::endl << std::endl;
	std::cout << "Print the 1st point (inizialized by 2-dimen coordinates):  " <<std::endl;
	a.print(std::cout);
	std::cout << "physical dimension:  " << a.dp() <<std::endl;
	std::cout << "search dimension:  " << a.dt() << std::endl;
	std::cout << "coordsize():  " << a.coordsize() << std::endl;
	std::cout << std::endl << std::endl;

	std::cout << "Print the 2nd point (inizialized by 2-dimen coordinates and Ids):  " <<std::endl;
	b.print(std::cout);
	std::cout << "physical dimension:  " << b.dp() <<std::endl;
	std::cout << "search dimension:  " << b.dt() << std::endl;
	std::cout << "coordsize():  " << b.coordsize() << std::endl;
	std::cout << std::endl << std::endl;

	std::cout << "Print the 3rd point (default constructor):  " <<std::endl;
	c.print(std::cout);
	std::cout << "physical dimension:  " << c.dp() << std::endl;
	std::cout << "search dimension:  " << c.dt() << std::endl;
	std::cout << "coordsize():  " << c.coordsize() << std::endl;
	std::cout << std::endl << std::endl;

	std::cout << "Print the 4th point (inizialized by 3-dimen coordinates):  " <<std::endl;
	d.print(std::cout);
	std::cout << "physical dimension:  " << d.dp() << std::endl;
	std::cout << "search dimension:  " << d.dt() << std::endl;
	std::cout << "coordsize():  " << d.coordsize() << std::endl;
	std::cout << std::endl << std::endl;

	std::cout << "Print the 5th point (inizialized by 3-dimen coordinates and Ids):  " <<std::endl;
	e.print(std::cout);
	std::cout << "physical dimension:  " << e.dp() << std::endl;
	std::cout << "search dimension:  " << e.dt() << std::endl;
	std::cout << "coordsize():  " << e.coordsize() << std::endl;
	std::cout << std::endl << std::endl;

	*n = 0;

	}
	catch(...)
	{
	std::cout <<"Exception caugth " << std::endl;
	*n = 1;
	}
}
}

extern "C"{
void Test_Element_c(int *n){
	try{
		////Point in 2-dimen
	Point a(1,2,1,1); //Point with Id, bcId and coordinates
	Point b(2,3,3,1); //Point with Id, bcId and coordinates
	Point c(3,4,2,2); //Point with Id, bcId and coordinates

	Point d(4,5,2.5,1.5); //Point with Id, bcId and coordinates
	Point e(5,6,1.5,1.5); //Point with Id, bcId and coordinates
	Point f(6,7,2,1); //Point with Id, bcId and coordinates

	Point test(2,1.5); //Point with coordinates
	std::vector<Point> vec{a,b,c};
	std::vector<Point> vec2{a,b,c,d,e,f};

	//case 1: Triangle in 2-dimen
	Element<3,2,2> Ta;		//default constructor
	Element<6,2,2> Tb;		//default constructor
	Element<3,2,2> Tc(1, vec);	//constructor from id and a vector of point
	Element<6,2,2> Td(2, vec2);	//constructor from id and a vector of point

	//Point in 3-dimen
	Point g(1,2,1,1,2); //Point with Id, bcId and coordinates
	Point h(2,3,3,1,3); //Point with Id, bcId and coordinates
	Point i(3,4,2,2,1); //Point with Id, bcId and coordinates
	Point j(4,5,4,3,2); //Point with Id, bcId and coordinates
	Point test2 (2.5,1.75,2);  //(2,1.5,2.5); //Point with coordinates

	std::vector<Point> vec3{g,h,i};
	std::vector<Point> vec4{g,h,i,j};

	//case 2: Triangle in 3-dimen
	Element<3,2,3> Te;		//default constructor
	Element<3,2,3> Tf(3, vec3);	//constructor from id and a vector of point

	//case 3: Tetrahedron in 3-dimen (there are 4 nodes for tetrahedron)
	Element<4,3,3> Tg(4, vec4);	//constructor from id and a vector of point


	std::cout << "Print the 1st triangle in 2 dimension (default constructor, 3 nodes):  " <<std::endl;
	Ta.print(std::cout);
	std::cout << "Print each node: " << std::endl;
	Ta[0].print(std::cout);
	Ta[1].print(std::cout);
	Ta[2].print(std::cout);
	std::cout << "physical dimension:  " << Ta.dp() <<std::endl;
	std::cout << "search dimension:  " << Ta.dt() <<std::endl;
	std::cout << "test: is point (2, 1.5) inside? (1 if true): " << Ta.isPointInside(test)<<std::endl;
	std::cout << "coordsize():  " << Ta.coordsize() << std::endl;
	std::cout << std::endl << std::endl;

	std::cout << "Print the 2nd triangle 2 dimension (default constructor, 6 nodes):  " <<std::endl;
	Tb.print(std::cout);
	std::cout << "Print each node: " << std::endl;
	Tb[0].print(std::cout);
	Tb[1].print(std::cout);
	Tb[2].print(std::cout);
	Tb[3].print(std::cout);
	Tb[4].print(std::cout);
	Tb[5].print(std::cout);
	std::cout << "physical dimension:  " << Tb.dp() <<std::endl;
	std::cout << "search dimension:  " << Tb.dt() <<std::endl;
	std::cout << "test: is point (2, 1.5) inside? (1 if true): " << Tb.isPointInside(test)<<std::endl;
	std::cout << "coordsize():  " << Tb.coordsize() << std::endl;
	std::cout << std::endl << std::endl;

	std::cout << "Print the 3rd triangle 2 dimension (construct by coordinates Points a(1,2,1,1), b(2,3,3,1),c(3,4,2,2), 3 nodes):  " <<std::endl;
	Tc.print(std::cout);
	std::cout << "Print each node: " << std::endl;
	Tc[0].print(std::cout);
	Tc[1].print(std::cout);
	Tc[2].print(std::cout);
	std::cout << "physical dimension:  " << Tc.dp() <<std::endl;
	std::cout << "search dimension:  " << Tc.dt() <<std::endl;
	std::cout << "test: is point (2, 1.5) inside? (1 if true): " << Tc.isPointInside(test)<<std::endl;
	std::cout << "coordsize():  " << Tc.coordsize() << std::endl;
	std::cout << std::endl << std::endl;

	std::cout << "Print the 4th triangle 2 dimension (construct by coordinates Points a(1,2,1,1), b(2,3,3,1), c(3,4,2,2), d(4,5,2.5,1.5), e(5,6,1.5,1.5), f(6,7,2,1), 6 nodes):  " <<std::endl;
	Td.print(std::cout);
	std::cout << "Print each node: " << std::endl;
	Td[0].print(std::cout);
	Td[1].print(std::cout);
	Td[2].print(std::cout);
	Td[3].print(std::cout);
	Td[4].print(std::cout);
	Td[5].print(std::cout);
	std::cout << "physical dimension:  " << Td.dp() <<std::endl;
	std::cout << "search dimension:  " << Td.dt() <<std::endl;
	std::cout << "test: is point (2, 1.5) inside? (1 if true): " << Td.isPointInside(test)<<std::endl;
	std::cout << "coordsize():  " << Td.coordsize() << std::endl;
	std::cout << std::endl << std::endl;

	std::cout << "Print the 5th triangle in 3 dimension (default constructor):  " <<std::endl;
	Te.print(std::cout);
	std::cout << "Print each node: " << std::endl;
	Te[0].print(std::cout);
	Te[1].print(std::cout);
	Te[2].print(std::cout);
	std::cout << "physical dimension:  " << Te.dp() <<std::endl;
	std::cout << "search dimension:  " << Te.dt() <<std::endl;
	std::cout << "test: is point (2,1.5,2.5) inside? (1 if true): " << Te.isPointInside(test2)<<std::endl;
	std::cout << "coordsize():  " << Te.coordsize() << std::endl;
	std::cout << std::endl << std::endl;

	std::cout << "Print the 6th triangle in 3 dimension (construct by coordinates Points):  " <<std::endl;
	Tf.print(std::cout);
	std::cout << "Print each node: " << std::endl;
	Tf[0].print(std::cout);
	Tf[1].print(std::cout);
	Tf[2].print(std::cout);
	std::cout << "physical dimension:  " << Tf.dp() <<std::endl;
	std::cout << "search dimension:  " << Tf.dt() <<std::endl;
	std::cout << "test: is point (2,1.5,2.5) inside? (1 if true): " << Tf.isPointInside(test2)<<std::endl;
	std::cout << "coordsize():  " << Tf.coordsize() << std::endl;
	std::cout << std::endl << std::endl;

	std::cout << "Print the 7th tetrahedron in 3 dimension (construct by coordinates Points):  " <<std::endl;
	Tg.print(std::cout);
	std::cout << "Print each node: " << std::endl;
	Tg[0].print(std::cout);
	Tg[1].print(std::cout);
	Tg[2].print(std::cout);
	Tg[3].print(std::cout);
	std::cout << "physical dimension:  " << Tg.dp() <<std::endl;
	std::cout << "search dimension:  " << Tg.dt() <<std::endl;
	std::cout << "test: is point (2,1.5,2.5) inside? (1 if true): " << Tg.isPointInside(test2)<<std::endl;
	std::cout << "coordsize():  " << Tg.coordsize() << std::endl;
	std::cout << std::endl << std::endl;
	*n = 0;
	}
	catch(...)
	{
	std::cout <<"Exception caugth " << std::endl;
	*n = 1;
	}
}
}


extern "C"{
void Test_Box_c(int *n){

	try{
	//Point in 2-dimen
	Point a(1,2,1,1); //Point with Id, bcId and coordinates
	Point b(2,3,3,1); //Point with Id, bcId and coordinates
	Point c(3,4,2,2); //Point with Id, bcId and coordinates

	Point d(4,5,2.5,1.5); //Point with Id, bcId and coordinates
	Point e(5,6,1.5,1.5); //Point with Id, bcId and coordinates
	Point f(6,7,2,1); //Point with Id, bcId and coordinates

	//Triangle in 2-dimen
	std::vector<Point> vec{a,b,c};
	std::vector<Point> vec2{a,b,c,d,e,f};
	Element<3,2,2> Tb(1, vec);
	Element<6,2,2> Td(2, vec2);


	//Point in 3-dimen
	Point g(1,2,1,1,2); //Point with Id, bcId and coordinates
	Point h(2,3,3,1,3); //Point with Id, bcId and coordinates
	Point i(3,4,2,2,1); //Point with Id, bcId and coordinates
	Point j(4,5,4,3,2); //Point with Id, bcId and coordinates

	//******* node validation to be tested
	// Point d(4,5,2.5,1.5); //Point with Id, bcId and coordinates
	// Point e(5,6,1.5,1.5); //Point with Id, bcId and coordinates
	// Point f(6,7,2,1); //Point with Id, bcId and coordinates

	//Triangle in 3-dimen
	std::vector<Point> vec3{g,h,i};
	Element<3,2,3> Ta(3, vec3);

	//Tetrahedron in 3-dimen
	std::vector<Point> vec4{g,h,i,j};
	Element<3,3,3> Tc(4, vec4);

	std::vector<Real> coord{1,1,4,3};
	std::vector<Real> coord2{1,1,1,4,3,5};
	std::vector<Real> coord3{2,2,5,4};
	std::vector<Real> coord4{2,2,2,5,4,6};


	//case 1: Triangle with 2-dimen Box
	Box<2> Ba;	//default constructor, 2D box
	Box<2> Bb(coord);	//construct by coordinates
	Box<2> Bc(Tb);	//construct from a Element<3,2,2>
	Box<2> Bd(Td);	//construct from a Element<6,2,2>
	Box<2> Be(Bb);	//copy constructor

	//case 2: Triangle with 3-dimen Box
	Box<3> Bf;	//default constructor, 3D box
	Box<3> Bg(coord2);	//construct by coordinates, 3D box
	Box<3> Bh(Ta);	//construct from a Element<3,2,3>

	//case 3: Tetrahedron with 3-dimen Box
	Box<3> Bi(Tc);	//construct from a Element<3,3,3>


	std::cout << "Print the 1st box (default constructor, NDIMP = 2):  " <<std::endl;
	Ba.print(std::cout);
	std::cout << "operator []: " << std::endl;
	std::cout << "(" << Ba[0] <<" , " << Ba[1] <<")" <<std::endl;
	std::cout << "(" << Ba[2] <<" , " << Ba[3] <<")" <<std::endl;
	std::cout << "physical dimension:  " << Ba.dp() <<std::endl;
	std::cout << "search dimension:  " << Ba.dt() <<std::endl;
	std::cout << "vector dimension:  " << Ba.coordsize() <<std::endl;
	std::cout << "change coordinates to {2,2,5,4}:  " << std::endl;
	Ba.set(coord3);
	Ba.print(std::cout);
	std::cout << std::endl << std::endl;

	std::cout << "Print the 2nd box (construct from a vector {1,1,4,3}, NDIMP = 2):  " <<std::endl;
	Bb.print(std::cout);
	std::cout << "operator []: " << std::endl;
	std::cout << "(" << Bb[0] <<" , " << Bb[1] <<")" <<std::endl;
	std::cout << "(" << Bb[2] <<" , " << Bb[3] <<")" <<std::endl;
	std::cout << "physical dimension:  " << Bb.dp() <<std::endl;
	std::cout << "search dimension:  " << Bb.dt() <<std::endl;
	std::cout << "vector dimension:  " << Bb.coordsize() <<std::endl;
	std::cout << "change coordinates to {2,2,5,4}:  " << std::endl;
	Bb.set(coord3);
	Bb.print(std::cout);
	std::cout << std::endl << std::endl;

	std::cout << "Print the 3rd box (construct from a 3 nodes triangle formed by (1,1),(3,1),(2,2), NDIMP = 2):  " <<std::endl;
	Bc.print(std::cout);
	std::cout << "operator []: " << std::endl;
	std::cout << "(" << Bc[0] <<" , " << Bc[1] <<")" <<std::endl;
	std::cout << "(" << Bc[2] <<" , " << Bc[3] <<")" <<std::endl;
	std::cout << "physical dimension:  " << Bc.dp() <<std::endl;
	std::cout << "search dimension:  " << Bc.dt() <<std::endl;
	std::cout << "vector dimension:  " << Bc.coordsize() <<std::endl;
	std::cout << "change coordinates to {2,2,5,4}:  " << std::endl;
	Bc.set(coord3);
	Bc.print(std::cout);
	std::cout << std::endl << std::endl;

	std::cout << "Print the 4th box (construct from 6 nodes triangle with vertex (1,1),(3,1),(2,2), NDIMP = 2):  " <<std::endl;
	Bd.print(std::cout);
	std::cout << "operator []: " << std::endl;
	std::cout << "(" << Bd[0] <<" , " << Bd[1] <<")" <<std::endl;
	std::cout << "(" << Bd[2] <<" , " << Bd[3] <<")" <<std::endl;
	std::cout << "physical dimension:  " << Bd.dp() <<std::endl;
	std::cout << "search dimension:  " << Bd.dt() <<std::endl;
	std::cout << "vector dimension:  " << Bd.coordsize() <<std::endl;
	std::cout << "change coordinates to {2,2,5,4}  " << std::endl;
	Bd.set(coord3);
	Bd.print(std::cout);
	std::cout << std::endl << std::endl;

	std::cout << "Print the 5th box (copy constructor from the 1st box {1,1,4,3}, NDIMP = 2):  " <<std::endl;
	Be.print(std::cout);
	std::cout << "operator []: " << std::endl;
	std::cout << "(" << Be[0] <<" , " << Be[1] <<")" <<std::endl;
	std::cout << "(" << Be[2] <<" , " << Be[3] <<")" <<std::endl;
	std::cout << "physical dimension:  " << Be.dp() <<std::endl;
	std::cout << "search dimension:  " << Be.dt() <<std::endl;
	std::cout << std::endl << std::endl;

	std::cout << "Print the 6th box (default constructor, NDIMP = 3):  " <<std::endl;
	Bf.print(std::cout);
	std::cout << "operator []: " << std::endl;
	std::cout << "(" << Bf[0] <<" , " << Bf[1] << " , " << Bf[2] <<")" <<std::endl;
	std::cout << "(" << Bf[3] <<" , " << Bf[4] << " , " << Bf[5] <<")" <<std::endl;
	std::cout << "physical dimension:  " << Bf.dp() <<std::endl;
	std::cout << "search dimension:  " << Bf.dt() <<std::endl;
	std::cout << "vector dimension:  " << Bf.coordsize() <<std::endl;
	std::cout << "change coordinates {2,2,2,5,4,6} " << std::endl;
	Bf.set(coord4);
	Bf.print(std::cout);
	std::cout << std::endl << std::endl;

	std::cout << "Print the 7th box (construct from a vector {1,1,1,4,3,5}, NDIMP = 3):  " <<std::endl;
	Bg.print(std::cout);
	std::cout << "operator []: " << std::endl;
	std::cout << "(" << Bg[0] <<" , " << Bg[1] << " , " << Bg[2] <<")" <<std::endl;
	std::cout << "(" << Bg[3] <<" , " << Bg[4] << " , " << Bg[5] <<")" <<std::endl;
	std::cout << "physical dimension:  " << Bg.dp() <<std::endl;
	std::cout << "search dimension:  " << Bg.dt() <<std::endl;
	std::cout << "vector dimension:  " << Bg.coordsize() <<std::endl;
	std::cout << "change coordinates to {2,2,2,5,4,6}  " << std::endl;
	Bg.set(coord4);
	Bg.print(std::cout);
	std::cout << std::endl << std::endl;

	std::cout << "Print the 8th box (construct from a 3 nodes triangle, NDIMP = 3):  " <<std::endl;
	Bh.print(std::cout);
	std::cout << "operator []: " << std::endl;
	std::cout << "(" << Bh[0] <<" , " << Bh[1] << " , " << Bh[2] <<")" <<std::endl;
	std::cout << "(" << Bh[3] <<" , " << Bh[4] << " , " << Bh[5] <<")" <<std::endl;
	std::cout << "physical dimension:  " << Bh.dp() <<std::endl;
	std::cout << "search dimension:  " << Bh.dt() <<std::endl;
	std::cout << "vector dimension:  " << Bh.coordsize() <<std::endl;
	std::cout << "change coordinates to {2,2,2,5,4,6}  " << std::endl;
	Bh.set(coord4);
	Bh.print(std::cout);
	std::cout << std::endl << std::endl;

	std::cout << "Print the 9th box (construct from a 4 nodes Tetrahedron, NDIMP = 3):  " <<std::endl;
	Bi.print(std::cout);
	std::cout << "operator []: " << std::endl;
	std::cout << "(" << Bi[0] <<" , " << Bi[1] << " , " << Bi[2] <<")" <<std::endl;
	std::cout << "(" << Bi[3] <<" , " << Bi[4] << " , " << Bi[5] <<")" <<std::endl;
	std::cout << "physical dimension:  " << Bi.dp() <<std::endl;
	std::cout << "search dimension:  " << Bi.dt() <<std::endl;
	std::cout << "vector dimension:  " << Bi.coordsize() <<std::endl;
	std::cout << "change coordinates to {2,2,2,5,4,6}  " << std::endl;
	Bi.set(coord4);
	Bi.print(std::cout);
	std::cout << std::endl << std::endl;

	*n = 0;
	}
	catch(...)
	{
	std::cout <<"Exception caugth " << std::endl;
	*n = 1;
	}
}
}

extern "C"{
void Test_Domain_c(int *n){

	try{
	std::vector<std::vector<Real>> coord;
	coord.resize(2);
	for (int i = 0; i < 2; i++) //2dimen
	{
		coord[i].resize(10);
		for(int j = 0; j < 10; j++)
			coord[i][j] = i+j;
	}
	// 0 1 2 3 4 5 6 7 8 9
	// 1 2 3 4 5 6 7 8 9 10

	std::vector<std::vector<Real>> coord2;
	coord2.resize(3);
	for (int i = 0; i < 3; i++) //3dimen
	{
		coord2[i].resize(10);
		for(int j = 0; j < 10; j++)
			coord2[i][j] = i+j;
	}
	// 0 1 2 3 4 5 6 7 8  9
	// 1 2 3 4 5 6 7 8 9  10
	// 2 3 4 5 6 7 8 9 10 11



	//default constructor
	Domain<Point> Pa;
	Domain<Element<3,2,2>> Ta;
	Domain<Box<2>> Ba;

	Domain<Element<3,2,3>> Tb;
	Domain<Box<3>> Bb;

	Domain<Element<3,3,3>> Tc;
	Domain<Box<3>> Bc;


	// 2-dimen coord with all type of Shape object
	Domain<Point> Pd(coord);
	Domain<Element<3,2,2>> Td(coord);
	Domain<Box<2>> Bd(coord);

	// 3-dimen coord with all type of Shape object
	Domain<Point> Pe(coord2);
	Domain<Element<3,2,3>> Te(coord2);
	Domain<Box<3>> Be(coord2);

	Domain<Element<3,3,3>> Tf(coord2);
	Domain<Box<3>> Bf(coord2);


	std::cout << "Print the 1st domain (constructor default, shape = Point, NDIMP = 3):  " <<std::endl;
	std::cout << Pa << std::endl;
	std::cout << "origin of the object: " << std::endl;
	for (int i=0; i<Pa.getoriginsize(); ++i)
		std::cout<< Pa.orig(i) <<" - ";
	std::cout<<std::endl;
	std::cout << "scaling factor:  " << std::endl;
	for (int i=0; i<Pa.getscalingsize(); ++i)
		std::cout<< Pa.scal(i) <<" - ";
	std::cout<<std::endl;
	std::cout << "print tolerance:  " << Pa.gettolerance() << std::endl;
	std::cout << "change tolerance to 0.01:  " << std::endl;
	Pa.settolerance(0.01);
	std::cout << "print new tolerance:  " << Pa.gettolerance() << std::endl;
	std::cout << "print minimum difference between coordinates:  " << Pa.getmindiff() << std::endl;
	std::cout << "change minimum difference between coordinates to 0.003:  " << std::endl;
	Pa.setmindiff(0.003);
	std::cout << "print new minimum difference between coordinates:  " << Pa.getmindiff() << std::endl;
	std::cout << std::endl << std::endl;



	std::cout << "Print the 2nd domain (constructor default, shape = Triangle, NDIMP = 2):  " <<std::endl;
	std::cout << Ta << std::endl;
	std::cout << "origin of the object: " << std::endl;
	for (int i=0; i<Ta.getoriginsize(); ++i)
		std::cout<< Ta.orig(i) <<" - ";
	std::cout<<std::endl;
	std::cout << "scaling factor:  " << std::endl;
	for (int i=0; i<Ta.getscalingsize(); ++i)
		std::cout<< Ta.scal(i) <<" - ";
	std::cout<<std::endl;
	std::cout << "print tolerance:  " << Ta.gettolerance() << std::endl;
	std::cout << "change tolerance to 0.01:  " << std::endl;
	Ta.settolerance(0.01);
	std::cout << "print new tolerance:  " << Ta.gettolerance() << std::endl;
	std::cout << "print minimum difference between coordinates:  " << Ta.getmindiff() << std::endl;
	std::cout << "change minimum difference between coordinates to 0.003:  " << std::endl;
	Ta.setmindiff(0.003);
	std::cout << "print new minimum difference between coordinates:  " << Ta.getmindiff() << std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 3rd domain (constructor default, shape = Box with Triangle, NDIMP = 2):  " <<std::endl;
	std::cout << Ba << std::endl;
	std::cout << "origin of the object: " << std::endl;
	for (int i=0; i<Ba.getoriginsize(); ++i)
		std::cout<< Ba.orig(i) <<" - ";
	std::cout<<std::endl;
	std::cout << "scaling factor:  " << std::endl;
	for (int i=0; i<Ba.getscalingsize(); ++i)
		std::cout<< Ba.scal(i) <<" - ";
	std::cout<<std::endl;
	std::cout << "print tolerance:  " << Ba.gettolerance() << std::endl;
	std::cout << "change tolerance to 0.01:  " << std::endl;
	Ba.settolerance(0.01);
	std::cout << "print new tolerance:  " << Ba.gettolerance() << std::endl;
	std::cout << "print minimum difference between coordinates:  " << Ba.getmindiff() << std::endl;
	std::cout << "change minimum difference between coordinates to 0.003:  " << std::endl;
	Ba.setmindiff(0.003);
	std::cout << "print new minimum difference between coordinates:  " << Ba.getmindiff() << std::endl;
	std::cout << std::endl << std::endl;



	std::cout << "Print the 4th domain (constructor default, shape = Triangle, NDIMP = 3):  " <<std::endl;
	std::cout << Tb << std::endl;
	std::cout << "origin of the object: " << std::endl;
	for (int i=0; i<Tb.getoriginsize(); ++i)
		std::cout<< Tb.orig(i) <<" - ";
	std::cout<<std::endl;
	std::cout << "scaling factor:  " << std::endl;
	for (int i=0; i<Tb.getscalingsize(); ++i)
		std::cout<< Tb.scal(i) <<" - ";
	std::cout<<std::endl;
	// std::cout << "print tolerance:  " << Tb.gettolerance() << std::endl;
	// std::cout << "change tolerance to 0.01:  " << std::endl;
	// Tb.settolerance(0.01);
	// std::cout << "print new tolerance:  " << Tb.gettolerance() << std::endl;
	// std::cout << "print minimum difference between coordinates:  " << Tb.getmindiff() << std::endl;
	// std::cout << "change minimum difference between coordinates to 0.003:  " << std::endl;
	// Tb.setmindiff(0.003);
	// std::cout << "print new minimum difference between coordinates:  " << Tb.getmindiff() << std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 5th domain (constructor default, shape = Box with Triangle , NDIMP = 3):  " <<std::endl;
	std::cout << Bb <<std::endl;
	std::cout << "origin of the object: " << std::endl;
	for (int i=0; i<Bb.getoriginsize(); ++i)
		std::cout<< Bb.orig(i) <<" - ";
	std::cout<<std::endl;
	std::cout << "scaling factor:  " << std::endl;
	for (int i=0; i<Bb.getscalingsize(); ++i)
		std::cout<< Bb.scal(i) <<" - ";
	std::cout<<std::endl;
	// std::cout << "print tolerance:  " << Bb.gettolerance() << std::endl;
	// std::cout << "change tolerance to 0.02:  " << std::endl;
	// Bb.settolerance(0.02);
	// std::cout << "print new tolerance:  " << Bb.gettolerance() << std::endl;
	// std::cout << "print minimum difference between coordinates:  " << Bb.getmindiff() << std::endl;
	// std::cout << "change minimum difference between coordinates to 0.005:  " << std::endl;
	// Bb.setmindiff(0.005);
	// std::cout << "print new minimum difference between coordinates:  " << Bb.getmindiff() << std::endl;
	std::cout << std::endl << std::endl;



	std::cout << "Print the 6th domain (constructor default, shape = Tetrahedron, NDIMP = 3):  " <<std::endl;
	std::cout << Tc << std::endl;
	std::cout << "origin of the object: " << std::endl;
	for (int i=0; i<Tc.getoriginsize(); ++i)
		std::cout<< Tc.orig(i) <<" - ";
	std::cout<<std::endl;
	std::cout << "scaling factor:  " << std::endl;
	for (int i=0; i<Tc.getscalingsize(); ++i)
		std::cout<< Tc.scal(i) <<" - ";
	std::cout<<std::endl;
	// std::cout << "print tolerance:  " << Tc.gettolerance() << std::endl;
	// std::cout << "change tolerance to 0.01:  " << std::endl;
	// Tc.settolerance(0.01);
	// std::cout << "print new tolerance:  " << Tc.gettolerance() << std::endl;
	// std::cout << "print minimum difference between coordinates:  " << Tc.getmindiff() << std::endl;
	// std::cout << "change minimum difference between coordinates to 0.003:  " << std::endl;
	// Tc.setmindiff(0.003);
	// std::cout << "print new minimum difference between coordinates:  " << Tc.getmindiff() << std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 7th domain (constructor default, shape = Box with Tetrahedron , NDIMP = 3):  " <<std::endl;
	std::cout << Bc <<std::endl;
	std::cout << "origin of the object: " << std::endl;
	for (int i=0; i<Bc.getoriginsize(); ++i)
		std::cout<< Bc.orig(i) <<" - ";
	std::cout<<std::endl;
	std::cout << "scaling factor:  " << std::endl;
	for (int i=0; i<Bc.getscalingsize(); ++i)
		std::cout<< Bc.scal(i) <<" - ";
	std::cout<<std::endl;
	// std::cout << "print tolerance:  " << Bc.gettolerance() << std::endl;
	// std::cout << "change tolerance to 0.02:  " << std::endl;
	// Bc.settolerance(0.02);
	// std::cout << "print new tolerance:  " << Bc.gettolerance() << std::endl;
	// std::cout << "print minimum difference between coordinates:  " << Bc.getmindiff() << std::endl;
	// std::cout << "change minimum difference between coordinates to 0.005:  " << std::endl;
	// Bc.setmindiff(0.005);
	// std::cout << "print new minimum difference between coordinates:  " << Bc.getmindiff() << std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 8th domain (constructor from coordinates (x1 from 0 to 9, x2 from 1 to 10), shape = Point, NDIMP = 2):  " <<std::endl;
	std::cout << Pd << std::endl;
	std::cout << "origin of the object: " << std::endl;
	for (int i=0; i<Pd.getoriginsize(); ++i)
		std::cout<< Pd.orig(i) <<" - ";
	std::cout<<std::endl;
	std::cout << "scaling factor:  " << std::endl;
	for (int i=0; i<Pd.getscalingsize(); ++i)
		std::cout<< Pd.scal(i) <<" - ";
	std::cout<<std::endl;
	// std::cout << "print tolerance:  " << Pd.gettolerance() << std::endl;
	// std::cout << "change tolerance to 0.01:  " << std::endl;
	// Pd.settolerance(0.01);
	// std::cout << "print new tolerance:  " << Pd.gettolerance() << std::endl;
	// std::cout << "print minimum difference between coordinates:  " << Pd.getmindiff() << std::endl;
	// std::cout << "change minimum difference between coordinates to 0.003:  " << std::endl;
	// Pd.setmindiff(0.003);
	// std::cout << "print new minimum difference between coordinates:  " << Pd.getmindiff() << std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 9th domain (constructor from coordinates (x1 from 0 to 9, x2 from 1 to 10), shape = Triangle, NDIMP = 2):  " <<std::endl;
	std::cout << Td << std::endl;
	std::cout << "origin of the object: " << std::endl;
	for (int i=0; i<Td.getoriginsize(); ++i)
		std::cout<< Td.orig(i) <<" - ";
	std::cout<<std::endl;
	std::cout << "scaling factor:  " << std::endl;
	for (int i=0; i<Td.getscalingsize(); ++i)
		std::cout<< Td.scal(i) <<" - ";
	std::cout<<std::endl;
	// std::cout << "print tolerance:  " << Td.gettolerance() << std::endl;
	// std::cout << "change tolerance to 0.01:  " << std::endl;
	// Td.settolerance(0.01);
	// std::cout << "print new tolerance:  " << Td.gettolerance() << std::endl;
	// std::cout << "print minimum difference between coordinates:  " << Td.getmindiff() << std::endl;
	// std::cout << "change minimum difference between coordinates to 0.003:  " << std::endl;
	// Td.setmindiff(0.003);
	// std::cout << "print new minimum difference between coordinates:  " << Td.getmindiff() << std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 10th domain (constructor from coordinates (x1 from 0 to 9, x2 from 1 to 10), shape = Box with Tiangle, NDIMP = 2):  " <<std::endl;
	std::cout << Bd << std::endl;
	std::cout << "origin of the object: " << std::endl;
	for (int i=0; i<Bd.getoriginsize(); ++i)
		std::cout<< Bd.orig(i) <<" - ";
	std::cout<<std::endl;
	std::cout << "scaling factor:  " << std::endl;
	for (int i=0; i<Bd.getscalingsize(); ++i)
		std::cout<< Bd.scal(i) <<" - ";
	std::cout<<std::endl;
	// std::cout << "print tolerance:  " << Bd.gettolerance() << std::endl;
	// std::cout << "change tolerance to 0.01:  " << std::endl;
	// Bd.settolerance(0.01);
	// std::cout << "print new tolerance:  " << Bd.gettolerance() << std::endl;
	// std::cout << "print minimum difference between coordinates:  " << Bd.getmindiff() << std::endl;
	// std::cout << "change minimum difference between coordinates to 0.003:  " << std::endl;
	// Bd.setmindiff(0.003);
	// std::cout << "print new minimum difference between coordinates:  " << Bd.getmindiff() << std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 11th domain (constructor from coordinates (x1 from 0 to 9, x2 from 1 to 10, x3 from 2 to 11), shape = Point, NDIMP = 3):  " <<std::endl;
	std::cout << Pe << std::endl;
	std::cout << "origin of the object: " << std::endl;
	for (int i=0; i<Pe.getoriginsize(); ++i)
		std::cout<< Pe.orig(i) <<" - ";
	std::cout<<std::endl;
	std::cout << "scaling factor:  " << std::endl;
	for (int i=0; i<Pe.getscalingsize(); ++i)
		std::cout<< Pe.scal(i) <<" - ";
	std::cout<<std::endl;
	// std::cout << "print tolerance:  " << Pe.gettolerance() << std::endl;
	// std::cout << "change tolerance to 0.01:  " << std::endl;
	// Pe.settolerance(0.01);
	// std::cout << "print new tolerance:  " << Pe.gettolerance() << std::endl;
	// std::cout << "print minimum difference between coordinates:  " << Pe.getmindiff() << std::endl;
	// std::cout << "change minimum difference between coordinates to 0.003:  " << std::endl;
	// Pe.setmindiff(0.003);
	// std::cout << "print new minimum difference between coordinates:  " << Pe.getmindiff() << std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 12th domain (constructor from coordinates (x1 from 0 to 9, x2 from 1 to 10, x3 from 2 to 11), shape = Triangle, NDIMP = 3):  " <<std::endl;
	std::cout << Te << std::endl;
	std::cout << "origin of the object: " << std::endl;
	for (int i=0; i<Te.getoriginsize(); ++i)
		std::cout<< Te.orig(i) <<" - ";
	std::cout<<std::endl;
	std::cout << "scaling factor:  " << std::endl;
	for (int i=0; i<Te.getscalingsize(); ++i)
		std::cout<< Te.scal(i) <<" - ";
	std::cout<<std::endl;
	// std::cout << "print tolerance:  " << Te.gettolerance() << std::endl;
	// std::cout << "change tolerance to 0.01:  " << std::endl;
	// Te.settolerance(0.01);
	// std::cout << "print new tolerance:  " << Te.gettolerance() << std::endl;
	// std::cout << "print minimum difference between coordinates:  " << Te.getmindiff() << std::endl;
	// std::cout << "change minimum difference between coordinates to 0.003:  " << std::endl;
	// Te.setmindiff(0.003);
	// std::cout << "print new minimum difference between coordinates:  " << Te.getmindiff() << std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 13rd domain (constructor from coordinates (x1 from 0 to 9, x2 from 1 to 10, x3 from 2 to 11), shape = Box with Tiangle, NDIMP = 3):  " <<std::endl;
	std::cout << Be << std::endl;
	std::cout << "origin of the object: " << std::endl;
	for (int i=0; i<Be.getoriginsize(); ++i)
		std::cout<< Be.orig(i) <<" - ";
	std::cout<<std::endl;
	std::cout << "scaling factor:  " << std::endl;
	for (int i=0; i<Be.getscalingsize(); ++i)
		std::cout<< Be.scal(i) <<" - ";
	std::cout<<std::endl;
	// std::cout << "print tolerance:  " << Be.gettolerance() << std::endl;
	// std::cout << "change tolerance to 0.01:  " << std::endl;
	// Be.settolerance(0.01);
	// std::cout << "print new tolerance:  " << Be.gettolerance() << std::endl;
	// std::cout << "print minimum difference between coordinates:  " << Be.getmindiff() << std::endl;
	// std::cout << "change minimum difference between coordinates to 0.003:  " << std::endl;
	// Be.setmindiff(0.003);
	// std::cout << "print new minimum difference between coordinates:  " << Be.getmindiff() << std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 14th domain (constructor from coordinates (x1 from 0 to 9, x2 from 1 to 10, x3 from 2 to 11), shape = Tetrahedron, NDIMP = 3):  " <<std::endl;
	std::cout << Tf << std::endl;
	std::cout << "origin of the object: " << std::endl;
	for (int i=0; i<Tf.getoriginsize(); ++i)
		std::cout<< Tf.orig(i) <<" - ";
	std::cout<<std::endl;
	std::cout << "scaling factor:  " << std::endl;
	for (int i=0; i<Tf.getscalingsize(); ++i)
		std::cout<< Tf.scal(i) <<" - ";
	std::cout<<std::endl;
	// std::cout << "print tolerance:  " << Tf.gettolerance() << std::endl;
	// std::cout << "change tolerance to 0.01:  " << std::endl;
	// Tf.settolerance(0.01);
	// std::cout << "print new tolerance:  " << Tf.gettolerance() << std::endl;
	// std::cout << "print minimum difference between coordinates:  " << Tf.getmindiff() << std::endl;
	// std::cout << "change minimum difference between coordinates to 0.003:  " << std::endl;
	// Tf.setmindiff(0.003);
	// std::cout << "print new minimum difference between coordinates:  " << Tf.getmindiff() << std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 15th domain (constructor from coordinates (x1 from 0 to 9, x2 from 1 to 10, x3 from 2 to 11), shape = Box with Tetrahedron, NDIMP = 3):  " <<std::endl;
	std::cout << Bf << std::endl;
	std::cout << "origin of the object: " << std::endl;
	for (int i=0; i<Bf.getoriginsize(); ++i)
		std::cout<< Bf.orig(i) <<" - ";
	std::cout<<std::endl;
	std::cout << "scaling factor:  " << std::endl;
	for (int i=0; i<Bf.getscalingsize(); ++i)
		std::cout<< Bf.scal(i) <<" - ";
	std::cout<<std::endl;
	// std::cout << "print tolerance:  " << Bf.gettolerance() << std::endl;
	// std::cout << "change tolerance to 0.01:  " << std::endl;
	// Bf.settolerance(0.01);
	// std::cout << "print new tolerance:  " << Bf.gettolerance() << std::endl;
	// std::cout << "print minimum difference between coordinates:  " << Bf.getmindiff() << std::endl;
	// std::cout << "change minimum difference between coordinates to 0.003:  " << std::endl;
	// Bf.setmindiff(0.003);
	// std::cout << "print new minimum difference between coordinates:  " << Bf.getmindiff() << std::endl;
	std::cout << std::endl << std::endl;

	*n = 0;
	}
	catch(...)
	{
	std::cout <<"Exception caugth " << std::endl;
	*n = 1;
	}

}
}


extern "C"{
void Test_TreeNode_c(int *n){

	try{

	Point a(1,2,1,1); //Point with Id, bcId and coordinates
	Point b(2,3,3,1); //Point with Id, bcId and coordinates
	Point c(3,4,2,2); //Point with Id, bcId and coordinates

	std::vector<Point> vec{a,b,c};

	Point e(1,1,1,1,1);
	Point f(1,1,3,1,1);
	Point g(1,1,2,1,3);
	Point h(1,1,4,1,7);
	std::vector<Point> vec2{e,f,g};
	std::vector<Point> vec3{e,f,g,h};

	Element<3,2,2> triangle(1, vec);
	Element<3,2,3> triangle2(2, vec2);
	Element<3,3,3> tetrahedron(3, vec3);

	std::vector<Real> coord{0,1,2,3}; //used below to set new coord
	// {100, 24, 35, 21}; ************doesn't validate xmin, ymin, xmax, ymax
	std::vector<Real> coord1 {1,1,4,3};
	// {4,3,1,1};  ************doesn't validate xmin, ymin, xmax, ymax
	std::vector<Real> coord2 {1,1,1,4,3,5};
	// {4,3,5,1,1,1}; ************doesn't validate xmin, ymin, xmax, ymax
	std::vector<Real> coord3{2,2,2,3,7,8}; //used below to set new coord
	// ************doesn't validate xmin, ymin, xmax, ymax
	Box<2> tmp;
	Box<3> tmp2;
	Box<2> ba(coord1);
	Box<3> bb(coord2);


	//default constructor
	TreeNode<Box<2>> treea; //1st
	TreeNode<Box<3>> treeb; //2nd

	TreeNode<Element<3,2,2>> treec(1, triangle); //3rd
	TreeNode<Element<3,2,3>> treed(2, triangle2); //4th
	TreeNode<Element<3,3,3>> treee(3, tetrahedron); //5th


	TreeNode<Box<2>> treef(2, ba); //6th
	TreeNode<Box<3>> treeg(3, bb); //7th


	std::cout << "Print the 1st TreeNode (constructor default, NDIMP = 2):  " <<std::endl;
	treea.print(std::cout);
	// std::cout << "father of the treenode: " << std::endl;
	// std::cout << treea.getfather() <<std::endl;
	// std::cout << "change father to 3: " << std::endl;
	// treea.setfather(3);
	// std::cout << treea.getfather() <<std::endl;
	std::cout << "left child of the treenode:  " << std::endl;
	std::cout << treea.getchild(0) <<std::endl;
	std::cout << "change left child of the treenode to 1:  " << std::endl;
	treea.setchild(0,1);
	std::cout << treea.getchild(0) <<std::endl;
	std::cout << "right child of the treenode:  " << std::endl;
	std::cout << treea.getchild(1) <<std::endl;
	std::cout << "change right child of the treenode to 5:  " << std::endl;
	treea.setchild(1,5);
	std::cout << treea.getchild(1) <<std::endl;
	std::cout << "Box coordinates of the treenode:  " << std::endl;
	std::cout << "Min_Point:  ( " << treea.getcoord(0) << " , " << treea.getcoord(1) << " )" <<std::endl;
	std::cout << "Max_Point:  ( " << treea.getcoord(2) << " , " << treea.getcoord(3) << " )" <<std::endl;
	std::cout << "change the coordinates to {0,1,2,3}:  " << std::endl;
	treea.setcoords(coord);
	std::cout << "Min_Point:  ( " << treea.getcoord(0) << " , " << treea.getcoord(1) << " )" <<std::endl;
	std::cout << "Max_Point:  ( " << treea.getcoord(2) << " , " << treea.getcoord(3) << " )" <<std::endl;
	std::cout << "id of the treenode:  " << std::endl;
	std::cout << treea.getid() <<std::endl;
	std::cout << "change id of the treenode to 10:  " << std::endl;
	treea.setid(10);
	std::cout << treea.getid() <<std::endl;
	std::cout<<std::endl;
	std::cout << "get Box:  " << std::endl;
	tmp = treea.getbox();
	tmp.print(std::cout);


	std::cout << "Print the 2nd TreeNode (constructor default, NDIMP = 3):  " <<std::endl;
	treeb.print(std::cout);
	// std::cout << "father of the treenode: " << std::endl;
	// std::cout << treeb.getfather() <<std::endl;
	// std::cout << "change father: " << std::endl;
	// treeb.setfather(3);
	// std::cout << treeb.getfather() <<std::endl;
	std::cout << "left child of the treenode:  " << std::endl;
	std::cout << treeb.getchild(0) <<std::endl;
	std::cout << "change left child of the treenode:  " << std::endl;
	treeb.setchild(0,1);
	std::cout << treeb.getchild(0) <<std::endl;
	std::cout << "right child of the treenode:  " << std::endl;
	std::cout << treeb.getchild(1) <<std::endl;
	std::cout << "change right child of the treenode:  " << std::endl;
	treeb.setchild(1,5);
	std::cout << treeb.getchild(1) <<std::endl;
	std::cout << "Box coordinates of the treenode:  " << std::endl;
	std::cout << "Min_Point:  ( " << treeb.getcoord(0) << " , " << treeb.getcoord(1) << " , " << treeb.getcoord(2)<< " )" <<std::endl;
	std::cout << "Max_Point:  ( " << treeb.getcoord(3) << " , " << treeb.getcoord(4)<< " , " << treeb.getcoord(5) << " )" <<std::endl;
	std::cout << "change the coordinates:  " << std::endl;
	treeb.setcoords(coord3);
	std::cout << "Min_Point:  ( " << treeb.getcoord(0) << " , " << treeb.getcoord(1) << " , " << treeb.getcoord(2)<< " )" <<std::endl;
	std::cout << "Max_Point:  ( " << treeb.getcoord(3) << " , " << treeb.getcoord(4)<< " , " << treeb.getcoord(5) << " )" <<std::endl;
	std::cout << "id of the treenode:  " << std::endl;
	std::cout << treeb.getid() <<std::endl;
	std::cout << "change id of the treenode:  " << std::endl;
	treeb.setid(10);
	std::cout << treeb.getid() <<std::endl;
	std::cout << "get Box:  " << std::endl;
	tmp2 = treeb.getbox();
	tmp2.print(std::cout);


	std::cout << "Print the 3rd TreeNode (constructor from Triangle, points a(1,1), b(3,1), c(2,2), NDIMP = 2):  " <<std::endl;
	treec.print(std::cout);
	// std::cout << "father of the treenode: " << std::endl;
	// std::cout << treec.getfather() <<std::endl;
	// std::cout << "change father: " << std::endl;
	// treec.setfather(3);
	// std::cout << treec.getfather() <<std::endl;
	std::cout << "left child of the treenode:  " << std::endl;
	std::cout << treec.getchild(0) <<std::endl;
	std::cout << "change left child of the treenode:  " << std::endl;
	treec.setchild(0,1);
	std::cout << treec.getchild(0) <<std::endl;
	std::cout << "right child of the treenode:  " << std::endl;
	std::cout << treec.getchild(1) <<std::endl;
	std::cout << "change right child of the treenode:  " << std::endl;
	treec.setchild(1,5);
	std::cout << treec.getchild(1) <<std::endl;
	std::cout << "Box coordinates of the treenode:  " << std::endl;
	std::cout << "Min_Point:  ( " << treec.getcoord(0) << " , " << treec.getcoord(1) << " )" <<std::endl;
	std::cout << "Max_Point:  ( " << treec.getcoord(2) << " , " << treec.getcoord(3) << " )" <<std::endl;
	std::cout << "change the coordinates:  " << std::endl;
	treec.setcoords(coord);
	std::cout << "Min_Point:  ( " << treec.getcoord(0) << " , " << treec.getcoord(1) << " )" <<std::endl;
	std::cout << "Max_Point:  ( " << treec.getcoord(2) << " , " << treec.getcoord(3) << " )" <<std::endl;
	std::cout << "id of the treenode:  " << std::endl;
	std::cout << treec.getid() <<std::endl;
	std::cout << "change id of the treenode:  " << std::endl;
	treec.setid(10);
	std::cout << treec.getid() <<std::endl;
	std::cout << "get Box:  " << std::endl;
	tmp = treec.getbox();
	tmp.print(std::cout);


	std::cout << "Print the 4th TreeNode (constructor from Triangle, points a(1,1,1), b(3,1,1), c(2,1,3), NDIMP = 3):  " <<std::endl;
	treed.print(std::cout);
	// std::cout << "father of the treenode: " << std::endl;
	// std::cout << treed.getfather() <<std::endl;
	// std::cout << "change father: " << std::endl;
	// treed.setfather(3);
	// std::cout << treed.getfather() <<std::endl;
	std::cout << "left child of the treenode:  " << std::endl;
	std::cout << treed.getchild(0) <<std::endl;
	std::cout << "change left child of the treenode:  " << std::endl;
	treed.setchild(0,1);
	std::cout << treed.getchild(0) <<std::endl;
	std::cout << "right child of the treenode:  " << std::endl;
	std::cout << treed.getchild(1) <<std::endl;
	std::cout << "change right child of the treenode:  " << std::endl;
	treed.setchild(1,5);
	std::cout << treed.getchild(1) <<std::endl;
	std::cout << "Box coordinates of the treenode:  " << std::endl;
	std::cout << "Min_Point:  ( " << treed.getcoord(0) << " , " << treed.getcoord(1) << " )" <<std::endl;
	std::cout << "Max_Point:  ( " << treed.getcoord(2) << " , " << treed.getcoord(3) << " )" <<std::endl;
	std::cout << "change the coordinates:  " << std::endl;
	treed.setcoords(coord);
	std::cout << "Min_Point:  ( " << treed.getcoord(0) << " , " << treed.getcoord(1) << " )" <<std::endl;
	std::cout << "Max_Point:  ( " << treed.getcoord(2) << " , " << treed.getcoord(3) << " )" <<std::endl;
	std::cout << "id of the treenode:  " << std::endl;
	std::cout << treed.getid() <<std::endl;
	std::cout << "change id of the treenode:  " << std::endl;
	treed.setid(10);
	std::cout << treed.getid() <<std::endl;
	std::cout << "get Box:  " << std::endl;
	tmp2 = treed.getbox();
	tmp2.print(std::cout);


	std::cout << "Print the 5th TreeNode (constructor from Tetrahedron, points a(1,1,1), b(3,1,1), c(2,1,3), d(4,1,7) NDIMP = 3):  " <<std::endl;
	treee.print(std::cout);
	// std::cout << "father of the treenode: " << std::endl;
	// std::cout << treee.getfather() <<std::endl;
	// std::cout << "change father: " << std::endl;
	// treee.setfather(3);
	// std::cout << treee.getfather() <<std::endl;
	std::cout << "left child of the treenode:  " << std::endl;
	std::cout << treee.getchild(0) <<std::endl;
	std::cout << "change left child of the treenode:  " << std::endl;
	treee.setchild(0,1);
	std::cout << treee.getchild(0) <<std::endl;
	std::cout << "right child of the treenode:  " << std::endl;
	std::cout << treee.getchild(1) <<std::endl;
	std::cout << "change right child of the treenode:  " << std::endl;
	treee.setchild(1,5);
	std::cout << treee.getchild(1) <<std::endl;
	std::cout << "Box coordinates of the treenode:  " << std::endl;
	std::cout << "Min_Point:  ( " << treee.getcoord(0) << " , " << treee.getcoord(1) << " )" <<std::endl;
	std::cout << "Max_Point:  ( " << treee.getcoord(2) << " , " << treee.getcoord(3) << " )" <<std::endl;
	std::cout << "change the coordinates:  " << std::endl;
	treee.setcoords(coord);
	std::cout << "Min_Point:  ( " << treee.getcoord(0) << " , " << treee.getcoord(1) << " )" <<std::endl;
	std::cout << "Max_Point:  ( " << treee.getcoord(2) << " , " << treee.getcoord(3) << " )" <<std::endl;
	std::cout << "id of the treenode:  " << std::endl;
	std::cout << treee.getid() <<std::endl;
	std::cout << "change id of the treenode:  " << std::endl;
	treee.setid(10);
	std::cout << treee.getid() <<std::endl;
	std::cout << "get Box:  " << std::endl;
	tmp2 = treee.getbox();
	tmp2.print(std::cout);


	std::cout << "Print the 6th TreeNode (constructor from Box, {1,1,4,3}, NDIMP = 2):  " <<std::endl;
	treef.print(std::cout);
	// std::cout << "father of the treenode: " << std::endl;
	// std::cout << treef.getfather() <<std::endl;
	// std::cout << "change father: " << std::endl;
	// treef.setfather(3);
	// std::cout << treef.getfather() <<std::endl;
	std::cout << "left child of the treenode:  " << std::endl;
	std::cout << treef.getchild(0) <<std::endl;
	std::cout << "change left child of the treenode:  " << std::endl;
	treef.setchild(0,1);
	std::cout << treef.getchild(0) <<std::endl;
	std::cout << "right child of the treenode:  " << std::endl;
	std::cout << treef.getchild(1) <<std::endl;
	std::cout << "change right child of the treenode:  " << std::endl;
	treef.setchild(1,5);
	std::cout << treef.getchild(1) <<std::endl;
	std::cout << "Box coordinates of the treenode:  " << std::endl;
	std::cout << "Min_Point:  ( " << treef.getcoord(0) << " , " << treef.getcoord(1) << " )" <<std::endl;
	std::cout << "Max_Point:  ( " << treef.getcoord(2) << " , " << treef.getcoord(3) << " )" <<std::endl;
	std::cout << "change the coordinates:  " << std::endl;
	treef.setcoords(coord);
	std::cout << "Min_Point:  ( " << treef.getcoord(0) << " , " << treef.getcoord(1) << " )" <<std::endl;
	std::cout << "Max_Point:  ( " << treef.getcoord(2) << " , " << treef.getcoord(3) << " )" <<std::endl;
	std::cout << "id of the treenode:  " << std::endl;
	std::cout << treef.getid() <<std::endl;
	std::cout << "change id of the treenode:  " << std::endl;
	treef.setid(10);
	std::cout << treef.getid() <<std::endl;
	std::cout << "get Box:  " << std::endl;
	tmp = treef.getbox();
	tmp.print(std::cout);



	std::cout << "Print the 7th TreeNode (constructor from Box {1,1,1,4,3,5}, NDIMP = 3):  " <<std::endl;
	treeg.print(std::cout);
	// std::cout << "father of the treenode: " << std::endl;
	// std::cout << treeg.getfather() <<std::endl;
	// std::cout << "change father: " << std::endl;
	// treeg.setfather(3);
	// std::cout << treeg.getfather() <<std::endl;
	std::cout << "left child of the treenode:  " << std::endl;
	std::cout << treeg.getchild(0) <<std::endl;
	std::cout << "change left child of the treenode:  " << std::endl;
	treeg.setchild(0,1);
	std::cout << treeg.getchild(0) <<std::endl;
	std::cout << "right child of the treenode:  " << std::endl;
	std::cout << treeg.getchild(1) <<std::endl;
	std::cout << "change right child of the treenode:  " << std::endl;
	treeg.setchild(1,5);
	std::cout << treeg.getchild(1) <<std::endl;
	std::cout << "Box coordinates of the treenode:  " << std::endl;
	std::cout << "Min_Point:  ( " << treeg.getcoord(0) << " , " << treeg.getcoord(1) << " , " << treeg.getcoord(2)<< " )" <<std::endl;
	std::cout << "Max_Point:  ( " << treeg.getcoord(3) << " , " << treeg.getcoord(4)<< " , " << treeg.getcoord(5) << " )" <<std::endl;
	std::cout << "change the coordinates to {2,2,2,3,7,8}:  " << std::endl;
	treeg.setcoords(coord3);
	std::cout << "Min_Point:  ( " << treeg.getcoord(0) << " , " << treeg.getcoord(1) << " , " << treeg.getcoord(2)<< " )" <<std::endl;
	std::cout << "Max_Point:  ( " << treeg.getcoord(3) << " , " << treeg.getcoord(4)<< " , " << treeg.getcoord(5) << " )" <<std::endl;
	std::cout << "id of the treenode:  " << std::endl;
	std::cout << treeg.getid() <<std::endl;
	std::cout << "change id of the treenode:  " << std::endl;
	treeg.setid(10);
	std::cout << treeg.getid() <<std::endl;
	std::cout << "get Box:  " << std::endl;
	tmp2 = treeg.getbox();
	tmp2.print(std::cout);

	*n = 0;
	}
	catch(...)
	{
	std::cout <<"Exception caugth " << std::endl;
	*n = 1;
	}
}
}

extern "C"{
void Test_TreeHeader_c(int *n){

	try{
	std::vector<std::vector<Real>> coord;
	coord.resize(2);
	for (int i = 0; i < 2; i++)
	{
		coord[i].resize(10);
		for(int j = 0; j < 10; j++)
			coord[i][j] = i+j;
	}
	// 0 1 2 3 4 5 6 7 8 9
	// 1 2 3 4 5 6 7 8 9 10

	std::vector<std::vector<Real>> coord2;
	coord2.resize(3);
	for (int i = 0; i < 3; i++)
	{
		coord2[i].resize(10);
		for(int j = 0; j < 10; j++)
			coord2[i][j] = i+j;
	}
	// 0 1 2 3 4 5 6 7 8  9
	// 1 2 3 4 5 6 7 8 9  10
	// 2 3 4 5 6 7 8 9 10 11

	Domain<Element<3,2,2>> mydomb(coord);
	Domain<Element<3,2,3>> mydomc(coord2);
	Domain<Element<3,3,3>> mydomd(coord2);

	Domain<Box<2>> mydome(coord);
	Domain<Box<3>> mydomf(coord2);

	//Default consturctor
	TreeHeader<Element<3,2,2>> THa; //1st
	TreeHeader<Element<3,2,3>> THb; //2nd
	TreeHeader<Element<3,3,3>> THc; //3rd
	TreeHeader<Box<2>> THd; //4th
	TreeHeader<Box<3>> THe; //5th

	//Default consturctor with createtreeheader function
	TreeHeader<Element<3,2,2>> THf = createtreeheader(1, mydomb); //6th
	TreeHeader<Element<3,2,3>> THg = createtreeheader(2, mydomc); //7th
	TreeHeader<Element<3,3,3>> THh = createtreeheader(3, mydomd); //8th
	TreeHeader<Box<2>> THi = createtreeheader(4, mydome); //9th
	TreeHeader<Box<3>> THj = createtreeheader(5, mydomf); //10th

	std::cout << "Print the 1st tree_header (constructor default, shape = Triangle, NDIMP=2):  " <<std::endl;
	std::cout << THa <<std::endl;
	std::cout << "number of tree memory locations: " << std::endl;
	std::cout  << THa.gettreeloc() <<std::endl;
	std::cout << "change number of tree memory locations to 3: " << std::endl;
	THa.settreeloc(3);
	std::cout  << THa.gettreeloc() <<std::endl;
	std::cout << "number of tree levels: " << std::endl;
	std::cout  << THa.gettreelev() <<std::endl;
	std::cout << "change number of tree levels to 4: " << std::endl;
	THa.settreelev(4);
	std::cout  << THa.gettreelev() <<std::endl;
	std::cout << "physical dimension: " << std::endl;
	std::cout  << THa.getndimp() <<std::endl;
	std::cout << "search dimension: " << std::endl;
	std::cout  << THa.getndimt() <<std::endl;
	std::cout << "number of location currently used: " << std::endl;
	std::cout  << THa.getnele() <<std::endl;
	std::cout << "change number of location currently used to 2: " << std::endl;
	THa.setnele(2);
	std::cout  << THa.getnele() <<std::endl;
	std::cout << "number of next available location: " << std::endl;
	std::cout  << THa.getiava() <<std::endl;
	std::cout << "change number of next available location to 3: " << std::endl;
	THa.setiava(3);
	std::cout  << THa.getiava() <<std::endl;
	std::cout << "number next available location in the tree free store: " << std::endl;
	std::cout  << THa.getiend() <<std::endl;
	std::cout << "change number next available location in the tree free store to 4: " << std::endl;
	THa.setiend(4);
	std::cout  << THa.getiend() <<std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 2nd tree_header (constructor default, shape = Triangle, NDIMP=3):  " <<std::endl;
	std::cout << THb <<std::endl;
	std::cout << "number of tree memory locations: " << std::endl;
	std::cout  << THb.gettreeloc() <<std::endl;
	std::cout << "change number of tree memory locations to 3: " << std::endl;
	THb.settreeloc(3);
	std::cout  << THb.gettreeloc() <<std::endl;
	std::cout << "number of tree levels: " << std::endl;
	std::cout  << THb.gettreelev() <<std::endl;
	std::cout << "change number of tree levels to 4: " << std::endl;
	THb.settreelev(4);
	std::cout  << THb.gettreelev() <<std::endl;
	std::cout << "physical dimension: " << std::endl;
	std::cout  << THb.getndimp() <<std::endl;
	std::cout << "search dimension: " << std::endl;
	std::cout  << THb.getndimt() <<std::endl;
	std::cout << "number of location currently used: " << std::endl;
	std::cout  << THb.getnele() <<std::endl;
	std::cout << "change number of location currently used to 2: " << std::endl;
	THb.setnele(2);
	std::cout  << THb.getnele() <<std::endl;
	std::cout << "number of next available location: " << std::endl;
	std::cout  << THb.getiava() <<std::endl;
	std::cout << "change number of next available location to 3: " << std::endl;
	THb.setiava(3);
	std::cout  << THb.getiava() <<std::endl;
	std::cout << "number next available location in the tree free store: " << std::endl;
	std::cout  << THb.getiend() <<std::endl;
	std::cout << "change number next available location in the tree free store to 4: " << std::endl;
	THb.setiend(4);
	std::cout  << THb.getiend() <<std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 3rd tree_header (constructor default, shape = Tetrahedron, NDIMP=3):  " <<std::endl;
	std::cout << THc <<std::endl;
	std::cout << "number of tree memory locations: " << std::endl;
	std::cout  << THc.gettreeloc() <<std::endl;
	std::cout << "change number of tree memory locations to 3: " << std::endl;
	THc.settreeloc(3);
	std::cout  << THc.gettreeloc() <<std::endl;
	std::cout << "number of tree levels: " << std::endl;
	std::cout  << THc.gettreelev() <<std::endl;
	std::cout << "change number of tree levels to 4: " << std::endl;
	THc.settreelev(4);
	std::cout  << THc.gettreelev() <<std::endl;
	std::cout << "physical dimension: " << std::endl;
	std::cout  << THc.getndimp() <<std::endl;
	std::cout << "search dimension: " << std::endl;
	std::cout  << THc.getndimt() <<std::endl;
	std::cout << "number of location currently used: " << std::endl;
	std::cout  << THc.getnele() <<std::endl;
	std::cout << "change number of location currently used to 2: " << std::endl;
	THc.setnele(2);
	std::cout  << THc.getnele() <<std::endl;
	std::cout << "number of next available location: " << std::endl;
	std::cout  << THc.getiava() <<std::endl;
	std::cout << "change number of next available location to 3: " << std::endl;
	THc.setiava(3);
	std::cout  << THc.getiava() <<std::endl;
	std::cout << "number next available location in the tree free store: " << std::endl;
	std::cout  << THc.getiend() <<std::endl;
	std::cout << "change number next available location in the tree free store to 4: " << std::endl;
	THc.setiend(4);
	std::cout  << THc.getiend() <<std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 4th tree_header (constructor default, shape = Box, NDIMP=2):  " <<std::endl;
	std::cout << THd <<std::endl;
	std::cout << "number of tree memory locations: " << std::endl;
	std::cout  << THd.gettreeloc() <<std::endl;
	std::cout << "change number of tree memory locations: " << std::endl;
	THd.settreeloc(3);
	std::cout  << THd.gettreeloc() <<std::endl;
	std::cout << "number of tree levels: " << std::endl;
	std::cout  << THd.gettreelev() <<std::endl;
	std::cout << "change number of tree levels: " << std::endl;
	THd.settreelev(4);
	std::cout  << THd.gettreelev() <<std::endl;
	std::cout << "physical dimension: " << std::endl;
	std::cout  << THd.getndimp() <<std::endl;
	std::cout << "search dimension: " << std::endl;
	std::cout  << THd.getndimt() <<std::endl;
	std::cout << "number of location currently used: " << std::endl;
	std::cout  << THd.getnele() <<std::endl;
	std::cout << "change number of location currently used: " << std::endl;
	THd.setnele(2);
	std::cout  << THd.getnele() <<std::endl;
	std::cout << "number of next available location: " << std::endl;
	std::cout  << THd.getiava() <<std::endl;
	std::cout << "change number of next available location: " << std::endl;
	THd.setiava(3);
	std::cout  << THd.getiava() <<std::endl;
	std::cout << "number next available location in the tree free store: " << std::endl;
	std::cout  << THd.getiend() <<std::endl;
	std::cout << "change number next available location in the tree free store: " << std::endl;
	THd.setiend(4);
	std::cout  << THd.getiend() <<std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 5th tree_header (constructor default, shape = Box, NDIMP=3):  " <<std::endl;
	std::cout << THe <<std::endl;
	std::cout << "number of tree memory locations: " << std::endl;
	std::cout  << THe.gettreeloc() <<std::endl;
	std::cout << "change number of tree memory locations: " << std::endl;
	THe.settreeloc(3);
	std::cout  << THe.gettreeloc() <<std::endl;
	std::cout << "number of tree levels: " << std::endl;
	std::cout  << THe.gettreelev() <<std::endl;
	std::cout << "change number of tree levels: " << std::endl;
	THe.settreelev(4);
	std::cout  << THe.gettreelev() <<std::endl;
	std::cout << "physical dimension: " << std::endl;
	std::cout  << THe.getndimp() <<std::endl;
	std::cout << "search dimension: " << std::endl;
	std::cout  << THe.getndimt() <<std::endl;
	std::cout << "number of location currently used: " << std::endl;
	std::cout  << THe.getnele() <<std::endl;
	std::cout << "change number of location currently used: " << std::endl;
	THe.setnele(2);
	std::cout  << THe.getnele() <<std::endl;
	std::cout << "number of next available location: " << std::endl;
	std::cout  << THe.getiava() <<std::endl;
	std::cout << "change number of next available location: " << std::endl;
	THe.setiava(3);
	std::cout  << THe.getiava() <<std::endl;
	std::cout << "number next available location in the tree free store: " << std::endl;
	std::cout  << THe.getiend() <<std::endl;
	std::cout << "change number next available location in the tree free store: " << std::endl;
	THe.setiend(4);
	std::cout  << THe.getiend() <<std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 6th tree_header (cconstructor from domain of shape = Triangle, NDIMP=2):  " <<std::endl;
	std::cout << THf <<std::endl;
	std::cout << "number of tree memory locations: " << std::endl;
	std::cout  << THf.gettreeloc() <<std::endl;
	std::cout << "change number of tree memory locations to 3: " << std::endl;
	THf.settreeloc(3);
	std::cout  << THf.gettreeloc() <<std::endl;
	std::cout << "number of tree levels: " << std::endl;
	std::cout  << THf.gettreelev() <<std::endl;
	std::cout << "change number of tree levels to 4: " << std::endl;
	THf.settreelev(4);
	std::cout  << THf.gettreelev() <<std::endl;
	std::cout << "physical dimension: " << std::endl;
	std::cout  << THf.getndimp() <<std::endl;
	std::cout << "search dimension: " << std::endl;
	std::cout  << THf.getndimt() <<std::endl;
	std::cout << "number of location currently used: " << std::endl;
	std::cout  << THf.getnele() <<std::endl;
	std::cout << "change number of location currently used to 2: " << std::endl;
	THf.setnele(2);
	std::cout  << THf.getnele() <<std::endl;
	std::cout << "number of next available location: " << std::endl;
	std::cout  << THf.getiava() <<std::endl;
	std::cout << "change number of next available location to 3: " << std::endl;
	THf.setiava(3);
	std::cout  << THf.getiava() <<std::endl;
	std::cout << "number next available location in the tree free store: " << std::endl;
	std::cout  << THf.getiend() <<std::endl;
	std::cout << "change number next available location in the tree free store to 4: " << std::endl;
	THf.setiend(4);
	std::cout  << THf.getiend() <<std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 7th tree_header (cconstructor from domain of shape = Triangle, NDIMP=3):  " <<std::endl;
	std::cout << THg <<std::endl;
	std::cout << "number of tree memory locations: " << std::endl;
	std::cout  << THg.gettreeloc() <<std::endl;
	std::cout << "change number of tree memory locations to 3: " << std::endl;
	THg.settreeloc(3);
	std::cout  << THg.gettreeloc() <<std::endl;
	std::cout << "number of tree levels: " << std::endl;
	std::cout  << THg.gettreelev() <<std::endl;
	std::cout << "change number of tree levels to 4: " << std::endl;
	THg.settreelev(4);
	std::cout  << THg.gettreelev() <<std::endl;
	std::cout << "physical dimension: " << std::endl;
	std::cout  << THg.getndimp() <<std::endl;
	std::cout << "search dimension: " << std::endl;
	std::cout  << THg.getndimt() <<std::endl;
	std::cout << "number of location currently used: " << std::endl;
	std::cout  << THg.getnele() <<std::endl;
	std::cout << "change number of location currently used to 2: " << std::endl;
	THg.setnele(2);
	std::cout  << THg.getnele() <<std::endl;
	std::cout << "number of next available location: " << std::endl;
	std::cout  << THg.getiava() <<std::endl;
	std::cout << "change number of next available location to 3: " << std::endl;
	THg.setiava(3);
	std::cout  << THg.getiava() <<std::endl;
	std::cout << "number next available location in the tree free store: " << std::endl;
	std::cout  << THg.getiend() <<std::endl;
	std::cout << "change number next available location in the tree free store to 4: " << std::endl;
	THg.setiend(4);
	std::cout  << THg.getiend() <<std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 8th tree_header (constructor from domain of shape = Tetrahedron, NDIMP=3):  " <<std::endl;
	std::cout << THh <<std::endl;
	std::cout << "number of tree memory locations: " << std::endl;
	std::cout  << THh.gettreeloc() <<std::endl;
	std::cout << "change number of tree memory locations to 3: " << std::endl;
	THh.settreeloc(3);
	std::cout  << THh.gettreeloc() <<std::endl;
	std::cout << "number of tree levels: " << std::endl;
	std::cout  << THh.gettreelev() <<std::endl;
	std::cout << "change number of tree levels to 4: " << std::endl;
	THh.settreelev(4);
	std::cout  << THh.gettreelev() <<std::endl;
	std::cout << "physical dimension: " << std::endl;
	std::cout  << THh.getndimp() <<std::endl;
	std::cout << "search dimension: " << std::endl;
	std::cout  << THh.getndimt() <<std::endl;
	std::cout << "number of location currently used: " << std::endl;
	std::cout  << THh.getnele() <<std::endl;
	std::cout << "change number of location currently used to 2: " << std::endl;
	THh.setnele(2);
	std::cout  << THh.getnele() <<std::endl;
	std::cout << "number of next available location: " << std::endl;
	std::cout  << THh.getiava() <<std::endl;
	std::cout << "change number of next available location to 3: " << std::endl;
	THh.setiava(3);
	std::cout  << THh.getiava() <<std::endl;
	std::cout << "number next available location in the tree free store: " << std::endl;
	std::cout  << THh.getiend() <<std::endl;
	std::cout << "change number next available location in the tree free store to 4: " << std::endl;
	THh.setiend(4);
	std::cout  << THh.getiend() <<std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 9th tree_header (constructor from domain of shape = Box, NDIMP=2):  " <<std::endl;
	std::cout << THi <<std::endl;
	std::cout << "number of tree memory locations: " << std::endl;
	std::cout  << THi.gettreeloc() <<std::endl;
	std::cout << "change number of tree memory locations: " << std::endl;
	THi.settreeloc(3);
	std::cout  << THi.gettreeloc() <<std::endl;
	std::cout << "number of tree levels: " << std::endl;
	std::cout  << THi.gettreelev() <<std::endl;
	std::cout << "change number of tree levels: " << std::endl;
	THi.settreelev(4);
	std::cout  << THi.gettreelev() <<std::endl;
	std::cout << "physical dimension: " << std::endl;
	std::cout  << THi.getndimp() <<std::endl;
	std::cout << "search dimension: " << std::endl;
	std::cout  << THi.getndimt() <<std::endl;
	std::cout << "number of location currently used: " << std::endl;
	std::cout  << THi.getnele() <<std::endl;
	std::cout << "change number of location currently used: " << std::endl;
	THi.setnele(2);
	std::cout  << THi.getnele() <<std::endl;
	std::cout << "number of next available location: " << std::endl;
	std::cout  << THi.getiava() <<std::endl;
	std::cout << "change number of next available location: " << std::endl;
	THi.setiava(3);
	std::cout  << THi.getiava() <<std::endl;
	std::cout << "number next available location in the tree free store: " << std::endl;
	std::cout  << THi.getiend() <<std::endl;
	std::cout << "change number next available location in the tree free store: " << std::endl;
	THi.setiend(4);
	std::cout  << THi.getiend() <<std::endl;
	std::cout << std::endl << std::endl;


	std::cout << "Print the 10th tree_header (constructor from domain of shape = Box, NDIMP=3):  " <<std::endl;
	std::cout << THj <<std::endl;
	std::cout << "number of tree memory locations: " << std::endl;
	std::cout  << THj.gettreeloc() <<std::endl;
	std::cout << "change number of tree memory locations: " << std::endl;
	THj.settreeloc(3);
	std::cout  << THj.gettreeloc() <<std::endl;
	std::cout << "number of tree levels: " << std::endl;
	std::cout  << THj.gettreelev() <<std::endl;
	std::cout << "change number of tree levels: " << std::endl;
	THj.settreelev(4);
	std::cout  << THj.gettreelev() <<std::endl;
	std::cout << "physical dimension: " << std::endl;
	std::cout  << THj.getndimp() <<std::endl;
	std::cout << "search dimension: " << std::endl;
	std::cout  << THj.getndimt() <<std::endl;
	std::cout << "number of location currently used: " << std::endl;
	std::cout  << THj.getnele() <<std::endl;
	std::cout << "change number of location currently used: " << std::endl;
	THj.setnele(2);
	std::cout  << THj.getnele() <<std::endl;
	std::cout << "number of next available location: " << std::endl;
	std::cout  << THj.getiava() <<std::endl;
	std::cout << "change number of next available location: " << std::endl;
	THj.setiava(3);
	std::cout  << THj.getiava() <<std::endl;
	std::cout << "number next available location in the tree free store: " << std::endl;
	std::cout  << THj.getiend() <<std::endl;
	std::cout << "change number next available location in the tree free store: " << std::endl;
	THj.setiend(4);
	std::cout  << THj.getiend() <<std::endl;
	std::cout << std::endl << std::endl;


	*n = 0;
	}
	catch(...)
	{
	std::cout <<"Exception caugth " << std::endl;
	*n = 1;
	}
}
}


extern "C"{
void Test_ADTree1_c(int *n){

	try{
	//costruisco una mesh semplice a mano salvata come quella che serve al costruttore
	/* the simple mesh is:
	 ___ ___
	|  /|  /|
	| / | / |
	|/__|/__|
	|  /|  /|
	| / | / |
	|/__|/__|

	nodes are numbered:
	6 7 8
	3 4 5
	0 1 2

	triangle are numbered:
	 ___ ___
	|4 /|6 /|
	| /5| /7|
	|/__|/__|
	|0 /|2 /|
	| /1| /3|
	|/__|/__|
	*/

	UInt num_nodes = 9;
	UInt num_triangle = 8;
	Real points[] = {0., 0.5, 1., 0.,  0.5, 1.,  0., 0.5, 1.,
					 0., 0.,  0., 0.5, 0.5, 0.5, 1., 1.,  1.}; //total 18
	UInt triangle[] = {0, 0, 1, 1, 3, 3, 4, 4,
					   4, 1, 5, 2, 7, 4, 8, 5,
					   3, 4, 4, 5, 6, 7, 7, 8}; //total 24

	std::vector<Real> coord = {0.3,0.3,0.3,0.6,0.6,0.3};
	std::vector<Real> coord2;
	Id id; //id2;
	std::vector<Real> region1 = {0.5,0.5,0.5,0.5};
	std::vector<Real> region2 = {1.,1.,1.,1.};
	std::vector<Real> region3 = {0.3,0.3, 0.3,0.3};
	std::vector<Real> region4 = {3,3, 3,3};
	std::set<int> found;
	int loc;


	//default tree_header
	// TreeHeader<Element<3,2,2>> THa;
	// TreeHeader<Element<3,2,3>> THb;
	// TreeHeader<Element<3,3,3>> THc;
	// TreeHeader<Box<2>> THd;
	// TreeHeader<Box<3>> THe;


	//default constructor
	// ADTree<Element<3,2,2>> ADTa; //1st
	// ADTree<Element<3,2,3>> ADTb; //2nd
	// ADTree<Element<3,3,3>> ADTc; //3rd
	// ADTree<Box<2>> ADTd; //4th
	// ADTree<Box<3>> ADTe; //5th

	//tree_header constructor
	// ADTree<Element<3,2,2>> ADTf(THa); //6th
	// ADTree<Element<3,2,3>> ADTg(THb); //7th
	// ADTree<Element<3,3,3>> ADTh(THc); //8th
	// ADTree<Box<2>> ADTi(THd); //9th
	// ADTree<Box<3>> ADTj(THe); //10th


	//*************only possible for triangle to have this kind of constructor
	//constructor from 2D/2.5D/3D mesh
	ADTree<Element<3,2,2>> ADTk(points, triangle, num_nodes, num_triangle); //11th
	// ADTree<Element<3,2,3>> ADTl(points, triangle, num_nodes, num_triangle); //12th
	// ADTree<Element<3,3,3>> ADTm(points, triangle, num_nodes, num_triangle); //13th
	///*******Can there be ADTree<Box<2>> or ADTree<Box<3>> as well??? I assume that
	//*********Shape will always be Element so I am using that logic at ADTree constructor




		std::cout << std::endl << std::endl;
	std::cout << "Print the 11th ADTree (construct from 2D mesh, Shape = Triangle, NDIMP = 2):  " <<std::endl;
	std::cout << ADTk <<std::endl;
	std::cout << "Print TreeHeader:  " <<std::endl;
	std::cout << ADTk.gettreeheader() <<std::endl;

	std::cout << std::endl;
	std::cout << "Print number of nodes (location) occupied:  " <<std::endl;
	std::cout << "Tree head + " << (ADTk.gettreeheader()).getnele() << " tree nodes" <<std::endl;
	std::cout << "Treenodes visualize as a vector:  " << std::endl;
	for (int i = 0; i < (ADTk.gettreeheader()).getnele()+1; i++) {
		(ADTk.gettreenode(i)).print(std::cout);
		//std::cout << "Shape id of father: " << ADTk.pointId(ADTk.gettreenode(i).getfather()) << std::endl;
		std::cout << "Shape id of left child: " << ADTk.pointId(ADTk.gettreenode(i).getchild(0)) << std::endl;
		std::cout << "Shape id of right child: " << ADTk.pointId(ADTk.gettreenode(i).getchild(1)) << std::endl;
	}

	std::cout << "visualize information about a node, for example node number 5: " << std::endl;
	ADTk.gettri(5, coord2, id);
	std::cout << "coordinate:  " ;
	for (size_t i = 0; i < coord2.size(); i++)
		std::cout << " - " << coord2[i] << " - ";
	std::cout << std::endl;
	std::cout << "Id :  " << id << std::endl;
	std::cout << "treenode:  " << std::endl;
	(ADTk.gettreenode(5)).print(std::cout);
	std::cout << "get box coord: " << std::endl;
	std::cout << "( " << ADTk.pointcoord(5+1, 0) << " , " << ADTk.pointcoord(5, 1) << " )" << std::endl;
	std::cout << "( " << ADTk.pointcoord(5+1, 2) << " , " << ADTk.pointcoord(5, 3) << " )" <<std::endl;
	std::cout << "get id: " << std::endl;
	std::cout <<  ADTk.pointId(5) <<std::endl <<std::endl;

	std::cout << "search point (0.5,0.5): " << ADTk.search(region1, found) << std::endl;
	for (std::set<int>::iterator i = found.begin(); i != found.end(); i++) {
  	 std::cout << " --- Location at: " << *i << std::endl;
	(ADTk.gettreenode(*i)).print(std::cout);
	}
	std::cout << std::endl << std::endl;

	std::cout << "search point (1,1): " << ADTk.search(region2, found) << std::endl;
	for (std::set<int>::iterator i = found.begin(); i != found.end(); i++) {
  	 std::cout << " --- Location at: " << *i << std::endl;
	(ADTk.gettreenode(*i)).print(std::cout);
	}
	std::cout << std::endl << std::endl;

	std::cout << "search point (0.3,0.3): " << ADTk.search(region3, found) << std::endl;
	for (std::set<int>::iterator i = found.begin(); i != found.end(); i++) {
  	 std::cout << " --- Location at: " << *i << std::endl;
	(ADTk.gettreenode(*i)).print(std::cout);
	}
	std::cout << std::endl << std::endl;

	std::cout << "search point (3,3): " << ADTk.search(region4, found) << std::endl;
	for (std::set<int>::iterator i = found.begin(); i != found.end(); i++) {
  	 std::cout << " --- Location at: " << *i << std::endl;
	(ADTk.gettreenode(*i)).print(std::cout);
	}
	std::cout << std::endl << std::endl;

	*n = 0;
	}
	catch(...)
	{
	std::cout <<"Exception caugth " << std::endl;
	*n = 1;
	}
}
}


extern "C"{
void Test_ADTree2_c(int *n){

	try{
	//costruisco una mesh semplice a mano salvata come quella che serve al costruttore
	/* the simple mesh is:
	 ___ ___ ___ ___
	|  /|  /|  /|  /|
	| / | / | / | / |
	|/__|/__|/__|/__|
	|  /|  /|  /|  /|
	| / | / | / | / |
	|/__|/__|/__|/__|
	|  /|  /|  /|  /|
	| / | / | / | / |
	|/__|/__|/__|/__|
	|  /|  /|  /|  /|
	| / | / | / | / |
	|/__|/__|/__|/__|


	nodes are numbered:
	20 21 22 23 24
	15 16 17 18 19
	10 11 12 13 14
	5  6  7  8  9
	0  1  2  3  4

	triangle are numbered:
	 ___ ___ ___ ___
	|24/|26/|28/|30/|
	| / | / | / | / |
	|/25|/27|/29|/31|
	|16/|18/|20/|22/|
	| / | / | / | / |
	|/17|/19|/21|/23|
	|8 /|10/|12/|14/|
	| / | / | / | / |
	|/_9|/11|/13|/15|
	|0 /|2 /|4 /|6 /|
	| / | / | / | / |
	|/_1|/_3|/_5|/_7|
*/
	UInt num_nodes = 25;
	UInt num_triangle = 32;
	Real points[] = {0., 0.25, 0.5, 0.75, 1., 0.,   0.25, 0.5,  0.75,  1.,   0.,  0.25, 0.5, 0.75, 1.,  0.,   0.25, 0.5,  0.75, 1.,   0., 0.25, 0.5, 0.75, 1.,
					 0., 0.,   0.,  0,    0,  0.25, 0.25, 0.25, 0.25,  0.25, 0.5, 0.5,  0.5, 0.5,  0.5, 0.75, 0.75, 0.75, 0.75, 0.75, 1,  1.,   1.,  1.,   1}; //total 25*2=50
	UInt triangle[] = {0, 0, 1, 1, 2, 2, 3, 3, 5,   5, 6,  6,  7, 7,  8,  8,  10, 10, 11, 11, 12, 12, 13, 13, 15, 15, 16, 16, 17, 17, 18, 18,
					   5, 1, 6, 2, 7, 3, 8, 4, 10,  6, 11, 7, 12, 8,  13, 9,  15, 11, 16, 12, 17, 13, 18, 14, 20, 16, 21, 17, 22, 18, 23, 19,
					   6, 6, 7, 7, 8, 8, 9, 9, 11, 11, 12, 12,13, 13, 14, 14, 16, 16, 17, 17, 18, 18, 19, 19, 21, 21, 22, 22, 23, 23, 24, 24}; //total 32*3=96

	std::vector<Real> coord = {0.3,0.3,0.3,0.6,0.6,0.3};
	std::vector<Real> coord2;
	Id id; //, id2;
	std::vector<Real> region1 = {0.5,0.5,0.5,0.5};
	std::vector<Real> region2 = {1.,1.,1.,1.};
	std::vector<Real> region3 = {0.3,0.3, 0.3,0.3};
	std::vector<Real> region4 = {3,3, 3,3};
	std::set<int> found;
	int loc;


	//default tree_header
	// TreeHeader<Element<3,2,2>> THa;
	// TreeHeader<Element<3,2,3>> THb;
	// TreeHeader<Element<3,3,3>> THc;
	// TreeHeader<Box<2>> THd;
	// TreeHeader<Box<3>> THe;


	//default constructor
	// ADTree<Element<3,2,2>> ADTa; //1st
	// ADTree<Element<3,2,3>> ADTb; //2nd
	// ADTree<Element<3,3,3>> ADTc; //3rd
	// ADTree<Box<2>> ADTd; //4th
	// ADTree<Box<3>> ADTe; //5th

	//tree_header constructor
	// ADTree<Element<3,2,2>> ADTf(THa); //6th
	// ADTree<Element<3,2,3>> ADTg(THb); //7th
	// ADTree<Element<3,3,3>> ADTh(THc); //8th
	// ADTree<Box<2>> ADTi(THd); //9th
	// ADTree<Box<3>> ADTj(THe); //10th


	//*************only possible for triangle to have this kind of constructor
	//constructor from 2D/2.5D/3D mesh
	ADTree<Element<3,2,2>> ADTk(points, triangle, num_nodes, num_triangle); //11th
	// ADTree<Element<3,2,3>> ADTl(points, triangle, num_nodes, num_triangle); //12th
	// ADTree<Element<3,3,3>> ADTm(points, triangle, num_nodes, num_triangle); //13th
	///*******Can there be ADTree<Box<2>> or ADTree<Box<3>> as well??? I assume that
	//*********Shape will always be Element so I am using that logic at ADTree constructor




	std::cout << std::endl << std::endl;
	std::cout << "Print the 11th ADTree (construct from 2D mesh, Shape = Triangle, NDIMP = 2):  " <<std::endl;
	std::cout << ADTk <<std::endl;
	std::cout << "Print TreeHeader:  " <<std::endl;
	std::cout << ADTk.gettreeheader() <<std::endl;

	std::cout << std::endl;
	std::cout << "Print number of nodes (location) occupied:  " <<std::endl;
	std::cout << "Tree head + " << (ADTk.gettreeheader()).getnele() << " tree nodes" <<std::endl;
	std::cout << "Treenodes visualize as a vector:  " << std::endl;
	for (int i = 0; i < (ADTk.gettreeheader()).getnele()+1; i++) {
		(ADTk.gettreenode(i)).print(std::cout);
		//std::cout << "Shape id of father: " << ADTk.pointId(ADTk.gettreenode(i).getfather()) << std::endl;
		std::cout << "Shape id of left child: " << ADTk.pointId(ADTk.gettreenode(i).getchild(0)) << std::endl;
		std::cout << "Shape id of right child: " << ADTk.pointId(ADTk.gettreenode(i).getchild(1)) << std::endl;
	}

	std::cout << "visualize information about a node, for example node number 5: " << std::endl;
	ADTk.gettri(5, coord2, id);
	std::cout << "coordinate:  " ;
	for (size_t i = 0; i < coord2.size(); i++)
		std::cout << " - " << coord2[i] << " - ";
	std::cout << std::endl;
	std::cout << "Id :  " << id << std::endl;
	std::cout << "treenode:  " << std::endl;
	(ADTk.gettreenode(5)).print(std::cout);
	std::cout << "get box coord: " << std::endl;
	std::cout << "( " << ADTk.pointcoord(5+1, 0) << " , " << ADTk.pointcoord(5, 1) << " )" << std::endl;
	std::cout << "( " << ADTk.pointcoord(5+1, 2) << " , " << ADTk.pointcoord(5, 3) << " )" <<std::endl;
	std::cout << "get id: " << std::endl;
	std::cout <<  ADTk.pointId(5) <<std::endl <<std::endl;

	std::cout << "search point (0.5,0.5): " << ADTk.search(region1, found) << std::endl;
	for (std::set<int>::iterator i = found.begin(); i != found.end(); i++) {
  	 std::cout << " --- Location at: " << *i << std::endl;
	(ADTk.gettreenode(*i)).print(std::cout);
	}
	std::cout << std::endl << std::endl;

	std::cout << "search point (1,1): " << ADTk.search(region2, found) << std::endl;
	for (std::set<int>::iterator i = found.begin(); i != found.end(); i++) {
  	 std::cout << " --- Location at: " << *i << std::endl;
	(ADTk.gettreenode(*i)).print(std::cout);
	}
	std::cout << std::endl << std::endl;

	std::cout << "search point (0.3,0.3): " << ADTk.search(region3, found) << std::endl;
	for (std::set<int>::iterator i = found.begin(); i != found.end(); i++) {
  	 std::cout << " --- Location at: " << *i << std::endl;
	(ADTk.gettreenode(*i)).print(std::cout);
	}
	std::cout << std::endl << std::endl;

	std::cout << "search point (3,3): " << ADTk.search(region4, found) << std::endl;
	for (std::set<int>::iterator i = found.begin(); i != found.end(); i++) {
  	 std::cout << " --- Location at: " << *i << std::endl;
	(ADTk.gettreenode(*i)).print(std::cout);
	}
	std::cout << std::endl << std::endl;

	*n = 0;
	}
	catch(...)
	{
	std::cout <<"Exception caugth " << std::endl;
	*n = 1;
	}
}
}


extern "C"{
void Test_MeshHandler_c(int *n){

	try{
	//costruisco una mesh semplice a mano salvata come quella che serve al costruttore
	Real points[] = {0., 0.5, 1.,  0., 0.5,  1.,  0., 0.5, 1.,
					 0., 0.,  0., 0.5, 0.5,  0.5, 1., 1.,  1.};
	UInt triangle[] = {0,0,1,1,3,3,4,4,
					   4,1,5,2,7,4,8,5,
					   3,4,4,5,6,7,7,8};
	UInt * edge;
	UInt * neighbors;
	UInt num_nodes = 9;
	UInt num_triangle = 8;
	UInt num_edges = 16;
	std::vector<Real> coord = {0.3,0.3,0.3,0.6,0.6,0.3};
	Point point1(0.5,0.5);
	Point point2(1.,1.);
	Point point3(0.8,0.3);
	Point point4(4.,4.);


/* the simple mesh is:
	 ___ ___
	|  /|  /|
	| / | / |
	|/__|/__|
	|  /|  /|
	| / | / |
	|/__|/__|

	nodes are numbered:
	6 7 8
	3 4 5
	0 1 2

	triangle are numbered:
	 ___ ___
	|4 /|6 /|
	| /5| /7|
	|/__|/__|
	|0 /|2 /|
	| /1| /3|
	|/__|/__|



*/
	MeshHandler<1,2,2> mesh(points, edge, triangle, neighbors, num_nodes, num_edges, num_triangle);
	Element<3,2,2> res;

	std::cout <<std::endl << std::endl;
	std::cout << "Element<3,2,2> mesh: " << std::endl;
	mesh.printPoints(std::cout);
	mesh.printElements(std::cout);
	std::cout << "find point 1: (0.5,0.5)   " << std::endl;
	res = mesh.findLocationTree(point1);
	res.print(std::cout);
	res[0].print(std::cout);
	std::cout << std::endl;
	res[1].print(std::cout);
	std::cout << std::endl;
	res[2].print(std::cout);
	std::cout << std::endl;
	std::cout <<std::endl << std::endl;
/*	std::cout << "find location naive: "<<std::endl;
	res = mesh.findLocationNaive(point1);
	res.print(std::cout);
	res[0].print(std::cout);
	std::cout << std::endl;
	res[1].print(std::cout);
	std::cout << std::endl;
	res[2].print(std::cout);
	std::cout << std::endl;
	std::cout <<std::endl << std::endl;
*/

	std::cout << "find point 1: (1.,1.)   " << std::endl;
	res = mesh.findLocationTree(point2);
	res.print(std::cout);
	res[0].print(std::cout);
	std::cout << std::endl;
	res[1].print(std::cout);
	std::cout << std::endl;
	res[2].print(std::cout);
	std::cout << std::endl;
	std::cout <<std::endl << std::endl;
/*	std::cout << "find location naive: "<<std::endl;
	res = mesh.findLocationNaive(point2);
	res.print(std::cout);
	res[0].print(std::cout);
	std::cout << std::endl;
	res[1].print(std::cout);
	std::cout << std::endl;
	res[2].print(std::cout);
	std::cout << std::endl;
	std::cout <<std::endl << std::endl;
*/
	std::cout << "find point 1: (0.8,0.3)   " << std::endl;
	res = mesh.findLocationTree(point3);
	res.print(std::cout);
	res[0].print(std::cout);
	std::cout << std::endl;
	res[1].print(std::cout);
	std::cout << std::endl;
	res[2].print(std::cout);
	std::cout << std::endl;
	std::cout <<std::endl << std::endl;
/*	std::cout << "find location naive: "<<std::endl;
	res = mesh.findLocationNaive(point3);
	res.print(std::cout);
	res[0].print(std::cout);
	std::cout << std::endl;
	res[1].print(std::cout);
	std::cout << std::endl;
	res[2].print(std::cout);
	std::cout << std::endl;
	std::cout <<std::endl << std::endl;
*/
	std::cout << "find point 1: (4.,4.)   " << std::endl;
	res = mesh.findLocationTree(point4);
	res.print(std::cout);
	res[0].print(std::cout);
	std::cout << std::endl;
	res[1].print(std::cout);
	std::cout << std::endl;
	res[2].print(std::cout);
	std::cout << std::endl;
	std::cout <<std::endl << std::endl;
/*	std::cout << "find location naive: "<<std::endl;
	res = mesh.findLocationNaive(point4);
	res.print(std::cout);
	res[0].print(std::cout);
	std::cout << std::endl;
	res[1].print(std::cout);
	std::cout << std::endl;
	res[2].print(std::cout);
	std::cout << std::endl;
	std::cout <<std::endl << std::endl;
*/
//	MeshHandler<2> mesh2(points, edge, triangle, neighbors, num_nodes, num_edges, num_triangle, 1);
//	Triangle<6> res2;

	delete[] edge;
	delete[] neighbors;
	*n = 0;
	}
	catch(...)
	{
	std::cout <<"Exception caugth " << std::endl;
	*n = 1;
	}
}
}

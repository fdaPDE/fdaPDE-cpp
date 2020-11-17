#ifndef __INTEGRATION_H__
#define __INTEGRATION_H__

#include <type_traits>

#include "../../FdaPDE.h"
#include "../../Mesh/Include/Point.h"

struct IntegratorTriangleP1{
	//Number of nodes
	static constexpr UInt NNODES = 1;
	static constexpr std::array<Real,NNODES> WEIGHTS{{1}};
	//Point locations (in barycentric coordinates)
	static constexpr std::array<Point<2>,NNODES> NODES{
		Point<2>({1./3,1./3})
	};
};

struct IntegratorTriangleP2{
	//Number of nodes
	static constexpr UInt NNODES = 3;
	static constexpr std::array<Real,NNODES> WEIGHTS{{1./3, 1./3, 1./3}};
	//Point locations (in barycentric coordinates)
	static constexpr std::array<Point<2>,NNODES> NODES{
		Point<2>({1./6,1./6}),
		Point<2>({2./3,1./6}),
		Point<2>({1./6,2./3})
	};
};

struct IntegratorTriangleP4{
	//Number of nodes
	static constexpr UInt NNODES = 6;
	static constexpr std::array<Real,NNODES> WEIGHTS{
		{0.223381589678011,
		 0.223381589678011,
		 0.223381589678011,
		 0.109951743655322,
		 0.109951743655322,
		 0.109951743655322}
	 };
	//Point locations (in barycentric coordinates)
	static constexpr std::array<Point<2>,NNODES> NODES{
		Point<2>({0.445948490915965,0.445948490915965}),
		Point<2>({0.445948490915965,0.108103018168070}),
		Point<2>({0.108103018168070,0.445948490915965}),
		Point<2>({0.091576213509771,0.091576213509771}),
		Point<2>({0.091576213509771,0.816847572980459}),
		Point<2>({0.816847572980459,0.091576213509771})
	};
};

struct IntegratorTriangleP5{
	//Number of nodes
	static constexpr UInt NNODES = 7;
	static constexpr std::array<Real,NNODES> WEIGHTS{
		{0.225000000000000,
		 0.125939180544827,
		 0.125939180544827,
		 0.125939180544827,
		 0.132394152788506,
		 0.132394152788506,
		 0.132394152788506}
	 };
	//Point locations (in barycentric coordinates)
	static constexpr std::array<Point<2>,NNODES> NODES{
		Point<2>({0.333333333333333,0.333333333333333}),
		Point<2>({0.101286507323456,0.101286507323456}),
		Point<2>({0.101286507323456,0.797426985353087}),
		Point<2>({0.797426985353087,0.101286507323456}),
		Point<2>({0.470142064105115,0.470142064105115}),
		Point<2>({0.470142064105115,0.059715871789770}),
		Point<2>({0.059715871789770,0.470142064105115}),
	};
};


struct IntegratorTetrahedronP1{
	//Number of nodes
	static constexpr UInt NNODES = 1;
	static constexpr std::array<Real,NNODES> WEIGHTS{{1}};
	//Point locations (in barycentric coordinates)
	static constexpr std::array<Point<3>,NNODES> NODES{
		Point<3>({1./4, 1./4, 1./4})
	};
};


struct IntegratorTetrahedronP2{
	//Number of nodes
	static constexpr UInt NNODES = 4;
	static constexpr std::array<Real,NNODES> WEIGHTS{{1./4, 1./4, 1./4, 1./4}};
	//Point locations (in barycentric coordinates)
	static constexpr std::array<Point<3>,NNODES> NODES{
		Point<3>({0.585410196624969,0.138196601125011,0.138196601125011}),
		Point<3>({0.138196601125011,0.138196601125011,0.138196601125011}),
		Point<3>({0.138196601125011,0.138196601125011,0.585410196624969}),
		Point<3>({0.138196601125011,0.585410196624969,0.138196601125011})
	};
};

struct IntegratorTetrahedronP4{
	//Number of nodes
	static constexpr UInt NNODES = 14;
	static constexpr std::array<Real,NNODES> WEIGHTS{
		{0.008501593297035,
		 0.008501593297035,
		 0.008501593297035,
		 0.008501593297035,
		 0.050615550377161,
		 0.050615550377161,
		 0.050615550377161,
		 0.050615550377161,
		 0.127255237550536,
		 0.127255237550536,
		 0.127255237550536,
		 0.127255237550536,
		 0.127255237550536,
		 0.127255237550536}
	};
	//Point locations (in barycentric coordinates)
	static constexpr std::array<Point<3>,NNODES> NODES{
		Point<3>({1./3, 1./3, 1./3}),
		Point<3>({	 0, 1./3, 1./3}),
		Point<3>({1./3,	   0, 1./3}),
		Point<3>({1./3, 1./3, 	 0}),
		Point<3>({0.076284841942834,0.076284841942834,0.076284841942834}),
		Point<3>({0.771145474171499,0.076284841942834,0.076284841942834}),
		Point<3>({0.076284841942834,0.771145474171499,0.076284841942834}),
		Point<3>({0.076284841942834,0.076284841942834,0.771145474171499}),
		Point<3>({0.405500312564648,0.094499687435353,0.094499687435353}),
		Point<3>({0.094499687435353,0.405500312564648,0.094499687435353}),
		Point<3>({0.094499687435353,0.094499687435353,0.405500312564648}),
		Point<3>({0.405500312564648,0.405500312564648,0.094499687435353}),
		Point<3>({0.405500312564648,0.094499687435353,0.405500312564648}),
		Point<3>({0.094499687435353,0.405500312564648,0.405500312564648})
	};
};

struct ElementIntegratorHelper{
	template<UInt NNODES, UInt mydim>
	using Integrator = typename std::conditional<mydim==2,
												typename std::conditional<NNODES==3, IntegratorTriangleP1, IntegratorTriangleP2>::type,
												typename std::conditional<NNODES==4, IntegratorTetrahedronP1, IntegratorTetrahedronP2>::type>::type;
};

struct SpaceIntegratorHelper{
	template<UInt ORDER, UInt mydim>
	using Integrator = typename std::conditional<mydim==2,
												typename std::conditional<ORDER==1, IntegratorTriangleP2, IntegratorTriangleP4>::type,
												typename std::conditional<ORDER==1, IntegratorTetrahedronP2, IntegratorTetrahedronP4>::type>::type;
};

struct IntegratorGaussP3{
	static constexpr UInt ORDER = 1;
	//Number of nodes
	static constexpr UInt NNODES = 2;

	static constexpr std::array<Real, NNODES> WEIGHTS{{1., 1.}};

	static constexpr std::array<Real, NNODES> NODES{{-0.577350269189626, 0.577350269189626}};
	};

struct IntegratorGaussP5{
	static constexpr UInt ORDER = 1;
	//Number of nodes
	static constexpr UInt NNODES = 3;
	//Point locations
	static constexpr std::array<Real, NNODES> WEIGHTS{{5./9, 8./9, 5./9}};
	
	static constexpr std::array<Real, NNODES> NODES{{-0.774596669241483, 0, 0.774596669241483}};
	};

// Gauss Legendre quadrature rules for triangles (n=3)
struct IntegratorGaussTriangle3{
	//Number of nodes
	static constexpr UInt NNODES = 9;
	// these sum up to 1/2
	static constexpr std::array<Real, NNODES> WEIGHTS{
		{0.068464377,
		 0.109543004,
		 0.068464377,
		 0.061728395,
		 0.098765432,
		 0.061728395,
		 0.008696116,
		 0.013913785,
		 0.008696116}
	 };
	static constexpr std::array<Point<2>, NNODES> NODES{
		Point<2>({0.112701665, 0.100000000}),
		Point<2>({0.112701665, 0.443649167}),
		Point<2>({0.112701665, 0.787298334}),
		Point<2>({0.500000000, 0.056350832}),
		Point<2>({0.500000000, 0.250000000}),
		Point<2>({0.500000000, 0.443649167}),
		Point<2>({0.887298334, 0.012701665}),
		Point<2>({0.887298334, 0.056350832}),
		Point<2>({0.887298334, 0.100000000})
	};
};

// Gauss Legendre quadrature rules for tetra (n=3)
struct IntegratorGaussTetra3{
	//Number of nodes
	static constexpr UInt NNODES = 27;
	//Point locations
	static constexpr std::array<Real, NNODES> WEIGHTS{
		{0.014972747367084,
		 0.014972747367084,
		 0.023956395787334,
		 0.001901788268649,
		 0.001901788268649,
		 0.003042861229838,
		 0.013499628508586,
		 0.013499628508586,
		 0.021599405613738,
		 0.000241558782106,
		 0.000241558782106,
		 0.000386494051369,
		 0.000030681988197,
		 0.000030681988197,
		 0.000049091181116,
		 0.000217792616242,
		 0.000217792616242,
		 0.000348468185988,
		 0.007607153074595,
		 0.007607153074595,
		 0.012171444919352,
		 0.000966235128423,
		 0.000966235128423,
		 0.001545976205477,
		 0.006858710562414,
		 0.006858710562414,
		 0.010973936899863}
	 };
	static constexpr std::array<Point<3>, NNODES> NODES{
		Point<3>({0.112701665379259, 0.100000000000000, 0.088729833462074}),
		Point<3>({0.112701665379259, 0.100000000000000, 0.698568501158667}),
		Point<3>({0.112701665379259, 0.100000000000000, 0.393649167310371}),
		Point<3>({0.112701665379259, 0.787298334620741, 0.011270166537926}),
		Point<3>({0.112701665379259, 0.787298334620741, 0.088729833462074}),
		Point<3>({0.112701665379259, 0.787298334620741, 0.050000000000000}),
		Point<3>({0.112701665379259, 0.443649167310371, 0.050000000000000}),
		Point<3>({0.112701665379259, 0.443649167310371, 0.393649167310371}),
		Point<3>({0.112701665379259, 0.443649167310371, 0.221824583655185}),
		Point<3>({0.887298334620741, 0.012701665379258, 0.011270166537926}),
		Point<3>({0.887298334620741, 0.012701665379258, 0.088729833462074}),
		Point<3>({0.887298334620741, 0.012701665379258, 0.050000000000000}),
		Point<3>({0.887298334620741, 0.100000000000000, 0.001431498841332}),
		Point<3>({0.887298334620741, 0.100000000000000, 0.011270166537926}),
		Point<3>({0.887298334620741, 0.100000000000000, 0.006350832689629}),
		Point<3>({0.887298334620741, 0.056350832689629, 0.006350832689629}),
		Point<3>({0.887298334620741, 0.056350832689629, 0.050000000000000}),
		Point<3>({0.887298334620741, 0.056350832689629, 0.028175416344815}),
		Point<3>({0.500000000000000, 0.056350832689629, 0.050000000000000}),
		Point<3>({0.500000000000000, 0.056350832689629, 0.393649167310371}),
		Point<3>({0.500000000000000, 0.056350832689629, 0.221824583655185}),
		Point<3>({0.500000000000000, 0.443649167310371, 0.006350832689629}),
		Point<3>({0.500000000000000, 0.443649167310371, 0.050000000000000}),
		Point<3>({0.500000000000000, 0.443649167310371, 0.028175416344815}),
		Point<3>({0.500000000000000, 0.250000000000000, 0.028175416344815}),
		Point<3>({0.500000000000000, 0.250000000000000, 0.221824583655185}),
		Point<3>({0.500000000000000, 0.250000000000000, 0.125000000000000})
	};
};


#endif

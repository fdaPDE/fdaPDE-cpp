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


#endif

#ifndef __INTEGRATION_HPP__
#define __INTEGRATION_HPP__

#include <type_traits>

#include "mesh_objects.h"

struct IntegratorTriangleP2{
	static constexpr UInt ORDER = 1;
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
	static constexpr UInt ORDER = 2;
	//Number of nodes
	static constexpr UInt NNODES = 6;
	static constexpr std::array<Real,NNODES> WEIGHTS{{0.223381589678011,0.223381589678011,0.223381589678011,0.109951743655322,0.109951743655322,0.109951743655322}};
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


struct IntegratorTetrahedronP2{
	static constexpr UInt ORDER = 1;
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
	static constexpr UInt ORDER = 2;
	//Number of nodes
	static constexpr UInt NNODES = 11;
	static constexpr std::array<Real,NNODES> WEIGHTS{{-0.078933333333333, 0.045733333333333, 0.045733333333333, 0.045733333333333, 0.045733333333333, 0.149333333333333, 0.149333333333333, 0.149333333333333, 0.149333333333333, 0.149333333333333, 0.149333333333333}};
	//Point locations (in barycentric coordinates)
	static constexpr std::array<Point<3>,NNODES> NODES{
		Point<3>({0.250000000000000,0.250000000000000,0.250000000000000}),
		Point<3>({0.071428571428571,0.071428571428571,0.071428571428571}),
		Point<3>({0.785714285714286,0.071428571428571,0.071428571428571}),
		Point<3>({0.071428571428571,0.785714285714286,0.071428571428571}),
		Point<3>({0.071428571428571,0.071428571428571,0.785714285714286}),
		Point<3>({0.399403576166799,0.100596423833201,0.100596423833201}),
		Point<3>({0.100596423833201,0.399403576166799,0.100596423833201}),
		Point<3>({0.100596423833201,0.100596423833201,0.399403576166799}),
		Point<3>({0.399403576166799,0.399403576166799,0.100596423833201}),
		Point<3>({0.399403576166799,0.100596423833201,0.399403576166799}),
		Point<3>({0.100596423833201,0.399403576166799,0.399403576166799})
	};
};

struct IntegratorHelper{
	template<UInt ORDER, UInt mydim>
	using Integrator = typename std::conditional<mydim==2,
												typename std::conditional<ORDER==1, IntegratorTriangleP2, IntegratorTriangleP4>::type,
												typename std::conditional<ORDER==1, IntegratorTetrahedronP2, IntegratorTetrahedronP4>::type>::type;
};


#endif

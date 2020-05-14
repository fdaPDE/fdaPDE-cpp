#ifndef __INTEGRATION_HPP__
#define __INTEGRATION_HPP__

#include "fdaPDE.h"
#include "mesh_objects.h"



class IntegratorTriangleP2{
	public:
	static const UInt ORDER = 1;
	//Number of nodes
	static const UInt NNODES = 3;
	//Point locations
	static const std::vector<Point> NODES;
	static const std::vector<Real> WEIGHTS;
	};

class IntegratorTriangleP4{
	public:
	static const UInt ORDER = 2;
	//Number of nodes
	static const UInt NNODES = 6;
	//Point locations
	static const std::vector<Point> NODES;
	static const std::vector<Real> WEIGHTS;
	};

class IntegratorTetrahedronP2{
	public:
	static const UInt ORDER = 1;
	//Number of nodes
	static const UInt NNODES = 4;
	//Point locations
	static const std::vector<Point> NODES;
	static const std::vector<Real> WEIGHTS;
	};

class IntegratorTetrahedronP1{
	public:
	static const UInt ORDER = 1;
	//Number of nodes
	static const UInt NNODES = 1;
	//Point locations
	static const std::vector<Point> NODES;
	static const std::vector<Real> WEIGHTS;
	};

	class IntegratorGaussP3{
	  public:
		static const UInt ORDER = 1;
		//Number of nodes
		static const UInt NNODES = 2;
		//Point locations
		static const std::vector<Real> NODES;
		static const std::vector<Real> WEIGHTS;
		};

	class IntegratorGaussP5{
	  public:
	  static const UInt ORDER = 1;
		//Number of nodes
		static const UInt NNODES = 3;
		//Point locations
		static const std::vector<Real> NODES;
		static const std::vector<Real> WEIGHTS;
		};

// Gauss Legendre quadrature rules for triangles (n=3)
class IntegratorGaussTriangle3{
	private:
		static const std::vector<Real> w;

	public:
		//Number of nodes
		static const UInt NNODES = 9;
		//Point locations
		static const std::vector<Point> NODES;
		static const VectorXr WEIGHTS;
};

// Gauss Legendre quadrature rules for tetra (n=3)
class IntegratorGaussTetra3
{
	private:
		static const std::vector<Real> w;

	public:
		//Number of nodes
		static const UInt NNODES = 27;
		//Point locations
		static const std::vector<Point> NODES;
		static const VectorXr WEIGHTS;
};



//#include "integration_imp.hpp"
#endif

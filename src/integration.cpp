#include "integration.h"

// Out of class definition is needed because these arrays are odr-used
// Not needed anymore from C++17
// Note: NNODES variables are not odr-used so no out of class definition is needed
constexpr std::array<Point<2>,IntegratorTriangleP2::NNODES> IntegratorTriangleP2::NODES;

constexpr std::array<Point<2>,IntegratorTriangleP4::NNODES> IntegratorTriangleP4::NODES;

constexpr std::array<Point<3>,IntegratorTetrahedronP2::NNODES> IntegratorTetrahedronP2::NODES;

constexpr std::array<Point<3>,IntegratorTetrahedronP4::NNODES> IntegratorTetrahedronP4::NODES;

constexpr std::array<Real,IntegratorTriangleP2::NNODES> IntegratorTriangleP2::WEIGHTS;

constexpr std::array<Real,IntegratorTriangleP4::NNODES> IntegratorTriangleP4::WEIGHTS;

constexpr std::array<Real,IntegratorTetrahedronP2::NNODES> IntegratorTetrahedronP2::WEIGHTS;

constexpr std::array<Real,IntegratorTetrahedronP4::NNODES> IntegratorTetrahedronP4::WEIGHTS;

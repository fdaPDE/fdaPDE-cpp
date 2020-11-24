#include "../Include/Integration.h"

// Out of class definition is needed because these arrays are odr-used
// Not needed anymore from C++17
// Note: NNODES variables are not odr-used so no out of class definition is needed
constexpr std::array<Point<2>,IntegratorTriangleP1::NNODES> IntegratorTriangleP1::NODES;

constexpr std::array<Point<2>,IntegratorTriangleP2::NNODES> IntegratorTriangleP2::NODES;

constexpr std::array<Point<2>,IntegratorTriangleP4::NNODES> IntegratorTriangleP4::NODES;

constexpr std::array<Point<2>,IntegratorTriangleP5::NNODES> IntegratorTriangleP5::NODES;

constexpr std::array<Point<3>,IntegratorTetrahedronP1::NNODES> IntegratorTetrahedronP1::NODES;

constexpr std::array<Point<3>,IntegratorTetrahedronP2::NNODES> IntegratorTetrahedronP2::NODES;

constexpr std::array<Point<3>,IntegratorTetrahedronP4::NNODES> IntegratorTetrahedronP4::NODES;

constexpr std::array<Real,IntegratorGaussP3::NNODES> IntegratorGaussP3::NODES;

constexpr std::array<Real,IntegratorGaussP5::NNODES> IntegratorGaussP5::NODES;



constexpr std::array<Real,IntegratorTriangleP1::NNODES> IntegratorTriangleP1::WEIGHTS;

constexpr std::array<Real,IntegratorTriangleP2::NNODES> IntegratorTriangleP2::WEIGHTS;

constexpr std::array<Real,IntegratorTriangleP4::NNODES> IntegratorTriangleP4::WEIGHTS;

constexpr std::array<Real,IntegratorTriangleP5::NNODES> IntegratorTriangleP5::WEIGHTS;

constexpr std::array<Real,IntegratorTetrahedronP1::NNODES> IntegratorTetrahedronP1::WEIGHTS;

constexpr std::array<Real,IntegratorTetrahedronP2::NNODES> IntegratorTetrahedronP2::WEIGHTS;

constexpr std::array<Real,IntegratorTetrahedronP4::NNODES> IntegratorTetrahedronP4::WEIGHTS;

constexpr std::array<Real,IntegratorGaussP3::NNODES> IntegratorGaussP3::WEIGHTS;

constexpr std::array<Real,IntegratorGaussP5::NNODES> IntegratorGaussP5::WEIGHTS;



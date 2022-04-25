#ifndef __FINITE_ELEMENT__
#define __FINITE_ELEMENT__

#include "../utils/CompileTime.h"
#include "../utils/Symbols.h"
#include "LagrangianBasis.h"

#include <array>
#include <cstddef>

// hard code the coordinates of nodes over an N-dimensional unit simplex for Lagrange interpolation.
// N dimension of the simplex, R order of the basis to be defined on it
template <unsigned int N, unsigned int R> struct ReferenceNodes;
template <unsigned int N, unsigned int R> using point_list = std::array<std::array<double, N>, R>;

template<> // 1D first order basis
struct ReferenceNodes<1,1>{
  static constexpr point_list<1,2> nodes  = {
    {{0}, {1}}
  };};
template<> // 1D second order basis
struct ReferenceNodes<1,2>{
  static constexpr point_list<1,3> nodes  = {
    {{0}, {0.5}, {1}}
  };};

template<> // 2D first order basis
struct ReferenceNodes<2,1>{
  static constexpr point_list<2,3> nodes  = {
    {{0, 0}, {1, 0}, {0, 1}}
  };};
template<> // 2D second order basis
struct ReferenceNodes<2,2>{
  static constexpr point_list<2,6> nodes  = {
    {{0, 0}, {1, 0}, {0, 1}, {0.5, 0.5}, {0, 0.5}, {0.5, 0}}
  };};

template<> // 3D first order basis
struct ReferenceNodes<3,1>{
  static constexpr point_list<3,4> nodes  = {
    {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}
  };};
template<> // 3D second order basis
struct ReferenceNodes<3,2>{
  static constexpr point_list<3,10> nodes = {
    {{0, 0, 0},   {1, 0, 0},   {0, 1, 0},     {0, 0, 1},     {0.5, 0.5, 0},
     {0, 0.5, 0}, {0.5, 0, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}, {0, 0, 0.5}}
  };};

// The refererence element defined on an N-dimensional simplex equipped with a
// space of polynomials of order R spanned by a Lagrangian basis
template <unsigned int N, unsigned int R>
class ReferenceElement{
private:
  // the set of nodes where the element basis is defined
  static constexpr std::array<std::array<double, N>, ct_binomial_coefficient(N+R, R)> nodes_ = ReferenceNodes<N, R>::nodes;
  LagrangianBasis<N, R> basis_;
public:
  // constructor
  ReferenceElement() : basis_(LagrangianBasis<N,R>(nodes_)) { };

  // getters
  LagrangianBasis<N, R> getBasis() const { return basis_; };
};


#endif // __FINITE_ELEMENT__

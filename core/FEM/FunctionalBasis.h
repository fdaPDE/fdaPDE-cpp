#ifndef __FUNCTIONAL_BASIS_H__
#define __FUNCTIONAL_BASIS_H__

#include "../utils/CompileTime.h"
#include "../utils/Symbols.h"
#include "../MESH/Element.h"
using fdaPDE::core::MESH::Element;

#include "LagrangianBasis.h"
#include <array>

// N space dimension, M local dimension of the element, R order of basis
template <unsigned int M, unsigned int N, unsigned int ORDER>
class FunctionalBasis {
private:
  // the set of nodes where the element basis is defined
  std::array<std::array<double, N>, ct_binomial_coefficient(N+ORDER, ORDER)> nodes_;
  LagrangianBasis<N, ORDER> basis_;
public:
  // construct from a set of nodes
  FunctionalBasis(const std::array<std::array<double, N>, ct_binomial_coefficient(N+ORDER, ORDER)>& nodes) : nodes_(nodes) {
    basis_ = LagrangianBasis<N, ORDER>(nodes_);
  }
  // construct from a mesh element
  FunctionalBasis(const Element<M, N>& e){
        // build functional basis over the current element e
    std::array<std::array<double, N>, ct_binomial_coefficient(N+ORDER, ORDER)> nodes{};
    for(size_t j = 0; j < e.getFESupport().size(); ++j) 
      nodes[j] = {e.getFESupport()[j].second[0], e.getFESupport()[j].second[1]};

    nodes_ = nodes;
    basis_ = LagrangianBasis<N, ORDER>(nodes_);
  }

  LagrangianBasis<N, ORDER> getBasis() const { return basis_; };      // getter
  // subscript operator to directly access basis elements
  MultivariatePolynomial<N, ORDER> operator[](size_t i) const { return basis_.getBasisElement(i); }
};

// hard code the coordinates of nodes over an N-dimensional unit simplex for Lagrange interpolation
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

// a functional basis built over the referece N-dimensional unit simplex
template <unsigned int M, unsigned int N, unsigned int ORDER>
struct ReferenceBasis : public FunctionalBasis<M, N, ORDER>{
  ReferenceBasis() : FunctionalBasis<M, N, ORDER>(ReferenceNodes<N, ORDER>::nodes){ };  
};

#endif // __FUNCTIONAL_BASIS_H__

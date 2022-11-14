#include <cstddef>
#include <gtest/gtest.h> // testing framework
#include <limits>
#include <string>
#include <type_traits>

#include "../../fdaPDE/core/utils/Symbols.h"
#include "../../fdaPDE/core/utils/CompileTime.h"
#include "../../fdaPDE/core/FEM/basis/LagrangianBasis.h"
#include "../../fdaPDE/core/utils/fields/VectorField.h"
using fdaPDE::core::FEM::LagrangianBasis;
using fdaPDE::core::FEM::ReferenceNodes;
using fdaPDE::core::FEM::point_list;
#include "../utils/Constants.h"
using fdaPDE::testing::DOUBLE_TOLERANCE;
#include "../utils/MeshLoader.h"
using fdaPDE::testing::MeshLoader;
using fdaPDE::testing::MESH_TYPE_LIST;

// a type representing a compile time pair of integer values
template <int i, int j> struct int_pair {
  static constexpr int first  = std::integral_constant<int, i>::value;
  static constexpr int second = std::integral_constant<int, j>::value;
};

template <typename E>
class LagrangianBasisTest : public ::testing::Test {
public:
  static constexpr int N = E::first;
  static constexpr int R = E::second;
  static constexpr std::size_t n_basis = ct_binomial_coefficient(N+R,R);

  // constructor
  LagrangianBasisTest() = default;
  
  // transforms the i-th point in a point_list to an SVector
  template <long unsigned int M, long unsigned int K>
  SVector<M> toSVector(const point_list<M, K> &pList, std::size_t i) const {
    SVector<M> result{};
    for(std::size_t j = 0; j < M; ++j){
      result[j] = pList[i][j];
    }
    return result;
  };
};

// pair format: <x,y> : x dimension of space, y order of mesh
using pairs = ::testing::Types<int_pair<1,1>, int_pair<2,1>, int_pair<3,1>   // order 1 elements (linear finite elements)
			      ,int_pair<1,2>, int_pair<2,2>, int_pair<3,2>>; // order 2 elements (quadratic finite elements)
TYPED_TEST_SUITE(LagrangianBasisTest, pairs);

// tests a Lagrangian basis can be successfully built over the reference unit simplex
TYPED_TEST(LagrangianBasisTest, DefineOverReferenceElement) {
  // create lagrangian basis over unit dimensional simplex
  LagrangianBasis<TestFixture::N, TestFixture::N, TestFixture::R> basis{};

  // expect correct number of basis functions
  EXPECT_EQ(basis.size(), TestFixture::n_basis);
  // check lagrangian property (each basis function is 1 in one and only one node and 0 in any other)
  for(const MultivariatePolynomial<TestFixture::N, TestFixture::R>& b : basis){
    std::size_t num_ones = 0, num_zeros = 0;
    for(std::size_t i = 0; i < TestFixture::n_basis; ++i){ // there are as many nodes as basis functions
      SVector<TestFixture::N> p = this->toSVector(ReferenceNodes<TestFixture::N, TestFixture::R>::nodes, i);
      if(std::abs(b(p) - 1) < DOUBLE_TOLERANCE){
	num_ones++;
      }else if(b(p) < DOUBLE_TOLERANCE){
	num_zeros++;
      }
    }
    EXPECT_EQ(num_ones + num_zeros, TestFixture::n_basis);
    EXPECT_EQ(num_ones, 1);
    EXPECT_EQ(num_zeros, TestFixture::n_basis - 1);
  }
}

// test linear elements behave correctly on reference element
TEST(LagrangianBasisTest, LinearElement) {
  // create finite linear elements over unit reference simplex
  LagrangianBasis<2,2,1> basis{};
  SVector<2> p(0,0); // define evaluation point

  // basis functions are defined in counterclockwise order starting from node (0,0)
  // for linear elements we get
  // (0,0) -> \Nabla phi_{0} = [-1 -1]
  // (1,0) -> \Nabla phi_{0} = [ 1  0]
  // (0,1) -> \Nabla phi_{0} = [ 0  1]
  std::vector<SVector<2>> gradients({SVector<2>(-1.0,-1.0), SVector<2>(1.0, 0.0), SVector<2>(0.0, 1.0)});

  // check gradient of each basis function equals the expected one
  for(std::size_t i = 0; i < basis.size(); ++i){
    // extract gradient of basis function
    VectorField<2> grad = basis[i].derive();
    for(std::size_t j = 0; j < 2; ++j)
      EXPECT_NEAR(grad(p)[j], gradients[i][j], DOUBLE_TOLERANCE);
  }
}

// fixture used to check if a lagrangian basis can be defined over a generic mesh (also manifold)
template <typename E>
struct LagrangianBasisMeshTest : public ::testing::Test {
  MeshLoader<E> meshLoader{}; // use default mesh
  static constexpr unsigned int M = MeshLoader<E>::M;
  static constexpr unsigned int N = MeshLoader<E>::N;

  // quantites specific to this test suite
  static constexpr bool is_manifold = MeshLoader<E>::manifold;
  static constexpr std::size_t n_basis = ct_binomial_coefficient(E::local_dimension + E::order, E::order);
};
TYPED_TEST_SUITE(LagrangianBasisMeshTest, MESH_TYPE_LIST);

// tests a Lagrangian basis can be built on any kind of mesh element
TYPED_TEST(LagrangianBasisMeshTest, DefineOverMeshElement) {
  // take element at random
  auto e = this->meshLoader.generateRandomElement();
  // build basis over the element (assume linear elements)
  LagrangianBasis<TestFixture::M, TestFixture::N, 1> basis(*e);
  // expect correct number of basis functions
  EXPECT_EQ(basis.size(), TestFixture::n_basis);
  // check lagrangian property (each basis function is 1 in one and only one node and 0 in any other)
  for(const MultivariatePolynomial<TestFixture::M, 1>& b : basis){
    std::size_t num_ones = 0, num_zeros = 0;
    for(std::size_t i = 0; i < TestFixture::n_basis; ++i){ // there are as many nodes as basis functions
      SVector<TestFixture::M> p;
      if constexpr(!TestFixture::is_manifold)
	p = e->coords()[i];
      else // in case of manifolds we must consider the projection of basis nodes onto the space spanned by the element
	p = e->spannedSpace().projectOnto(e->coords()[i]);
      
      if(std::abs(b(p) - 1) < DOUBLE_TOLERANCE){
	num_ones++;
      }else if(b(p) < DOUBLE_TOLERANCE){
	num_zeros++;
      }
    }
    EXPECT_EQ(num_ones + num_zeros, TestFixture::n_basis);
    EXPECT_EQ(num_ones, 1);
    EXPECT_EQ(num_zeros, TestFixture::n_basis - 1);
  }
  // if the basis function represents a linear function built over an element, its evaluation at the element's midpoint must be
  // equal to 1/n_basis. This check is fundamental to assess if what a LagrangianBasis is building behaves as expected
  EXPECT_NEAR(basis[0](e->spannedSpace().projectOnto(e->midPoint())), 1.0/TestFixture::n_basis, DOUBLE_TOLERANCE);
}

// test linear elements behave correctly on generic mesh elements
TEST(LagrangianBasisMeshTest, LinearElement){
  // load sample mesh
  MeshLoader<Mesh2D<>> CShaped("c_shaped");
  auto e = CShaped.mesh.element(175); // reference element for this test
  
  // define a linear finite element basis over e
  LagrangianBasis<2,2,1> basis(*e);
  SVector<2> p = e->coords()[0]; // define evaluation point
  // define vector of expected gradients, in the order with which basis are looped
  std::vector<SVector<2>> gradients({
      SVector<2>(-5.2557081783567776, -3.5888585000106943),
      SVector<2>( 6.2494499783110298, -1.2028086954015513),
      SVector<2>(-0.9937417999542519,  4.7916671954122458)
    });

  // check gradient of each basis function equals the expected one
  for(std::size_t i = 0; i < basis.size(); ++i){
    // extract gradient of basis function
    VectorField<2> grad = basis[i].derive();
    for(std::size_t j = 0; j < 2; ++j)
      EXPECT_NEAR(grad(p)[j], gradients[i][j], DOUBLE_TOLERANCE);
  }

  // the same can be obtained using the barycentric matrix of e and finite element basis defined over the reference element
  Eigen::Matrix<double, 2, 2> invJ = e->invBarycentricMatrix().transpose();
  LagrangianBasis<2,2,1> refBasis{};

  for(std::size_t i = 0; i < refBasis.size(); ++i){
    VectorField<2, 2> grad = invJ * refBasis[i].derive();
    for(std::size_t j = 0; j < 2; ++j)
      EXPECT_NEAR(grad(p)[j], gradients[i][j], DOUBLE_TOLERANCE);
  }
}

#include <cstddef>
#include <gtest/gtest.h> // testing framework
#include <gtest/internal/gtest-type-util.h>
#include <limits>
#include <string>
#include <type_traits>
#include <random>

#include "../fdaPDE/core/utils/Symbols.h"
#include "../fdaPDE/core/FEM/basis/LagrangianBasis.h"
#include "../fdaPDE/core/utils/CompileTime.h"
using fdaPDE::core::FEM::LagrangianBasis;
using fdaPDE::core::FEM::ReferenceNodes;
using fdaPDE::core::FEM::point_list;
#include "../fdaPDE/core/FEM/basis/MultivariatePolynomial.h"
#include "../fdaPDE/core/MESH/Mesh.h"
using fdaPDE::core::MESH::Mesh2D;
using fdaPDE::core::MESH::Mesh3D;
using fdaPDE::core::MESH::NetworkMesh;
using fdaPDE::core::MESH::SurfaceMesh;
using fdaPDE::core::MESH::is_manifold;

// a type representing a pair of integer values
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
  double tolerance = 20*std::numeric_limits<double>::epsilon();

  LagrangianBasisTest() = default;
  
  // transforms the i-th point in a point_list to an SVector
  template <long unsigned int M, long unsigned int K>
  SVector<M> toSVector(const point_list<M, K> &pList, std::size_t i) const;
};

// transforms a point_list in an SVector
template <typename E>
template <long unsigned int M, long unsigned int K>
SVector<M> LagrangianBasisTest<E>::toSVector(const point_list<M, K> &pList, std::size_t i) const{
  SVector<M> result{};
  for(std::size_t j = 0; j < M; ++j){
    result[j] = pList[i][j];
  }
  return result;
}

using pairs = ::testing::Types<int_pair<1,1>, int_pair<2,1>, int_pair<3,1>   // order 1 elements (linear finite elements)
			      ,int_pair<1,2>, int_pair<2,2>, int_pair<3,2>>; // order 2 elements (quadratic finite elements)
TYPED_TEST_SUITE(LagrangianBasisTest, pairs);

// tests a Lagrangian basis can be successfully built over the reference unit simplex
TYPED_TEST(LagrangianBasisTest, ReferenceElement) {
  // create lagrangian basis over unit dimensional simplex
  LagrangianBasis<TestFixture::N, TestFixture::N, TestFixture::R> basis;

  // expect correct number of basis functions
  EXPECT_EQ(basis.size(), TestFixture::n_basis);
  // check lagrangian property (each basis function is 1 in one and only one node and 0 in any other)
  for(const MultivariatePolynomial<TestFixture::N, TestFixture::R>& b : basis){
    double sum = 0;
    for(std::size_t i = 0; i < TestFixture::n_basis; ++i){ // there are as many nodes as basis functions
      sum += b(this->toSVector(ReferenceNodes<TestFixture::N, TestFixture::R>::nodes, i));
    }
    EXPECT_NEAR(sum, 1.0, this->tolerance);
  }
}

// fixture used to check if a lagrangian basis can be defined over a generic mesh (also manifold)
template <typename E>
class LagrangianBasisMeshTest : public ::testing::Test {
private:
  // RNG stuffs
  std::default_random_engine rng{};

public:
  E m; // mesh
  static constexpr unsigned int M = E::local_dimension;
  static constexpr unsigned int N = E::embedding_dimension;
  static constexpr std::size_t n_basis = ct_binomial_coefficient(M+1,1);
  double tolerance = std::pow(0.1, 13);
  
  // load mesh from .csv files
  LagrangianBasisMeshTest() {
    std::string dim = (E::manifold) ? std::to_string(M) + ".5" : std::to_string(M);
    // compute file names
    std::string point    = "data/m" + dim + "D_points.csv";
    std::string edges    = "data/m" + dim + "D_edges.csv";
    std::string elements = "data/m" + dim + "D_elements.csv";
    std::string neigh    = "data/m" + dim + "D_neigh.csv";
    std::string boundary = "data/m" + dim + "D_boundary.csv";
    // initialize test objects
    m = E(point, edges, elements, neigh, boundary);
    std::random_device rd;
    rng = std::default_random_engine(rd());
  };

  // generate element at random inside mesh m
  std::shared_ptr<Element<E::local_dimension, E::embedding_dimension>> generateRandomElement();
};

template <typename E>
std::shared_ptr<Element<E::local_dimension, E::embedding_dimension>>
LagrangianBasisMeshTest<E>::generateRandomElement() {
  std::uniform_int_distribution<int> randomID(0, this->m.elements()-1);
  int ID = randomID(rng);  
  return m.element(ID);
}

using meshList = ::testing::Types<Mesh2D<>, SurfaceMesh<>, Mesh3D<>, NetworkMesh<>>;
TYPED_TEST_SUITE(LagrangianBasisMeshTest, meshList);

// tests a Lagrangian basis can be built on any kind of mesh element
TYPED_TEST(LagrangianBasisMeshTest, DefineOverElement) {
  // take element at random
  auto e = this->generateRandomElement();
  // build basis over the element
  LagrangianBasis<TestFixture::M, TestFixture::N, 1> basis(*e);
  // expect correct number of basis functions
  EXPECT_EQ(basis.size(), TestFixture::n_basis);
  // check lagrangian property (each basis function is 1 in one and only one node and 0 in any other)
  for(const MultivariatePolynomial<TestFixture::M, 1>& b : basis){
    // a possible way to do so is to sum the value of the basis function at all nodes and check at the end it equals 1
    double sum = 0;
    for(std::size_t i = 0; i < TestFixture::n_basis; ++i){ // there are as many nodes as basis functions
      if constexpr(!is_manifold<TestFixture::M, TestFixture::N>::value){
	sum += b(e->coords()[i]);
      }else{
	// in case of manifolds we must consider the projection of basis nodes onto the space spanned by the element
	sum += b(e->spannedSpace().projectOnto(e->coords()[i]));
      }
    }
    EXPECT_NEAR(sum, 1.0, this->tolerance);
  }
}

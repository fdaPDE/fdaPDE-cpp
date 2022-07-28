#include <gtest/gtest.h> // testing framework
#include <limits>
#include <memory>
#include <random>
#include <utility>

#include "../fdaPDE/core/utils/Symbols.h"
#include "../fdaPDE/core/MESH/Mesh.h"
using fdaPDE::core::MESH::NetworkMesh;
using fdaPDE::core::MESH::Mesh2D;
using fdaPDE::core::MESH::Mesh3D;
using fdaPDE::core::MESH::SurfaceMesh;
#include "../fdaPDE/core/MESH/Element.h"
using fdaPDE::core::MESH::Element;
#include "../fdaPDE/core/MESH/engines/BruteForce/BruteForce.h"
using fdaPDE::core::MESH::BruteForce;
#include "../fdaPDE/core/MESH/engines/BarycentricWalk/BarycentricWalk.h"
using fdaPDE::core::MESH::BarycentricWalk;
#include "../fdaPDE/core/MESH/engines/AlternatingDigitalTree/ADT.h"
using fdaPDE::core::MESH::ADT;
#include "../fdaPDE/core/utils/IO/CSVReader.h"

// correctness of contained information in an Element object has been asserted in MeshTest.cpp test suite, please refer to it if you believe
// that some information stored in an Element object is not correct (is indeed responsibility of Mesh.h to build the single Element instance).
// Here we assume informations packed in an instance of Element class are correct and test for correctness of operations on elements.

// sample meshes used for testing: see m*D_*.csv data series in test/data/ folder to inspect raw informations
//     * 2D:   unit disk with:   441  2D points, 780   elements, 1220  edges. /test/data/m2D_*.csv
//     * 3D:   unit sphere with: 4193 3D points, 22200 elements, 45380 faces. /test/data/m3D_*.csv
//     * 2.5D: manifold with:    340  3D points, 616   elements, 956   edges. /test/data/m2.5D_*.csv
//     * 1.5D: graph with:       204  2D points, 559   elements               /test/data/m1.5D_*.csv

// test fixture. ADT and bruteforce can work on this fixture (barycentric walk is not able to handle manifolds)
template <typename E>
class SearchEngineTest : public ::testing::Test {
public:
  E m; // mesh
  static constexpr unsigned int M = E::local_dimension;
  static constexpr unsigned int N = E::embedding_dimension;

  // RNG stuffs
  uint32_t seed = time(NULL);
  std::default_random_engine rng;
  
  // load mesh from .csv files
  SearchEngineTest() {
    std::string dim = (E::manifold) ? std::to_string(M) + ".5" : std::to_string(M);
    // compute file names
    std::string point    = "data/m" + dim + "D_points.csv";
    std::string edges    = "data/m" + dim + "D_edges.csv";
    std::string elements = "data/m" + dim + "D_elements.csv";
    std::string neigh    = "data/m" + dim + "D_neigh.csv";
    std::string boundary = "data/m" + dim + "D_boundary.csv";
    // initialize test object
    m = E(point, edges, elements, neigh, boundary);
    rng = std::default_random_engine(seed);
  };

  // generate element at random inside mesh m
  std::shared_ptr<Element<E::local_dimension, E::embedding_dimension>> generateRandomElement();
  
  // generate point at random inside element e
  SVector<E::embedding_dimension> generateRandomPoint(
      const std::shared_ptr<Element<E::local_dimension, E::embedding_dimension>>
        &e);

  // generates a random set of pairs <element ID, point> where point is an element contained in element ID
  std::vector<std::pair<std::size_t, SVector<E::embedding_dimension>>> generateTestSet();
};

template <typename E>
std::shared_ptr<Element<E::local_dimension, E::embedding_dimension>>
SearchEngineTest<E>::generateRandomElement() {
  std::uniform_int_distribution<int> randomID(0, this->m.elements()-1);
  int ID = randomID(rng);  
  return m.element(ID);
}

template <typename E>
SVector<E::embedding_dimension> SearchEngineTest<E>::generateRandomPoint(
    const std::shared_ptr<Element<E::local_dimension, E::embedding_dimension>>
        &e) {
  std::uniform_real_distribution<double> T(0,1);
  // let t, s, u ~ U(0,1) and P1, P2, P3, P4 a set of points, observe that:
  //     * if P1 and P2 are the vertices of a linear element, p = t*P1 + (1-t)*P2 lies into it for any t ~ U(0,1)
  //     * if P1, P2, P3 are vertices of a triangle, the point P = (1-t)P1 + t((1-s)P2 + sP3) is in the triangle
  //       for any choice of t, s ~ U(0,1)
  //     * if P1, P2, P3, P4 are vertices of a tetrahedron, then letting Q = (1-t)P1 + t((1-s)P2 + sP3) and
  //       P = (1-u)P4 + uQ, P belongs to the tetrahedron for any choice of t, s, u ~ U(0,1)
  double t = T(rng);
  SVector<N> p = t*e->coords()[0] + (1-t)*e->coords()[1];
  for(std::size_t j = 1; j < M; ++j){
    t = T(rng);
    p = (1-t)*e->coords()[1+j] + t*p;
  }
  return p;
}

template <typename E>
std::vector<std::pair<std::size_t, SVector<E::embedding_dimension>>>
SearchEngineTest<E>::generateTestSet() {
  // draw 0.1*m.elements() random pairs <ID, point> on the mesh
  std::vector<std::pair<std::size_t, SVector<E::embedding_dimension>>> result{};
  std::size_t n = 0.1*m.elements();
  result.resize(n);

  for(std::size_t i = 0; i < n; ++i){
    auto e = this->generateRandomElement();
    SVector<E::embedding_dimension> p = this->generateRandomPoint(e);
    result[i] = std::make_pair(e->ID(), p);
  }

  return result;
}

using meshList = ::testing::Types<Mesh2D<>, SurfaceMesh<>, Mesh3D<>, NetworkMesh<>>;
TYPED_TEST_SUITE(SearchEngineTest, meshList);

// in the following a test is passed if **all** the queries are correctly satisfied

TYPED_TEST(SearchEngineTest, BruteForce) {
  // build search engine
  BruteForce<TestFixture::M, TestFixture::N> engine(this->m);
  // build test set
  std::vector<std::pair<std::size_t, SVector<TestFixture::N>>> testSet = this->generateTestSet();
  // test all queries in test set
  for(auto query : testSet){
    auto e = engine.search(query.second);
    EXPECT_TRUE(e->ID() == query.first);
  }
}

TYPED_TEST(SearchEngineTest, AlternatingDigitalTree) {
  // build search engine
  ADT<TestFixture::M, TestFixture::N> engine(this->m);
  // build test set
  std::vector<std::pair<std::size_t, SVector<TestFixture::N>>> testSet = this->generateTestSet();
  // test all queries in test set
  for(auto query : testSet){
    auto e = engine.search(query.second);
    EXPECT_TRUE(e->ID() == query.first);
  }
}

// barycentric walk cannot be applied to manifold mesh, filter out manifold cases at compile time
TYPED_TEST(SearchEngineTest, BarycentricWalkTest) {
  if constexpr(TestFixture::N == TestFixture::M){
    BarycentricWalk<TestFixture::M, TestFixture::N> engine(this->m);
    // build test set
    std::vector<std::pair<std::size_t, SVector<TestFixture::N>>> testSet = this->generateTestSet();
    // test all queries in test set
    for(auto query : testSet){
      auto e = engine.search(query.second);
      EXPECT_TRUE(e->ID() == query.first);
    }    
  }else{
    // nothing to do in manifold cases here.
    SUCCEED();
  }
}


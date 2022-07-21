#include <gtest/gtest.h> // testing framework
#include <limits>
#include <memory>
#include <random>

#include "../fdaPDE/core/utils/Symbols.h"
#include "../fdaPDE/core/MESH/Mesh.h"
using fdaPDE::core::MESH::NetworkMesh;
using fdaPDE::core::MESH::Mesh2D;
using fdaPDE::core::MESH::Mesh3D;
using fdaPDE::core::MESH::SurfaceMesh;
#include "../fdaPDE/core/MESH/Element.h"
using fdaPDE::core::MESH::Element;
#include "../fdaPDE/core/utils/IO/CSVReader.h"

// correctness of contained information in an Element object has been asserted in MeshTest.cpp test suite, please refer to it if you believe
// that some information stored in an Element object is not correct (is indeed responsibility of Mesh.h to build the single Element instance).
// Here we assume informations packed in an instance of Element class are correct and test for correctness of operations on elements.

// sample meshes used for testing: see m*D_*.csv data series in test/data/ folder to inspect raw informations
//     * 2D:   unit disk with:   441  2D points, 780   elements, 1220  edges. /test/data/m2D_*.csv
//     * 3D:   unit sphere with: 4193 3D points, 22200 elements, 45380 faces. /test/data/m3D_*.csv
//     * 2.5D: manifold with:    340  3D points, 616   elements, 956   edges. /test/data/m2.5D_*.csv
//     * 1.5D: graph with:       204  2D points, 559   elements               /test/data/m1.5D_*.csv

// test fixture
template <typename E>
class ElementTest : public ::testing::Test {
public:
  E m; // mesh
  static constexpr unsigned int M = E::local_dimension;
  static constexpr unsigned int N = E::embedding_dimension;

  // RNG stuffs
  uint32_t seed = time(NULL);
  std::default_random_engine rng;
  
  // load mesh from .csv files
  ElementTest() {
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
};

template <typename E>
std::shared_ptr<Element<E::local_dimension, E::embedding_dimension>>
ElementTest<E>::generateRandomElement() {
  std::uniform_int_distribution<int> randomID(0, this->m.elements()-1);
  int ID = randomID(rng);  
  return m.element(ID);
}

template <typename E>
SVector<E::embedding_dimension> ElementTest<E>::generateRandomPoint(
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

using meshList = ::testing::Types<Mesh2D, SurfaceMesh, Mesh3D, NetworkMesh>;
TYPED_TEST_SUITE(ElementTest, meshList);

// check computation of barycentric coordinates is coherent with well known properties
TYPED_TEST(ElementTest, BarycentricCoordinates) {
  // generate a random element, a random point inside it and compute its barycentric coordinates
  auto e = this->generateRandomElement();
  SVector<TestFixture::N> p = this->generateRandomPoint(e);
  SVector<TestFixture::M + 1> q = e->toBarycentricCoords(p); // compute barycentric coordinates of point p

  // the barycentric coordinates of a point inside the element sums to 1
  EXPECT_DOUBLE_EQ(q.sum(), 1);
  // the barycentric coordinates of a point inside the element are all positive
  EXPECT_TRUE((q.array() >= 0).all());
  // a point outside the element has at least one negative barycentric coordinate
  auto f = this->generateRandomElement();
  while(e->ID() == f->ID()){
    f = this->generateRandomElement();
  }
  EXPECT_FALSE((e->toBarycentricCoords(f->midPoint()).array() > 0).all()) << f->midPoint() << "\n----" << e->toBarycentricCoords(f->midPoint()) << "\n-----" << f->coords()[0] << "\n~~~~~~" << f->coords()[1] << "\nAAA " << f->ID() << "\nBBB " << e->ID() << "\n~~~~~~~" << e->coords()[0] << "\n --" << e->coords()[1];
  
  // the barycentric coordinates of the mid point of an element are all equal to (1+M)^{-1} (M is the dimension of the space
  // where the element belongs)
  SVector<TestFixture::M + 1> expected = SVector<TestFixture::M + 1>::Constant(1.0/(1+TestFixture::M));
  EXPECT_TRUE((e->toBarycentricCoords(e->midPoint()) - expected).norm() < 50*std::numeric_limits<double>::epsilon());

  // a vertex has all its barycentric coordinates equal to 0 except a single one
  for(std::size_t i = 0; i < e->coords().size(); ++i){
    SVector<TestFixture::N> node = e->coords()[i];
    q = e->toBarycentricCoords(node);
    EXPECT_TRUE(((q.array()-1).abs() < 50*std::numeric_limits<double>::epsilon()).count() == 1 &&
		(q.array() < 50*std::numeric_limits<double>::epsilon()).count() == TestFixture::M);
  }
}

// test midpoint is correctly computed
TYPED_TEST(ElementTest, MidPoint) {
  // generate random element
  auto e = this->generateRandomElement();
  SVector<TestFixture::N> m = e->midPoint();
  SVector<TestFixture::M + 1> b = e->toBarycentricCoords(m);
  // the midpoint of an element is strictly inside it <-> its barycentric coordinates are all strictly positive
  EXPECT_TRUE((b.array() > 0).all());
  // the barycentric coordinates of the midpoint are all approximately equal (the midpoint is the center of mass)
  for(std::size_t i = 0; i < TestFixture::M; ++i){
    double x = b[i];
    for(std::size_t j = i+1; j < TestFixture::M + 1; ++j){
      double y = b[j];
      EXPECT_NEAR(x, y, 50*std::numeric_limits<double>::epsilon());
    }
  }
}

// chcek .contains is able to correctly evaluate when a point is contained or not in an element
TYPED_TEST(ElementTest, CanAssertIfPointIsInside) {
  // draw some element at random
  auto e = this->generateRandomElement();
  // expect the mid point of the element is contained in the element itself
  EXPECT_TRUE(e->contains(e->midPoint()));
  
  // generate random points inside the element and check they are all contained into it
  for(std::size_t i = 0; i < 100; ++i){
    SVector<TestFixture::N> p = this->generateRandomPoint(e);
    EXPECT_TRUE(e->contains(p));
  }

  // expect the midpoint of a different element is not contained in e
  auto f = this->generateRandomElement();
  while(e->ID() == f->ID()){
    f = this->generateRandomElement();
  }
  EXPECT_FALSE(e->contains(f->midPoint()));
}

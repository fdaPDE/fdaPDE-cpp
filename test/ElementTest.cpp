#include <cstddef>
#include <gtest/gtest-typed-test.h>
#include <gtest/gtest.h> // testing framework
#include <random>
#include <string>
#include <string_view>
#include <type_traits>
#include <unistd.h>
#include <unordered_set>
#include <vector>

#include "../fdaPDE/core/utils/Symbols.h"
#include "../fdaPDE/core/MESH/Mesh.h"
using fdaPDE::core::MESH::NetworkMesh;
using fdaPDE::core::MESH::Mesh2D;
using fdaPDE::core::MESH::Mesh3D;
using fdaPDE::core::MESH::SurfaceMesh;
using fdaPDE::core::MESH::is_linear_network;
#include "../fdaPDE/core/MESH/Element.h"
using fdaPDE::core::MESH::Element;
#include "../fdaPDE/core/utils/IO/CSVReader.h"

// quality of contained information in an Element object has been asserted in MeshTest.cpp test suite, please refer to it if you believe
// that some information stored in an Element object is correct (is indeed responsibility of Mesh.h to build the single Element instance).
// here we assume informations packed in an instance of Element class is correct.

// tests are performed against a single element of the mesh for any kind of mesh

template <typename E>
class ElementTest : public ::testing::Test {
public:
  E m; // mesh
  static constexpr unsigned int M = E::local_dimension;
  static constexpr unsigned int N = E::embedding_dimension;
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
  };
};

using meshList = ::testing::Types<Mesh2D, SurfaceMesh, Mesh3D, NetworkMesh>;
TYPED_TEST_SUITE(ElementTest, meshList);

TYPED_TEST(ElementTest, contains) {
  // prepare RNG
  std::default_random_engine rng;
  std::uniform_int_distribution<int> randomID(0, this->m.elements()-1);
  // draw some element at random
  int ID1 = randomID(rng);
  auto e = this->m.element(ID1);
  // expect the mid point of the element is contained in the element itself
  EXPECT_TRUE(e->contains(e->midPoint()));
  
  // now generate random points inside the element and check they are all contained into it
  std::uniform_real_distribution<double> T(0,1);
  for(std::size_t i = 0; i < 100; ++i){
    // let t, s, u ~ U(0,1) and P1, P2, P3, P4 a set of points, observe that:
    //     * if P1 and P2 are the vertices of a linear element, p = t*P1 + (1-t)*P2 lies into it for any t ~ U(0,1)
    //     * if P1, P2, P3 are vertices of a triangle, the point P = (1-t)P1 + t((1-s)P2 + sP3) is in the triangle
    //       for any choice of t, s ~ U(0,1)
    //     * if P1, P2, P3, P4 are vertices of a tetrahedron, then letting Q = (1-t)P1 + t((1-s)P2 + sP3) and
    //       P = (1-u)P4 + uQ, P belongs to the tetrahedron for any choice of t, s, u ~ U(0,1)
    double t = T(rng);
    SVector<TestFixture::N> p = t*e->coords()[0] + (1-t)*e->coords()[1];
    for(std::size_t j = 1; j < TestFixture::M; ++j){
      t = T(rng);
      p = (1-t)*e->coords()[1+j] + t*p;
    }
    EXPECT_TRUE(e->contains(p));
  }

  // expect that the midpoint of a different element is not contained in e
  int ID2 = randomID(rng);
  while(ID1 == ID2) ID2 = randomID(rng);
  auto f = this->m.element(ID2);
  EXPECT_FALSE(e->contains(f->midPoint()));
}

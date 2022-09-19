#include <gtest/gtest.h> // testing framework
#include <memory>

#include "../fdaPDE/core/utils/Symbols.h"
#include "../fdaPDE/core/MESH/Element.h"
using fdaPDE::core::MESH::Element;
#include "../fdaPDE/core/NLA/VectorSpace.h"
using fdaPDE::core::NLA::VectorSpace;

#include "utils/MeshLoader.h"
using fdaPDE::testing::MeshLoader;
using fdaPDE::testing::standard_mesh_selector;
using fdaPDE::testing::MESH_TYPE_LIST;
#include "utils/Constants.h"
using fdaPDE::testing::DOUBLE_TOLERANCE;

// test fixture
template <typename E>
class ElementTest : public ::testing::Test {
public:
  MeshLoader<E> meshWrapper;
  static constexpr unsigned int M = MeshLoader<E>::M;
  static constexpr unsigned int N = MeshLoader<E>::N;
  
  ElementTest() : meshWrapper(standard_mesh_selector(E::local_dimension, E::embedding_dimension)) {};
};
TYPED_TEST_SUITE(ElementTest, MESH_TYPE_LIST);

// check computation of barycentric coordinates is coherent with well known properties
TYPED_TEST(ElementTest, BarycentricCoordinates) {
  // generate a random element, a random point inside it and compute its barycentric coordinates
  auto e = this->meshWrapper.generateRandomElement();
  SVector<TestFixture::N> p = this->meshWrapper.generateRandomPoint(e);
  SVector<TestFixture::M + 1> q = e->toBarycentricCoords(p); // compute barycentric coordinates of point p

  // the barycentric coordinates of a point inside the element sums to 1
  EXPECT_DOUBLE_EQ(q.sum(), 1);
  // the barycentric coordinates of a point inside the element are all positive
  EXPECT_TRUE((q.array() >= 0).all());

  // a point outside the element has at least one negative barycentric coordinate
  auto f = this->meshWrapper.generateRandomElement();
  while(e->ID() == f->ID()) f = this->meshWrapper.generateRandomElement();
  
  if constexpr(TestFixture::N == TestFixture::M){
    EXPECT_FALSE((e->toBarycentricCoords(f->midPoint()).array() > 0).all());
  }else{
    // for manifolds we have to consider if the point x belongs to the space spanned by the element
    VectorSpace<TestFixture::M, TestFixture::N> vs = e->spannedSpace();
    EXPECT_FALSE(vs.distance(f->midPoint()) < DOUBLE_TOLERANCE && (e->toBarycentricCoords(f->midPoint()).array() > 0).all());
  }
  
  // the barycentric coordinates of the mid point of an element are all equal to (1+M)^{-1} (M is the dimension of the space
  // where the element belongs)
  SVector<TestFixture::M + 1> expected = SVector<TestFixture::M + 1>::Constant(1.0/(1+TestFixture::M));
  EXPECT_TRUE((e->toBarycentricCoords(e->midPoint()) - expected).norm() < DOUBLE_TOLERANCE);

  // a vertex has all its barycentric coordinates equal to 0 except a single one
  for(std::size_t i = 0; i < e->coords().size(); ++i){
    SVector<TestFixture::N> node = e->coords()[i];
    q = e->toBarycentricCoords(node);
    EXPECT_TRUE(((q.array()-1).abs() < DOUBLE_TOLERANCE).count() == 1 &&
		(q.array() < DOUBLE_TOLERANCE).count() == TestFixture::M);
  }
}

// test midpoint is correctly computed
TYPED_TEST(ElementTest, MidPoint) {
  // generate random element
  auto e = this->meshWrapper.generateRandomElement();
  SVector<TestFixture::N> m = e->midPoint();
  SVector<TestFixture::M + 1> b = e->toBarycentricCoords(m);
  // the midpoint of an element is strictly inside it <-> its barycentric coordinates are all strictly positive
  EXPECT_TRUE((b.array() > 0).all());
  // the barycentric coordinates of the midpoint are all approximately equal (the midpoint is the center of mass)
  for(std::size_t i = 0; i < TestFixture::M; ++i){
    double x = b[i];
    for(std::size_t j = i+1; j < TestFixture::M + 1; ++j){
      double y = b[j];
      EXPECT_NEAR(x, y, DOUBLE_TOLERANCE);
    }
  }
}

// chcek .contains is able to correctly evaluate when a point is contained or not in an element
TYPED_TEST(ElementTest, CanAssertIfPointIsInside) {
  // draw some element at random
  auto e = this->meshWrapper.generateRandomElement();
  // expect the mid point of the element is contained in the element itself
  EXPECT_TRUE(e->contains(e->midPoint()));
  
  // generate random points inside the element and check they are all contained into it
  for(std::size_t i = 0; i < 100; ++i){
    SVector<TestFixture::N> p = this->meshWrapper.generateRandomPoint(e);
    EXPECT_TRUE(e->contains(p));
  }

  // expect the midpoint of a different element is not contained in e
  auto f = this->meshWrapper.generateRandomElement();
  while(e->ID() == f->ID()) f = this->meshWrapper.generateRandomElement();
  EXPECT_FALSE(e->contains(f->midPoint()));
}

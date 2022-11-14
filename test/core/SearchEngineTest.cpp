#include <gtest/gtest.h> // testing framework
#include <memory>

#include "../../fdaPDE/core/utils/Symbols.h"
#include "../../fdaPDE/core/MESH/Element.h"
using fdaPDE::core::MESH::Element;
#include "../../fdaPDE/core/MESH/engines/BruteForce/BruteForce.h"
#include "../../fdaPDE/core/MESH/engines/BarycentricWalk/BarycentricWalk.h"
#include "../../fdaPDE/core/MESH/engines/AlternatingDigitalTree/ADT.h"
using fdaPDE::core::MESH::BruteForce;
using fdaPDE::core::MESH::BarycentricWalk;
using fdaPDE::core::MESH::ADT;

#include "../utils/MeshLoader.h"
using fdaPDE::testing::MeshLoader;
using fdaPDE::testing::MESH_TYPE_LIST;

// test fixture. ADT and bruteforce can work on this fixture (barycentric walk is not able to handle manifolds)
template <typename E>
struct SearchEngineTest : public ::testing::Test {
  MeshLoader<E> meshLoader{}; // use default mesh
  static constexpr unsigned int M = MeshLoader<E>::M;
  static constexpr unsigned int N = MeshLoader<E>::N;
  static constexpr unsigned int R = MeshLoader<E>::R;
};
TYPED_TEST_SUITE(SearchEngineTest, MESH_TYPE_LIST);

// in the following a test is passed if **all** the queries are correctly satisfied
TYPED_TEST(SearchEngineTest, BruteForce) {
  // build search engine
  BruteForce<TestFixture::M, TestFixture::N, TestFixture::R> engine(this->meshLoader.mesh);
  // build test set
  std::vector<std::pair<std::size_t, SVector<TestFixture::N>>> testSet
    = this->meshLoader.sample(100);
  // test all queries in test set
  std::size_t matches = 0;
  for(auto query : testSet){
    auto e = engine.search(query.second);
    if(e != nullptr && e->ID() == query.first)
      matches++;
  }
  EXPECT_EQ(matches, 100);
}

TYPED_TEST(SearchEngineTest, AlternatingDigitalTree) {
  // build search engine
  ADT<TestFixture::M, TestFixture::N, TestFixture::R> engine(this->meshLoader.mesh);
  // build test set
  std::vector<std::pair<std::size_t, SVector<TestFixture::N>>> testSet
    = this->meshLoader.sample(100);
  // test all queries in test set
  std::size_t matches = 0;
  for(auto query : testSet){
    auto e = engine.search(query.second);
    if(e != nullptr && e->ID() == query.first){
      matches++;
    }
  }
  EXPECT_EQ(matches, 100);
}

// barycentric walk cannot be applied to manifold mesh, filter out manifold cases at compile time
TYPED_TEST(SearchEngineTest, BarycentricWalkTest) {
  if constexpr(TestFixture::N == TestFixture::M){
    BarycentricWalk<TestFixture::M, TestFixture::N, TestFixture::R> engine(this->meshLoader.mesh);
    // build test set
    std::vector<std::pair<std::size_t, SVector<TestFixture::N>>> testSet
      = this->meshLoader.sample(100);
    // test all queries in test set
    std::size_t matches = 0;
    for(auto query : testSet){
      auto e = engine.search(query.second);
      if(e != nullptr && e->ID() == query.first)
	matches++;
    }
    EXPECT_EQ(matches, 100);
  }else{
    // nothing to do in manifold cases here.
    SUCCEED();
  }
}

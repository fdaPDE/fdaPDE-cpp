#include <gtest/gtest.h> // testing framework
#include <vector>
#include <unordered_set>
#include <memory>
#include <random>

#include "../../../fdaPDE/core/utils/Symbols.h"
#include "../../../fdaPDE/core/MESH/Element.h"
using fdaPDE::core::MESH::Element;
#include "../../utils/MeshLoader.h"
using fdaPDE::testing::MeshLoader;
using fdaPDE::testing::standard_mesh_selector;
using fdaPDE::testing::MESH_TYPE_LIST;

// test suite for testing both non-manifold meshes (2D/3D) and manifold mesh (2.5D/1.5D)
template <typename E>
class MeshTest : public ::testing::Test {
public:
  MeshLoader<E> meshWrapper;
  static constexpr unsigned int M = MeshLoader<E>::M;
  static constexpr unsigned int N = MeshLoader<E>::N;
  
  MeshTest() : meshWrapper(standard_mesh_selector(E::local_dimension, E::embedding_dimension)) {};
};
TYPED_TEST_SUITE(MeshTest, MESH_TYPE_LIST);

// check that mesh correctly load in memory the mesh information
TYPED_TEST(MeshTest, DimensionsAreCorrect) {
  // check that at least dimensions are correct
  std::size_t ne = this->meshWrapper.elementsCSV.rows();
  std::size_t nn = this->meshWrapper.pointsCSV.rows();
  
  EXPECT_EQ(this->meshWrapper.mesh_.elements(), ne);
  EXPECT_EQ(this->meshWrapper.mesh_.nodes(),    nn);
}

// check points' coordinate embedded in an element are loaded correctly
TYPED_TEST(MeshTest, PointCoordinatesAreLoadedCorrectly) {
  for(std::size_t i = 0; i < this->meshWrapper.mesh_.elements(); ++i){
    // request element with ID i
    std::unique_ptr<Element<TestFixture::M, TestFixture::N>> e = this->meshWrapper.mesh_.element(i);
    
    // check coordinates stored in element built from Mesh object match raw informations
    int j = 0;
    for(int nodeID : this->meshWrapper.elementsCSV.row(i)){
      // recall that raw files haven't the index correction of one unit, need to subtract 1 from nodeID
      std::vector<double> coords = this->meshWrapper.pointsCSV.row(nodeID-1);
      SVector<TestFixture::N> ePoint = e->coords()[j];
      
      for(std::size_t idx = 0; idx < TestFixture::N; ++idx)
	EXPECT_DOUBLE_EQ(coords[idx], ePoint[idx]);
      j++;
    }
  }
}

// check neighboring identifiers embedded in an element are loaded correctly
TYPED_TEST(MeshTest, NeighboringStructureIsLoadedCorrectly) {
  for(std::size_t i = 0; i < this->meshWrapper.mesh_.elements(); ++i){
    // request element with ID i
    std::unique_ptr<Element<TestFixture::M, TestFixture::N>> e = this->meshWrapper.mesh_.element(i);
    // request data from raw file
    std::vector<int> neigh = this->meshWrapper.neighCSV.row(i);
    // take neighboring information packed inside the element built from Mesh
    auto eNeigh = e->neighbors();

    // check that all claimed neighbors are indeed so
    for(int n : neigh){
      // need to subtract 1 from n for index alignment
      auto search_it = std::find(eNeigh.begin(), eNeigh.end(), n-1);
      EXPECT_TRUE(search_it != eNeigh.end());
      eNeigh.erase(search_it); 
    }
    // at the end we expect there are no more neighbors in e than the ones stored in the raw data file
    EXPECT_TRUE(eNeigh.empty());
  }
}

// performs some checks on the mesh topology, e.g. checks that stated neighbors shares exactly M points
TYPED_TEST(MeshTest, MeshTopologyChecks) {
  // cycle over all mesh elements
  for(std::size_t i = 0; i < this->meshWrapper.mesh_.nodes(); ++i){
    // request the element with ID i
    std::unique_ptr<Element<TestFixture::M, TestFixture::N>> e = this->meshWrapper.mesh_.element(i);
    // check that neighboing elements have always M points in common, so that if we consider two elements as
    // neighbors they are geometrically so (there is a face in common, which is what we expect from neighboring relation)
    for(int neighID : e->neighbors()){
      if(!e->isOnBoundary()){
	// request neighboring element from mesh
	std::unique_ptr<Element<TestFixture::M, TestFixture::N>> n = this->meshWrapper.mesh_.element(neighID);
	// take nodes of both elements
	std::array<SVector<TestFixture::N>, TestFixture::M + 1> eList = e->coords(), nList = n->coords();
	// check that the points in common between the two are exactly M
        std:size_t matches = 0;
	for(SVector<TestFixture::N> p : eList){
	  if(std::find(nList.begin(), nList.end(), p) != nList.end())
	    matches++;
	  }
	EXPECT_TRUE(matches == TestFixture::M);
      }else{
	// check that Mesh is choerent with the fact that e is on bonudary, i.e. at least one vertex of e is detected as boundary point
	bool element_on_boundary = false;
	auto boundaryNodes = e->boundaryNodes();
	for(auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it){
	  if(this->meshWrapper.mesh_.isOnBoundary(it->first)){ // mesh detects this point is a boundary point
	    element_on_boundary = true;
	  }
	}
	EXPECT_TRUE(element_on_boundary);
      }
    }
  }
}

// check the range for loop scans the whole mesh element by element
TYPED_TEST(MeshTest, RangeForLoop) {
  // prepare set with all indexes of IDs to touch
  std::unordered_set<int> meshIDs{};
  for(int i = 0; i < this->meshWrapper.mesh_.elements(); ++i)
    meshIDs.insert(i);

  // range-for over all elements removing the element's ID from the above set when the element is visited
  for(const auto& e : this->meshWrapper.mesh_){
    // check element ID still present in the IDs set (ID not visisted by means of a different element)
    EXPECT_TRUE(meshIDs.find(e->ID()) != meshIDs.end());
    meshIDs.erase(e->ID());
  }
  // check that no ID is left in the initial set
  EXPECT_TRUE(meshIDs.empty());
}

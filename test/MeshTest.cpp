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

// test suite for testing both non-manifold meshes (2D/3D) and manifold mesh (2.5D/1.5D)

// sample meshes used for testing: see m*D_*.csv data series in test/data/ folder to inspect raw informations
//     * 2D:   unit disk with:   441 2D points, 780  elements, 1220 edges. /test/data/m2D_*.csv
//     * 3D:   unit sphere with: 587 3D points, 2775 elements, 5795 faces. /test/data/m3D_*.csv
//     * 2.5D: manifold with:    340 3D points, 616  elements, 956  edges. /test/data/m2.5D_*.csv
//     * 1.5D: graph with:       204 2D points, 559  elements              /test/data/m1.5D_*.csv
// observe that 1.5D testing requires special logic. Specialized code is compiled to test the neighboring structure stored by 
// LinearNetworkMesh objects

template <typename E>
class MeshTest : public ::testing::Test {
public:
  E m; // mesh
  static constexpr unsigned int M = E::local_dimension;
  static constexpr unsigned int N = E::embedding_dimension;
  
  CSVReader reader{}; // csv parser
  // raw files
  CSVFile<double> pointsData;
  CSVFile<int> elementsData;
  // cope with different storage strategy adopted by linear network meshes
  typename std::conditional<
    !is_linear_network<M, N>::value, CSVFile<int>, CSVSparseFile<int>
    >::type neighboringData;

  // load mesh from .csv files
  MeshTest() {
    std::string dim = (E::manifold) ? std::to_string(M) + ".5" : std::to_string(M);
    // compute file names
    std::string point    = "data/m" + dim + "D_points.csv";
    std::string edges    = "data/m" + dim + "D_edges.csv";
    std::string elements = "data/m" + dim + "D_elements.csv";
    std::string neigh    = "data/m" + dim + "D_neigh.csv";
    std::string boundary = "data/m" + dim + "D_boundary.csv";
    // initialize test objects
    m = E(point, edges, elements, neigh, boundary);
    pointsData = reader.parseFile<double>(point);
    elementsData = reader.parseFile<int>(elements);
    // proper parse the neighboring information in case of 1.5D meshes
    if constexpr(!is_linear_network<M, N>::value)  
      neighboringData = reader.parseFile<int>(neigh);
    else
      neighboringData = reader.parseSparseFile<int>(neigh);
  };

  std::vector<SVector<E::embedding_dimension>> collectPointsFromFile(std::size_t elementID) const;
  std::vector<int> collectNeighborsFromFile(std::size_t elementID) const;
  
};

template <typename E>
std::vector<SVector<E::embedding_dimension>> MeshTest<E>::collectPointsFromFile(std::size_t elementID) const {
  // take raw informations from file
  std::vector<SVector<N>> rawPointSet; // the set of vertices' coordinates of e directly coming from raw file
  auto rawElementFile = elementsData.getRawParsedFile();
  for(auto it = rawElementFile.begin(); it != rawElementFile.end(); ++it){	
    // C++ starts counting from 0 -> the element with index i > 0 in mesh has index i-1 in parsed .csv
    std::size_t rawPointID = it->second[elementID] - 1; // subtract 1 to realign indexes!
    // observe that Mesh automatically subtracts 1 to all indexes at construction time, so that this realignment is never required.
    // Here we are doing bare-metal testing to stress the Mesh module hence we should not trust what is coming from Mesh
    SVector<N> rawPoint{};
    for(std::size_t coord = 0; coord < N; ++coord){
      std::string columnName = "V" + std::to_string(coord+1); // column name as stored in .csv files
      rawPoint[coord] = pointsData.getRawParsedFile()[columnName][rawPointID];
    }
    rawPointSet.push_back(rawPoint);
  }
  return rawPointSet;
}

template <typename E>
std::vector<int> MeshTest<E>::collectNeighborsFromFile(std::size_t elementID) const {
  // take raw informations from file
  std::vector<int> rawNeighboringSet{};
  auto rawNeighboring_file = neighboringData.getRawParsedFile();
    
  // indexes of neighbors' element built from Mesh are all contained in the indexes set built from raw information?
  for(auto it = rawNeighboring_file.begin(); it != rawNeighboring_file.end(); ++it){
    if constexpr(!is_linear_network<M, N>::value){
      int neighborID = it->second[elementID] - 1; // subtract 1 to realign indexes!
      rawNeighboringSet.push_back(neighborID);
    }else{
      for(auto neighID = it->second[elementID].begin(); neighID != it->second[elementID].end(); ++neighID){
	int neighborID = *neighID - 1; // subtract 1 to realign indexes!
	rawNeighboringSet.push_back(neighborID);
      }
    }
  }
  return rawNeighboringSet;
}

using meshList = ::testing::Types<Mesh2D<>, SurfaceMesh<>, Mesh3D<>, NetworkMesh<>>;
TYPED_TEST_SUITE(MeshTest, meshList);

// check that mesh correctly load in memory the mesh information
TYPED_TEST(MeshTest, DimensionsAreCorrect) {
  // check that at least dimensions are correct
  std::size_t ne = this->elementsData.rows();
  std::size_t nn = this->pointsData.rows();
  
  EXPECT_EQ(this->m.elements(), ne);
  EXPECT_EQ(this->m.nodes(),    nn);  
}

// check points' coordinate embedded in an element are loaded correctly
TYPED_TEST(MeshTest, PointCoordinatesArePackedCorrectly) {
  // prepare RNG
  std::default_random_engine rng;
  std::uniform_int_distribution<int> randomID(0, this->m.elements()-1);
  // draw some elements at random and test if contained information reflects raw information
  for(std::size_t i = 0; i < 0.1*this->m.nodes(); ++i){
    std::size_t elementID = randomID(rng); // draw an ID at random
    
    // take raw informations from file
    std::vector<SVector<TestFixture::N>> rawPointSet = this->collectPointsFromFile(elementID);
    
    // test MESH: request the element with elementID ID...
    std::shared_ptr<Element<TestFixture::M, TestFixture::N>> e = this->m.element(elementID);
    // and check coordinates coming from the element built from Mesh match raw informations
    EXPECT_TRUE(e->coords().size() == TestFixture::M + 1);
    for(std::size_t j = 0; j < e->coords().size(); ++j){
      SVector<TestFixture::N> ePoint = e->coords()[j];
      EXPECT_TRUE(std::find(rawPointSet.begin(), rawPointSet.end(), ePoint) != rawPointSet.end());
    }
  }
}

// check neighboring identifiers embedded in an element are loaded correctly
TYPED_TEST(MeshTest, NeighboringInformationsArePackedCorrectly) {
  // prepare RNG
  std::default_random_engine rng;
  std::uniform_int_distribution<int> randomID(0, this->m.elements()-1);
  // draw some elements at random and test if contained information reflects raw information
  for(std::size_t i = 0; i < 0.1*this->m.nodes(); ++i){
    std::size_t elementID = randomID(rng); // draw an ID at random

    // take raw informations from file
    std::vector<int> rawNeighboringSet = this->collectNeighborsFromFile(elementID);

    // test MESH: request the element with elementID ID...
    std::shared_ptr<Element<TestFixture::M, TestFixture::N>> e = this->m.element(elementID);
    // take neighboring information packed inside the element built from Mesh
    auto neighbors = e->neighbors();
    // and check that all claimed neighbors are indeed so
    for(int n : rawNeighboringSet){
      auto search_it = std::find(neighbors.begin(), neighbors.end(), n);
      EXPECT_TRUE(search_it != neighbors.end());
      neighbors.erase(search_it); 
    }
    // at the end we expect there are no more neighbors in e than the ones stored in the raw data file
    EXPECT_TRUE(neighbors.empty());
  }
}

// performs some checks on the mesh topology, e.g. checks that stated neighbors shares exactly M points
TYPED_TEST(MeshTest, MeshTopologyChecks) {
  // prepare RNG
  std::default_random_engine rng;
  std::uniform_int_distribution<int> randomID(0, this->m.elements()-1);

  for(std::size_t i = 0; i < 0.1*this->m.nodes(); ++i){
    std::size_t elementID = randomID(rng); // draw an ID at random
    // request the element with that ID
    std::shared_ptr<Element<TestFixture::M, TestFixture::N>> e = this->m.element(elementID);
    // check that neighboing elements have always M points in common, this tests that if we consider two elements as
    // neighbors they are geometrically linked (there is a face in common, which is what we expect from neighboring relation)
    for(int neighID : e->neighbors()){
      if(!e->isOnBoundary()){
	// request neighboring element from mesh
	std::shared_ptr<Element<TestFixture::M, TestFixture::N>> n = this->m.element(neighID);
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
	  if(this->m.isOnBoundary(it->first)){ // mesh detects this point is a boundary point
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
  for(int i = 0; i < this->m.elements(); ++i)
    meshIDs.insert(i);

  // range-for over all elements removing the element's ID from the above set when the element is visited
  for(const auto& e : this->m){
    // check element ID still present in the IDs set (ID not visisted by means of a different element)
    EXPECT_TRUE(meshIDs.find(e->ID()) != meshIDs.end());
    meshIDs.erase(e->ID());
  }
  // check that no ID is left in the initial set
  EXPECT_TRUE(meshIDs.empty());
}

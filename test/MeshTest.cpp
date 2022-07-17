#include <cstddef>
#include <gtest/gtest.h> // testing framework
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

#include "../fdaPDE/core/utils/Symbols.h"
#include "../fdaPDE/core/MESH/Mesh.h"
using fdaPDE::core::MESH::Mesh2D;
#include "../fdaPDE/core/MESH/Element.h"
using fdaPDE::core::MESH::Element;
#include "../fdaPDE/core/MESH/CSVReader.h"

// test suite for testing non-manifold meshes (2D/3D)
// mesh used for testing is a triangulated unit disk maden by 441 2D points (nodes), 780 elements and 1220 edges.
// see cirlce_*.csv file in test_data/ folder to inspect raw information

class Mesh2DTest : public ::testing::Test {
protected:
  Mesh2D m; // loaded mesh
  CSVReader reader{}; // csv parser
  // raw files
  CSVFile<double> pointFile;
  CSVFile<int> elementFile, neighborFile;
  // load mesh from .csv files
  Mesh2DTest() : m(Mesh2D("circle_points.csv", "circle_edges.csv", "circle_triangles.csv", "circle_neighbors.csv", "circle_boundary_markers.csv")),
		 pointFile(reader.parseFile<double>("circle_points.csv")),
		 elementFile(reader.parseFile<int>("circle_triangles.csv")),
		 neighborFile(reader.parseFile<int>("circle_neighbors.csv")) {};
};

// check that mesh correctly load in memory the mesh information
TEST_F(Mesh2DTest, DimensionsAreCorrect) {
  // check that at least dimensions are correct
  EXPECT_EQ(m.getNumberOfElements(), 780);
  EXPECT_EQ(m.getNumberOfNodes(), 441);  
}

TEST_F(Mesh2DTest, PointCoordinatesAreLoadedCorrectly) {
  // prepare RNG
  std::default_random_engine rng;
  std::uniform_int_distribution<int> randomID(0, m.getNumberOfElements()-1);
  // draw some elements at random and test if contained information reflects raw information
  for(std::size_t i = 0; i < 0.1*m.getNumberOfNodes(); ++i){
    std::size_t elementID = randomID(rng); // draw an ID at random
    // request the element with that ID
    std::shared_ptr<Element<2,2>> e = m.requestElementById(elementID);

    // check vertices coordinates are loaded correctly
    std::vector<SVector<2>> rawPointSet; // the set of vertices' coordinates of e directly coming from raw file
    auto rawElementFile = elementFile.getRawParsedFile();
    for(auto it = rawElementFile.begin(); it != rawElementFile.end(); ++it){	
      // C++ starts counting from 0 -> the element with index i > 0 in mesh has index i-1 in fdaPDE internals
      std::size_t rawPointID = it->second[elementID] - 1; // subtract 1 to realign indexes!
      // observe that Mesh automatically subtracts 1 to all indexes at construction time, so that this realignment is
      // never required. Here we are doing bare-metal testing to strees the Mesh module hence we should not trust what is coming from Mesh
      double x = pointFile.getRawParsedFile()["V1"][rawPointID];
      double y = pointFile.getRawParsedFile()["V2"][rawPointID];
      SVector<2> rawPoint(x,y);
      rawPointSet.push_back(rawPoint);
    }
    // check coordinates coming from the element built from Mesh match raw informations
    for(std::size_t j = 0; j < e->getCoords().size(); ++j){
      SVector<2> ePoint = e->getCoords()[j];
      EXPECT_TRUE(std::find(rawPointSet.begin(), rawPointSet.end(), ePoint) != rawPointSet.end());
    }
  }
}

TEST_F(Mesh2DTest, NeighboringInformationsAreLoadedCorrectly) {
  // prepare RNG
  std::default_random_engine rng;
  std::uniform_int_distribution<int> randomID(0, m.getNumberOfElements()-1);
  // draw some elements at random and test if contained information reflects raw information
  for(std::size_t i = 0; i < 0.1*m.getNumberOfNodes(); ++i){
    std::size_t elementID = randomID(rng); // draw an ID at random
    // request the element with that ID
    std::shared_ptr<Element<2,2>> e = m.requestElementById(elementID);

    // check neighboring information are loaded correctly
    std::vector<int> rawNeighSet;
    auto rawNeighFile = neighborFile.getRawParsedFile();
    auto eNeighbors   = e->getNeighbors();
    // indexes of neighbors' element built from Mesh are all contained in the indexes set built from raw information?
    for(auto it = rawNeighFile.begin(); it != rawNeighFile.end(); ++it){
      int neighID = it->second[elementID] - 1; // subtract 1 to realign indexes!
      EXPECT_TRUE(std::find(eNeighbors.begin(), eNeighbors.end(), neighID) != eNeighbors.end());
    }
  }
}

TEST_F(Mesh2DTest, MeshTopologyChecks) {
  // prepare RNG
  std::default_random_engine rng;
  std::uniform_int_distribution<int> randomID(0, m.getNumberOfElements()-1);

  for(std::size_t i = 0; i < 0.1*m.getNumberOfNodes(); ++i){
    std::size_t elementID = randomID(rng); // draw an ID at random
    // request the element with that ID
    std::shared_ptr<Element<2,2>> e = m.requestElementById(elementID);
    // check that neighboing elements have always 2 points in common, this tests that if we consider two elements as
    // neighbors they are geometrically linked (there is a face in common, which is what we expect from neighboring relation)
    for(int neighID : e->getNeighbors()){
      if(neighID > 0){ 	// by convention, a negative value in the neighboring table means the element is on boundary
	std::shared_ptr<Element<2,2>> n = m.requestElementById(neighID);
	std::array<SVector<2>, 3> pList = e->getCoords(), nList = n->getCoords();
        std:size_t matches = 0;
	for(SVector<2> p : pList){
	  if(std::find(nList.begin(), nList.end(), p) != nList.end())
	    matches++;
	}
	EXPECT_TRUE(matches == 2); // exactly two vertices in common
      }else{ // we are considering the boundary of the domain
	EXPECT_TRUE(e->isOnBoundary());
	// check that Mesh is choerent with this information, i.e. at least one vertex of e is detected as boundary point
	bool element_on_boundary = false;
	for(auto it = e->getBoundaryMarkers().begin(); it != e->getBoundaryMarkers().end(); ++it){
	  if(m.isOnBoundary(it->first)) // mesh detects this point is a boundary point
	    element_on_boundary = true;
	}
	EXPECT_TRUE(element_on_boundary);
      }
    }
  }
}

// check the range for loop scans the whole mesh element by element
TEST_F(Mesh2DTest, RangeForLoop) {
  // prepare set with all indexes of IDs to touch
  std::unordered_set<int> meshIDs{};
  for(int i = 0; i < m.getNumberOfElements(); ++i)
    meshIDs.insert(i);

  // range-for over all elements removing the element's ID from the above set when the element is visited
  for(const auto& e : m){
    // check element ID still present in the IDs set (ID not visisted by means of a different element)
    EXPECT_TRUE(meshIDs.find(e->getID()) != meshIDs.end());
    meshIDs.erase(e->getID());
  }
  // check that no ID is left in the initial set
  EXPECT_TRUE(meshIDs.empty());
}

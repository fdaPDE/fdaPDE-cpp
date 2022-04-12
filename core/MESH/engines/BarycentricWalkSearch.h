#ifndef __BARYCENTRIC_WALK_SEARCH_H__
#define __BARYCENTRIC_WALK_SEARCH_H__

#include "Mesh.h"
#include <random>

typedef std::mt19937 RNG;  // the Mersenne Twister with a popular choice of parameters

template <unsigned int N, unsigned int M>
class BarycentricWalkSearch{

private:
  const Mesh<N,M>& mesh_;

  // specific members of the engine
  uint32_t seed;
  RNG rng;
  std::uniform_int_distribution<uint32_t> uniform_int;
  
public:
  BarycentricWalkSearch(const Mesh<N,M>& mesh) : mesh_(mesh) {
    seed = time(NULL);  // seed for RNG
    rng  = RNG(seed);   // define RNG
    
    // define uniform distribution over the ID space
    uniform_int = std::uniform_int_distribution<uint32_t>(0, mesh_.getNumberOfElements()-1);
  }

  // applies a barycentric walk strategy to search for the element containing a given point
  std::shared_ptr<Element<N, M>> search(const SVector<N>& point);
};
  
// applies a barycentric walk search
template <unsigned int N, unsigned int M>
std::shared_ptr<Element<N, M>> BarycentricWalkSearch<N, M>::search(const SVector<N>& point){

  // start from an element at random
  std::shared_ptr<Element<N,M>> element = mesh_.requestElementById(uniform_int(rng)); 

  if(element->contains(point)){
    return element;
  }
  
  while(!element->contains(point)){
    // compute barycantric coordinates with respect to the element
    SVector<N+1> baryCoord = element->computeBarycentricCoordinates(point);
  
    // Pick the vertices corresponding to the n highest coordinates, and move into the adjacent element that
    // shares those vertices. This is equivalent to find the minimum baricentric coordinate and move to
    // the element adjacent to the face opposite to this point
    unsigned int minBaryCoordIndex;
    baryCoord.minCoeff(&minBaryCoordIndex);

    // by construction barycentric coordinate at position i is relative to vertex i of the mesh element
    // we can move to the next element in O(1) exploiting the memory representation of mesh in memory
    unsigned int nextID = element->getNeighbors()[minBaryCoordIndex];
    element = mesh_.requestElementById(nextID);
  }
  
  return element;
}

#endif // __BARYCENTRIC_WALK_SEARCH_H__

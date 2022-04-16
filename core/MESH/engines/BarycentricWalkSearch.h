#ifndef __BARYCENTRIC_WALK_SEARCH_H__
#define __BARYCENTRIC_WALK_SEARCH_H__

#include "Mesh.h"
#include <random>

typedef std::mt19937 RNG;  // the Mersenne Twister with a popular choice of parameters

template <unsigned int M, unsigned int N>
class BarycentricWalkSearch{

private:
  const Mesh<M,N>& mesh_;

  // a random number generator is used as part of the initialization procedure of the walking algorithm
  uint32_t seed;
  RNG rng;
  std::uniform_int_distribution<uint32_t> uniform_int;
  
public:
  BarycentricWalkSearch(const Mesh<M,N>& mesh) : mesh_(mesh) {
    seed = time(NULL);  // seed for RNG
    rng  = RNG(seed);   // define RNG
    
    // define uniform distribution over the ID space
    uniform_int = std::uniform_int_distribution<uint32_t>(0, mesh_.getNumberOfElements()-1);
  }

  // applies a barycentric walk strategy to search for the element containing a given point
  std::shared_ptr<Element<M, N>> search(const SVector<N>& point);
};

#include "BarycentricWalkSearch.tpp"
  
#endif // __BARYCENTRIC_WALK_SEARCH_H__

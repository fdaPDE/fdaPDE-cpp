#ifndef __BRUTEFORCE_SEARCH_H__
#define __BRUTEFORCE_SEARCH_H__

#include "Mesh.h"
#include "Element.h"
#include <memory>

template <unsigned int M, unsigned int N>
class BruteForceSearch{

private:
  Mesh<M,N>& mesh_;

public:
  // constructor
  BruteForceSearch(Mesh<M,N>& mesh) : mesh_(mesh) {} ;

  // applies a brute force strategy to search for the element containing a given point
  std::shared_ptr<Element<M, N>> search(const SVector<N>& point);
  
};

#include "BruteForceSearch.tpp"

#endif // __BRUTEFORCE_SEARCH_H__

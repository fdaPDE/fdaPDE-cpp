#ifndef __BRUTEFORCE_SEARCH_H__
#define __BRUTEFORCE_SEARCH_H__

#include "Mesh.h"
#include "Element.h"
#include <memory>

template <unsigned int N, unsigned int M>
class BruteForceSearch{

private:
  Mesh<N,M>& mesh_;

public:
  BruteForceSearch(Mesh<N,M>& mesh) : mesh_(mesh) {} ;

  // applies a brute force strategy to search for the element containing a given point
  std::shared_ptr<Element<N, M>> search(const SVector<N>& point);
  
};

template <unsigned int N, unsigned int M>
std::shared_ptr<Element<N, M>> BruteForceSearch<N, M>::search(const SVector<N>& point) {
  
  for(std::shared_ptr<Element<N,M>> element : mesh_){
    if(element->contains(point))
      return element;
  }

  // no element in mesh found
  return std::shared_ptr<Element<N,M>>();
}

#endif // __BRUTEFORCE_SEARCH_H__

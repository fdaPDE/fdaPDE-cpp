#ifndef __GRID_H__
#define __GRID_H__

#include <array>
#include <tuple>
#include "../../utils/Symbols.h"
#include "../../utils/fields/ScalarField.h"
#include "../extensions/Extension.h"

namespace fdaPDE{
namespace core{
namespace OPT{
  // optimize a given scalar field over an N-dimensional grid of equidistant points
  template <unsigned int N>
  class GridOptimizer{

  private:
    std::array<std::pair<double,double>, N> domainLimits; /* an array of pairs where each pair indicate inferior limit
							     and superior limit of the grid along that dimension */
    unsigned int spaceDimension;                          // dimension of the space where the grid is embedded
    std::array<double, N> steps;                          /* an array of double representing the increment step
							     along each dimension (element at position i refers to dimension i) */
  public:
    // constructor
    GridOptimizer(std::array<std::pair<double,double>, N> domain_, std::array<double, N> steps_)
      : spaceDimension(N), domainLimits(domain_), steps(steps_) {}
  
    // optimization routine
    template <typename... Args>
    std::pair<SVector<N>, double> findMinimum(const ScalarField<N>& objective_, const Args&... args);

    std::string description = "Grid optimization over " + spaceDimension + "D space";
  };

#include "Grid.tpp"
}}}

#endif // __GRID_H__

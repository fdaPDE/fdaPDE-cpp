#ifndef __GRID_OPTIMIZER_H__
#define __GRID_OPTIMIZER_H__

#include <array>
#include <tuple>
#include "Function.h"
#include "Optimizer.h"

using std::array;
using std::pair;

// optimize a given scalar field over an N-dimensional grid of equidistant points
template <unsigned int N>
class GridOptimizer : public Optimizer<N>{

private:
  // an array of pairs where each pair indicate inferior limit and superior limit of the grid along that dimension
  array<pair<double,double>, N> domainLimits;
  // dimension of the space where the grid is embedded
  unsigned int spaceDimension;
  // an array of double representing the increment step along each dimension (element at position i refers to dimension i)
  array<double, N> steps;

  const ScalarField<N>& objective;   // objective function to optimize
  
public:
  // constructor
  GridOptimizer(array<pair<double,double>, N> domain_, array<double, N> steps_, const ScalarField<N>& objective_) :
    spaceDimension(N), domainLimits(domain_), steps(steps_), objective(objective_) {}
  // optimization routine
  std::pair<SVector<N>, double> findMinimum() override;
};

#include "GridOptimizer.tpp"

#endif // __GRID_OPTIMIZER_H__

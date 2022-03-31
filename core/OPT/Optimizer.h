#ifndef __OPTIMIZER_H__
#define __OPTIMIZER_H__

#include "Utils.h"
#include <tuple>

// abstract class representing a generic optimization method
template <unsigned int N>
class Optimizer {
  virtual std::pair<SVector<N>, double> findMinimum() = 0;
};

#endif // __OPTIMIZER_H__

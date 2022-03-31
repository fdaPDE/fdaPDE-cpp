#ifndef __ITERATIVE_OPTIMIZER_H__
#define __ITERATIVE_OPTIMIZER_H__

#include <unordered_map>
#include <list>
#include <functional>
#include <string>
#include "Optimizer.h"

// Base class for iterative optimizers.
// Extending this class gives the opportunity for the user
// of the class you are extending to customize the flow of
// execution of the optimization. Please, observe that customization is intended
// for the final user of the algorithm, not for who is writing the algorithm itself
template <unsigned int N>
class IterativeOptimizer : public Optimizer<N> {
protected:
  // use this map to store temporary data used by the controller during the optimization.
  // this structure is returned in output as the result of the optimization process
  std::unordered_map<std::string, std::list<double>> controllerData;
public:
  // controller injected callables, default does nothing
  virtual void init()          { return; };
  virtual void preStep()       { return; };
  virtual void postStep()      { return; };
  virtual bool stopCondition() { return false; };
  virtual void finalize()      { return; };
};

#endif // __ITERATIVE_OPTIMIZER_H__

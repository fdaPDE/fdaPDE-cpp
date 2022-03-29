#ifndef __CUSTOMIZABLE_H__
#define __CUSTOMIZABLE_H__

#include <unordered_map>
#include <list>
#include <functional>
#include <string>
#include "Optimizer.h"

// Base class for customizable optimizers. Extending this class gives the
// opportunity for the user of the class you are extending to customize the flow of
// execution of the optimization. Please, observe that customization is intended
// for the final user of the algorithm, not for who is writing the algorithm
// itself
template <unsigned int N>
class CustomizableOptimizer : public Optimizer<N> {
protected:
  // use this map to store temporary data used by the controller during the optimization.
  // this structure is returned in output as the result of the optimization process
  std::unordered_map<std::string, std::list<double>> controllerData;
public:
  // controller injected callables, default does nothing (about 40 microseconds overhead on average)
  std::function<void(void)> init           = [](void) -> void { return; };
  std::function<void(void)> beginIteration = [](void) -> void { return; };
  std::function<void(void)> endIteration   = [](void) -> void { return; };
  std::function<bool(void)> stopCondition  = [](void) -> bool { return false; };
  std::function<void(void)> finalize       = [](void) -> void { return; };

  // getters and setters for internal data structure
  std::unordered_map<std::string, std::list<double>> getControllerData() { return controllerData; };
  void setControllerData(std::string key, double value) { controllerData[key].push_back(value); };
  void setControllerData(std::string key, std::list<double> value) { controllerData[key] = value; };
};

#endif // __CUSTOMIZABLE_H__

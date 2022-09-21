#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <limits>

// a set of usefull constants
namespace fdaPDE{
namespace testing{

  // the treshold under which two doubles are considered equal
  const double DOUBLE_TOLERANCE = 10*std::numeric_limits<double>::epsilon(); // approx 2*10^-15
  const double MACHINE_EPSILON  = std::numeric_limits<double>::epsilon();
  
}}

#endif // __CONSTANTS_H__

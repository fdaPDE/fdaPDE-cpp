#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <limits>

// a set of usefull constants
namespace fdaPDE{
namespace testing{

  // the treshold under which two doubles are considered equal
  const double MODEL_TOLERANCE  = std::pow(0.1, 12);
  const double DOUBLE_TOLERANCE = 50*std::numeric_limits<double>::epsilon(); // approx 10^-14
  const double MACHINE_EPSILON  = std::numeric_limits<double>::epsilon();    // approx 2.22*10^-16

  // hardcoded value of pi
  constexpr double pi = 3.14159265358979323846;
}}

#endif // __CONSTANTS_H__

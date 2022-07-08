#include <functional>
#include <iostream>

#include "../src/core/utils/fields/ScalarField.h"
#include "utils/Symbols.h"
using fdaPDE::core::ScalarField;

#include <gtest/gtest.h>

// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {

  std::function<double(SVector<2>)> fun = [](SVector<2> x) -> double {
    return x[0] + x[1];
  };

  ScalarField<2> field(fun);
  SVector<2> p(1,1);
  
  // Expect equality.
  EXPECT_EQ(field(p), 2);
}

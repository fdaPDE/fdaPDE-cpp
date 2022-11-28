#include <cstddef>
#include <deque>
#include <gtest/gtest.h> // testing framework

#include "../../fdaPDE/core/utils/Symbols.h"
#include "../../fdaPDE/models/regression/splines/SplineBasis.h"
using fdaPDE::models::SplineBasis;
#include "../utils/Constants.h"
using fdaPDE::testing::DOUBLE_TOLERANCE;
#include "../utils/MeshLoader.h"
using fdaPDE::testing::MeshLoader;
using fdaPDE::testing::MESH_TYPE_LIST;

// test definition of spline basis over a closed interval [0,1]
TEST(SplineBasisTest, Definition) {
  // define vector of equidistant knots on unit interval [0,1]
  DVector<double> knots;
  knots.resize(11);
  std::size_t i = 0;
  for(double x = 0; x <= 1; x += 0.1, ++i) knots[i] = x;
  
  // define cubic B-spline basis over [0,1]
  SplineBasis<3> basis(knots);

  // evaluate spline number 0 over [0,1]
  std::vector<double> result;
  result.reserve(101);
  for(double x = 0; x <= 1.01; x += 0.01){ 
    result.push_back(basis[0](SVector<1>(x)));
  }
  // set expected results
  std::deque<double> expected = {
    1.0000000000000000, 0.7290000000000003, 0.5120000000000001, 0.3430000000000001, 0.2160000000000001,
    0.1250000000000000, 0.0640000000000000, 0.0270000000000000, 0.0080000000000000, 0.0010000000000000};
  for(std::size_t j = expected.size(); j < result.size(); ++j) expected.push_back(0.0);
  // check results within double tolerance
  for(std::size_t j = 0; j < result.size(); ++j){
    EXPECT_NEAR(result[j], expected[j], DOUBLE_TOLERANCE);
  }
  expected.clear(); // prepare for next spline
  result.clear();

  // test spline number 5
  for(double x = 0; x <= 1.01; x += 0.01){ 
    result.push_back(basis[5](SVector<1>(x)));
  }
  // set expected results
  expected = {
    0.0001666666666667, 0.0013333333333333, 0.0045000000000000, 0.0106666666666667, 0.0208333333333334, 
    0.0360000000000001, 0.0571666666666668, 0.0853333333333335, 0.1215000000000003, 0.1666666666666670, 
    0.2211666666666672, 0.2826666666666673, 0.3481666666666673, 0.4146666666666674, 0.4791666666666674, 
    0.5386666666666674, 0.5901666666666674, 0.6306666666666672, 0.6571666666666669, 0.6666666666666666, 
    0.6571666666666662, 0.6306666666666659, 0.5901666666666657, 0.5386666666666655, 0.4791666666666653, 
    0.4146666666666651, 0.3481666666666651, 0.2826666666666650, 0.2211666666666650, 0.1666666666666655, 
    0.1214999999999990, 0.0853333333333325, 0.0571666666666660, 0.0359999999999995, 0.0208333333333330, 
    0.0106666666666664, 0.0044999999999999, 0.0013333333333333, 0.0001666666666667};
  for(double x = 0; x < 0.21; x += 0.01) expected.push_front(0.0);
  for(std::size_t j = expected.size(); j < result.size(); ++j) expected.push_back(0.0);
  // check results within double tolerance
  for(std::size_t j = 0; j < result.size(); ++j){
    EXPECT_NEAR(result[j], expected[j], DOUBLE_TOLERANCE);
  }
  expected.clear(); // prepare for next spline
  result.clear();

  // test spline number 12
  for(double x = 0; x <= 1.01; x += 0.01){ 
    result.push_back(basis[12](SVector<1>(x)));
  }  
  expected = {
    0.0010000000000002, 0.0080000000000008, 0.0270000000000019, 0.0640000000000034, 0.1250000000000054,
    0.2160000000000079, 0.3430000000000109, 0.5120000000000144, 0.7290000000000184, 0.0000000000000000};
  for(double x = 0; x < 0.91; x += 0.01) expected.push_front(0.0);
  // check results within double tolerance
  for(std::size_t j = 0; j < result.size(); ++j){
    EXPECT_NEAR(result[j], expected[j], DOUBLE_TOLERANCE);
  }
}

// test definition of spline basis over a closed interval [0,1]
TEST(SplineBasisTest, SecondDerivative) {
  // define vector of equidistant knots on unit interval [0,1]
  DVector<double> knots;
  knots.resize(11);
  std::size_t i = 0;
  for(double x = 0; x <= 1; x += 0.1, ++i) knots[i] = x;
  
  // define cubic B-spline basis over [0,1]
  SplineBasis<3> basis(knots);

  // evaluate spline number 0 over [0,1]
  std::vector<double> result;
  result.reserve(101);
  for(double x = 0; x <= 1.01; x += 0.01){ 
    result.push_back(basis[0].derive<2>()(SVector<1>(x)));
  }
  // set expected results
  std::deque<double> expected = {
    600.0000000000000000, 540.0000000000001137, 480.0000000000000000, 420.0000000000000568, 360.0000000000000568,
    300.0000000000000000, 240.0000000000000000, 180.0000000000000000, 120.0000000000000284,  60.0000000000000568,
      0.0000000000000833};
  for(std::size_t j = expected.size(); j < result.size(); ++j) expected.push_back(0.0);
  // check results within double tolerance
  for(std::size_t j = 0; j < result.size(); ++j){
    EXPECT_NEAR(result[j], expected[j], DOUBLE_TOLERANCE);
  }
  expected.clear(); // prepare for next spline
  result.clear();

  // test spline number 5
  for(double x = 0; x <= 1.01; x += 0.01){ 
    result.push_back(basis[5].derive<2>()(SVector<1>(x)));
  }
  // set expected results
  expected = {
      0.0000000000000278,  10.0000000000000338,  20.0000000000000355,  30.0000000000000426,  40.0000000000000426,
     50.0000000000000142,  60.0000000000000355,  70.0000000000000284,  80.0000000000000284,  90.0000000000000568,
     99.9999999999998295,  69.9999999999997868,  39.9999999999997442,   9.9999999999997158, -20.0000000000003268,
    -50.0000000000003695, -80.0000000000003979,-110.0000000000004263,-140.0000000000004847,-170.0000000000004832,
   -199.9999999999995737,-169.9999999999994884,-139.9999999999994900,-109.9999999999994031, -79.9999999999993747,
    -49.9999999999993392, -19.9999999999992966,  10.0000000000007390,  40.0000000000007816,  70.0000000000008142,
     99.9999999999998010,  89.9999999999997868,  79.9999999999997726,  69.9999999999997584,  59.9999999999997513,
     49.9999999999997442,  39.9999999999997229,  29.9999999999997087,  19.9999999999996945,   9.9999999999996803};
  for(double x = 0; x < 0.2; x += 0.01) expected.push_front(0.0);
  for(std::size_t j = expected.size(); j < result.size(); ++j) expected.push_back(0.0);
  // check results within double tolerance
  for(std::size_t j = 0; j < result.size(); ++j){
    EXPECT_NEAR(result[j], expected[j], DOUBLE_TOLERANCE);
  }
  expected.clear(); // prepare for next spline
  result.clear();

  // test spline number 12
  for(double x = 0; x <= 1.01; x += 0.01){ 
    result.push_back(basis[12].derive<2>()(SVector<1>(x)));
  }  
  expected = {
      0.0000000000039968,  60.0000000000040785, 120.0000000000041638, 180.0000000000042633, 240.0000000000043201,
    300.0000000000044338, 360.0000000000044906, 420.0000000000046043, 480.0000000000047180, 540.0000000000047748,
      0.0000000000000000};
  for(double x = 0; x < 0.90; x += 0.01) expected.push_front(0.0);
  // check results within double tolerance
  for(std::size_t j = 0; j < result.size(); ++j){
    EXPECT_NEAR(result[j], expected[j], DOUBLE_TOLERANCE);
  }
}

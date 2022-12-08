#include <cstddef>
#include <deque>
#include <gtest/gtest.h> // testing framework
#include <unsupported/Eigen/SparseExtra>

#include "../../fdaPDE/core/utils/Symbols.h"
#include "../../fdaPDE/models/space_time/SplineBasis.h"
using fdaPDE::models::SplineBasis;
#include "../utils/Constants.h"
using fdaPDE::testing::DOUBLE_TOLERANCE;
#include "../utils/MeshLoader.h"
using fdaPDE::testing::MeshLoader;
using fdaPDE::testing::MESH_TYPE_LIST;
#include "../utils/Utils.h"
using fdaPDE::testing::almost_equal;

// test definition of spline basis over a closed interval [0,1]
TEST(SplineBasisTest, Definition) {
  // define vector of equidistant knots on unit interval [0,1]
  DVector<double> knots;
  knots.resize(11);std::size_t i = 0;
  for(double x = 0; x <= 1; x += 0.1, ++i) knots[i] = x;
  
  // define cubic B-spline basis over [0,1]
  SplineBasis<3> basis(knots);

  // load expected results. file misc/spline.mtx contains the evaluation of a correct definition of a Bspline
  // basis over the interval [0,1] defined on knots vector and evaluated at points placed at a distance of 0.01 each
  SpMatrix<double> expected;
  Eigen::loadMarket(expected, "data/misc/spline.mtx");
  
  for(std::size_t i = 0; i < basis.size(); ++i){
    // evaluate i-th spline over [0,1]
    std::vector<double> result;
    result.reserve(101);
    for(double x = 0; x <= 1.01; x += 0.01){ 
      result.push_back(basis[i](SVector<1>(x)));
    }
    // check results within double tolerance
    for(std::size_t j = 0; j < result.size(); ++j){
      EXPECT_TRUE(almost_equal(result[j], expected.coeff(j,i)));
    }    
  }  
}

// test definition of spline basis over a closed interval [0,1]
TEST(SplineBasisTest, SecondDerivative) {
  // define vector of equidistant knots on unit interval [0,1]
  DVector<double> knots;
  knots.resize(11); std::size_t i = 0;
  for(double x = 0; x <= 1; x += 0.1, ++i) knots[i] = x;
  
  // define cubic B-spline basis over [0,1]
  SplineBasis<3> basis(knots);

  // load expected results. file misc/spline_der.mtx contains the evaluation of the second derivative of a correct definition
  // of a Bspline basis over the interval [0,1] defined on knots vector and evaluated at points placed at a distance of 0.01 each
  SpMatrix<double> expected;
  Eigen::loadMarket(expected, "data/misc/spline_der.mtx");
  
  for(std::size_t i = 0; i < basis.size(); ++i){
    // evaluate i-th spline over [0,1]
    std::vector<double> result;
    result.reserve(101);
    for(double x = 0; x <= 1.01; x += 0.01){ 
      result.push_back(basis[i].derive<2>()(SVector<1>(x)));
    }
    // check results within double tolerance
    for(std::size_t j = 0; j < result.size(); ++j){
      EXPECT_TRUE(almost_equal(result[j], expected.coeff(j,i)));
    }    
  }  
}

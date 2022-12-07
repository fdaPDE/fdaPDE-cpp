#include <complex>
#include <cstddef>
#include <gtest/gtest.h> // testing framework
#include <unsupported/Eigen/SparseExtra>

#include "../fdaPDE/core/utils/Symbols.h"
#include "../fdaPDE/models/space_time/SpaceTime.h"
using fdaPDE::models::Phi;
using fdaPDE::models::TimeAssembler;
using fdaPDE::models::TimeMass;
using fdaPDE::models::TimePenalty;

#include "../utils/Utils.h"
using fdaPDE::testing::spLInfinityNorm;
#include "../utils/Constants.h"
using fdaPDE::testing::DOUBLE_TOLERANCE;

// this test suite is for testing the correctness of the assembled matrices required for space-time models

TEST(SpaceTime, CubicSplineEvaluationMatrix) {
  // define time domain
  DVector<double> time;
  time.resize(11);
  std::size_t time_i = 0;
  for(double x = 0; x <= 2; x+=0.2, ++time_i){
    time[time_i] = x;
  }

  // compute \Phi matrix: [\Phi]_{ij} = \phi_i(t_j) using a cubic spline basis
  SpMatrix<double> computedPhi = Phi<SplineBasis<3>>(time);

  // load expected data from file
  SpMatrix<double> expectedPhi;
  Eigen::loadMarket(expectedPhi, "data/misc/Phi.mtx");

  // check for equality under DOUBLE_TOLERANCE
  EXPECT_TRUE( spLInfinityNorm(expectedPhi, computedPhi) < DOUBLE_TOLERANCE);
}

TEST(SpaceTime, CubicSplineTimeMassMatrix) {
  // define time domain
  DVector<double> time;
  time.resize(11);
  std::size_t time_i = 0;
  for(double x = 0; x <= 2; x+=0.2, ++time_i){
    time[time_i] = x;
  }

  // compute P0 time mass matrix using a cubic B-spline basis: [P0]_{ij} = \int_{time} \phi_i*\phi_j
  TimeAssembler<SplineBasis<3>> assembler(time);
  SpMatrix<double> computedP0 = assembler.assemble(TimeMass<SplineBasis<3>>());

  // load expected data from file
  SpMatrix<double> expectedP0;
  Eigen::loadMarket(expectedP0, "data/misc/P0.mtx");

  // check for equality under DOUBLE_TOLERANCE
  EXPECT_TRUE( spLInfinityNorm(expectedP0, computedP0) < DOUBLE_TOLERANCE);
}

TEST(SpaceTime, CubicSplineTimePenaltyMatrix) {
  // define time domain
  DVector<double> time;
  time.resize(11);
  std::size_t time_i = 0;
  for(double x = 0; x <= 2; x+=0.2, ++time_i){
    time[time_i] = x;
  }

  // compute Pt time penalization matrix using a cubic B-spline basis: [Pt]_{ij} = \int_{time} (\phi_i)_tt*(\phi_j)_tt
  TimeAssembler<SplineBasis<3>> assembler(time);
  SpMatrix<double> computedPt = assembler.assemble(TimePenalty<SplineBasis<3>>());

  // load expected data from file
  SpMatrix<double> expectedPt;
  Eigen::loadMarket(expectedPt, "data/misc/Pt.mtx");

  // check for equality under DOUBLE_TOLERANCE
  EXPECT_TRUE( spLInfinityNorm(expectedPt, computedPt) < DOUBLE_TOLERANCE );
}

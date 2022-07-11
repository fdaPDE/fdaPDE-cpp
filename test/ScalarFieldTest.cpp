#include <functional>
#include <vector>
#include <gtest/gtest.h> // testing framework

#include "../src/core/utils/fields/ScalarField.h"
using fdaPDE::core::ScalarField;
using fdaPDE::core::DifferentiableScalarField;
using fdaPDE::core::TwiceDifferentiableScalarField;
#include "utils/Symbols.h"
#include "../src/core/utils/fields/VectorField.h"
using fdaPDE::core::VectorField;

// check if ScalarField class wraps lambda correctly
TEST(ScalarFieldTest, ScalarFieldWrapsCorrectly) {
  // define field expression
  auto fieldExpr = [](SVector<2> x) -> double { // e^x + x^2*y*log(y)
    return std::exp(x[0]) + std::pow(x[0],2)*x[1]*std::log(x[1]);
  };

  // build the ScalarField object
  ScalarField<2> field(fieldExpr);
  // test if the ScalarField wraps correctly the lambda
  SVector<2> p(1,1);
  double trueResult = std::exp(1);
  // expect equality
  EXPECT_EQ(field(p), trueResult);

  // initialize directly with a lambda
  ScalarField<2> lambda_field([](SVector<2> x) -> double { return std::pow(x[0], 2) + x[1]; });
  trueResult = 2;
  EXPECT_EQ(lambda_field(p), trueResult);
}

// checks if ScalarField approximates correctly its analytical gradient
TEST(ScalarFieldTest, GradientApproximation) {
  auto fieldExpr = [](SVector<2> x) -> double { // [e^(2x+y)]/x
      return std::exp(2*x[0]+x[1])/x[0];
    };
  // build the ScalarField object
  ScalarField<2> field(fieldExpr);
  // get the gradient function
  VectorField<2> grad = field.derive(0.01);
  // define evaluation point and true gradient vector
  SVector<2> p(1,1);
  SVector<2> trueGradient = SVector<2>(std::exp(3), std::exp(3));
  // test if the obtained gradient is approximated correctly (truncation error is of order O(h^2), h is the step size)
  EXPECT_TRUE((grad(p) - trueGradient).squaredNorm() < std::pow(0.01, 2));

  // access to gradient approximation directly, without passing to a VectorField object
  SVector<2> approxGradient = field.approxGradient(p, 0.01);
  EXPECT_TRUE((approxGradient - trueGradient).squaredNorm() < std::pow(0.01, 2));
}

// checks if ScalarField approximates correctly its analytical hessian
TEST(ScalarFieldTest, HessianApproximation) {
  auto fieldExpr = [](SVector<2> x) -> double { // e^x + x^2*y*log(y)
      return std::exp(x[0]) + std::pow(x[0],2)*x[1]*std::log(x[1]);
    };
  // build the ScalarField object
  ScalarField<2> field(fieldExpr);
  // get the hessian function
  std::function<SMatrix<2>(SVector<2>)> hess = field.deriveTwice(0.01);
  // define evaluation point and true hessian matrix
  SVector<2> p(1,1);
  SMatrix<2> trueHessian = SMatrix<2>{{std::exp(1), 2},
				      {2,           1}};
  // test if the obtained gradient is approximated correctly (truncation error is of order O(h^2), h is the step size)
  EXPECT_TRUE((hess(p) - trueHessian).squaredNorm() < std::pow(0.01, 2));
}

// check if ScalarField class handles a discontinuous function correctly
TEST(ScalarFieldTest, DiscontinuousField) {
  // define field expression
  auto fieldExpr = [](SVector<2> x) -> double { // I_{x>0, y>0}
    if(x[0] > 0 && x[1] > 0)
      return 1;
    return 0;
  };
    
  // build the ScalarField object
  ScalarField<2> field(fieldExpr);
  // test if the ScalarField wraps correctly the lambda
  SVector<2> p(1,1);
  double trueResult = 1;
  
  // expect equality
  EXPECT_EQ(field(p), trueResult);

  VectorField<2> grad = field.derive(0.01);
  // define true gradient vector
  SVector<2> trueGradient = SVector<2>(0,0);
  // test if the obtained gradient is approximated correctly (truncation error is of order O(h^2), h is the step size)
  EXPECT_TRUE((grad(p) - trueGradient).squaredNorm() < std::pow(0.01, 2));

  // get the hessian function
  std::function<SMatrix<2>(SVector<2>)> hess = field.deriveTwice(0.01);
  SMatrix<2> trueHessian = SMatrix<2>{{0, 0},
				      {0, 0}};
  // test if the obtained gradient is approximated correctly (truncation error is of order O(h^2), h is the step size)
  EXPECT_TRUE((hess(p) - trueHessian).squaredNorm() < std::pow(0.01, 2));  
}

// check expression template mechanism for ScalarField
TEST(ScalarFieldTest, ExprTemplate) {
  // define two scalar fields sf1 and sf2
  auto fieldExpr1 = [](SVector<2> x) -> double { // x^3 + y
    return std::pow(x[0], 3) + x[1];
  };
  ScalarField<2> sf1(fieldExpr1);
  auto fieldExpr2 = [](SVector<2> x) -> double { // [e^(2x+y)]/x
    return std::exp(2*x[0]+x[1])/x[0];
  };
  ScalarField<2> sf2(fieldExpr2);
  // define evaluation point
  SVector<2> p(1,1);

  // evaluate fields at point p
  double x1 = sf1(p);
  double x2 = sf2(p);
  
  // build various expressions and test for equality
  auto sf3 = sf1 + sf2;
  ASSERT_DOUBLE_EQ(sf3(p), (x1 + x2));
  auto sf4 = sf1 - sf2;
  ASSERT_DOUBLE_EQ(sf4(p), (x1 - x2));
  auto sf5 = sf1 * sf2;
  ASSERT_DOUBLE_EQ(sf5(p), (x1 * x2));
  auto sf6 = sf1 / sf2;
  ASSERT_DOUBLE_EQ(sf6(p), (x1 / x2));
  // combine many previous expressions and scalar values
  auto sf7 = sf1 + sf2/sf3 - sf2*2 + 5;
  ASSERT_DOUBLE_EQ(sf7(p), (x1 + x2/(x1 + x2) - x2*2 + 5));

  // define a discontinuous field and combine with a polynomial field via a summation
  std::function<double(SVector<2>)> poly = [](SVector<2> x) -> double {
    return x[0] + x[1];
  };
  ScalarField<2> polyField(poly);
  // this evaluates 2 in (1,1)
  std::function<double(SVector<2>)> step = [](SVector<2> x) -> double {
    if(x[0] > 0 && x[1] > 0)
      return 1;
    return std::pow(x[0], 2);
  };
  ScalarField<2> stepField(step);
  // this evaluates 1 in (1,1)
  
  auto sf8 = polyField + stepField;
  ASSERT_DOUBLE_EQ(sf8(p), 3);
}

// use case for wrapping the result of a field expression in a valid ScalarField
TEST(ScalarFieldTest, ExprCanBeWrapped) {
  // define a field expression
  auto fieldExpr1 = [](SVector<2> x) -> double { // x^3 + y
    return std::pow(x[0], 3) + x[1];
  };
  ScalarField<2> sf1(fieldExpr1);
  auto fieldExpr2 = [](SVector<2> x) -> double { // [e^(2x+y)]/x
    return std::exp(2*x[0]+x[1])/x[0];
  };
  ScalarField<2> sf2(fieldExpr2);

  auto sf3 = sf1 + sf2/sf1; // x^3 + y + (e^(2x+y))/(x^4 + xy)

  // wrap sf3 in a scalar field
  auto fieldExpr3 = [=](SVector<2> x) -> double {
    return sf3(x);
  };
  ScalarField<2> sf4(fieldExpr3);

  // check expression is wrapped correctly
  SVector<2> p(1,1);
  EXPECT_EQ(sf4(p), sf1(p) + sf2(p)/sf1(p));

  VectorField<2> grad = sf4.derive(0.01);
  // exact gradient evaluated at p
  SVector<2> trueGradient = SVector<2>(3 - 0.25*std::exp(3), 1 + 0.25*std::exp(3));
  EXPECT_TRUE((grad(p) - trueGradient).squaredNorm() < std::pow(0.01, 2));
}

// checks if DEF_FIELD_UNARY_OPERATOR and DEF_FIELD_UNARY_FUNCTOR defines valid operators
TEST(ScalarFieldTest, UnaryExpr) {
  // define a scalar field
  auto fieldExpr1 = [](SVector<2> x) -> double { // x^3 + y
    return std::pow(x[0], 3) + x[1];
  };
  ScalarField<2> sf1(fieldExpr1);
  // define a new scalar field as the application of sin to sf1
  ScalarField<2> sf2 = sin(sf1);
  // define evaluation point
  SVector<2> p(1,1);
  // field sf1 evaluated in (1,1) gives 2, sf2 expects equal to sin(2) ~ 0,03489949670250097165
  EXPECT_DOUBLE_EQ(std::sin(2), sf2(p));

  // apply sin function to a field expression
  auto sf3 = sin(2*sf1 + sf2) + sf1;
  // evaluate it at p. expected value: sin(2*2 + sin(2)) + 2
  EXPECT_DOUBLE_EQ(std::sin(2*2 + std::sin(2)) + 2, sf3(p));
}

// checks if a ScalarField can be built from a field expression
TEST(ScalarFieldTest, FieldExprImplicitConversion) {
  // define scalar field object
  auto fieldExpr = [](SVector<2> x) -> double { // x^3 + y
    return std::pow(x[0], 3) + x[1];
  };
  ScalarField<2> sf1(fieldExpr);
  // define field expression
  auto expr = 2*log(sf1) + sf1;
  // convert expr to a scalar field
  ScalarField<2> sf2(expr);

  SVector<2> p(1,1);
  EXPECT_EQ(expr(p), sf2(p));
}

// checks DifferentiableScalarField calls analytical gradient on .derive() call
TEST(ScalarFieldTest, DifferentiableField) {
  // define analytical gradient expression
  std::function<double(SVector<2>)> baseField = [](SVector<2> x) -> double { // x^3 + y
    return std::pow(x[0], 3) + x[1];
  };
  std::function<SVector<2>(SVector<2>)> grad = [](SVector<2> x) -> SVector<2> {
    return SVector<2>(3*std::pow(x[0], 2), 1);
  };
  // wrap all in a differentiable scalar field
  ScalarField sf(baseField);

  // let DifferentiableScalarField to build the gradient VectorField (inefficient but more abstract)
  DifferentiableScalarField df(baseField, grad);
  // or build explicitly the gradient field using list initialization
  std::function<double(SVector<2>)> dx = [](SVector<2> x) -> double { return 3*std::pow(x[0], 2); };
  std::function<double(SVector<2>)> dy = [](SVector<2> x) -> double { return 1; };
  DifferentiableScalarField df_(baseField, {dx, dy});
    
  // define evaluation point
  SVector<2> p(1,1);

  // check both initializations define the same object
  EXPECT_TRUE(df(p) == df_(p));
  
  SVector<2> trueGradient = SVector<2>(3,1);
  // check derive() calls exact gradient
  ASSERT_FALSE(sf.derive()(p) == df.derive()(p));
  // check exact derivative is actually called
  ASSERT_TRUE(df.derive()(p) == trueGradient);
}

// checks TwiceDifferentiableScalarField calls analytical hessian on .deriveTwice() call
TEST(ScalarFieldTest, TwiceDifferentiableField) {
  // define analytical gradient and hessian expression
  std::function<double(SVector<2>)> baseField = [](SVector<2> x) -> double { // x^3 + y
    return std::pow(x[0], 3) + x[1];
  };
  std::function<SVector<2>(SVector<2>)> grad  = [](SVector<2> x) -> SVector<2> {
    return SVector<2>(3*std::pow(x[0], 2), 1);
  };
  std::function<SMatrix<2>(SVector<2>)> hess  = [](SVector<2> x) -> SMatrix<2> {
    return SMatrix<2>{{6*x[0], 0},
		      {0,      0}};
  };
  
  // wrap all in a twice differentiable scalar field
  ScalarField<2> sf(baseField);
  TwiceDifferentiableScalarField<2> tdf(baseField, grad, hess);
  // define evaluation point
  SVector<2> p(1,1);
  SMatrix<2> trueHessian = SMatrix<2>{{6,0},
				      {0,0}};
  // check deriveTwice() calls exact hessian
  ASSERT_FALSE(sf.deriveTwice()(p) == tdf.deriveTwice()(p));
  // check exact hessian is actually called
  ASSERT_TRUE(tdf.deriveTwice()(p) == trueHessian);
}

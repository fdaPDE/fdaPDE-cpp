#include <array>
#include <cmath>
#include <functional>
#include <gtest/gtest.h> // testing framework

#include "../fdaPDE/core/utils/fields/ScalarField.h"
using fdaPDE::core::ScalarField;
#include "../fdaPDE/core/utils/Symbols.h"
#include "../fdaPDE/core/utils/fields/VectorField.h"
using fdaPDE::core::VectorField;
#include "../fdaPDE/core/utils/fields/DotProduct.h"
using fdaPDE::core::DotProduct;
#include "../fdaPDE/core/utils/fields/Divergence.h"
using fdaPDE::core::Divergence;

// test different constructors of VectorField
TEST(VectorFieldTest, VectorFieldWrapsCorrectly) {
  // define vector field using a vectorial lambda expression
  std::function<SVector<2>(SVector<2>)> exprField = [](SVector<2> x) -> SVector<2> {
    return SVector<2>(std::exp(x[0]*x[1]) + 2, (std::pow(x[0], 2) + x[1])/x[1]); // e^{xy} + 2; (x^2 + y)/y
  };
  // wrap it in a VectorField
  VectorField<2> field1(exprField);
  // define evaluation point
  SVector<2> p(1,1);
  // field evaluates in p: (e+2, 2)
  SVector<2> trueEvaluation(std::exp(1)+2, 2);
  for(std::size_t i = 0; i < 2; ++i) EXPECT_DOUBLE_EQ(field1(p)[i], trueEvaluation[i]);
  
  // define vector field using list initialization
  std::function<double(SVector<2>)> x_comp = [](SVector<2> x) -> double {
    return std::exp(x[0]*x[1]) + 2;
  };
  std::function<double(SVector<2>)> y_comp = [](SVector<2> x) -> double {
    return (std::pow(x[0], 2) + x[1])/x[1];
  };
  VectorField<2> field2({x_comp, y_comp});
  for(std::size_t i = 0; i < 2; ++i) EXPECT_DOUBLE_EQ(field2(p)[i], trueEvaluation[i]);
  
  // define vector field explicitly declaring an array of lambdas
  std::array<std::function<double(SVector<2>)>, 2> comp_array = {x_comp, y_comp};
  VectorField<2> field3(comp_array);
  for(std::size_t i = 0; i < 2; ++i) EXPECT_DOUBLE_EQ(field3(p)[i], trueEvaluation[i]);

  // if all asserts are verified by transitivity all the definitions give origin to the same object
}

// check if field obtained by subscript is the same as the one supplied in input
TEST(VectorFieldTest, ConstSubscript) {
  // define vector field
  typedef std::function<double(SVector<2>)> field_component;
  field_component x_comp = [](SVector<2> x) -> double { return x[0] + 2; };
  field_component y_comp = [](SVector<2> x) -> double { return x[0]*x[1]/4; };
  // wrap it in a VectorField
  VectorField<2> field1({x_comp, y_comp});

  // extract first dimension using const subscript
  ScalarField<2> x_extracted = field1[0];
  // evaluation point
  SVector<2> p(1,1);
  EXPECT_EQ(x_comp(p), x_extracted(p));
}

// check if a vector field can be assigned using the subscript operator and any FieldExpr on the rhs
TEST(VectorFieldTest, AssingBySubscript) {
  // define an empty vector field
  VectorField<2> v_field;
  // define a lambda expression valid to be wrapped by a ScalarField
  std::function<double(SVector<2>)> scalarField = [](SVector<2> x) -> double {
    return std::log(x[0] + x[1]) - std::pow(x[0],2)*x[1] + 5;
  };
  ScalarField<2> s_field(scalarField);
  // set the field coordinates
  v_field[0] = scalarField;   // converting constructor ScalarField(std::function<double(SVector<2>)>) called
  v_field[1] = 2*s_field + 1; // built directly from a field expression

  // evaluation point
  SVector<2> p(1,1);
  // evaluation of field at p: (ln(2) + 4, 2ln(2) + 9)
  SVector<2> trueEvaluation(std::log(2) + 4, 2*std::log(2) + 9);
  for(std::size_t i = 0; i < 2; ++i){
    EXPECT_DOUBLE_EQ(v_field(p)[i], trueEvaluation[i]);
  }
}

// checks if the inner product between two fields is correct
TEST(VectorFieldTest, InnerProductVField) {
  // define two vector fields
  VectorField<2> v_field1;
  v_field1[0] = [](SVector<2> x) -> double { return std::pow(x[0], 2) + 1; };
  v_field1[1] = [](SVector<2> x) -> double { return std::exp(x[0]*x[1]);   };

  VectorField<2> v_field2;
  v_field2[0] = v_field1[0] + 2;
  v_field2[1] = 2*v_field1[0]*v_field1[1];

  // build a vector field as the dot product of the previous two fields
  DotProduct dotProduct = v_field1.dot(v_field2);
  // evaluation point
  SVector<2> p(1,1);
  // dot product functor evaluates same as dot product of evaluated fields
  EXPECT_DOUBLE_EQ(dotProduct(p), v_field1(p).dot(v_field2(p)));

  // dot product evaluates correctly
  // dotProduct encodes the scalar field of equation: (x^2 + 1)*(x*2 + 3) + 2*e^{xy}*(x^2 + 1)*e^{xy}
  // dotProduct evaluates at p: (2)*(4) + 2*e^{2}*2 = 8 + 4*e^{2}
  EXPECT_DOUBLE_EQ(dotProduct(p), 8 + 4*std::exp(2));
}

// checks if the inner product between a vector field and an SVector is correct
TEST(VectorFieldTest, InnerProductSVector) {
  // vector field definition
  VectorField<2> field;
  field[0] = [](SVector<2> x) -> double { return std::pow(x[0], 2) + 1; };
  field[1] = [](SVector<2> x) -> double { return std::exp(x[0]*x[1]);   };
  // define an SVector
  SVector<2> coeff(5,2);
  // perform dot product
  ScalarField<2> dotProduct = field.dot(coeff); // 5*(x^2 + 1) + 2*e^{xy}
  // evaluation point
  SVector<2> p(1,1);
  EXPECT_DOUBLE_EQ(10 + 2*std::exp(1), dotProduct(p));
}

// checks if the matrix * VectorField product returns a valid VectorField
TEST(VectorFieldTest, MatrixProduct) {
  // vector field definition
  std::function<SVector<2>(SVector<2>)> exprField = [](SVector<2> x) -> SVector<2> {
    return SVector<2>(std::exp(x[0]*x[1]) + 2, (std::pow(x[0], 2) + x[1])/x[1]); // e^{xy} + 2; (x^2 + y)/y
  };
  VectorField<2> field(exprField);

  // define a coefficient matrix M
  SMatrix<2> M;
  M << 1, 2, 3, 4;
  // obtain product field M*field
  VectorField<2> productField = M*field; // (e^{xy} + 2 + 2*(x^2 + y)/y; 3*e^{xy} + 6 + 4*(x^2 + y)/y)
  // evaluation point
  SVector<2> p(1,1);
  SVector<2> trueEvaluation(std::exp(1) + 6, 3*std::exp(1) + 14);
  for(std::size_t i = 0; i < 2; ++i){
    EXPECT_DOUBLE_EQ(productField(p)[i], trueEvaluation[i]);
  }
}

// check expression template mechanism for VectorField
TEST(VectorFieldTest, FieldExpr) {
  // create vector field from assignment operator
  VectorField<2> vf1;
  vf1[0] = [](SVector<2> x) -> double { return std::pow(x[0], 2) + 1; }; // x^2 + 1
  vf1[1] = [](SVector<2> x) -> double { return std::exp(x[0]*x[1]);   }; // e^{xy}

  // define vector field using a vectorial lambda expression
  std::function<SVector<2>(SVector<2>)> exprField = [](SVector<2> x) -> SVector<2> {
    return SVector<2>(std::exp(x[0]*x[1]) + 2, (std::pow(x[0], 2) + x[1])/x[1]); // e^{xy} + 2; (x^2 + y)/y
  };
  VectorField<2> vf2(exprField);

  // evaluation point
  SVector<2> p(1,1); 
  SVector<2> vf1_eval(2, std::exp(1));     // vf1 -> (2, e)
  SVector<2> vf2_eval(std::exp(1) + 2, 2); // vf2 -> (e + 2, 2)
  
  // build various expressions and test for equality
  auto vf3 = vf1 + vf2;
  for(std::size_t i = 0; i < 2; ++i)
    EXPECT_DOUBLE_EQ(vf3(p)[i], (vf1_eval + vf2_eval)[i]);
  auto vf4 = vf1 - vf2; 
  for(std::size_t i = 0; i < 2; ++i)
    EXPECT_DOUBLE_EQ(vf4(p)[i], (vf1_eval - vf2_eval)[i]);

  SMatrix<2> M;
  M << 1, 2, 3, 4;
  auto vf5 = M*vf1 - vf2 + vf3;
  for(std::size_t i = 0; i < 2; ++i)
    EXPECT_DOUBLE_EQ(vf5(p)[i], (M*vf1_eval + vf1_eval)[i]);

  // the following are ScalarField expressions coming from VectorField operations
  auto vf6 = vf1.dot(vf2) + vf3[0];
  EXPECT_DOUBLE_EQ(vf6(p), vf1_eval[0]*vf2_eval[0] + vf1_eval[1]*vf2_eval[1] + (vf1_eval[0] + vf2_eval[0]));
}

// check a vectorial expression can be wrapped in a VectorField
TEST(VextorFieldTest, ExprToField) {
  // define some field expression
  VectorField<2> vf1;
  vf1[0] = [](SVector<2> x) -> double { return std::pow(x[0], 2) + 1; }; // x^2 + 1
  vf1[1] = [](SVector<2> x) -> double { return std::exp(x[0]*x[1]);   }; // e^{xy}

  std::function<SVector<2>(SVector<2>)> exprField = [](SVector<2> x) -> SVector<2> {
    return SVector<2>(std::exp(x[0]*x[1]) + 2, (std::pow(x[0], 2) + x[1])/x[1]); // e^{xy} + 2; (x^2 + y)/y
  };
  VectorField<2> vf2(exprField);
  auto vf3 = vf1 + vf2;
  // convert the expression in a valid vector field
  VectorField<2> vf4(vf3);
  // evaluation point
  SVector<2> p(1,1); 
  SVector<2> vf1_eval(2, std::exp(1));     // vf1 -> (2, e)
  SVector<2> vf2_eval(std::exp(1) + 2, 2); // vf2 -> (e + 2, 2)

  for(std::size_t i = 0; i < 2; ++i)
    EXPECT_DOUBLE_EQ(vf4(p)[i], (vf1_eval + vf2_eval)[i]);  
}

// check divergence of VectorField
TEST(VectorFieldTest, Divergence) {
  // define VectorField
  std::function<SVector<2>(SVector<2>)> fieldExpr = [](SVector<2> x) -> SVector<2> {
    return SVector<2>(std::pow(x[0], 3)*x[1], std::log(x[1])); // (x^3*y, ln(y))
  };
  VectorField<2> field(fieldExpr);
  // get divergence operator of the field
  Divergence fieldDivergece(field);
  // evaluation point
  SVector<2> p(1,1);
  // divergence evaluates in p: 3*x^2*y + 1/y -> 3 + 1 = 4
  // precision is bounded by the step used in the central difference formula for the approximation of derivatives
  // moreover approximation errors accumulate in the summation dFx/dx + dFy/dy defining the divergence operator (this is the reason
  // of the 2* in the error bound used in the next assertion)
  EXPECT_TRUE(std::abs(fieldDivergece(p) - 4) < 2*std::pow(0.001, 2));

  // divergence can be used in a field expression
  auto sf1 = fieldDivergece + 2*field[0];
  EXPECT_TRUE(std::abs(sf1(p) - (4 + 2)) < 2*std::pow(0.001, 2));
}

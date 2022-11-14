#include <cstddef>
#include <gtest/gtest-typed-test.h>
#include <gtest/gtest.h> // testing framework

#include "../../fdaPDE/core/utils/Symbols.h"
#include "../../fdaPDE/core/MESH/Element.h"
using fdaPDE::core::MESH::Element;
#include "../../fdaPDE/core/NLA/VectorSpace.h"
using fdaPDE::core::NLA::VectorSpace;

#include "../utils/MeshLoader.h"
using fdaPDE::testing::MESH_TYPE_LIST;
using fdaPDE::testing::MeshLoader;
#include "../utils/Constants.h"
using fdaPDE::testing::DOUBLE_TOLERANCE;
using fdaPDE::testing::MACHINE_EPSILON;

#include "../../fdaPDE/core/utils/CompileTime.h"
#include "../../fdaPDE/core/FEM/integration/Integrator.h"
using fdaPDE::core::FEM::Integrator;
#include "../../fdaPDE/core/FEM/basis/LagrangianBasis.h"
using fdaPDE::core::FEM::LagrangianBasis;
#include "../../fdaPDE/core/FEM/operators/Identity.h"
#include "../../fdaPDE/core/FEM/operators/Gradient.h"
#include "../../fdaPDE/core/FEM/operators/Laplacian.h"
using fdaPDE::core::FEM::Identity;
using fdaPDE::core::FEM::Gradient;
using fdaPDE::core::FEM::Laplacian;

// test fixture
template <typename E>
struct IntegratorTest : public ::testing::Test {
  MeshLoader<E> meshLoader{}; // use default mesh
  static constexpr unsigned int M = MeshLoader<E>::M;
  static constexpr unsigned int N = MeshLoader<E>::N;
};
TYPED_TEST_SUITE(IntegratorTest, MESH_TYPE_LIST);

// tests if the integration of the constant field 1 over an element equals its measure
TYPED_TEST(IntegratorTest, ElementMeasure){
  // generate random element from mesh
  auto e = this->meshLoader.generateRandomElement();
  Integrator<TestFixture::M, 1> integrator; // define integrator
  // the integral of the constant field 1 over the mesh element equals its measure
  std::function<double(SVector<TestFixture::N>)> f = [](SVector<TestFixture::N> x) -> double { return 1; };
  EXPECT_NEAR(e->measure(), integrator.integrate(*e, f), DOUBLE_TOLERANCE);
}

// test if linear fields can be integrated over mesh elements. In particular a closed formula for
// the volume of a truncated prism defined over an element e having height h1, h2, ..., hm at the m vertices is known as
//     e.measure()*(h1 + h2 + ... hm)/m
TYPED_TEST(IntegratorTest, LinearFieldsAreIntegratedCorrectly){
    // generate random element from mesh
    auto e = this->meshLoader.generateRandomElement();
    Integrator<TestFixture::M, 1> integrator; // define integrator
    // a linear function over an element e defines a truncated prism over e
    std::function<double(SVector<TestFixture::N>)> f = [](SVector<TestFixture::N> x) -> double {
      return x[0] + x[1];
    };
    // compute volume of truncated rectangular prism: 1/(M+1)*V*(h1 + h2 + ... + hM), where V is the element's measure
    double h = 0;
    for(auto p : *e)
      h += f(p);
    double measure = e->measure()*h/(TestFixture::M+1);
    // test for equality
    EXPECT_NEAR(measure, integrator.integrate(*e, f), DOUBLE_TOLERANCE);
}

// test if is possible to integrate a field over the entire mesh
TEST(IntegratorTest, CanIntegrateFieldOverMesh) {
  // load sample mesh
  MeshLoader<Mesh2D<>> CShaped("unit_square");
  Integrator<2,1> integrator{};
  // define field to integrate
  std::function<double(SVector<2>)> f = [](SVector<2> x) -> double { return 1; };
  EXPECT_NEAR(1, integrator.integrate(CShaped.mesh, f), DOUBLE_TOLERANCE);
}

// test if all integrator tables produce the same result (this proves weights and quadrature nodes are correct)
template <typename E>
class IntegratorTablesTest : public ::testing::Test {};




// TYPED_TEST(IntegratorTest, IdentityOperator) {
//   // generate random element from mesh
//   auto e = this->generateRandomElement();
//   Integrator<TestFixture::M> integrator;
//   // define basis over the element
//   LagrangianBasis<TestFixture::M, TestFixture::N, 1> elementBasis(*e);
//   // define identity operator
//   auto form = Identity();

//   // wrap the M-dimesional function into an N-dimensional one
//   std::function<double(SVector<TestFixture::N>)> f = [elementBasis, e](SVector<TestFixture::N> x) -> double {
//     // project point on space spanned by the element
//     SVector<TestFixture::M> p = e->spannedSpace().projectOnto(x);
//     return (elementBasis[0]*elementBasis[0])(p);
//   };
//   // differential operators assume to be integrated over a basis defined on the reference element.
//   LagrangianBasis<TestFixture::M, TestFixture::N, 1> referenceBasis;
//   EXPECT_NEAR(integrator.integrate(referenceBasis, *e, 0, 0, form), integrator.integrate(*e, f), this->tolerance);
// }

// TYPED_TEST(IntegratorTest, LaplacianOperator) {
//   // generate random element from mesh
//   auto e = this->generateRandomElement();
//   Integrator<TestFixture::M> integrator;
//   // define basis over the element
//   LagrangianBasis<TestFixture::M, TestFixture::N, 1> elementBasis(*e);
//   // define laplacian operator
//   auto form = Laplacian();
  
//   // wrap the M-dimesional function into an N-dimensional one
//   std::function<double(SVector<TestFixture::N>)> f = [elementBasis, e](SVector<TestFixture::N> x) -> double {
//     // project point on space spanned by the element
//     SVector<TestFixture::M> p = e->spannedSpace().projectOnto(x);
//     return (elementBasis[1].derive().dot(elementBasis[0].derive()))(p);
//   };
//   // differential operators assume to be integrated over a basis defined on the reference element.
//   LagrangianBasis<TestFixture::M, TestFixture::N, 1> referenceBasis{};
//   EXPECT_NEAR(integrator.integrate(referenceBasis, *e, 1, 0, form), integrator.integrate(*e, f), this->tolerance);
// }

// TYPED_TEST(IntegratorTest, GradientOperator) {
//   if constexpr(TestFixture::M == TestFixture::N){ // tested only for non-manifold
//     // generate random element from mesh
//     auto e = this->generateRandomElement();
//     Integrator<TestFixture::M> integrator;
//     // define basis over the element
//     LagrangianBasis<TestFixture::M, TestFixture::N, 1> elementBasis(*e);
//     // define gradient operator
//     SVector<TestFixture::N> b = SVector<TestFixture::N>::Ones();
//     auto form = dot(b, Gradient());

//     // wrap the M-dimesional function into an N-dimensional one (this is useless in the non manifold case)
//     std::function<double(SVector<TestFixture::N>)> f = [elementBasis, e](SVector<TestFixture::N> x) -> double {
//       // project point on space spanned by the element
//       SVector<TestFixture::M> p = e->spannedSpace().projectOnto(x);
//       return (elementBasis[0]*elementBasis[0].derive().dot(SVector<TestFixture::M>::Ones()))(p);
//     };
    
//     // differential operators assume to be integrated over a basis defined on the reference element.
//     LagrangianBasis<TestFixture::M, TestFixture::N, 1> referenceBasis{};
//     EXPECT_NEAR(integrator.integrate(referenceBasis, *e, 0, 0, form), integrator.integrate(*e, f), this->tolerance);
//   }else{
//     // skip test for manifold case: NEED FIX!
//     GTEST_SKIP();
//   }
// }

#include <gtest/gtest-typed-test.h>
#include <gtest/gtest.h> // testing framework
#include <limits>
#include <string>
#include <type_traits>
#include <random>
#include <Eigen/Dense>

#include "../fdaPDE/core/MESH/Element.h"
#include "core/utils/Symbols.h"
using fdaPDE::core::MESH::Element;
#include "../fdaPDE/core/MESH/Mesh.h"
using fdaPDE::core::MESH::Mesh2D;
using fdaPDE::core::MESH::Mesh3D;
using fdaPDE::core::MESH::NetworkMesh;
using fdaPDE::core::MESH::SurfaceMesh;
#include "../fdaPDE/core/utils/CompileTime.h"
#include "../fdaPDE/core/FEM/integration/Integrator.h"

template <typename E>
class IntegratorTest : public ::testing::Test {
private:
  // RNG stuffs
  std::default_random_engine rng{};

public:
  E m; // mesh
  static constexpr unsigned int M = E::local_dimension;
  static constexpr unsigned int N = E::embedding_dimension;
  static constexpr std::size_t n_basis = ct_binomial_coefficient(M+1,1);
  double tolerance = std::pow(0.1, 15);
  
  // load mesh from .csv files
  IntegratorTest() {
    std::string dim = (E::manifold) ? std::to_string(M) + ".5" : std::to_string(M);
    // compute file names
    std::string point    = "data/m" + dim + "D_points.csv";
    std::string edges    = "data/m" + dim + "D_edges.csv";
    std::string elements = "data/m" + dim + "D_elements.csv";
    std::string neigh    = "data/m" + dim + "D_neigh.csv";
    std::string boundary = "data/m" + dim + "D_boundary.csv";
    // initialize test objects
    m = E(point, edges, elements, neigh, boundary);
    std::random_device rd;
    rng = std::default_random_engine(rd());
  };

  // generate element at random inside mesh m
  std::shared_ptr<Element<E::local_dimension, E::embedding_dimension>> generateRandomElement();
};

template <typename E>
std::shared_ptr<Element<E::local_dimension, E::embedding_dimension>>
IntegratorTest<E>::generateRandomElement() {
  std::uniform_int_distribution<int> randomID(0, this->m.elements()-1);
  int ID = randomID(rng);  
  return m.element(ID);
}

using MeshList = ::testing::Types<Mesh2D<>, SurfaceMesh<>, Mesh3D<>, NetworkMesh<>>;
TYPED_TEST_SUITE(IntegratorTest, MeshList);

TYPED_TEST(IntegratorTest, ElementMeasure){
  // generate random element from mesh
  auto e = this->generateRandomElement();
  Integrator<TestFixture::N> integrator; // define integrator
  // the integral of the constant field 1 over the mesh element equals its measure
  std::function<double(SVector<TestFixture::N>)> f = [](SVector<TestFixture::N> x) -> double { return 1; };
  EXPECT_NEAR(e->measure(), integrator.integrate(*e, f), this->tolerance);
}

TYPED_TEST(IntegratorTest, LinearField){
  if constexpr(TestFixture::N == TestFixture::M){
    // generate random element from mesh
    auto e = this->generateRandomElement();
    Integrator<TestFixture::N> integrator; // define integrator
    std::function<double(SVector<TestFixture::N>)> f = [](SVector<TestFixture::N> x) -> double {
      return x[0] + x[1];
    };
    // compute volume of truncated rectangular prism: 1/(M+1)*V*(h1 + h2 + ... + hM)
    double h = 0;
    for(auto p : *e) h += f(p);
    double measure = e->measure()*h/(TestFixture::M+1);
    EXPECT_NEAR(measure, integrator.integrate(*e, f), this->tolerance);
  }
  // if constexpr(TestFixture::N == 2 && TestFixture::M == 1){
  //   // generate random element from mesh
  //   auto e = this->generateRandomElement();
  //   Integrator<TestFixture::N> integrator; // define integrator
  //   std::function<double(SVector<TestFixture::N>)> f = [](SVector<TestFixture::N> x) -> double {
  //     return x[0] + x[1];
  //   };
  //   // compute volume of trapezoid
  //   double h = 0;
  //   for(auto p : *e){ h += f(p); std::cout << f(p) << "\n----\n" << p << "\n....." << std::endl; }
  //   double measure = e->measure()*h/2;
  //   EXPECT_NEAR(measure, integrator.integrate(*e, f), this->tolerance);
  // }

}

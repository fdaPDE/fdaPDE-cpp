#ifndef __ASSEMBLER_H__
#define __ASSEMBLER_H__

#include <Eigen/Sparse>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <Eigen/src/SparseCore/SparseUtil.h>
#include <cstddef>
#include <limits>
#include <memory>
#include <ostream>

#include "../MESH/Mesh.h"
#include "../MESH/Element.h"
#include "Integrator.h"
#include "FiniteElement.h"
#include "MultivariatePolynomial.h"
#include "../utils/fields/VectorField.h"
#include "../utils/fields/ScalarField.h"

using fdaPDE::core::VectorField;
using fdaPDE::core::ScalarField;
using fdaPDE::core::MESH::Element;
using fdaPDE::core::MESH::Mesh;

template <unsigned int M, unsigned int N, unsigned int ORDER>
class Assembler {

private:
  Mesh<M, N>& mesh_;
  constexpr static unsigned n_basis = ct_binomial_coefficient(N+ORDER, ORDER);
  
public:
  Assembler(Mesh<M, N>& mesh) : mesh_(mesh) {};

  // assemble stiffness matrix
  Eigen::SparseMatrix<double> assemble();

  // assemble forcing vector
  Eigen::Matrix<double, Eigen::Dynamic, 1> forcingTerm(const ScalarField<N>& f);
};

// stiff matrix for laplacian operator
template <unsigned int M, unsigned int N, unsigned int ORDER>
Eigen::SparseMatrix<double> Assembler<M,N,ORDER>::assemble() {

  std::vector<Eigen::Triplet<double>> tripletList;  // store triplets (node_i, node_j, integral_value)
  Eigen::SparseMatrix<double> stiffnessMatrix;      // stiffness matrix is sparse due to the local support of basis functions
  Integrator<2,6> integrator;                       // quadrature rule to approximate integrals
  ReferenceElement<N, ORDER> fe;                    // reference element where functional information is defined

  // properly allocate memory to avoid reallocations
  tripletList.reserve(n_basis*n_basis*mesh_.getNumberOfElements());
  stiffnessMatrix.resize(mesh_.getNumberOfNodes(), mesh_.getNumberOfNodes());

  // cycle over all mesh elements
  for(std::shared_ptr<Element<M,N>> e : mesh_){
    // consider all pair of nodes
    for(size_t i = 0; i < n_basis; ++i){
      for(size_t j = 0; j < n_basis; ++j){
	// get basis elements
	VectorField<N> NablaPhi_i = fe.getBasis().getBasisElement(i).gradient();
	VectorField<N> NablaPhi_j = fe.getBasis().getBasisElement(j).gradient();

	// callable object ready to be integrated
	auto bilinear_form = NablaPhi_i.dot(NablaPhi_j); // for laplacian: \nabla phi_i * \nabla * phi_j
	
	// integrate the functional over the reference element
	double value = integrator.integrate(*e, bilinear_form);
	
	// From Eigen doucmentation: A triplet is a tuple (i,j,value) defining a non-zero element. The input list of triplets
	// does not have to be sorted, and can contains duplicated elements. In any case, the result is a sorted and compressed
	// sparse matrix where the duplicates have been summed up.
	// Linearity of the integral is implicitly used during matrix construction! (there is no explicit sum in this function)
	tripletList.emplace_back(e->getFESupport()[i].first, e->getFESupport()[j].first, value);
      }
    }
  }
  
  // stiff matrix assembled
  stiffnessMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
  stiffnessMatrix.prune(std::numeric_limits<double>::epsilon() * 10); // remove almost zero entries

  return stiffnessMatrix;
};

template <unsigned int M, unsigned int N, unsigned int ORDER>
Eigen::Matrix<double, Eigen::Dynamic, 1> Assembler<M,N,ORDER>::forcingTerm(const ScalarField<N>& f) {

  Eigen::Matrix<double, Eigen::Dynamic, 1> result{};
  result.resize(mesh_.getNumberOfNodes(), 1); // there are as many basis functions as number of nodes in the mesh
  result.fill(0);                            // init result vector to zero
  
  Integrator<2,6> integrator;

  // define the reference element
  ReferenceElement<N, ORDER> fe;

  for(std::shared_ptr<Element<M,N>> e : mesh_){    
    for(size_t i = 0; i < n_basis; ++i){
	auto phi_i = fe.getBasis().getBasisElement(i);	// basis function
	auto functional = f*phi_i;	                // functional to integrate

	// perform integration
	double value = integrator.integrate(*e, functional);
	
	// store result exploiting additiviy of the integral
        result[e->getFESupport()[i].first] += value;
    }
  }
  return result;
}


#endif // __ASSEMBLER_H__

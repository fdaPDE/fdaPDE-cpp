#ifndef __ASSEMBLER_H__
#define __ASSEMBLER_H__

#include <Eigen/Sparse>
#include <cstddef>
#include <limits>
#include <memory>

#include "../utils/fields/VectorField.h"
#include "../utils/fields/ScalarField.h"
#include "../MESH/Mesh.h"
#include "../MESH/Element.h"

using fdaPDE::core::VectorField;
using fdaPDE::core::ScalarField;
using fdaPDE::core::MESH::Element;
using fdaPDE::core::MESH::Mesh;

// FEM module includes
#include "integration/Integrator.h"
#include "FunctionalBasis.h"
#include "MultivariatePolynomial.h"

template <unsigned int M, unsigned int N, unsigned int ORDER>
class Assembler {
private:
  constexpr static unsigned n_basis = ct_binomial_coefficient(N+ORDER, ORDER);
  Mesh<M, N>& mesh_;                                // mesh
  Integrator<2,6> integrator{};                     // quadrature rule to approximate integrals
  ReferenceBasis<M, N, ORDER> referenceBasis{};     // functional basis over reference N-dimensional unit simplex
  
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

  // properly preallocate memory to avoid reallocations
  tripletList.reserve(n_basis*n_basis*mesh_.getNumberOfElements());
  stiffnessMatrix.resize(mesh_.getNumberOfNodes(), mesh_.getNumberOfNodes());

  // cycle over all mesh elements
  for(std::shared_ptr<Element<M,N>> e : mesh_){
    // consider all pair of nodes
    for(size_t i = 0; i < n_basis; ++i){
      for(size_t j = 0; j < n_basis; ++j){
	// get a pair of basis functions over the reference element
	VectorField<N> NablaPhi_i = referenceBasis[i].gradient();
	VectorField<N> NablaPhi_j = referenceBasis[j].gradient();

	// for laplacian: \nabla phi_i * \nabla * phi_j
        InnerProduct<N> bilinear_form = NablaPhi_i.dot(NablaPhi_j);
	double value = integrator.integrate(*e, bilinear_form); // integrate

	// From Eigen doucmentation: The input list of triplets does not have to be sorted, and can contains duplicated elements.
	// In any case, the result is a sorted and compressed sparse matrix where the duplicates have been summed up.
	// Linearity of the integral is implicitly used during matrix construction by eigen!
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
  result.fill(0);                             // init result vector to zero

  // build forcing vector
  for(std::shared_ptr<Element<M,N>> e : mesh_){

    // build functional basis over the current element e
    FunctionalBasis<M, N, ORDER> basis(*e);
    
    // integrate on each node
    for(size_t i = 0; i < n_basis; ++i){      
	auto phi_i = basis[i];	      // basis function
	auto functional = f*phi_i;    // functional to integrate
	
	// perform integration and store result exploiting additiviy of the integral
        result[e->getFESupport()[i].first] += integrator.integrate(*e, functional);
    }
  }
  return result;
}

#endif // __ASSEMBLER_H__

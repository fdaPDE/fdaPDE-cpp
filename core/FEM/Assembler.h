#ifndef __ASSEMBLER_H__
#define __ASSEMBLER_H__

#include <Eigen/Sparse>
#include <Eigen/src/Core/util/Constants.h>
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
#include "operators/BilinearFormExpressions.h"

template <unsigned int M, unsigned int N, unsigned int ORDER>
class Assembler {
private:
  constexpr static unsigned n_basis = ct_binomial_coefficient(N+ORDER, ORDER);
  Mesh<M, N>& mesh_;                                // mesh
  Integrator<2U,6U> integrator{};                   // quadrature rule to approximate integrals
  ReferenceBasis<M, N, ORDER> referenceBasis{};     // functional basis over reference N-dimensional unit simplex
  
public:
  Assembler(Mesh<M, N>& mesh) : mesh_(mesh) {};

  // assemble stiffness matrix
  template <typename E>
  Eigen::SparseMatrix<double> assemble(const E& bilinearForm);
  // assemble forcing vector
  Eigen::Matrix<double, Eigen::Dynamic, 1> forcingTerm(const Eigen::Matrix<double, Eigen::Dynamic, 1>& f);
};

// stiff matrix for laplacian operator
template <unsigned int M, unsigned int N, unsigned int ORDER>
template <typename E>
Eigen::SparseMatrix<double> Assembler<M,N,ORDER>::assemble(const E& bilinearForm) {

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
	// any integral computation for the construction of the stiffness matrix is performed on the reference element
	double value = integrator.integrate(referenceBasis, *e, i, j, bilinearForm);
	
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

  // impose homogeneous boundary condition to remove not necessary degrees of freedom
  // (otherwise the corresponding linear system is undetermined system!)
  for(size_t i = 0; i < stiffnessMatrix.rows(); ++i){
    if(mesh_.isOnBoundary(i)){ 
      stiffnessMatrix.row(i) *= 0;              // zero all entries of this row
      stiffnessMatrix.coeffRef(i,i) = 1;        // set diagonal element to 1 to impose equation u_j = b_j
    }
  }

  // return just half of the discretization matrix if the form is symmetric (lower triangular part)
  if (bilinearForm.isSymmetric())
    return stiffnessMatrix.selfadjointView<Eigen::Lower>();
  else
    return stiffnessMatrix;
};

template <unsigned int M, unsigned int N, unsigned int ORDER>
Eigen::Matrix<double, Eigen::Dynamic, 1> Assembler<M,N,ORDER>::forcingTerm(const Eigen::Matrix<double, Eigen::Dynamic, 1>& f) {

  Eigen::Matrix<double, Eigen::Dynamic, 1> result{};
  result.resize(mesh_.getNumberOfNodes(), 1); // there are as many basis functions as number of nodes in the mesh
  result.fill(0);                             // init result vector to zero

  // build forcing vector
  for(std::shared_ptr<Element<M,N>> e : mesh_){

    // build functional basis over the current element e
    FunctionalBasis<M, N, ORDER> basis(*e);

    // integrate on each node
    for(size_t i = 0; i < n_basis; ++i){
      if(e->getBoundaryMarkers()[i].second == 0){ // skip computation if node is a boundary node
	auto phi_i = basis[i];	                  // basis function
	// f[e->getID()] is the value of the discretized forcing field (given as datum) over the current element
	auto functional = f[e->getID()]*phi_i;    // functional to integrate
	
	// perform integration and store result exploiting additiviy of the integral
        result[e->getFESupport()[i].first] += integrator.integrate(*e, functional);
      }
      else result[e->getFESupport()[i].first] += 0;
    }
  }
  return result;
}

#endif // __ASSEMBLER_H__

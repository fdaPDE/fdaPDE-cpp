#ifndef __ASSEMBLER_H__
#define __ASSEMBLER_H__

#include <Eigen/Sparse>
#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <Eigen/src/SparseCore/SparseUtil.h>
#include <cstddef>
#include <memory>

#include "../MESH/Mesh.h"
#include "../MESH/Element.h"
#include "Integrator.h"
#include "IntegratorTables.h"
#include "FiniteElement.h"
#include "MultivariatePolynomial.h"

using fdaPDE::core::MESH::Element;
using fdaPDE::core::MESH::Mesh;

template <unsigned int M, unsigned int N, unsigned int ORDER>
class Assembler {

public:
  void assemble();
  
};

// stiff matrix for laplacian operator
template <unsigned int M, unsigned int N, unsigned int ORDER>
void Assembler<M,N,ORDER>::assemble() {
  std::vector<Eigen::Triplet<double>> tripletList;
  Eigen::SparseMatrix<double> m;

  Mesh<M,N> mesh;
  Integrator<IntegratorTable<M,N>> integrator(T2D6P); // non funziona tanto, cioè non mi piace come si usa manco un po

  // the bilinear form induced by the laplacian is symmetric, therefore compute and store just half of the values
  
  // fill m in some way...
  for(std::shared_ptr<Element<M,N>> e : mesh){
    
    // move the mesh element over the reference element where the functional information is stored (on the reference we have a
    // basis for a space of functions)
    ReferenceElement<N, ORDER> fe(e);
    
    for(size_t i = 0; i < BASIS; ++i){
      for(size_t j = 0; j < BASIS; ++j){
	// get basis elements
	auto phi_i = fe.getBasis().getBasisElement(i).gradient();
	auto phi_j = fe.getBasis().getBasisElement(j).gradient();

	// callable object ready to be integrated
	auto functional = phi_i.dot(phi_j); // for laplacian: \nabla phi_i * \nabla * phi_j

	// integrate the functional over the reference element
	double value = integrator.integrate(functional); 

	/*From Eigen doucmentation: A triplet is a tuple (i,j,value) defining a non-zero element. The input list of triplets
	  does not have to be sorted, and can contains duplicated elements. In any case, the result is a sorted and compressed
	  sparse matrix where the duplicates have been summed up.
	  Linearity of the integral is implicitly used during matrix construction! (there is no explicit sum in this code)
	 */
	tripletList.emplace_back(e->getFESupport()[i].first, e->getFESupport()[j].first, value);
      }
    }        
  }

  // stiff matrix assembled
  m.setFromTriplets(tripletList.begin(), tripletList.end());
};


#endif // __ASSEMBLER_H__

#ifndef __INTERNALS_H__
#define __INTERNALS_H__

#include "../core/utils/Symbols.h"
#include "../core/FEM/PDE.h"
using fdaPDE::core::FEM::PDE;
#include <memory>

namespace fdaPDE{
namespace regression{
namespace internal{

  // compute n x N \Psi matrix.
  template <typename M>
  std::shared_ptr<SpMatrix<double>> psi(M& model) {
    // preallocate space for Psi matrix
    std::shared_ptr<SpMatrix<double>> psi = std::make_shared<SpMatrix<double>>();
    // detect number of observations and number of nodes from model object
    unsigned int n = model.obs();
    unsigned int N = model.loc();    
    psi->resize(n, N);    
    // triplet list to fill sparse matrix
    std::list<Eigen::Triplet<double>> tripletList;

    // data observed at mesh nodes
    if(model.dataAtNodes()){
      // if data locations are equal to mesh nodes then \Psi is the identity matrix
      if(n == N){
	// (\psi_i(p_i) = 1 and \psi_i(p_j) = 0 \forall i \neq j)
	for(std::size_t i = 0; i < n; ++i){
	  tripletList.push_back(Eigen::Triplet<double>(i, i, 1.0));
	}
      }else{ // if data locations are a subset of the mesh nodes then \Psi will have some of its columns set to zero
	// (\psi_i(p_j) = 0 \forall j \in {1, ..., n} such that no data is observed at location i)
	for(std::size_t i = 0; i < n; ++i){
	  tripletList.push_back(Eigen::Triplet<double>(i, model.z_idx()(i,0), 1.0));
	}
      }
    }else{ // general case in which locations are provided as plain coordinates to the model
      // cycle over all locations
      for(std::size_t i = 0; i < model.locations().rows(); ++i){
	SVector<M::embedding_dimension> p_i(model.locations().row(i));
	// search element containing the point
	auto e = model.searchEngine()->search(p_i);
	// update \Psi matrix
	for(std::size_t j = 0; j < model.pde()->basis()[e->ID()]->size(); ++j){
	  std::size_t h = model.pde()->basis()[e->ID()]->nodes()[j]; // column index of \Psi matrix
	  // extract \phi_h from basis
	  auto phi_h = (*model.pde()->basis()[e->ID()])[j];
	  // evaluate \phi_h(p_i) (value of the basis function centered in mesh node h and evaluated in point p_i)
	  tripletList.push_back(Eigen::Triplet<double>(i, h, phi_h(p_i)));
	}
      }
    }
    // finalize construction
    psi->setFromTriplets(tripletList.begin(), tripletList.end());
    psi->makeCompressed();

    return psi;
  }
  
  // an utility to set dirichlet boundary conditions to the regression problem system. Note that this is different from setting
  // dirichlet conditions in the FEM system
  template <unsigned int M, unsigned int N, unsigned int R, typename E>
  void setDirichletBC(const PDE<M,N,R,E>& pde, SpMatrix<double>& A, DVector<double>& b){
    // get number of rows of block in matrix A
    std::size_t n = A.rows()/2;

    for(std::size_t i = 0; i < n; ++i){
      if(pde.domain().isOnBoundary(i)){
	A.row(i) *= 0;       // zero all entries of this row
	A.coeffRef(i,i) = 1; // set diagonal element to 1 to impose equation u_j = b_j

	A.row(i+n) *= 0;
	A.coeffRef(i+n,i+n) = 1;

        // boundaryDatum is a pair (nodeID, boundary value)
	double boundaryDatum = pde.boundaryData().empty() ? 0 : pde.boundaryData().at(i)[0];
	b.coeffRef(i,0) = boundaryDatum; // impose boundary value
	b.coeffRef(i+n, 0) = 0;
      }
    }
    return;
  }
  
}}}

#endif // __INTERNALS_H__

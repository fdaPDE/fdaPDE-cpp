#ifndef __INTERNALS_H__
#define __INTERNALS_H__

#include "../core/utils/Symbols.h"
#include "../core/FEM/PDE.h"
#include <Eigen/src/Core/util/Meta.h>
using fdaPDE::core::FEM::PDE;
#include <memory>

namespace fdaPDE{
namespace regression{
namespace internal{

  // type of possible sampling strategies
  enum SamplingStrategy{ GeostatisticalAtNodes, GeostatisticalAtLocations, Areal };

  // a Triplet type (almost identical with respect to Eigen::Triplet<T>) but allowing for non const access to stored value
  template <typename T>
  class Triplet{
  private:
    Eigen::Index row_, col_;
    T value_;
  public:
    Triplet() = default;
    Triplet(const Eigen::Index& row, const Eigen::Index& col, const T& value)
      : row_(row), col_(col), value_(value) {};
    
    const Eigen::Index& row() const { return row_; }
    const Eigen::Index& col() const { return col_; }
    const T& value() const { return value_; }
    T& value() { return value_; } // allow for modifications of stored value, this not allowed by Eigen::Triplet
  };
  
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
    std::vector<Triplet<double>> tripletList;
    tripletList.reserve(n*N);
    
    switch(model.sampling()){
    case SamplingStrategy::GeostatisticalAtNodes:
      // * if data locations are equal to mesh nodes then \Psi is the identity matrix.
      //     (\psi_i(p_i) = 1 and \psi_i(p_j) = 0 \forall i \neq j)
      // * if data are observed in a subset of the spatial locations then \Psi will have some of its columns set to zero
      //     (\psi_i(p_j) = 0 \forall j \in {1, ..., n} such that no data is observed at location i)
      // This is automatically handled by looking at the observation index
      for(std::size_t i = 0; i < n; ++i){
	tripletList.push_back(Triplet<double>(i, model.z_loc()(i,0), 1.0));
      }
      break;
    case SamplingStrategy::GeostatisticalAtLocations:
      // general case in which locations are provided as plain coordinates to the model. In this case \Psi has no particular structure
      for(std::size_t i = 0; i < model.locations().rows(); ++i){ // cycle over all locations
	SVector<M::embedding_dimension> p_i(model.locations().row(i));
	// search element containing the point
	auto e = model.searchEngine()->search(p_i);
	// update \Psi matrix
	for(std::size_t j = 0; j < model.pde()->basis()[e->ID()].size(); ++j){
	  std::size_t h = e->nodeIDs()[j]; // column index of \Psi matrix
	  // extract \phi_h from basis
	  auto phi_h = model.pde()->basis()[e->ID()][j];
	  // evaluate \phi_h(p_i) (value of the basis function centered in mesh node h and evaluated in point p_i)
	  tripletList.push_back(Triplet<double>(i, h, phi_h(p_i)));
	}
      }
      break;
    case SamplingStrategy::Areal:
      // incidence matrix is a dxM sparse matrix (d number of subdomains, M number of elements) where D_{ij} = 1 \iff element j belongs
      // to subdomain i. In this case matrix \Psi is such that [\Psi]_{ij} = \int_{D_i} \psi_j = \sum_{e \in D_i} \int_{e} \psi_j
      // std::size_t tail = 0;
      // for(std::size_t k = 0; k < model.subdomains().outerSize(); ++k){
      // 	std::size_t head = 0;
      // 	double Di = 0; // measure of subdomain D_i
      // 	for(SpMatrix<int>::InnerIterator it(model.subdomains(),k); it; ++it){
      // 	  // get element with this ID
      // 	  auto e = model.pde()->domain().element(it.col());
      // 	  // compute \int_e \phi_h \forall \phi_h defined on e
      // 	  std::size_t j = 0;
      // 	  for(const auto& phi : *(model.pde()->basis()[e->ID()])){
      // 	    std::size_t h = model.pde()->basis()->nodes()[j]; // node where \phi is centered
      // 	    // evaluate \int_e \phi_h and insert in tripletList. summation is implicitly resolved by Eigen::setFromTriplets
      // 	    tripletList.push_back(Triplet<double>(it.row(), h, model.pde()->integrator().integrate(*e, phi)));
      // 	    head++, j++; // increment counters
      // 	  }
      // 	  Di += e->measure(); // update measure of subdomain D_i
      // 	}
      // 	// divide each \int_{D_i} \psi_j by the measure of subdomain D_i
      // 	for(std::size_t j = 0; j < head; ++j){
      // 	  tripletList[tail + j].value() /= Di;
      // 	}
      // 	tail += head;
      // }
      break;
    }

    // finalize construction
    psi->setFromTriplets(tripletList.begin(), tripletList.end());
    psi->makeCompressed();
    return psi;
  }
  
  // an utility to set dirichlet boundary conditions to the regression problem system. Note that this is different from setting
  // dirichlet conditions in the FEM system
  template <unsigned int M, unsigned int N, unsigned int R, typename E, typename F, typename B, typename I, typename S>
  void setDirichletBC(const PDE<M,N,R,E,F,B,I,S>& pde, SpMatrix<double>& A, DVector<double>& b){
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

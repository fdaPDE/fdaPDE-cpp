#ifndef __SAMPLING_DESIGN_H__
#define __SAMPLING_DESIGN_H__

#include "../core/utils/Symbols.h"
#include "ModelBase.h"
#include "ModelTraits.h"
using fdaPDE::models::is_space_time;

namespace fdaPDE{
namespace models{

  // base classes for the implenetation of the different sampling designs. Here is computed the matrix of spatial basis
  // evaluations \Psi = [\Psi]_{ij} = \psi_i(p_j) whose construction strongly depends on the type of sampling in space

  template <typename Model, Sampling S> class SamplingDesign {};
  
  // data sampled at mesh nodes
  template <typename Model>
  class SamplingDesign<Model, GeoStatMeshNodes> {
  private:
    // recover model object
    inline const Model& model() const { return static_cast<const Model&>(*this); }
  public:
    // constructor
    SamplingDesign() = default;
    SamplingDesign(const DMatrix<double>&) {};
    // init sampling data structures
    void init_sampling() {
      // preallocate space for Psi matrix
      std::size_t n = model().domain().nodes();
      std::size_t N = model().n_basis();    
      Psi_.resize(n, N);    
      // triplet list to fill sparse matrix
      std::vector<fdaPDE::Triplet<double>> tripletList;
      tripletList.reserve(n*N);

      if(!model().data().hasBlock(INDEXES_BLK)){
	// if data locations are equal to mesh nodes then \Psi is the identity matrix.
	// \psi_i(p_i) = 1 and \psi_i(p_j) = 0 \forall i \neq j
	for(std::size_t i = 0; i < n; ++i)
	  tripletList.emplace_back(i, i, 1.0);
      }else{
	// if data are observed in a subset of the spatial locations then \Psi will have some of its columns set to zero
	// \psi_i(p_j) = 0 \forall j \in {1, ..., n} such that no data is observed at location i
	// ** this branch is here for supporting cross validation of models with GeoStatMeshNodes sampling **
	for(std::size_t i = 0; i < n; ++i)
	  tripletList.emplace_back(i, model().idx()(i,0), 1.0); 
      }
      // finalize construction
      Psi_.setFromTriplets(tripletList.begin(), tripletList.end());
      Psi_.makeCompressed();
    }
    
    // getters
    SpMatrix<double> Psi_{}; // n x N matrix \Psi = [\psi_{ij}] = \psi_j(p_i) of spatial basis evaluation at data locations p_i
    auto PsiTD() const { return model().Psi().transpose(); }
    std::size_t n_locs() const { return model().domain().nodes(); }
    DMatrix<double> locs() const { return DMatrix<double>(); }
  };

  // data sampled at general locations p_1, p_2, ... p_n
  template <typename Model>
  class SamplingDesign<Model, GeoStatLocations> {
  private:
    DMatrix<double> locs_;   // matrix of spatial locations p_1, p2_, ... p_n

    // recover model object
    inline Model& model() { return static_cast<Model&>(*this); }
    inline const Model& model() const { return static_cast<const Model&>(*this); }
  public:   
    // constructor
    SamplingDesign(const DMatrix<double>& locs) : locs_(locs) {};
    // init sampling data structures
    void init_sampling() {    
      // preallocate space for Psi matrix
      std::size_t n = locs_.rows();
      std::size_t N = model().n_basis();    
      Psi_.resize(n, N);    
      // triplet list to fill sparse matrix
      std::vector<fdaPDE::Triplet<double>> tripletList;
      tripletList.reserve(n*N);

      auto gse = model().gse(); // geometric search engine
      // cycle over all locations
      for(std::size_t i = 0; i < locs_.rows(); ++i){ 
	SVector<model_traits<Model>::PDE::local_dimension> p_i(locs_.row(i));
	// search element containing the point
	auto e = gse.search(p_i);
	// update \Psi matrix
	for(std::size_t j = 0; j < model().pde().basis()[e->ID()].size(); ++j){
	  std::size_t h = e->nodeIDs()[j]; // column index of \Psi matrix
	  // extract \phi_h from basis
	  auto psi_h = model().pde().basis()[e->ID()][j];
	  // evaluate \phi_h(p_i) (value of the basis function centered in mesh node h and evaluated in point p_i)
	  tripletList.emplace_back(i, h, psi_h(p_i));
	}
      }
      // finalize construction
      Psi_.setFromTriplets(tripletList.begin(), tripletList.end());
      Psi_.makeCompressed();
    };

    // getters
    SpMatrix<double> Psi_{}; // n x N matrix \Psi = [\psi_{ij}] = \psi_j(p_i) of spatial basis evaluation at data locations p_i
    auto PsiTD() const { return model().Psi().transpose(); }
    std::size_t n_locs() const { return locs_.rows(); }
    const DMatrix<double>& locs() const { return locs_; }
  };

  // data sampled at subdomains D_1, D_2, ... D_d
  template <typename Model>
  class SamplingDesign<Model, Areal> {
  private:
    DMatrix<int> subdomains_; // incidence matrix D = [D]_{ij} = 1 \iff element j belongs to subdomain i.
    DiagMatrix<double> D_;    // diagonal matrix of subdomains' measures
    
    // recover model object
    inline const Model& model() const { return static_cast<const Model&>(*this); }
  public:   
    // constructor
    SamplingDesign(const DMatrix<int>& subdomains) : subdomains_(subdomains) {};
    // init sampling data structures
    void init_sampling() {
      // preallocate space for Psi matrix
      std::size_t n = subdomains_.rows();
      std::size_t N = model().n_basis();    
      Psi_.resize(n, N);    
      // triplet list to fill sparse matrix
      std::vector<fdaPDE::Triplet<double>> tripletList;
      tripletList.reserve(n*N);

      DVector<double> D; // store measure of subdomains, this will be ported to a diagonal matrix at the end
      D.resize(subdomains_.rows());      
      // start construction of \Psi matrix
      std::size_t tail = 0;
      for(std::size_t k = 0; k < subdomains_.rows(); ++k){
	std::size_t head = 0;
	double Di = 0; // measure of subdomain D_i
	for(std::size_t l = 0; l < subdomains_.cols(); ++l){
	  if(subdomains_(k,l) == 1){ // element with ID l belongs to k-th subdomain
	    // get element with this ID
	    auto e = model().domain().element(l);
	    // compute \int_e \phi_h \forall \phi_h defined on e
	    std::size_t j = 0;
	    for(const auto& phi : model().pde().basis()[e->ID()]){
	      std::size_t h = phi.node(); // if we write the finite element as \phi_h, this is h
	      // evaluate \int_e \phi_h and insert in tripletList. summation is implicitly resolved by Eigen::setFromTriplets
	      tripletList.emplace_back(k, h, model().pde().integrator().integrate(*e, phi));
	      head++, j++; // increment counters
	    }
	    Di += e->measure(); // update measure of subdomain D_i
	  }
	}
	// divide each \int_{D_i} \psi_j by the measure of subdomain D_i
	for(std::size_t j = 0; j < head; ++j){
	  tripletList[tail + j].value() /= Di;
	}
	D[k] = Di; // store measure of subdomain
	tail += head;
      }
      // here we must be carefull of the type of model (space-only or space-time) we are handling
      if constexpr(is_space_time<Model>::value){
	// store I_m \kron D
	std::size_t m = model().time_domain().rows();
	std::size_t n = n_locs();
	DVector<double> IkronD(n*m);
	for(std::size_t i = 0; i < m; ++i) IkronD.segment(i*n, n) = D;
	// compute and store result
	D_ = IkronD.asDiagonal();
      }else{
	// for space-only problems store diagonal matrix D_ = diag(D_1, D_2, ... ,D_d) as it is
	D_ = D.asDiagonal();
      }
      // finalize construction
      Psi_.setFromTriplets(tripletList.begin(), tripletList.end());
      Psi_.makeCompressed();
    };

    // getters
    SpMatrix<double> Psi_{};  // n x N matrix \Psi = [\Psi]_{ij} = \int_{D_i} \psi_j = \sum_{e \in D_i} \int_{e} \psi_j
    auto PsiTD() const { return model().Psi().transpose()*D_; }
    std::size_t n_locs() const { return subdomains_.rows(); }
    const DiagMatrix<double>& D() const { return D_; }
    const DMatrix<int>& locs() const { return subdomains_; }
  };  
  
}}

#endif // __SAMPLING_DESIGN_H__

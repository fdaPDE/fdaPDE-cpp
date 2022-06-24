#ifndef __PDE_H__
#define __PDE_H__

#include "../utils/Symbols.h"
#include "../MESH/Mesh.h"
#include "operators/BilinearFormTraits.h"
#include <unordered_map>
#include <vector>
using fdaPDE::core::MESH::Mesh;
#include "Assembler.h"
#include "operators/Identity.h"

// top level class to describe a partial differential equation. PDE objects are used by solvers to obtain a solution to the problem

// N and M are the problem dimensions: these informations are strictly releated to the mesh used for
// discretizing the domain. In particular N is the dimension of the problem domain, M is the local dimension of
// the manifold describing the domain. Refer to Mesh documentation for more details
template <unsigned int M,  // dimension of the mesh embedding space
	  unsigned int N,  // local dimension of the mesh
	  typename E>      // type for the BilinearFormExpr
class PDE{
private:
  const Mesh<M,N>& domain_;           // problem domain
  E bilinearForm_;                    // the differential operator of the problem in its weak formulation
  std::vector<DVector> forcingData_;  // forcing data, a vector is used to handle space-time problems
  DVector initialCondition_;          // initial condition, used in space-time problems only
  
  // memorize boundary data in a sparse structure, by storing the index of the boundary node and the relative boundary value.
  // a vector is used as mapped type to handle space-time problems
  std::unordered_map<unsigned, std::vector<double>> boundaryData_;
  
public:
  // expose space dimensions to PDE solver
  static constexpr unsigned M_ = M;
  static constexpr unsigned N_ = N;
  
  // space-only specific methods
  PDE(const Mesh<M,N>& domain, E bilinearForm, const DVector& forcingData) : domain_(domain), bilinearForm_(bilinearForm) {
    forcingData_.push_back(forcingData);
  };
  
  void setDirichletBC(const DVector& data);
  
  // space-time specific methods
  PDE(const Mesh<M,N>& domain, E bilinearForm, const std::vector<DVector>& forcingData) :
    domain_(domain), bilinearForm_(bilinearForm), forcingData_(forcingData) { };

  // vector at index j gives the values at the domain boundary at time step j
  void setDirichletBC(const std::vector<DVector>& data);
  
  void setInitialCondition(const DVector& data) { initialCondition_ = data; };

  //void setNeumannBC();
  
  // getters
  const Mesh<M, N>& getDomain() const { return domain_; }
  E getBilinearForm() const { return bilinearForm_; }
  const std::unordered_map<unsigned, std::vector<double>>& getBoundaryData() const { return boundaryData_; };
  std::vector<DVector> getForcingData() const { return forcingData_; }
  DVector getInitialCondition() const { return initialCondition_; }
};

// argument deduction rule for PDE object
template <unsigned int M, unsigned int N, typename E>
PDE(const Mesh<M,N>& domain, E bilinearForm, const DVector& forcingData) -> PDE<M, N, decltype(bilinearForm)>;

// store in the format (boundaryID, boundaryValue) the dirichlet boundary conditions
// we assume that data is a vector storing at index j the boundary value of mesh node having ID equal to j.
template <unsigned int M, unsigned int N, typename E>
void PDE<M, N, E>::setDirichletBC(const DVector& data){
  // stop if this is a space-time PDE
  static_assert(!is_parabolic<decltype(bilinearForm_)>::value,
		"You are applying a space-only method on a space-time PDE. Aborting");
  
  for(size_t j = 0; j < domain_.getNumberOfNodes(); ++j){
    // if j is a node on the domain boundary store the pair (node ID - boundary value)
    if(domain_.isOnBoundary(j)){
      boundaryData_[j].push_back(data[j]); // O(1) complexity
    }
  }
  return;
}

// store in the format (boundaryID, { ... }) the dirichlet boundary conditions, where { ... } is the time series of the
// data at boundary for boundary node boundaryID
template <unsigned int M, unsigned int N, typename E>
void PDE<M, N, E>::setDirichletBC(const std::vector<DVector>& data){
  // stop if this is a space-only PDE
  static_assert(is_parabolic<decltype(bilinearForm_)>::value,
		"You are applying a space-time method on a space-only PDE. Aborting");

  for(size_t j = 0; j < domain_.getNumberOfNodes(); ++j){
    // if j is a node on the domain boundary store the pair (node ID - boundary value)
    if(domain_.isOnBoundary(j)){
      std::vector<double> boundaryTimeSeries;
      for(const DVector& d : data){
	boundaryTimeSeries.push_back(d[j]);
      }
      
      boundaryData_[j] = boundaryTimeSeries; // O(1) complexity
    }
  }
  return;
}

#endif // __PDE_H__

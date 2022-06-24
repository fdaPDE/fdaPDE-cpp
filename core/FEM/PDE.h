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
  const Mesh<M,N>& domain_;          // problem domain
  E bilinearForm_;                   // the differential operator of the problem in weak formulation
  DVector forcingData_; // forcing data, a vector of vectors is used to handle space-time problems

  // memorize boundary data in a sparse structure, by storing the index of the boundary node and the relative boundary value.
  // a vector is used as mapped type to handle space-time problems
  std::unordered_map<unsigned, std::vector<double>> boundaryData_;
  
public:
  // expose space dimensions to PDE solver
  static constexpr unsigned M_ = M;
  static constexpr unsigned N_ = N;
  
  // construct the PDE object
  PDE(const Mesh<M,N>& domain, E bilinearForm, const DVector& forcingData) :
    domain_(domain), bilinearForm_(bilinearForm), forcingData_(forcingData) {}; 
  
  void setDirichletBC(const DVector& data); // specify dirichlet boundary conditions
  void setNeumannBC();

  // getters
  const Mesh<M, N>& getDomain() const { return domain_; }
  E getBilinearForm() const { return bilinearForm_; }
  const std::unordered_map<unsigned, std::vector<double>>& getBoundaryData() const { return boundaryData_; };
  DVector getForcingData() const { return forcingData_; }
  
};

// argument deduction rule for PDE object
template <unsigned int M, unsigned int N, typename E>
PDE(const Mesh<M,N>& domain, E bilinearForm, const DVector& forcingData) -> PDE<M, N, decltype(bilinearForm)>;

// store in the format (boundaryID, boundaryValue) dirichlet boundary conditions
// we assume that data is a vector storing at index j the boundary value of mesh node having ID equal to j.
template <unsigned int M, unsigned int N, typename E>
void PDE<M, N, E>::setDirichletBC(const DVector& data){
  for(size_t j = 0; j < domain_.getNumberOfNodes(); ++j){
    // if j is a node on the domain boundary store the pair (node ID - boundary value)
    if(domain_.isOnBoundary(j)){
      boundaryData_[j].push_back(data[j]); // O(1) complexity
    }
  }
  return;
}

#endif // __PDE_H__

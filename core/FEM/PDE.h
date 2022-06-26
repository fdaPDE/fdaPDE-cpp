#ifndef __PDE_H__
#define __PDE_H__

#include "../utils/Symbols.h"
#include "../MESH/Mesh.h"
#include "operators/BilinearFormTraits.h"
#include <cstddef>
#include <type_traits>
#include <unordered_map>
#include <vector>
using fdaPDE::core::MESH::Mesh;
#include "Assembler.h"
#include "operators/Identity.h"

// forward declarations
template <typename B, typename I> class FEMStandardSpaceTimeSolver;
template <typename B, typename I> class FEMStandardSpaceSolver;

// trait to select the space-only or the space-time variant of the PDE standard solver
template <typename E, typename B, typename I> struct pde_standard_solver_selector {
  using type = typename std::conditional<is_parabolic<E>::value,
					 FEMStandardSpaceTimeSolver<B, I>,
					 FEMStandardSpaceSolver<B, I>>::type;
};

// top level class to describe a partial differential equation. PDE objects are used by solvers to obtain a solution to the problem

// N and M are the problem dimensions: these informations are strictly releated to the mesh used for
// discretizing the domain. In particular N is the dimension of the problem domain, M is the local dimension of
// the manifold describing the domain. Refer to Mesh documentation for more details
template <unsigned int M,  // dimension of the mesh embedding space
	  unsigned int N,  // local dimension of the mesh
	  typename E>      // type for the BilinearFormExpr
class PDE{
private:
  const Mesh<M,N>& domain_;     // problem domain
  E bilinearForm_;              // the differential operator of the problem in its weak formulation
  DMatrix forcingData_{};         // forcing data, a vector is used to handle space-time problems
  DVector initialCondition_{};    // initial condition, used in space-time problems only
  
  // memorize boundary data in a sparse structure, by storing the index of the boundary node and the relative boundary value.
  // a vector is used as mapped type to handle space-time problems
  std::unordered_map<unsigned, DVector> boundaryData_{};

  DMatrix solution_;
public:
  // expose space dimensions to PDE solver
  static constexpr unsigned M_ = M;
  static constexpr unsigned N_ = N;

  typedef std::unordered_map<unsigned, DVector> boundary_data;
  
  // constructor, a DMatrix is accepted as forcingData to handle also space-time problems
  PDE(const Mesh<M,N>& domain, E bilinearForm, const DMatrix& forcingData) :
    domain_(domain), bilinearForm_(bilinearForm), forcingData_(forcingData) {};

  // vector at index j gives the values at the domain boundary at time step j
  void setDirichletBC(const DMatrix& data);
  //void setNeumannBC();
  
  void setInitialCondition(const DVector& data) { initialCondition_ = data; };

  // solves the differential equation using the supplied solver.
  template <typename B, typename I, typename SOLVER = typename pde_standard_solver_selector<E, B, I>::type, typename... Args >
  void solve(const B& base, const I& integrator, Args... args);
  
  // getters
  const Mesh<M, N>& getDomain()          const { return domain_; }
  E getBilinearForm()                    const { return bilinearForm_; }
  const boundary_data& getBoundaryData() const { return boundaryData_; };
  DMatrix getForcingData()               const { return forcingData_; }
  DVector getInitialCondition()          const { return initialCondition_; }

  DMatrix getSolution() const { return solution_; }

  void setSolution(const DMatrix& solution) { solution_ = solution; }
};

// argument deduction rule for PDE object
template <unsigned int M, unsigned int N, typename E>
PDE(const Mesh<M,N>& domain, E bilinearForm, const DVector& forcingData) -> PDE<M, N, decltype(bilinearForm)>;

// store in the format (boundaryID, { ... }) the dirichlet boundary conditions, where { ... } is the time series of the
// data at boundary for boundary node boundaryID
template <unsigned int M, unsigned int N, typename E>
void PDE<M, N, E>::setDirichletBC(const DMatrix& data){
 for(size_t j = 0; j < domain_.getNumberOfNodes(); ++j){
    // if j is a node on the domain boundary store the pair (node ID - boundary value)
    if(domain_.isOnBoundary(j)){
      boundaryData_[j] = data.row(j); // O(1) complexity
    }
  }
  return;
}

template <unsigned int M, unsigned int N, typename E>
template <typename B, typename I, typename SOLVER, typename... Args>
void PDE<M, N, E>::solve(const B &base, const I &integrator, Args... args) {
  // define solver and call solve method on it
  SOLVER solver(base, integrator);
  solver.solve(*this, args...);

  return;
}


#endif // __PDE_H__

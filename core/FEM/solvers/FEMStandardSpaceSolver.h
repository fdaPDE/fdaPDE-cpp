#ifndef __FEM_STANDARD_SPACE_SOLVER_H__
#define __FEM_STANDARD_SPACE_SOLVER_H__

#include "solvers/FEMBaseSolver.h"

template <typename B, typename I>
class FEMStandardSpaceSolver {
private:
  const I& integrator_; // integrator used to approximate integrals
  const B& basis_;      // basis used for functional domain discretization
  
  // informations produced by FEM. Filled by .solve() method
  DVector solution_;                // vector of coefficients of the approximate solution written in terms of the chosen basis
  Eigen::SparseMatrix<double> R1_;  // result of the discretization of the bilinear form, also known as R1 (stiff matrix)
  DVector forcingVector_;           // right-hand side of the linear system giving the FEM solution

  // informations not directly related to PDE solution but usefull for higher layers of the architecture
  Eigen::SparseMatrix<double> R0_;
 
  // some informations about solution error

  // initializes internal FEM solver status
  template <unsigned int M, unsigned int N, typename E> 
  void init(const PDE<M, N, E>& pde);
  
public:
  // constructor
  FEMStandardSpaceSolver(const B& basis, const I& integrator) : basis_(basis), integrator_(integrator) {};
  // flag used to notify is something was wrong during computation of solution
  bool success = true;

  // solve PDE
  template <unsigned int M, unsigned int N, typename E> 
  void solve(const PDE<M, N, E>& pde);

  // getters
  DVector getSolution() const { return solution_; }
};

// fill all internal data structures required by FEM to solve the problem. These
// operations constitute the core of FEM and should be independent on any specific solver (both space or space-time)
template <typename B, typename I>
template <unsigned int M, unsigned int N, typename E> 
void FEMStandardSpaceSolver<B, I>::init(const PDE<M, N, E>& pde) {
  Assembler<M, N, B, I> assembler(pde.getDomain(), basis_, integrator_); // create assembler object
  R1_ = assembler.assemble(pde.getBilinearForm());       // fill discretization matrix for current operator
  // SparseQR solver needs its matrix in compressed form (see Eigen documentation for details)
  R1_.makeCompressed();
  
  // differs between space and space-time solvers
  forcingVector_ = assembler.forcingTerm(pde.getForcingData()); // fill discretization of rhs for FEM linear system

  // impose boundary conditions
  for(const auto& boundaryDatum : pde.getBoundaryData()){
    // boundaryDatum is a pair (nodeID, boundary value)
    
    // To impose a Dirichlet boundary condition means to introduce an equation of the kind u_j = b_j where j is the index
    // of the boundary node and b_j is the boundary value we want to impose on this node. This actually removes one degree
    // of freedom from the system. We do so by zeroing out the j-th row of the stiff matrix and set the corresponding
    // diagonal element to 1
    R1_.row(boundaryDatum.first) *= 0;                             // zero all entries of this row
    R1_.coeffRef(boundaryDatum.first, boundaryDatum.first) = 1;    // set diagonal element to 1 to impose equation u_j = b_j
    forcingVector_[boundaryDatum.first] = boundaryDatum.second[0]; // impose boundary value
  }

  // R0_ is a mass matrix ([R0]_{ij} = \int_{\Omega} \phi_i \phi_j). This quantity can be obtained by computing
  // the discretization of the Identity() operator
  R0_ = assembler.assemble(Identity());
  return;
}

template <typename B, typename I>
template <unsigned int M, unsigned int N, typename E> 
void FEMStandardSpaceSolver<B, I>::solve(const PDE<M, N, E>& pde){
  init(pde); // init solver for this PDE
  
  // define eigen system solver, use QR decomposition.
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
  solver.compute(R1_);
  
  // stop if something was wrong...
  if(solver.info()!=Eigen::Success) {
    success = false;
    return;
  }
  
  // solve FEM linear system: discretizationMatrix_*solution_ = forcingVector_;
  solution_ = solver.solve(forcingVector_);
  return;
}

#endif // __SPACE_SOLVER_H__

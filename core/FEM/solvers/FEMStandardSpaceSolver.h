#ifndef __FEM_STANDARD_SPACE_SOLVER_H__
#define __FEM_STANDARD_SPACE_SOLVER_H__

#include "../Assembler.h"
#include "../../utils/Symbols.h"
#include "../PDE.h"

template <typename B, typename I>
class FEMStandardSpaceSolver {
private:
  const I& integrator_; // integrator used to approximate integrals
  const B& basis_;      // basis used as approximation of the infinite dimensional space where the PDE solution is searched
  
  DVector solution_;                // vector of coefficients of the approximate solution written in terms of the chosen basis
  DVector forcingVector_;           // right-hand side of the linear system giving the FEM solution
  Eigen::SparseMatrix<double> R1_;  // result of the discretization of the bilinear form, also known as R1.
  Eigen::SparseMatrix<double> R0_;  // mass matrix, needed by components in higher levels of the architecture, known as R0.
 
  // some informations about solution error

  // initializes internal FEM solver status
  template <unsigned int M, unsigned int N, typename E> 
  void init(const PDE<M, N, E>& pde);
  
public:
  // constructor
  FEMStandardSpaceSolver(const B& basis, const I& integrator) : basis_(basis), integrator_(integrator) {};
  // flag used to notify is something was wrong during computation of solution
  bool success = true;

  // solves the PDE using the classical FEM approach: compute stiffness matrix using some finite element basis R1_ and forcing
  // vector b, then solves the linear system R1_*u = b where u is the searched PDE approximation
  template <unsigned int M, unsigned int N, typename E> 
  void solve(const PDE<M, N, E>& pde);

  // getters
  DVector getSolution() const { return solution_; }
  DVector getForce() const { return forcingVector_; }
  Eigen::SparseMatrix<double> getR1() const { return R1_; }
  Eigen::SparseMatrix<double> getR0() const { return R0_; }
};

// fill internal data structures required by FEM to solve the problem.
template <typename B, typename I>
template <unsigned int M, unsigned int N, typename E> 
void FEMStandardSpaceSolver<B, I>::init(const PDE<M, N, E>& pde) {
  Assembler<M, N, B, I> assembler(pde.getDomain(), basis_, integrator_); // create assembler object
  R1_ = assembler.assemble(pde.getBilinearForm());       // fill discretization matrix for current operator
  // SparseQR solver needs its matrix in compressed form (see Eigen documentation for details)
  R1_.makeCompressed();

  forcingVector_ = assembler.forcingTerm(pde.getForcingData()); // fill discretization of rhs for FEM linear system

  // impose boundary conditions
  for(std::size_t i = 0; i < pde.getDomain().getNumberOfNodes(); ++i){
    if(pde.getDomain().isOnBoundary(i)){
      // boundaryDatum is a pair (nodeID, boundary value)
      double boundaryDatum = pde.getBoundaryData().at(i)[0];

      // To impose a Dirichlet boundary condition means to introduce an equation of the kind u_j = b_j where j is the index
      // of the boundary node and b_j is the boundary value we want to impose on this node. This actually removes one degree
      // of freedom from the system. We do so by zeroing out the j-th row of the stiff matrix and set the corresponding
      // diagonal element to 1
      R1_.row(i) *= 0;                   // zero all entries of this row
      R1_.coeffRef(i, i) = 1;            // set diagonal element to 1 to impose equation u_j = b_j
      forcingVector_[i] = boundaryDatum; // impose boundary value
    }
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

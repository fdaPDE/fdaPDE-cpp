#ifndef __FEM_BASE_SOLVER_H__
#define __FEM_BASE_SOLVER_H__

template <typename B, typename I>
class FEMBaseSolver{
protected:
  const I& integrator_;             // integrator used to approximate integrals
  const B& basis_;                  // basis used as approximation of the infinite dimensional space where the PDE solution is searched
  
  DMatrix solution_;                // vector of coefficients of the approximate solution written in terms of the chosen basis
  DMatrix forcingVector_;           // right-hand side of the linear system giving the FEM solution
  Eigen::SparseMatrix<double> R1_;  // result of the discretization of the bilinear form, also known as R1.
  Eigen::SparseMatrix<double> R0_;  // mass matrix, needed by components in higher levels of the architecture, known as R0.

  // some informations about solution error

  // initializes internal FEM solver status
  template <unsigned int M, unsigned int N, typename E> 
  void init(const PDE<M, N, E>& pde);
  
public:
  FEMBaseSolver(const B& basis, const I& integrator) : basis_(basis), integrator_(integrator) {};
    // flag used to notify is something was wrong during computation of solution
  bool success = true;

  // getters
  DMatrix getSolution() const { return solution_; }
  DMatrix getForce() const { return forcingVector_; }
  Eigen::SparseMatrix<double> getR1() const { return R1_; }
  Eigen::SparseMatrix<double> getR0() const { return R0_; }
};

// fill internal data structures required by FEM to solve the problem.
template <typename B, typename I>
template <unsigned int M, unsigned int N, typename E> 
void FEMBaseSolver<B, I>::init(const PDE<M, N, E>& pde) {
  Assembler<M, N, B, I> assembler(pde.getDomain(), basis_, integrator_); // create assembler object
  R1_ = assembler.assemble(pde.getBilinearForm());       // fill discretization matrix for current operator
  // SparseQR solver needs its matrix in compressed form (see Eigen documentation for details)
  R1_.makeCompressed();
  forcingVector_.resize(pde.getDomain().getNumberOfNodes(), pde.getForcingData().cols());

  for(std::size_t i = 0; i < pde.getForcingData().cols(); ++i){
    forcingVector_.col(i) = assembler.forcingTerm(pde.getForcingData().col(i)); // fill discretization of rhs for FEM linear system
  }

  // R0_ is a mass matrix ([R0]_{ij} = \int_{\Omega} \phi_i \phi_j). This quantity can be obtained by computing
  // the discretization of the Identity() operator
  R0_ = assembler.assemble(Identity());
  return;
}

#endif // __FEM_BASE_SOLVER_H__

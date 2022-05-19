#ifndef __PDE_H__
#define __PDE_H__

#include "../utils/Symbols.h"
#include "../MESH/Mesh.h"
using fdaPDE::core::MESH::Mesh;
#include "Assembler.h"
#include "operators/Identity.h"

// top level class to describe a partial differential equation

// N and M are the problem dimensions: these informations are strictly releated to the mesh used for
// discretizing the domain. In particular N is the dimension of the problem domain, M is the local dimension of
// the manifold describing the domain. Refer to Mesh documentation for more details
template <unsigned int M,     // dimension of the mesh embedding space
	  unsigned int N,     // local dimension of the mesh
	  unsigned int ORDER, // order of FEM approximation basis
	  typename E,         // type for the BilinearFormExpr
	  typename I>         // integrator type
class PDE{
private:
  const I& integrator_;                                   // integrator used for integral approximations
  Mesh<M,N>& domain_;                                     // problem domain
  E bilinearForm_;                                        // the differential operator of the problem
  Eigen::Matrix<double, Eigen::Dynamic, 1> forcingData_;  // forcing data
  Eigen::Matrix<double, Eigen::Dynamic, 1> boundaryData_; // boundary conditions

  // informations stored inside a PDE object. Filled by .solve() method
  Eigen::Matrix<double, Eigen::Dynamic, 1> solution_;      // vector of coefficients of the approximate solution written in terms of the chosen basis
  Eigen::SparseMatrix<double> R1_;                         // result of the discretization of the bilinear form, also known as R1
  Eigen::Matrix<double, Eigen::Dynamic, 1> forcingVector_; // right-hand side of the linear system giving the FEM solution

  // informations not directly related to PDE solution but usefull for higher layers of the architecture
  Eigen::SparseMatrix<double> R0_;
  // basis function... matrix big Phi here

  // some informations about solution error
  
  void init(); // assemble system matrices and internal data structures
public:
  bool success = true; // flag used to notify is something was wrong during computation of solution

  // construct a PDE object, homogeneous boundary conditions implicitly assumed
  PDE(Mesh<M,N>& domain, E bilinearForm, const Eigen::Matrix<double, Eigen::Dynamic, 1>& forcingData, const I& integrator) :
    domain_(domain), bilinearForm_(bilinearForm), forcingData_(forcingData), integrator_(integrator) { init(); }; 
  
  void setDirichletBC(const Eigen::Matrix<double, Eigen::Dynamic, 1>& boundaryData); // impose dirichlet boundary conditions
  void setNeumannBC();

  void solve(); // solves the PDE
 
  // getters
  Eigen::Matrix<double, Eigen::Dynamic, 1> getSolution() const { return solution_; }
};

// fill all internal data structures relative to the problem. Automatically called by constructor
template <unsigned int M, unsigned int N, unsigned int ORDER, typename E, typename I>
void PDE<M, N, ORDER, E, I>::init() {
  Assembler<M,N,ORDER> assembler(domain_);              // create assembler object 
  R1_ = assembler.assemble(bilinearForm_);              // fill discretization matrix for current operator
  // SparseQR solver needs its matrix in compressed form (see Eigen documentation for details)
  R1_.makeCompressed();
  forcingVector_ = assembler.forcingTerm(forcingData_); // fill discretization of rhs for FEM linear system

  // R0_ is a mass matrix ([R0]_{ij} = \int_{\Omega} \phi_i \phi_j). This quantity can be obtained by computing
  // the discretization of the Identity() operator
  R0_ = assembler.assemble(Identity());
  
  return;
}

template <unsigned int M, unsigned int N, unsigned int ORDER, typename E, typename I>
void PDE<M, N, ORDER, E, I>::solve() {
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

template <unsigned int M, unsigned int N, unsigned int ORDER, typename E, typename I>
void PDE<M, N, ORDER, E, I>::setDirichletBC(const Eigen::Matrix<double, Eigen::Dynamic, 1>& boundaryData){
  // impose dirichet conditions in R1_
  for(size_t j = 0; j < domain_.getNumberOfNodes(); ++j){
    // if this node is a boundary node, change the corresponding row in the stiff matrix
    // To impose a Dirichlet boundary condition means to introduce an equation of the kind u_j = b_j where j is the index
    // of the boundary node and b_j is the boundary value we want to impose on this node. This actually removes one degree
    // of freedom from the system. We do so by zeroing out the j-th row of the stiff matrix and set the corresponding
    // diagonal element to 1
    if(domain_.isOnBoundary(j)){
      R1_.row(j) *= 0;                     // zero all entries of this row
      R1_.coeffRef(j,j) = 1;               // set diagonal element to 1 to impose equation u_j = b_j
      forcingVector_[j] = boundaryData[j]; // impose boundary value
    }
  }
  return;
}


#endif // __PDE_H__

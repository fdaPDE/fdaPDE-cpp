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
  
  // serve anche l'informazione sull'integratore da utilizzare... mettine uno di default ma comunque bisogna lasciare opportunità di settarlo dall'esterno
  
public:
  bool success = true; // flag used to notify is something was wrong during computation of solution

  // construct a PDE object, homogeneous boundary conditions implicitly assumed
  PDE(Mesh<M,N>& domain, E bilinearForm, const Eigen::Matrix<double, Eigen::Dynamic, 1>& forcingData, const I& integrator) :
    domain_(domain), bilinearForm_(bilinearForm), forcingData_(forcingData), integrator_(integrator) {}; 
  
  void assemble(); // assemble FEM matrices without solving the system
  void setDirichletConditions(); // impose boundary conditions
  void setNeumannConditions();

  void solve(); // solve the system providing the solution vector

  // getters
  Eigen::Matrix<double, Eigen::Dynamic, 1> getSolution() const { return solution_; }
};

// solve the differential problem and fill all internal data structures
template <unsigned int M, unsigned int N, unsigned int ORDER, typename E, typename I>
void PDE<M, N, ORDER, E, I>::solve() {
  Assembler<M,N,ORDER> assembler(domain_);              // create assembler object 
  R1_ = assembler.assemble(bilinearForm_);              // fill discretization matrix for current operator
  forcingVector_ = assembler.forcingTerm(forcingData_); // fill discretization of rhs for FEM linear system

  // SparseQR solver needs its matrix in compressed form (see Eigen documentation for details)
  R1_.makeCompressed();
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
  solver.compute(R1_);

  // stop if something was wrong...
  if(solver.info()!=Eigen::Success) {
    success = false;
    return;
  }

  // solve FEM linear system: discretizationMatrix_*solution_ = forcingVector_;
  solution_ = solver.solve(forcingVector_);

  // R0_ is a mass matrix ([R0]_{ij} = \int_{\Omega} \phi_i \phi_j). This quantity can be obtained by computing
  // the discretization of the Identity() operator
  R0_ = assembler.assemble(Identity());
  
  return;
}




#endif // __PDE_H__

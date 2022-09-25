#ifndef __INTERNALS_H__
#define __INTERNALS_H__

#include "../core/utils/Symbols.h"
#include "../core/FEM/PDE.h"
#include <Eigen/src/Core/Matrix.h>
#include <memory>
#include <tuple>
using fdaPDE::core::FEM::PDE;

namespace fdaPDE{
namespace regression{
namespace internal{

  // compute Psi matrix assuming locations equal to mesh's nodes (there is 1 only at mesh nodes and 0 elsewhere due to support of lagrangian basis)
  // in general it is not diagonal!
  template <unsigned int M, unsigned int N, unsigned int R, typename E>
  std::unique_ptr<SpMatrix<double>> psi(const PDE<M,N,R,E>& pde) {
    // preallocate space for Psi matrix
    std::unique_ptr<SpMatrix<double>> psi = std::make_unique<SpMatrix<double>>();
    unsigned int locations = pde.domain().nodes();
    unsigned int nbasis = pde.domain().nodes();
    psi->resize(locations, nbasis);
  
    // fill matrix
    std::list<Eigen::Triplet<double>> tripletList;  
    for(std::size_t i = 0; i < locations; ++i){
      tripletList.push_back(Eigen::Triplet<double>(i, i, 1.0));
    }
  
    psi->setFromTriplets(tripletList.begin(), tripletList.end());
    psi->makeCompressed();
    return psi;
  }
  
  // an efficient way to perform a left multiplication by Q. The following method is an implementation of the algorithm
  //  given the design matrix W and x
  //    compute, factorize and store Y = W^T*W
  //    compute v = W^T*x
  //    solve Yz = v
  //    return x - Wz = Qx
  // it is required to having assigned a design matrix W to the model before calling this method
  template <typename M>
  DMatrix<double> lmbQ(const M& model, const DMatrix<double>& x){
    DMatrix<double> v = model.W()->transpose()*x; // W^T*x
    DMatrix<double> z = model.invWTW().solve(v);  // (W^T*W)^{-1}*W^T*x
    // compute x - W*z = x - (W*(W^T*W)^{-1}*W^T)*x = (I - H)*x = Q*x
    return x - (*model.W())*z;
  }


  // an utility to set dirichlet boundary conditions to the regression problem system. Note that this is different from setting
  // dirichlet conditions in the FEM system
  template <unsigned int M, unsigned int N, unsigned int R, typename E>
  void setDirichletBC(const PDE<M,N,R,E>& pde, SpMatrix<double>& A, DVector<double>& b){
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

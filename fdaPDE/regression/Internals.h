#ifndef __INTERNALS_H__
#define __INTERNALS_H__

#include "../core/utils/Symbols.h"
#include "../core/FEM/PDE.h"
#include <Eigen/src/Core/Matrix.h>
#include <memory>
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
  
    // fill psi matrix
    std::list<Eigen::Triplet<double>> tripletList;  
    for(std::size_t i = 0; i < locations; ++i){
      tripletList.push_back(Eigen::Triplet<double>(i, i, 1.0));
    }
  
    psi->setFromTriplets(tripletList.begin(), tripletList.end());
    psi->makeCompressed();
    return psi;
  }

  std::unique_ptr<DMatrix<double>> H(const DMatrix<double>& W){
    // compute transpose of W once here
    DMatrix<double> Wt = W.transpose();
    // H = W*(W*W^T)^{-1}*W^T
    std::unique_ptr<DMatrix<double>> H = std::make_unique<DMatrix<double>>(W*(Wt*W).ldlt().solve(Wt));
    return H;
  }

  std::unique_ptr<DMatrix<double>> Q(const DMatrix<double>& H){
    DMatrix<double> I = DMatrix<double>::Identity(H.rows(), H.cols());
    // Q = I - H
    return std::make_unique<DMatrix<double>>(I - H);
  }
  
}}}

#endif // __INTERNALS_H__

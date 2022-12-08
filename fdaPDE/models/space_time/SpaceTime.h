#ifndef __SPACE_TIME_H__
#define __SPACE_TIME_H__

#include "../../core/utils/Symbols.h"
#include "SplineBasis.h"
using fdaPDE::models::SplineBasis;
#include "../../core/FEM/integration/IntegratorTables.h"
using fdaPDE::core::FEM::IntegratorTable;
using fdaPDE::core::FEM::GaussLegendre;
#include <type_traits>

namespace fdaPDE {
namespace models {

  // general class for the assembly of matrix \Phi = [\Phi]_{ij} = \phi_i(t_j),
  // being \phi_1, \phi_2, ... \phi_M a basis system of type B
  template <typename B> class PhiAssembler;

  // specialization of \Phi for a spline basis of order R
  template <unsigned int R> struct PhiAssembler<SplineBasis<R>> {
    typedef SplineBasis<R> B;
    static SpMatrix<double> compute(const DVector<double>& time_);
  };
  
  template <unsigned int R>
  SpMatrix<double> PhiAssembler<SplineBasis<R>>::compute(const DVector<double>& time_) {
    // create time basis
    B basis(time_);
    // define Phi matrix dimensions
    int m = time_.rows();
    int M = basis.size();
    // resize result matrix
    SpMatrix<double> Phi;
    Phi.resize(m, M);
    
    // triplet list to fill sparse matrix \Phi
    std::vector<fdaPDE::Triplet<double>> tripletList;
    tripletList.reserve(m*M);

    for(int i = 0; i < M; ++i){
      // exploit local support of splines
      for(int j = ((int)(i-B::order) > 0 ? (int)(i-B::order) : 0); j <= std::min(m-1, i+1); ++j){
	tripletList.emplace_back(j, i, basis[i](SVector<1>(time_[j])));
      }
    }
    // finalize construction
    Phi.setFromTriplets(tripletList.begin(), tripletList.end());
    Phi.makeCompressed();    
    return Phi;
  }

  // convenient function template for the construction of \Phi
  template <typename B>
  SpMatrix<double> Phi(const DVector<double>& time) {
    return PhiAssembler<B>::compute(time);
  }
  
}}

#endif // __SPACE_TIME_H__

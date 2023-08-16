// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __EXACT_EDF_H__
#define __EXACT_EDF_H__

#include <fdaPDE/utils.h>
#include "../models/regression/regression_base.h"
using fdapde::models::is_regression_model;

namespace fdapde {
namespace calibration{

  // Evaluates exactly the trace of matrix S = \Psi*T^{-1}*\Psi^T*Q. Uses the cyclic property of the trace
  // operator: Tr[S] = Tr[\Psi*T^{-1}*\Psi^T*Q] = Tr[Q*\Psi*T^{-1}*\Psi^T]
  template <typename Model>
  class ExactEDF {
    static_assert(is_regression_model<Model>::value);
  private:
    Model& m_;

    // S = Q*\Psi*T^{-1}*\Psi^T
    const DMatrix<double>& S() {
      // compute \Psi^T*D*Q (take into account of areal sampling)
      // if(m_.hasCovariates())
      // 	E_ = m_.PsiTD()*m_.Q();
      // else E_ = m_.PsiTD();

      // factorize matrix T
      invT_ = m_.T().partialPivLu();
      DMatrix<double> E_ = m_.PsiTD(); // need to cast to dense for PartialPivLU::solve()
      S_ = m_.lmbQ(m_.Psi()*invT_.solve(E_));
      
      // V_ = invT_.solve(E_); // V = invT*E = T^{-1}*\Psi^T*Q
      // S_ = m_.Psi()*V_; // S = \Psi*V
      return S_;
    };    

  public:
    // let public access since those might be required for the computation of GCV derivatives
    Eigen::PartialPivLU<DMatrix<double>> invT_{}; // T^{-1}
    DMatrix<double> E_{}; // \Psi^T*Q
    DMatrix<double> V_{}; // T^{-1}*\Psi^T*Q
    DMatrix<double> S_{}; // \Psi*T^{-1}*\Psi^T*Q = \Psi*V_

    // constructor
    explicit ExactEDF(Model& m) : m_(m) {};
    double compute() { return S().trace(); } // computes Tr[S]
  };

}}
#endif // __EXACT_EDF_H__

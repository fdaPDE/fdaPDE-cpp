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

namespace fdapde {
namespace models {

// Evaluates exactly the trace of matrix S = \Psi*T^{-1}*\Psi^T*Q. Uses the cyclic property of the trace
// operator: Tr[S] = Tr[\Psi*T^{-1}*\Psi^T*Q] = Tr[Q*\Psi*T^{-1}*\Psi^T]
template <typename Model> class ExactEDF {
   private:
    Model& m_;
    // computes smoothing matrix S = Q*\Psi*T^{-1}*\Psi^T
    const DMatrix<double>& S() {
        // factorize matrix T
        invT_ = m_.T().partialPivLu();
        DMatrix<double> E_ = m_.PsiTD();    // need to cast to dense for PartialPivLU::solve()
        S_ = m_.lmbQ(m_.Psi() * invT_.solve(E_));   // \Psi*T^{-1}*\Psi^T*Q
        return S_;
    };
   public:
    Eigen::PartialPivLU<DMatrix<double>> invT_ {};   // T^{-1}
    DMatrix<double> S_ {};                           // \Psi*T^{-1}*\Psi^T*Q = \Psi*V_
  
    explicit ExactEDF(Model& m) : m_(m) {};
    double compute() { return S().trace(); }   // computes Tr[S]
};

}   // namespace models
}   // namespace fdapde

#endif   // __EXACT_EDF_H__

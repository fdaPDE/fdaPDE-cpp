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

#ifndef __SPECKMAN_H__
#define __SPECKMAN_H__

#include <fdaPDE/linear_algebra.h>
#include <fdaPDE/utils.h>
using fdapde::core::FSPAI;
using fdapde::core::lump;
using fdapde::core::is_empty;

#include "../model_macros.h"
#include "../model_traits.h"
#include "../model_base.h"
#include "srpde.h"
#include "strpde.h"
#include "exact_edf.h"
#include "stochastic_edf.h"
#include "inference_base.h"
#include "inference.h"


namespace fdapde {
namespace models {

template <typename Model, typename Strategy> class Speckman: public InferenceBase<Model> {

   private:
      struct ExactInverse{
         DMatrix<double> compute(Model m){
            return inverse(m.E());       
         }
      };

      struct NonExactInverse{
         SpMatrix<double> compute(Model m){
            SpMatrix<double> invE_ = Base::invE_approx(m);
            return Base::invE_approx(m);
         }
      };

      DMatrix<double> Lambda_ {};

   public: 
      using Base = InferenceBase<Model>;
      using Base::m_;
      using Base::V_;
      using Base::beta_;
      using Solver = typename std::conditional<std::is_same<Strategy, exact>::value, ExactInverse, NonExactInverse>::type;
      Solver s_; 


      // constructors
      Speckman() = default;                   // deafult constructor
      Speckman(const Model& m): Base(m) {};     // constructor    

      // return Lambda_^2
      DMatrix<double> Lambda() {
         int n = m_.n_obs();
         DMatrix<double> pax =  s_.compute(m_);
         DMatrix<double> Lambda = DMatrix<double>::Identity(n, n) - m_.Psi() * pax * m_.PsiTD();
         DMatrix<double> Lambda_squared = Lambda * Lambda;
         return Lambda_squared;
      }

      void beta() override{
         if(is_empty(Lambda_)){
               Lambda_ = Lambda();
         }
         DMatrix<double> W = m_.X();
         DMatrix<double> invWtW = inverse(W.transpose() * Lambda_ * (W));      
         beta_ = invWtW * W.transpose() * Lambda_ * (m_.y());         
      }

      void V() override{
         if(is_empty(Lambda_)){
               Lambda_ = Lambda();
         }
         DMatrix<double> W = m_.X();
         DMatrix<double> invWtW = inverse(W.transpose() * Lambda_ * (W));
         DVector<double> eps_ = (m_.y() - m_.fitted());
         DVector<double> Res2 = eps_.array() * eps_.array();            
         // resize the variance-covariance matrix
         V_.resize(m_.q(), m_.q());                   
         DMatrix<double> W_t = W.transpose();           
         DMatrix<double> diag = Res2.asDiagonal();
         V_ = invWtW * (W_t) * Lambda_ * Res2.asDiagonal() * Lambda_ * (W) * invWtW;      
      }
 
};


} // namespace models
} // namespace fdapde

#endif   // __SPECKMAN_H__
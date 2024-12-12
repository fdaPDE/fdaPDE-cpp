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

#ifndef __INFERENCE_BASE_H__
#define __INFERENCE_BASE_H__

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
#include "inference.h"

#include <cmath>
#include <random>

//serve per salvare le matrici 
#include <unsupported/Eigen/SparseExtra> 
#include <fstream>



namespace fdapde {
namespace models {

template <typename Model> class InferenceBase{

   protected:

      Model m_;
      DMatrix<double> V_ {};
      DMatrix<double> C_ {};            // inference matrix C (p x q) matrix  
      DVector<double> beta_ {};        // sol of srpde ( q x 1 ) matrix          
      DVector<double> beta0_ {};        // inference hypothesis H0 (p x 1) matrix
      DVector<double> f0_ {};          // inference hypothesis H0
      double alpha_ = 0;                // level of the confidence intervals for beta
      double alpha_f_ = 0.05;              // level of confidence intervals for f 
      DVector<int> locations_f_ {};   // indexes of the subset of locations if locations are exctracted from existing ones


   public: 
     
      InferenceBase() = default;                   // deafult constructor
      InferenceBase(const Model& m): m_(m) {};     // constructor    

      virtual void beta() {};

      virtual void V() = 0;

      virtual DMatrix<double> computeCI(CIType type){ 

         fdapde_assert(!is_empty(C_));  

         if(alpha_ == 0.0) {
            setAlpha(0.05);         // default value 5%
         }
         if(is_empty(V_)){
            V();
         }
         if(is_empty(beta_)){
            beta();
         }
         
         int p = C_.rows();
         int size = std::min(C_.rows(), V_.rows());
         DVector<double> diagon(size);
         for (int i = 0; i < C_.rows(); ++i) {
            DVector<double> ci = C_.row(i);
            diagon[i] = ci.transpose() * V_ * ci;
         }

         DVector<double> lowerBound(size);
         DVector<double> upperBound(size);

         if(type == simultaneous){ 
            // SIMULTANEOUS
            double quantile = inverseChiSquaredCDF(1 - alpha_, p);
            lowerBound = (C_ * beta_).array() - (quantile * diagon.array()).sqrt();
            upperBound = (C_ * beta_).array() + (quantile * diagon.array()).sqrt();         
         }
         else if (type == bonferroni){
            // BONFERRONI
            double quantile = normal_standard_quantile(1 - alpha_/(2 * p));            
            lowerBound = (C_ * beta_).array() - quantile * (diagon.array()).sqrt();
            upperBound = (C_ * beta_).array() + quantile * (diagon.array()).sqrt();
         }
         else if (type == one_at_the_time){
            // ONE AT THE TIME
            double quantile = normal_standard_quantile(1 - alpha_/2);            
            lowerBound = (C_ * beta_).array() - quantile * (diagon.array()).sqrt();
            upperBound = (C_ * beta_).array() + quantile * (diagon.array()).sqrt();
         }
         else{
            // inserire errore: nome intervallo non valido
            return DMatrix<double>::Zero(1, 1);
         }

         DMatrix<double> CIMatrix(p, 2);      //matrix of confidence intervals
         CIMatrix.col(0) = lowerBound;
         CIMatrix.col(1) = upperBound;
         return CIMatrix;
      }

      virtual DVector<double> p_value(CIType type){

         fdapde_assert(!is_empty(C_));      
         if(is_empty(beta0_)){
            if(is_empty(beta_)){
               beta();
            }
            setBeta0(DVector<double>::Zero(beta_.size())); 
         }
         if(is_empty(beta_)){
            beta();
         }
         if(is_empty(V_)){
            V();
         }
         int p = C_.rows();
         DVector<double> statistics(p);         
         if(type == simultaneous){
            // SIMULTANEOUS 
            DVector<double> diff = C_ * beta_ - beta0_;          
            DMatrix<double> Sigma = C_ * V_ * C_.transpose();
            DMatrix<double> Sigmadec_ = inverse(Sigma);
            double stat = diff.adjoint() * Sigmadec_ * diff;            
            statistics.resize(p);
            double pvalue = chi_squared_cdf(stat, p);
            if(pvalue < 0){ 
               statistics(0) = 1;
            }
            if(pvalue > 1){
               statistics(0) = 0;
            }
            else{
               statistics(0) = 1 - pvalue;
            }
            for(int i = 1; i < C_.rows(); i++){
               statistics(i) = 10e20;
            }
            return statistics; 
         }

         else if (type == one_at_the_time){
            //  ONE AT THE TIME 
            int p = C_.rows();
            statistics.resize(p);
            for(int i = 0; i < p; i++){
               DVector<double> col = C_.row(i);
               double diff = col.adjoint() * beta_ - beta0_[i];
               double sigma = col.adjoint() * V_ *col;
               double stat = diff/std::sqrt(sigma);
               double pvalue = 2 * gaussian_cdf(-std::abs(stat), 0, 1);
               if(pvalue < 0){
                  statistics(i) = 0;
               }
               if(pvalue > 1){
                  statistics(i) = 1;
               }
               else{
                  statistics(i) = pvalue;
               }
            }
            return statistics;
         }
         else{
            //inserire messaggio di errore 
            return DVector<double>::Zero(1);
         }
      }

      // return the sparse approx of E^{-1}
      static SpMatrix<double> invE_approx(const Model& m){
       /* SpMatrix<double> decR0_ = lump(m.R0()); 
        DiagMatrix<double> invR0_(decR0_.rows());
        invR0_.setZero(); 
        for (int i = 0; i < decR0_.rows(); ++i) {
            double diagElement = decR0_.diagonal()[i];  
            invR0_.diagonal()[i] = 1.0 / diagElement; 
        }*/
        //DMatrix<double> Et_ = m.PsiTD()* m.Psi()+ m.lambda_D() * m.R1().transpose() * invR0_ * m.R1();

        //applico FSPAI su Atilde
        int alpha = 10;    // Numero di aggiornamenti del pattern di sparsità per ogni colonna di A (perform alpha steps of approximate inverse update along column k)
        int beta = 10;      // Numero di indici da aggiungere al pattern di sparsità di Lk per ogni passo di aggiornamento
        double epsilon = 0.005; // Soglia di tolleranza per l'aggiornamento del pattern di sparsità (the best improvement is higher than accetable treshold)
            //questi sono quelli trovati nella libreria vecchia 
            //std::string tol_Inverse     = "0.005";  oppure 0.05                     // Controls the quality of approximation, default 0.005 
            //std::string max_Step_Col    = "20";     oppure 10                     // Max number of improvement steps per columns
            // std::string max_New_Nz      = "20";     oppure 10                     // Max number of new nonzero candidates per step
            //Et_ should be stored as a sparse matrix 

       // questo serve se voglio invertire R0 cojn fspai           
        FSPAI fspai_R0(m.R0());
        fspai_R0.compute(m.R0(),alpha, beta, epsilon);
        SpMatrix<double> invR0_ = fspai_R0.inverse();  

        DMatrix<double> Et_ = m.PsiTD()* m.Psi()+ m.lambda_D() * m.R1().transpose() * invR0_ * m.R1();


        
        SpMatrix<double> Et_sparse = Et_.sparseView();
          
        //Eigen::saveMarket(Et_sparse, "Edainvertire.mtx");  
        FSPAI fspai_E(Et_sparse);
        fspai_E.compute(Et_sparse,alpha, beta, epsilon);
        SpMatrix<double> invE_ = fspai_E.inverse();
        //Eigen::saveMarket(invE_, "inversaE2.mtx");  
        //SpMatrix<double> precondE= fspai_E.getL();
        //Eigen::saveMarket(precondE, "precondE.mtx");
        
        /*
        //CALCOLO DELL'INVERSA PER CONFRONTI CON LIBRERIA ORIGINALE DI FSPAI 
        //SpMatrix<double> Edainvertire;
        //Eigen::loadMarket(Edainvertire, "Edainvertire.mtx");
        //std::cout<<"righe di Edainveritre: "<<Edainvertire.rows()<<std::endl;

        FSPAI fspai_Edainvertire(Edainvertire);
        fspai_Edainvertire.compute(alpha, beta, epsilon);
        SpMatrix<double> precond_nostro = fspai_Edainvertire.getL();
        Eigen::saveMarket(precond_nostro, "precond_nostro.mtx");  
*/




        //SpMatrix<double> risultatoFSPAI;
        //Eigen::loadMarket(risultatoFSPAI, "risultatoFSPAI.mtx");
        //std::cout<<"righe di fspai"<<risultatoFSPAI.rows()<<std::endl;
        //return risultatoFSPAI;

        return invE_;  
      }

     
      // setter for matrix of combination of coefficients C
      void setC(const DMatrix<double>& C) {
         C_ = C;
      }

      // setter for alpha
      void setAlpha(double alpha){
         fdapde_assert(0 <= alpha && alpha <= 1);      // throw an exception if condition is not met  
         if( 0 <= alpha && alpha <= 1) {
            alpha_ = alpha;
         }
      }

      // setter for alpha f
      void setAlpha_f(double alpha){
         fdapde_assert(0 <= alpha&& alpha <= 1);      // throw an exception if condition is not met  
         if( 0 <= alpha && alpha <= 1) {
            alpha_f_ = alpha;
         }
      }      

      // setter for beta0_ 
      void setBeta0(DVector<double> beta0){
         beta0_ = beta0;
      }

      // setter for f0_
      void setf0(DVector<double> f0){ 
         f0_ = f0;
      }

      void setLocationsF(DVector<int> locs){
         locations_f_ = locs;
      }


};

} // namespace models
} // namespace fdapde

#endif   // __INFERENCE_BASE_H__
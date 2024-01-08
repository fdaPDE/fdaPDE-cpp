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


#include <fdaPDE/core.h>
#include <gtest/gtest.h>   // testing framework

#include <cstddef>
using fdapde::core::advection;
using fdapde::core::diffusion;
using fdapde::core::fem_order;
using fdapde::core::FEM;
using fdapde::core::Grid;
using fdapde::core::Mesh; 
using fdapde::core::laplacian;
using fdapde::core::bilaplacian;
using fdapde::core::SPLINE;
using fdapde::core::spline_order;
using fdapde::core::MatrixDataWrapper;
using fdapde::core::PDE;
using fdapde::core::VectorDataWrapper;

#include "../../fdaPDE/models/sampling_design.h"
#include "../../fdaPDE/models/regression/qsrpde.h"
using fdapde::models::QSRPDE;
using fdapde::models::SpaceTimeSeparable;
using fdapde::models::SpaceTimeParabolic;

#include "../../fdaPDE/models/regression/gcv.h"
using fdapde::models::ExactEDF;
using fdapde::models::GCV;
using fdapde::models::StochasticEDF;
using fdapde::models::Sampling;

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_mtx;
using fdapde::testing::read_csv;


// for time and memory performances
#include <chrono>
#include <iomanip>
using namespace std::chrono;
#include <unistd.h>
#include <fstream>


// /* test 1 SQRPDE - Time
//    domain:       unit square [0,1] x [0,1] (coarse)
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
//  */
// TEST(gcv_sqrpde_time_test, laplacian_nonparametric_samplingatnodes_spacetimeseparable_gridexact) {
  
  
//   // Parameters 
//   const std::string TestNumber = "1"; 
  
//   std::vector<double> alphas = {0.1, 0.5, 0.9}; 

//   // number of simulations 
//   unsigned int n_sim = 5; 

//   // Choose missing strategy and proportion
//   std::string missing = "_missing"; 
//   // std::string missing = ""; 

//   const std::string p_string = "/p_50";

//   // Marco
//   // std::string path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/PACS_project_shared"; 
//   // Ilenia 
//   std::string path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared"; 

//   std::string path_test = path + "/models/space_time/Test_" + TestNumber + p_string;


//   // define time domain  
//   DVector<double> time_mesh;
//   time_mesh.resize(6);
//   std::size_t i = 0;
//   for(double x = 0; x <= 2; x+=0.4, ++i) time_mesh[i] = x;
  
//   // define spatial domain 
//   MeshLoader<Mesh2D> domain("unit_square_coarse");
//   // define regularizing PDE    
//   auto L = -laplacian<FEM>();
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.rows(), 1);
//   PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
//   for(double alpha : alphas){

//     unsigned int alpha_int = alpha*100; 

//     std::cout << "------------------------------------------alpha=" << std::to_string(alpha_int) << "%-------------------------------------------------" << std::endl; 

  
//     // define statistical model
//     SQRPDE<decltype(problem), SpaceTimeSeparable, GeoStatMeshNodes, MonolithicSolver> model(problem, time_mesh, alpha);


//     for(unsigned int sim = 1; sim <= n_sim; ++sim){ 

//       // import data from files
//       DMatrix<double> y = read_csv<double>(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/y" + missing + ".csv");

//       // set model data
//       BlockFrame<double, int> df;
//       df.stack(OBSERVATIONS_BLK, y);
//       model.set_data(df);
//       model.init(); 

//       // define GCV function and grid of \lambda_D values
//       GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
//       ScalarField<2, decltype(GCV)> obj(GCV);
//       std::vector<SVector<2>> lambdas_d_t;

//       if(std::abs(alpha - 0.5) < 0.5){
//       for (double x_s = -5.0; x_s <= -3.9; x_s +=0.20){
//         for (double x_t = -7.0; x_t <= -3.9; x_t +=1.00) lambdas_d_t.push_back(SVector<2>(std::pow(10, x_s), std::pow(10,x_t)));
//       }
//       }
//       if(std::abs(alpha - 0.1) < 0.5){
//         for (double x_s = -5.8; x_s <= -4.5; x_s +=0.20){
//           for (double x_t = -7.0; x_t <= -4.9; x_t +=1.00) lambdas_d_t.push_back(SVector<2>(std::pow(10, x_s), std::pow(10,x_t)));
//       }
//       }
//       else{
//         for (double x_s = -6.2; x_s <= -4.9; x_s +=0.20){
//           for (double x_t = -7.0; x_t <= -3.9; x_t +=1.00) lambdas_d_t.push_back(SVector<2>(std::pow(10, x_s), std::pow(10,x_t)));
//       }
//       }


//       // optimize GCV
//       Grid<2> opt;
//       opt.optimize(obj, lambdas_d_t);
//       SVector<2> best_lambda = opt.optimum();
      
//       std::cout << "Best lambda = " << best_lambda << std::endl ; 

//       // Save Lambda opt
//       std::ofstream fileLambdaoptS(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/GCV/Exact/LambdaS" + missing + ".csv");
//       if (fileLambdaoptS.is_open()){
//       fileLambdaoptS << std::setprecision(16) << best_lambda[0];
//       fileLambdaoptS.close();
//       }

//       std::ofstream fileLambdaoptT(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/GCV/Exact/LambdaT" + missing + ".csv");
//       if (fileLambdaoptT.is_open()){
//       fileLambdaoptT << std::setprecision(16) << best_lambda[1];
//       fileLambdaoptT.close();
//       }
           
      
//       // Save GCV scores
//       std::ofstream fileGCV_scores(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/GCV/Exact/GCV_scores" + missing + ".csv");
//       for(std::size_t i = 0; i < GCV.values().size(); ++i) 
//         fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.values()[i]) << "\n" ; 

//       fileGCV_scores.close(); 

//     }

//   }


// }


// /* test 3 
//    domain:            c-shaped
//    space sampling:    locations != nodes
//    time sampling:     locations = nodes
//    penalization:      simple laplacian
//    missing:           no
//    covariates:        yes
//    BC:                no
//    order FE:          1
//    GCV optimization:  grid exact
//    time penalization: separable (mass penalization)
//  */
// TEST(sqrpde_time_test, laplacian_semiparametric_samplingatlocations_timelocations_separable_monolithic) {

//     // Marco 
//     std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/space_time/Test_2"; 
//     //   // Ilenia 
//     //   std::string path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/space_time/Test_2"; 

//     std::vector<double> alphas = {0.1, 0.5, 0.9}; 

//     // define temporal domain
//     DVector<double> time_mesh;
//     time_mesh.resize(5);
//     for (std::size_t i = 0; i < 5; ++i) time_mesh[i] = (fdapde::testing::pi / 4) * i;
//     // define spatial domain and regularizing PDE
//     MeshLoader<Mesh2D> domain("c_shaped");
//     // import data from files
//     DMatrix<double> space_locs = read_csv<double>(R_path + "/locs.csv");
//     DMatrix<double> X = read_csv<double>(R_path + "/X.csv");

//     // define regularizing PDE
//     auto L = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.rows(), 1);
//     PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

//     // lambdas sequence 
//     std::vector<SVector<2>> lambdas_d_t;
//     for(double x = -6.0; x <= -3.0; x +=0.5) lambdas_d_t.push_back(SVector<2>(std::pow(10,x), std::pow(10,x)));

//     unsigned int n_sim = 10; 
//     for(unsigned int sim = 1; sim <= n_sim; ++sim){
//       std::cout << "---------------------------Simulation #" << sim << "--------------------------" << std::endl; 
//       for(double alpha : alphas){

//         unsigned int alpha_int = alpha*100; 
//         std::string alpha_string = std::to_string(alpha_int); 

//         SQRPDE<decltype(problem), SpaceTimeSeparable, GeoStatLocations, MonolithicSolver> model(problem, time_mesh, alpha);

//         // load data from .csv files
//         DMatrix<double> y = read_csv<double>(R_path + "/sim_" + std::to_string(sim) + "/y.csv"); 

//         // set model's data
//         model.set_spatial_locations(space_locs);
        
//         BlockFrame<double, int> df;
//         df.stack(OBSERVATIONS_BLK, y);
//         df.stack(DESIGN_MATRIX_BLK, X);
//         model.set_data(df);
//         model.init();

//         // define GCV function and grid of \lambda_D values
//         GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
//         ScalarField<2, decltype(GCV)> obj(GCV);  
//         // optimize GCV
//         Grid<2> opt;
//         opt.optimize(obj, lambdas_d_t);
//         SVector<2> best_lambda = opt.optimum();

//         // Save Lambda opt
//         std::ofstream fileLambdaoptS(R_path + "/sim_" + std::to_string(sim) + "/alpha_" + alpha_string  + "/lambda_opt.csv");
//         if (fileLambdaoptS.is_open()){
//           fileLambdaoptS << std::setprecision(16) << best_lambda[0];
//           fileLambdaoptS.close();
//         }
//         // Save GCV scores
//         std::ofstream fileGCV_scores(R_path + "/sim_" + std::to_string(sim) + "/alpha_" + alpha_string  + "/gcv_scores.csv");
//         for(std::size_t i = 0; i < GCV.values().size(); ++i) 
//           fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.values()[i]) << "\n" ; 

//         fileGCV_scores.close();

//       }
//     }
// }


  
/* test 5 
   domain:       c-shaped
   space sampling: locations != nodes
   time sampling:  locations != nodes
   penalization: simple laplacian
   missing:      yes
   covariates:   no
   BC:           no
   order FE:     1
   GCV optimization: grid exact
   time penalization: separable (mass penalization)
 */
TEST(sqrpde_time_test, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

    // Marco 
    std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/space_time/Test_3"; 
    //   // Ilenia 
    //   std::string path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/models/space_time/Test_3"; 

    std::vector<double> alphas = {0.1, 0.5, 0.9}; 
    //std::string data_type = "all";   // all d
    std::vector<std::string> data_types = {"all"}; 
    std::string p_string = "50";   

    // define temporal domain
    unsigned int M = 3; 
    std::string M_string = std::to_string(M);
    double tf = fdapde::testing::pi;   // final time 
    Mesh<1, 1> time_mesh(0, tf, M-1);
    // define spatial domain and regularizing PDE
    MeshLoader<Mesh2D> domain("c_shaped_adj");
    // MeshLoader<Mesh2D> domain("c_shaped_504");   // mesh fine 

    // import locs from files
    DMatrix<double> space_locs = read_csv<double>(R_path + "/space_locs.csv");
    DMatrix<double> time_locs = read_csv<double>(R_path + "/time_locs.csv");

    // define regularizing PDE in space 
    auto Ld = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
    PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

    // define regularizing PDE in time
    auto Lt = -bilaplacian<SPLINE>();
    PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

    // lambdas sequence 
    std::vector<DVector<double>> lambdas_d_t;

    unsigned int n_sim = 5; 
    std::string eps_string = ""; 

    for(auto data_type : data_types){
      if(data_type == "all")
        std::cout << "--------------------------------------ALL DATA----------------------------" << std::endl; 
      else 
        std::cout << "---------------------------------------MISSING DATA----------------------------" << std::endl;

      std::vector<DVector<double>> lambdas10_d_t;
      std::vector<DVector<double>> lambdas50_d_t;
      std::vector<DVector<double>> lambdas90_d_t;
      if(data_type == "all"){
        // 10% 
        for(double xs = -3.6; xs <= -1.8; xs +=0.05)
          for(double xt = -7.0; xt <= -6.0; xt +=1.0) 
            lambdas10_d_t.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));
        // 50% 
        for(double xs = -3.2; xs <= -1.8; xs +=0.05)
          for(double xt = -7.0; xt <= -6.0; xt +=1.0) 
            lambdas50_d_t.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));
        // 90%
        for(double xs = -3.8; xs <= -1.8; xs +=0.05)
          for(double xt = -7.0; xt <= -6.0; xt +=1.0) 
            lambdas90_d_t.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));
      }
      if(data_type == "d"){
        // 10% 
        for(double xs = -5.0; xs <= -3.0; xs +=0.2)
          for(double xt = -9.0; xt <= -7.0; xt +=0.5) 
            lambdas10_d_t.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));
        // 50% 
        for(double xs = -6.0; xs <= -4.0; xs +=0.2)
          for(double xt = -11.0; xt <= -7.0; xt +=0.5) 
            lambdas50_d_t.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));
        // 90%
        for(double xs = -5.0; xs <= -2.0; xs +=0.2)
          for(double xt = -9.0; xt <= -7.0; xt +=0.5) 
            lambdas90_d_t.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));
      }
   
      // simulations 
      for(unsigned int sim = 1; sim <= n_sim; ++sim){
        std::cout << "---------------Simulation #" << sim << "--------------" << std::endl; 
        for(double alpha : alphas){

          unsigned int alpha_int = alpha*100; 
          std::string alpha_string = std::to_string(alpha_int);
          std::cout << "--------alpha=" << alpha_string << "%" << std::endl;  

          QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);

          // load data from .csv files
          DMatrix<double> y; 
          if(data_type == "all")
            y = read_csv<double>(R_path + "/simulations/all/sim_" + std::to_string(sim) + "/y_all.csv");
          else
            y = read_csv<double>(R_path + "/simulations/miss_strategy_" + data_type + "/p_" + p_string + "/sim_" + std::to_string(sim) + "/y.csv");

          if(almost_equal(alpha, 0.1)){
            lambdas_d_t = lambdas10_d_t; 
          }  
          if(almost_equal(alpha, 0.5)){
            lambdas_d_t = lambdas50_d_t; 
          }  
          if(almost_equal(alpha, 0.9)){
            lambdas_d_t = lambdas90_d_t; 
          }  

          // set model's data
          model.set_spatial_locations(space_locs);
          model.set_temporal_locations(time_locs);
          
          BlockFrame<double, int> df;
          df.stack(OBSERVATIONS_BLK, y);
          model.set_data(df);
          model.init();

          // define GCV function and grid of \lambda_D values
          auto GCV = model.gcv<ExactEDF>();
          // optimize GCV
          Grid<fdapde::Dynamic> opt;
          opt.optimize(GCV, lambdas_d_t);
          SVector<2> best_lambda = opt.optimum();

          std::string solutions_path; 
          if(data_type == "all")
            solutions_path = R_path + "/simulations/all/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + "/M_" + M_string; 
          else
            solutions_path = R_path + "/simulations/miss_strategy_" + data_type + "/p_" + p_string + "/sim_" + std::to_string(sim) + "/alpha_" + alpha_string;

          // Save Lambda opt
          std::ofstream fileLambdaoptS(solutions_path + "/lambda_s_opt.csv");
          if(fileLambdaoptS.is_open()){
            fileLambdaoptS << std::setprecision(16) << best_lambda[0];
            fileLambdaoptS.close();
          }
          std::ofstream fileLambdaoptT(solutions_path + "/lambda_t_opt.csv");
          if (fileLambdaoptT.is_open()){
            fileLambdaoptT << std::setprecision(16) << best_lambda[1];
            fileLambdaoptT.close();
          }
          // Save GCV scores
          std::ofstream fileGCV_scores(solutions_path + "/gcv_scores.csv");
          for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
            fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.gcvs()[i]) << "\n" ; 

          fileGCV_scores.close();

        }
      }
      }

}
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
using fdapde::core::PDE;

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


// Test ufficiale 

// Test NON PARAMETRICO
/* test 3 
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
TEST(sqrpde_time_test_gcv, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

    // Marco 
    std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/space_time/Test_3"; 
    // // Ilenia 
    // std::string R_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/models/space_time/Test_3"; 

    std::vector<double> alphas =  {0.1}; // , 0.5, 0.9};  
    std::vector<std::string> data_types = {"d"};  // "all" per il test senza missing
    // std::string p_string = "70";  // 50
    std::vector<std::string> p_string_vec = {"50"};
    // std::string lambda_selection_type = "gcv_smooth_eps1e-1";   //  +0
    std::vector<std::string> lambda_selection_type_vec = {"gcv_smooth_eps1e-1"}; // "gcv_smooth_eps1e-3", "gcv_smooth_eps1e-1", "gcv_smooth_eps1e+0"
    bool check_for_merge = true; 
    bool no_init = false;
    bool exactGCV = false; 
    bool no_switch = false;
    bool artificial_weights = true; 

    bool save = false; 

    // for stochasticEDF
    std::size_t seed = 438172;
    unsigned int MC_run = 100; 

    // define temporal domain
    double tf = fdapde::testing::pi;   // final time 
  
    unsigned int M = 7;              // CAMBIA !
    std::string path_mesh = "/M_7";  // ATT : da cambiare
    Mesh<1, 1> time_mesh(0, tf, M-1);     // t0, tf, #subintervals   

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

    unsigned int n_sim = 1; 
    std::string eps_string = ""; 

    for(auto data_type : data_types){
      if(data_type == "all")
        std::cout << "--------------------------------------ALL DATA----------------------------" << std::endl; 
      else 
        std::cout << "---------------------------------------MISSING DATA----------------------------" << std::endl;


      // lambdas sequence 
      std::vector<DVector<double>> lambdas_d_t; std::vector<double> lambdas_d; std::vector<double> lambdas_t;
      std::vector<double> lambdas10_d; std::vector<double> lambdas10_t; std::vector<DVector<double>> lambdas10_d_t;
      std::vector<double> lambdas50_d; std::vector<double> lambdas50_t; std::vector<DVector<double>> lambdas50_d_t; 
      std::vector<double> lambdas90_d; std::vector<double> lambdas90_t; std::vector<DVector<double>> lambdas90_d_t; 


      if(data_type == "all"){

        // 10% 
        {
            for(double xs = -4.0; xs <= -1.6; xs += 0.05)
                lambdas10_d.push_back(std::pow(10,xs));

            for(double xt = -5.0; xt <= -4.0; xt += 2.0)
                lambdas10_t.push_back(std::pow(10,xt));    

            for(auto i = 0; i < lambdas10_d.size(); ++i)
                for(auto j = 0; j < lambdas10_t.size(); ++j) 
                    lambdas10_d_t.push_back(SVector<2>(lambdas10_d[i], lambdas10_t[j]));
        }
        
        // 50% 
        {
            for(double xs = -4.0; xs <= -1.1; xs += 0.05)
                lambdas50_d.push_back(std::pow(10,xs));

            for(double xt = -5.0; xt <= -4.0; xt += 2.0)
                lambdas50_t.push_back(std::pow(10,xt)); 

            for(auto i = 0; i < lambdas50_d.size(); ++i)
                for(auto j = 0; j < lambdas50_t.size(); ++j) 
                    lambdas50_d_t.push_back(SVector<2>(lambdas50_d[i], lambdas50_t[j]));
        }

        // 90% 
        {
            for(double xs = -4.2; xs <= -1.4; xs += 0.05)
                lambdas90_d.push_back(std::pow(10,xs));

            for(double xt = -5.0; xt <= -4.0; xt += 2.0)
                lambdas90_t.push_back(std::pow(10,xt));    

            for(auto i = 0; i < lambdas90_d.size(); ++i)
                for(auto j = 0; j < lambdas90_t.size(); ++j) 
                    lambdas90_d_t.push_back(SVector<2>(lambdas90_d[i], lambdas90_t[j]));
        }
 
      }
      if(data_type == "a"){
        
        // 10% 
        {
            for(double xs = -4.0; xs <= -1.6; xs += 1.0)
                lambdas10_d.push_back(std::pow(10,xs));

            for(double xt = -5.0; xt <= -4.0; xt += 2.0)
                lambdas10_t.push_back(std::pow(10,xt));    

            for(auto i = 0; i < lambdas10_d.size(); ++i)
                for(auto j = 0; j < lambdas10_t.size(); ++j) 
                    lambdas10_d_t.push_back(SVector<2>(lambdas10_d[i], lambdas10_t[j]));
        }
        
        // 50% 
        {
            for(double xs = -4.0; xs <= -1.6; xs += 0.05)
                lambdas50_d.push_back(std::pow(10,xs));

            for(double xt = -5.0; xt <= -4.0; xt += 4.0)
                lambdas50_t.push_back(std::pow(10,xt)); 

            // lambdas50_t.push_back(std::pow(10,-5)); 
            // lambdas50_t.push_back(std::pow(10,-1)); 
            // lambdas50_t.push_back(std::pow(10,+1));    

            for(auto i = 0; i < lambdas50_d.size(); ++i)
                for(auto j = 0; j < lambdas50_t.size(); ++j) 
                    lambdas50_d_t.push_back(SVector<2>(lambdas50_d[i], lambdas50_t[j]));
        }

        // 90% 
        {
            for(double xs = -4.2; xs <= -1.6; xs += 0.05)
                lambdas90_d.push_back(std::pow(10,xs));

            for(double xt = -5.0; xt <= -4.0; xt += 2.0)
                lambdas90_t.push_back(std::pow(10,xt));    

            for(auto i = 0; i < lambdas90_d.size(); ++i)
                for(auto j = 0; j < lambdas90_t.size(); ++j) 
                    lambdas90_d_t.push_back(SVector<2>(lambdas90_d[i], lambdas90_t[j]));
        }
      }
        
      if(data_type == "d"){
        
        // 10% 
        {
            for(double xs = -4.0; xs <= -1.6; xs += 0.05)
                lambdas10_d.push_back(std::pow(10,xs));

            for(double xt = -5.0; xt <= -4.0; xt += 2.0)
                lambdas10_t.push_back(std::pow(10,xt));    

            for(auto i = 0; i < lambdas10_d.size(); ++i)
                for(auto j = 0; j < lambdas10_t.size(); ++j) 
                    lambdas10_d_t.push_back(SVector<2>(lambdas10_d[i], lambdas10_t[j]));
        }
        
        // 50% 
        {
            for(double xs = -4.0; xs <= -1.1; xs += 0.05)
                lambdas50_d.push_back(std::pow(10,xs));

            for(double xt = -5.0; xt <= -4.0; xt += 2.0)
                lambdas50_t.push_back(std::pow(10,xt)); 

            for(auto i = 0; i < lambdas50_d.size(); ++i)
                for(auto j = 0; j < lambdas50_t.size(); ++j) 
                    lambdas50_d_t.push_back(SVector<2>(lambdas50_d[i], lambdas50_t[j]));
        }

        // 90% 
        {
            for(double xs = -4.2; xs <= -1.4; xs += 0.05)
                lambdas90_d.push_back(std::pow(10,xs));

            for(double xt = -5.0; xt <= -4.0; xt += 2.0)
                lambdas90_t.push_back(std::pow(10,xt));    

            for(auto i = 0; i < lambdas90_d.size(); ++i)
                for(auto j = 0; j < lambdas90_t.size(); ++j) 
                    lambdas90_d_t.push_back(SVector<2>(lambdas90_d[i], lambdas90_t[j]));
        }
      }
   
      // simulations 
      for(auto p_string : p_string_vec){
        for(auto lambda_selection_type : lambda_selection_type_vec){

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
                y = read_csv<double>(R_path + "/simulations/all_new/sim_" + std::to_string(sim) + "/y_all.csv");  // ATT: MESSO NEW
              else
                y = read_csv<double>(R_path + "/simulations/miss_strategy_" + data_type + "/p_" + p_string + "/sim_" + std::to_string(sim) + "/y.csv");

              if(alpha_string == "10"){
                    lambdas_d = lambdas10_d; 
                    lambdas_t = lambdas10_t;
                    lambdas_d_t = lambdas10_d_t;
                }
              if(alpha_string == "50"){
                    lambdas_d = lambdas50_d; 
                    lambdas_t = lambdas50_t;
                    lambdas_d_t = lambdas50_d_t;
                }  
              if(alpha_string == "90"){
                    lambdas_d = lambdas90_d; 
                    lambdas_t = lambdas90_t;
                    lambdas_d_t = lambdas90_d_t;
                }

              // set model's data
              model.set_spatial_locations(space_locs);
              model.set_temporal_locations(time_locs);

              if(lambda_selection_type == "gcv_smooth_eps1e-3")
                model.set_eps_power(-3.0);
              if(lambda_selection_type == "gcv_smooth_eps1e-2")
                model.set_eps_power(-2.0);
              if(lambda_selection_type == "gcv_smooth_eps1e-1.5")
                model.set_eps_power(-1.5);
              if(lambda_selection_type == "gcv_smooth_eps1e-1")
                model.set_eps_power(-1.0);
              if(lambda_selection_type == "gcv_smooth_eps1e-0.5")
                model.set_eps_power(-0.5);
              if(lambda_selection_type == "gcv_smooth_eps1e+0")
                model.set_eps_power(0.0);
              if(lambda_selection_type == "gcv_smooth_eps1e+0.5")
                model.set_eps_power(0.5);
              if(lambda_selection_type == "gcv_smooth_eps1e+1")
                model.set_eps_power(1.0);

              std::string solutions_path; 
              if(data_type == "all")
                solutions_path = R_path + "/simulations/all_new/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + path_mesh + "/" + lambda_selection_type; 
              else
                solutions_path = R_path + "/simulations/miss_strategy_" + data_type + "/p_" + p_string + "/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + 
                                path_mesh + "/" + lambda_selection_type;

              std::cout << "Solution path: " << solutions_path << std::endl ; 

              BlockFrame<double, int> df;
              df.stack(OBSERVATIONS_BLK, y);
              model.set_data(df);
              model.init();

              // // exact
              // auto GCV = model.gcv<ExactEDF>();

              // stochastic
              auto GCV = model.gcv<StochasticEDF>(MC_run, seed);

              // optimize GCV
              Grid<fdapde::Dynamic> opt;
              opt.optimize(GCV, lambdas_d_t);
              SVector<2> best_lambda = opt.optimum();

              // if(check_for_merge){
              //   if(no_init){
              //     solutions_path += "/check_merge_no_init";
              //   } else{
              //     if(exactGCV){
              //       solutions_path += "/check_merge_exactGCV"; 
              //     }
              //     if(no_switch){
              //       solutions_path += "/check_merge_no_switch"; 
              //     }
              //     if(!no_switch && !exactGCV){
              //       solutions_path += "/check_merge"; 
              //     }
                    
              //   }
                
              // }

              if(check_for_merge){
                if(artificial_weights){
                  solutions_path += "/check_merge_weights";
                } else{
                  solutions_path += "/check_merge";
                }
                 
              }

              if(save){
                // Save lambda sequence 
                std::ofstream fileLambda_S_Seq(solutions_path + "/lambdas_S_seq.csv");
                for(std::size_t i = 0; i < lambdas_d.size(); ++i) 
                    fileLambda_S_Seq << std::setprecision(16) << lambdas_d[i] << "\n"; 
                fileLambda_S_Seq.close();

                std::ofstream fileLambda_T_Seq(solutions_path + "/lambdas_T_seq.csv");
                for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
                    fileLambda_T_Seq << std::setprecision(16) << lambdas_t[i] << "\n"; 
                fileLambda_T_Seq.close();

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
      }

    }

}





// // Test SEMI PARAMETRICO
// /* test 3 
//    domain:       c-shaped
//    space sampling: locations != nodes
//    time sampling:  locations != nodes
//    penalization: simple laplacian
//    missing:      no
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
//    time penalization: separable (mass penalization)
//  */
// TEST(sqrpde_time_test_gcv_semiparam, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

//     // Marco 
//     // std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/space_time/Test_3"; 
//     // Ilenia 
//     std::string R_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/models/space_time/Test_3"; 

//     std::vector<double> alphas = {0.1};   // , 0.5, 0.9};  
//     std::vector<std::string> data_types = {"all"};  // "all" per il test senza missing
//     // std::string p_string = "70";  // 50
//     std::vector<std::string> p_string_vec = {"50"};
//     // std::string lambda_selection_type = "gcv_smooth_eps1e-1";   //  +0
//     std::vector<std::string> lambda_selection_type_vec = {"gcv_smooth_eps1e-1"}; // "gcv_smooth_eps1e-3", "gcv_smooth_eps1e-1", "gcv_smooth_eps1e+0"

//     // for stochasticEDF
//     std::size_t seed = 438172;
//     unsigned int MC_run = 100; 

//     // define temporal domain
//     double tf = fdapde::testing::pi;   // final time 
  
//     unsigned int M = 7;              // CAMBIA !
//     std::string path_mesh = "/M_7";  // ATT : da cambiare
//     Mesh<1, 1> time_mesh(0, tf, M-1);     // t0, tf, #subintervals   

//     // define spatial domain and regularizing PDE
//     MeshLoader<Mesh2D> domain("c_shaped_adj");
//     // MeshLoader<Mesh2D> domain("c_shaped_504");   // mesh fine 

//     // import locs from files
//     DMatrix<double> space_locs = read_csv<double>(R_path + "/space_locs.csv");
//     DMatrix<double> time_locs = read_csv<double>(R_path + "/time_locs.csv");

//     // import covariate (they do not depend on the sim)
//     DMatrix<double> X = read_csv<double>(R_path + "/simulations/all_covariate/X.csv");
//     std::cout << "Dimension X = ( " << X.rows() << " , " << X.cols() << " )" << std::endl ; 

//     // define regularizing PDE in space 
//     auto Ld = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
//     PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

//     // define regularizing PDE in time
//     auto Lt = -bilaplacian<SPLINE>();
//     PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

//     // lambdas sequence 
//     std::vector<DVector<double>> lambdas_d_t;

//     unsigned int n_sim = 1; 
//     std::string eps_string = ""; 

//     for(auto data_type : data_types){
//       if(data_type == "all")
//         std::cout << "--------------------------------------ALL DATA----------------------------" << std::endl; 
//       else 
//         std::cout << "---------------------------------------MISSING DATA----------------------------" << std::endl;


//       // lambdas sequence 
//       std::vector<DVector<double>> lambdas_d_t; std::vector<double> lambdas_d; std::vector<double> lambdas_t;
//       std::vector<double> lambdas10_d; std::vector<double> lambdas10_t; std::vector<DVector<double>> lambdas10_d_t;
//       std::vector<double> lambdas50_d; std::vector<double> lambdas50_t; std::vector<DVector<double>> lambdas50_d_t; 
//       std::vector<double> lambdas90_d; std::vector<double> lambdas90_t; std::vector<DVector<double>> lambdas90_d_t; 


//       if(data_type == "all"){
//                 // 10% 
//         {
//             for(double xs = -4.5; xs <= -1.5; xs += 0.05)
//                 lambdas10_d.push_back(std::pow(10,xs));

//             for(double xt = -5.0; xt <= -4.0; xt += 2.0)
//                 lambdas10_t.push_back(std::pow(10,xt));    

//             for(auto i = 0; i < lambdas10_d.size(); ++i)
//                 for(auto j = 0; j < lambdas10_t.size(); ++j) 
//                     lambdas10_d_t.push_back(SVector<2>(lambdas10_d[i], lambdas10_t[j]));
//         }
        
//         // 50% 
//         {
//             for(double xs = -4.5; xs <= -1.5; xs += 0.05)
//                 lambdas50_d.push_back(std::pow(10,xs));

//             for(double xt = -5.0; xt <= -4.0; xt += 2.0)
//                 lambdas50_t.push_back(std::pow(10,xt)); 

//             for(auto i = 0; i < lambdas50_d.size(); ++i)
//                 for(auto j = 0; j < lambdas50_t.size(); ++j) 
//                     lambdas50_d_t.push_back(SVector<2>(lambdas50_d[i], lambdas50_t[j]));
//         }

//         // 90% 
//         {
//             for(double xs = -4.5; xs <= -1.5; xs += 0.05)
//                 lambdas90_d.push_back(std::pow(10,xs));

//             for(double xt = -5.0; xt <= -4.0; xt += 2.0)
//                 lambdas90_t.push_back(std::pow(10,xt));    

//             for(auto i = 0; i < lambdas90_d.size(); ++i)
//                 for(auto j = 0; j < lambdas90_t.size(); ++j) 
//                     lambdas90_d_t.push_back(SVector<2>(lambdas90_d[i], lambdas90_t[j]));
//         }
//       }
//       if(data_type == "a"){
        
//         // 10% 
//         {
//             for(double xs = -4.0; xs <= -1.6; xs += 0.05)
//                 lambdas10_d.push_back(std::pow(10,xs));

//             for(double xt = -5.0; xt <= -4.0; xt += 2.0)
//                 lambdas10_t.push_back(std::pow(10,xt));    

//             for(auto i = 0; i < lambdas10_d.size(); ++i)
//                 for(auto j = 0; j < lambdas10_t.size(); ++j) 
//                     lambdas10_d_t.push_back(SVector<2>(lambdas10_d[i], lambdas10_t[j]));
//         }
        
//         // 50% 
//         {
//             for(double xs = -4.0; xs <= -1.6; xs += 0.05)
//                 lambdas50_d.push_back(std::pow(10,xs));

//             for(double xt = -5.0; xt <= -4.0; xt += 4.0)
//                 lambdas50_t.push_back(std::pow(10,xt)); 

//             for(auto i = 0; i < lambdas50_d.size(); ++i)
//                 for(auto j = 0; j < lambdas50_t.size(); ++j) 
//                     lambdas50_d_t.push_back(SVector<2>(lambdas50_d[i], lambdas50_t[j]));
//         }

//         // 90% 
//         {
//             for(double xs = -4.2; xs <= -1.6; xs += 0.05)
//                 lambdas90_d.push_back(std::pow(10,xs));

//             for(double xt = -5.0; xt <= -4.0; xt += 2.0)
//                 lambdas90_t.push_back(std::pow(10,xt));    

//             for(auto i = 0; i < lambdas90_d.size(); ++i)
//                 for(auto j = 0; j < lambdas90_t.size(); ++j) 
//                     lambdas90_d_t.push_back(SVector<2>(lambdas90_d[i], lambdas90_t[j]));
//         }
//       }  
//       if(data_type == "d"){
        
//         // 10% 
//         {
//             for(double xs = -4.0; xs <= -1.6; xs += 0.05)
//                 lambdas10_d.push_back(std::pow(10,xs));

//             for(double xt = -5.0; xt <= -4.0; xt += 2.0)
//                 lambdas10_t.push_back(std::pow(10,xt));    

//             for(auto i = 0; i < lambdas10_d.size(); ++i)
//                 for(auto j = 0; j < lambdas10_t.size(); ++j) 
//                     lambdas10_d_t.push_back(SVector<2>(lambdas10_d[i], lambdas10_t[j]));
//         }
        
//         // 50% 
//         {
//             for(double xs = -4.0; xs <= -1.1; xs += 0.05)
//                 lambdas50_d.push_back(std::pow(10,xs));

//             for(double xt = -5.0; xt <= -4.0; xt += 2.0)
//                 lambdas50_t.push_back(std::pow(10,xt)); 

//             for(auto i = 0; i < lambdas50_d.size(); ++i)
//                 for(auto j = 0; j < lambdas50_t.size(); ++j) 
//                     lambdas50_d_t.push_back(SVector<2>(lambdas50_d[i], lambdas50_t[j]));
//         }

//         // 90% 
//         {
//             for(double xs = -4.2; xs <= -1.4; xs += 0.05)
//                 lambdas90_d.push_back(std::pow(10,xs));

//             for(double xt = -5.0; xt <= -4.0; xt += 2.0)
//                 lambdas90_t.push_back(std::pow(10,xt));    

//             for(auto i = 0; i < lambdas90_d.size(); ++i)
//                 for(auto j = 0; j < lambdas90_t.size(); ++j) 
//                     lambdas90_d_t.push_back(SVector<2>(lambdas90_d[i], lambdas90_t[j]));
//         }
//       }
   
//       // simulations 
//       for(auto p_string : p_string_vec){
//         for(auto lambda_selection_type : lambda_selection_type_vec){

//           for(unsigned int sim = 1; sim <= n_sim; ++sim){
//           std::cout << "---------------Simulation #" << sim << "--------------" << std::endl; 

//             for(double alpha : alphas){
//           unsigned int alpha_int = alpha*100; 
//           std::string alpha_string = std::to_string(alpha_int);
//           std::cout << "--------alpha=" << alpha_string << "%" << std::endl;
//           QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);

//           // load data from .csv files
//           DMatrix<double> y; 
//           if(data_type == "all")
//             y = read_csv<double>(R_path + "/simulations/all_covariate/sim_" + std::to_string(sim) + "/y.csv");
//           else
//             y = read_csv<double>(R_path + "/simulations/covariate_miss_strategy_" + data_type + "/p_" + p_string + "/sim_" + std::to_string(sim) + "/y.csv");

//           if(alpha_string == "10"){
//                 lambdas_d = lambdas10_d; 
//                 lambdas_t = lambdas10_t;
//                 lambdas_d_t = lambdas10_d_t;
//             }
//           if(alpha_string == "50"){
//                 lambdas_d = lambdas50_d; 
//                 lambdas_t = lambdas50_t;
//                 lambdas_d_t = lambdas50_d_t;
//             }  
//           if(alpha_string == "90"){
//                 lambdas_d = lambdas90_d; 
//                 lambdas_t = lambdas90_t;
//                 lambdas_d_t = lambdas90_d_t;
//             }

//           // set model's data
//           model.set_spatial_locations(space_locs);
//           model.set_temporal_locations(time_locs);

//           model.set_exact_gcv(lambda_selection_type == "gcv");
//           if(lambda_selection_type == "gcv_smooth_eps1e-3")
//             model.set_eps_power(-3.0);
//           if(lambda_selection_type == "gcv_smooth_eps1e-2")
//             model.set_eps_power(-2.0);
//           if(lambda_selection_type == "gcv_smooth_eps1e-1.5")
//             model.set_eps_power(-1.5);
//           if(lambda_selection_type == "gcv_smooth_eps1e-1")
//             model.set_eps_power(-1.0);
//           if(lambda_selection_type == "gcv_smooth_eps1e-0.5")
//             model.set_eps_power(-0.5);
//           if(lambda_selection_type == "gcv_smooth_eps1e+0")
//             model.set_eps_power(0.0);
//           if(lambda_selection_type == "gcv_smooth_eps1e+0.5")
//             model.set_eps_power(0.5);
//           if(lambda_selection_type == "gcv_smooth_eps1e+1")
//             model.set_eps_power(1.0);

//           std::string solutions_path; 
//           if(data_type == "all")
//             solutions_path = R_path + "/simulations/all_covariate/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + path_mesh + "/" + lambda_selection_type; 
//           else
//             solutions_path = R_path + "/simulations/covariate_miss_strategy_" + data_type + "/p_" + p_string + "/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + 
//                              path_mesh + "/" + lambda_selection_type ;

//           std::cout << "Solution path: " << solutions_path << std::endl ; 

//           BlockFrame<double, int> df;
//           df.stack(OBSERVATIONS_BLK, y);
//           df.insert(DESIGN_MATRIX_BLK, X);
//           model.set_data(df);
//           model.init();

//           // // exact
//           // auto GCV = model.gcv<ExactEDF>();

//           // stochastic
//           auto GCV = model.gcv<StochasticEDF>(MC_run, seed);

//           // optimize GCV
//           Grid<fdapde::Dynamic> opt;
//           opt.optimize(GCV, lambdas_d_t);
//           SVector<2> best_lambda = opt.optimum();

//           // Save lambda sequence 
//           std::ofstream fileLambda_S_Seq(solutions_path + "/lambdas_S_seq.csv");
//           for(std::size_t i = 0; i < lambdas_d.size(); ++i) 
//               fileLambda_S_Seq << std::setprecision(16) << lambdas_d[i] << "\n"; 
//           fileLambda_S_Seq.close();

//           std::ofstream fileLambda_T_Seq(solutions_path + "/lambdas_T_seq.csv");
//           for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
//               fileLambda_T_Seq << std::setprecision(16) << lambdas_t[i] << "\n"; 
//           fileLambda_T_Seq.close();

//           // Save Lambda opt
//           std::ofstream fileLambdaoptS(solutions_path + "/lambda_s_opt.csv");
//           if(fileLambdaoptS.is_open()){
//             fileLambdaoptS << std::setprecision(16) << best_lambda[0];
//             fileLambdaoptS.close();
//           }
//           std::ofstream fileLambdaoptT(solutions_path + "/lambda_t_opt.csv");
//           if (fileLambdaoptT.is_open()){
//             fileLambdaoptT << std::setprecision(16) << best_lambda[1];
//             fileLambdaoptT.close();
//           }
//           // Save GCV scores
//           std::ofstream fileGCV_scores(solutions_path + "/gcv_scores.csv");
//           for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
//             fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.gcvs()[i]) << "\n" ; 

//           fileGCV_scores.close();
//         }
  
//           }
      
//         }
//       }

//     }

// }






// ------------ TEST PER PALU 

// // test 1
// //  domain:       c-shaped
// //  space sampling: locations != nodes
// //  time sampling:  locations != nodes
// //  penalization: simple laplacian
// //  missing:      yes
// //  covariates:   no
// //  BC:           no
// //  order FE:     1
// //  GCV optimization: grid exact
// //  time penalization: separable (mass penalization)
// TEST(gcv_qstrpde_test, laplacian_nonparametric_samplingatlocations_gridexact) {
  
//     std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/space_time/Test_3/per_pull_request"; 
  
//     // define domain
//     MeshLoader<Mesh2D> domain("c_shaped_adj");
//     unsigned int M = 3;  
//     Mesh<1, 1> time_mesh(0, fdapde::testing::pi, M-1);     

//     // import data from files
//     DMatrix<double> space_locs = read_csv<double>("../data/models/gcv/2D_test9/locs.csv");
//     DMatrix<double> time_locs = read_csv<double>("../data/models/gcv/2D_test9/time_locations.csv"); 
//     DMatrix<double> y = read_csv<double>("../data/models/gcv/2D_test9/y.csv");

//     // define regularizing PDE in space 
//     auto Ld = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
//     PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

//     // define regularizing PDE in time
//     auto Lt = -bilaplacian<SPLINE>();
//     PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

//     // define model
//     double alpha = 0.5;
//     QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);
//     model.set_spatial_locations(space_locs);
//     model.set_temporal_locations(time_locs);
//     // set model's data
//     BlockFrame<double, int> df;
//     df.stack(OBSERVATIONS_BLK, y);
//     model.set_data(df);
//     model.init();
//     // define GCV function and grid of \lambda_D values
//     auto GCV = model.gcv<ExactEDF>();
//     std::vector<DVector<double>> lambdas_d_t;
//     for(double xs = -4.0; xs <= -2.0; xs +=1.0)
//       for(double xt = -7.0; xt <= -5.0; xt +=1.0) 
//         lambdas_d_t.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));

//     // optimize GCV
//     Grid<fdapde::Dynamic> opt;
//     opt.optimize(GCV, lambdas_d_t);
//     // test correctness
//     // EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/models/gcv/2D_test9/edfs.mtx"));
//     // EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/models/gcv/2D_test9/gcvs.mtx"));

//     std::ofstream fileEDF_scores(R_path + "/2D_test1/edfs.mtx");
//     for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
//       fileEDF_scores << std::setprecision(16) << GCV.edfs()[i] << "\n" ; 
//     fileEDF_scores.close();

//     std::ofstream fileGCV_scores(R_path + "/2D_test1/gcvs.mtx");
//     for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
//       fileGCV_scores << std::setprecision(16) << GCV.gcvs()[i] << "\n" ; 
//     fileGCV_scores.close();
  
// }

// // test 2
// //  domain:       c-shaped
// //  space sampling: locations != nodes
// //  time sampling:  locations != nodes
// //  penalization: simple laplacian
// //  missing:      yes
// //  covariates:   no
// //  BC:           no
// //  order FE:     1
// //  GCV optimization: grid stochastic
// //  time penalization: separable (mass penalization)
// TEST(gcv_qstrpde_test, laplacian_nonparametric_samplingatlocations_gridstochastic) {

//     std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/space_time/Test_3/per_pull_request"; 
  
//     // define domain
//     MeshLoader<Mesh2D> domain("c_shaped_adj");
//     unsigned int M = 3;  
//     Mesh<1, 1> time_mesh(0, fdapde::testing::pi, M-1);
//     // import data from files
//     DMatrix<double> space_locs = read_csv<double>("../data/models/gcv/2D_test10/locs.csv");
//     DMatrix<double> time_locs = read_csv<double>("../data/models/gcv/2D_test10/time_locations.csv"); 
//     DMatrix<double> y = read_csv<double>("../data/models/gcv/2D_test10/y.csv");

//     // define regularizing PDE in space 
//     auto Ld = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
//     PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

//     // define regularizing PDE in time
//     auto Lt = -bilaplacian<SPLINE>();
//     PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

//     // define model
//     double alpha = 0.5;
//     QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);
//     model.set_spatial_locations(space_locs);
//     model.set_temporal_locations(time_locs);
//     // set model's data
//     BlockFrame<double, int> df;
//     df.stack(OBSERVATIONS_BLK, y);
//     model.set_data(df);
//     model.init();
//     // define GCV function and grid of \lambda_D values
//     std::size_t seed = 66546513;
//     auto GCV = model.gcv<StochasticEDF>(100, seed);
//     std::vector<DVector<double>> lambdas_d_t;
//     for(double xs = -4.0; xs <= -2.0; xs +=1.0)
//       for(double xt = -7.0; xt <= -5.0; xt +=1.0)
//         lambdas_d_t.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));

//     // optimize GCV
//     Grid<fdapde::Dynamic> opt;
//     opt.optimize(GCV, lambdas_d_t);

//     // test correctness
//     // EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/models/gcv/2D_test10/edfs.mtx"));
//     // EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/models/gcv/2D_test10/gcvs.mtx"));

//     std::ofstream fileEDF_scores(R_path + "/2D_test2/edfs.mtx");
//     for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
//       fileEDF_scores << std::setprecision(16) << GCV.edfs()[i] << "\n" ; 
//     fileEDF_scores.close();

//     std::ofstream fileGCV_scores(R_path + "/2D_test2/gcvs.mtx");
//     for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
//       fileGCV_scores << std::setprecision(16) << GCV.gcvs()[i] << "\n" ; 
//     fileGCV_scores.close();
// }
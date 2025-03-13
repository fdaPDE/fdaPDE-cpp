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
using fdapde::core::dt;
using fdapde::core::FEM;
using fdapde::core::fem_order;
using fdapde::core::bilaplacian;
using fdapde::core::laplacian;
using fdapde::core::PDE;
using fdapde::core::Triangulation;
using fdapde::core::SPLINE;
using fdapde::core::spline_order;

#include "../../fdaPDE/models/sampling_design.h"
#include "../../fdaPDE/models/regression/qsrpde.h"

using fdapde::models::QSRPDE;
using fdapde::models::SpaceTimeSeparable;
using fdapde::models::SpaceTimeParabolic;
using fdapde::models::SpaceOnly;
using fdapde::models::Sampling;

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_mtx;
using fdapde::testing::read_csv;


// /* test 1 
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    time penalization: separable (mass penalization)
//  */
// TEST(sqrpde_time_test, laplacian_nonparametric_samplingatnodes_separable_monolithic) {

//   // Parameters 
//   const std::string TestNumber = "1"; 

//   std::vector<double> alphas = {0.1, 0.5, 0.9};  

//   // number of simulations
//   unsigned int n_sim = 1;

//   // Choose missing strategy and proportion
//   std::string missing = "_missing"; 
//   // std::string missing = ""; 

//   const std::string p_string = "/p_50";
  
//   // Marco
//   // std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/PACS_project_shared"; 
//   // Ilenia 
//   std::string path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared"; 

//   std::string path_test = path + "/models/space_time/Test_" + TestNumber + p_string ; 


//   // define time domain 
//   DVector<double> time_mesh;
//   time_mesh.resize(6);
//   std::size_t i = 0;
//   for(double x = 0; x <= 2; x+=0.4, ++i) time_mesh[i] = x;
  
//   // define spatial domain 
//   MeshLoader<Mesh2D> domain("unit_square_coarse");
//   // define regularizing PDE    
//   auto L = -laplacian<FEM>();
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3 * time_mesh.rows(), 1);
//   PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
  
//   // define statistical model
//   double lambda_D ; // smoothing in space
//   double lambda_T ; // smoothing in time

//   for(double alpha : alphas){

//     unsigned int alpha_int = alpha*100; 

//     std::cout << "------------------------------------------alpha=" << std::to_string(alpha_int) << "%-------------------------------------------------" << std::endl; 

//     SQRPDE<decltype(problem), SpaceTimeSeparable, GeoStatMeshNodes, MonolithicSolver> model(problem, time_mesh, alpha);

//     for(unsigned int sim = 1; sim <= n_sim; ++sim){
    
//     // Read lambda
//     std::ifstream fileLambdaS(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/GCV/Exact/LambdaS" + missing + ".csv");
//     if (fileLambdaS.is_open()){
//       fileLambdaS >> lambda_D; 
//       fileLambdaS.close();
//     }

//     std::ifstream fileLambdaT(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/GCV/Exact/LambdaT" + missing + ".csv");
//     if (fileLambdaT.is_open()){
//       fileLambdaT >> lambda_T; 
//       fileLambdaT.close();
//     }

//     model.set_lambda_D(lambda_D);
//     model.set_lambda_T(lambda_T);

//     // import data from files
//     DMatrix<double> y = read_csv<double>(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/y" + missing + ".csv");

//     // set model data
//     BlockFrame<double, int> df;
//     df.stack(OBSERVATIONS_BLK, y);
//     model.set_data(df);
    
//     // solve smoothing problem
//     model.init();
//     model.solve();

//     // Save C++ solution 
//     DMatrix<double> computedF = model.f();
//     const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//     // std::ofstream filef(path_test + "/alpha_" + alpha_string + "/fCpp.csv");
//     std::ofstream filef(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/fCpp" + missing + "_new.csv");
//     if (filef.is_open()){
//           filef << computedF.format(CSVFormatf);
//           filef.close();
//     }

//     }

//   }

// }


// // test 2
// //    domain:         c-shaped
// //    space sampling: locations != nodes
// //    time sampling:  locations = nodes
// //    missing data:   no
// //    penalization:   simple laplacian
// //    covariates:     yes
// //    BC:             no
// //    order FE:       1
// //    time penalization: separable (mass penalization)
// TEST(sqrpde_time_test, laplacian_semiparametric_samplingatlocations_timelocations_separable_monolithic) {

//     // Marco 
//     std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/space_time/Test_2"; 
//     //   // Ilenia 
//     //   std::string path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/models/space_time/Test_2"; 

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
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3 * time_mesh.rows(), 1);
//     PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

//     unsigned int n_sim = 10; 
//     for(unsigned int sim = 1; sim <= n_sim; ++sim){
//       std::cout << "---------------------------Simulation #" << sim << "--------------------------" << std::endl; 
//       for(double alpha : alphas){

//         unsigned int alpha_int = alpha*100; 
//         std::string alpha_string = std::to_string(alpha_int); 

//         // load data from .csv files
//         DMatrix<double> y = read_csv<double>(R_path + "/sim_" + std::to_string(sim) + "/y.csv");

//         // optima lambdas 
//         double lambda_D = read_csv<double>(R_path + "/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + "/lambda_opt.csv")(0,0); 
//         double lambda_T = lambda_D;

//         SQRPDE<decltype(problem), SpaceTimeSeparable, GeoStatLocations, MonolithicSolver> model(problem, time_mesh, alpha);
//         model.set_lambda_D(lambda_D);
//         model.set_lambda_T(lambda_T);
//         model.set_spatial_locations(space_locs);
//         // set model's data
//         BlockFrame<double, int> df;
//         df.stack(OBSERVATIONS_BLK, y);
//         df.stack(OBSERVATIONS_BLK, X);
//         model.set_data(df);
//         // solve smoothing problem
//         model.init();
//         model.solve();

//         // Save C++ solution 
//         DMatrix<double> computedF = model.f();
//         const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filef(R_path + "/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + "/f.csv");
//         if (filef.is_open()){
//           filef << computedF.format(CSVFormatf);
//           filef.close();
//         }

//         DMatrix<double> computedFn = model.Psi()*model.f();
//         const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filefn(R_path + "/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + "/fn.csv");
//         if (filefn.is_open()){
//           filefn << computedFn.format(CSVFormatfn);
//           filefn.close();
//         }

//         DMatrix<double> computedbeta = model.beta();
//         const static Eigen::IOFormat CSVFormatbeta(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filebeta(R_path + "/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + "/beta.csv");
//         if (filebeta.is_open()){
//           filebeta << computedbeta.format(CSVFormatbeta);
//           filebeta.close();
//         } 

//       }
//     }
// }


// // Test ufficiale 

// TEST NON-PARAMETRICO
// test 3   
//    domain:         c-shaped
//    space sampling: locations != nodes
//    time sampling:  locations != nodes
//    missing data:   yes
//    penalization:   simple laplacian
//    covariates:     no
//    BC:             no
//    order FE:       1
//    time penalization: separable (mass penalization)
TEST(sqrpde_time_test, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

    // Marco 
    std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/space_time/Test_3"; 
    // Ilenia 
    // std::string R_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/models/space_time/Test_3"; 

    std::vector<double> alphas = {0.1}; // , 0.5, 0.9};   // , 0.9}; 
    std::vector<std::string> data_types = {"all"}; 
    std::string p_string = "50";   
    std::string lambda_selection_type = "gcv_smooth_eps1e-1";   // gcv gcv_smooth manual 
    bool check_for_merge = true; 

    // define temporal domain
    unsigned int M = 7; 
    std::string M_string = std::to_string(M);
    std::string path_mesh = "/M_7";

    double tf = fdapde::testing::pi;   // final time 
    Triangulation<1, 1> time_mesh(0, tf, M-1);     // t0, tf, #subintervals 

    // define spatial domain and regularizing PDE
    MeshLoader<Triangulation<2,2>> domain("c_shaped_adj");

    // import locs from files
    DMatrix<double> space_locs = read_csv<double>(R_path + "/space_locs.csv");
    DMatrix<double> time_locs = read_csv<double>(R_path + "/time_locs.csv");

    // define regularizing PDE in space 
    auto Ld = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3 * time_mesh.n_nodes(), 1);
    PDE<Triangulation<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

    // define regularizing PDE in time
    auto Lt = -bilaplacian<SPLINE>();
    PDE<Triangulation<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

    unsigned int n_sim = 10; 

    // // lambdas_d_t for RMSE
    // std::vector<SVector<2>> lambdas_try; 
    // for(double xs = -4.0; xs <= -1.5; xs +=0.5)
    //   for(double xt = -7.0; xt <= -6.0; xt +=1.0) 
    //     lambdas_try.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));

    for(auto data_type : data_types){
      if(data_type == "all")
        std::cout << "---------------------------------------------ALL DATA----------------------------" << std::endl; 
      else 
        std::cout << "---------------------------------------------MISSING DATA----------------------------" << std::endl; 
      

      for(unsigned int sim = 1; sim <= n_sim; ++sim){
        
          std::cout << "---------------------------Simulation #" << sim << "--------------------------" << std::endl; 
          for(double alpha : alphas){

            // std::size_t count_l = 1; 
            // for(SVector<2> l : lambdas_try){
            //   const std::string lambda_string = std::to_string(count_l); 

              unsigned int alpha_int = alpha*100; 
              std::string alpha_string = std::to_string(alpha_int); 

              // load data from .csv files
              DMatrix<double> y; 
              if(data_type == "all")
                y = read_csv<double>(R_path + "/simulations/all_new/sim_" + std::to_string(sim) + "/y_all.csv");  // ATT: messo all_new
              else{
                y = read_csv<double>(R_path + "/simulations/miss_strategy_" + data_type + "/p_" + p_string + "/sim_" + std::to_string(sim) + "/y.csv");
              }

              std::string solutions_path; 
              if(data_type == "all")
                solutions_path = R_path + "/simulations/all_new/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + path_mesh + "/" + lambda_selection_type ; // ATT: messo all_new
              else
                solutions_path = R_path + "/simulations/miss_strategy_" + data_type + "/p_" + p_string + "/sim_" + std::to_string(sim) + "/alpha_" + alpha_string +
                                 path_mesh + "/" + lambda_selection_type ;

              std::cout << "Solution path: " << solutions_path << std::endl ; 


              // optima lambdas 
              double lambda_D;  
              double lambda_T;  
              std::ifstream fileLambdaS_gcv(solutions_path + "/lambda_s_opt.csv");
              if(fileLambdaS_gcv.is_open()){
                fileLambdaS_gcv >> lambda_D; 
                fileLambdaS_gcv.close();
              }
              std::ifstream fileLambdaT(solutions_path + "/lambda_t_opt.csv");
              if(fileLambdaT.is_open()){
                fileLambdaT >> lambda_T; 
                fileLambdaT.close();
              }
             
             std::cout << "lambdaS = " << lambda_D << std::endl ;  
             std::cout << "lambdaT = " << lambda_T << std::endl ;  

              QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);
              model.set_lambda_D(lambda_D);
              model.set_lambda_T(lambda_T);
              model.set_spatial_locations(space_locs);
              model.set_temporal_locations(time_locs);

              // set model's data
              BlockFrame<double, int> df;
              df.stack(OBSERVATIONS_BLK, y);
              model.set_data(df);
              // solve smoothing problem
              model.init();
              model.solve();

              if(check_for_merge){
                solutions_path += "/check_merge"; 
              }

              // Save C++ solution 
              DMatrix<double> computedF = model.f();
              const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
              std::ofstream filef(solutions_path + "/f.csv");
              if (filef.is_open()){
                filef << computedF.format(CSVFormatf);
                filef.close();
              }

              DMatrix<double> computedFn = model.Psi()*model.f();
              const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
              std::ofstream filefn(solutions_path + "/fn.csv");
              if (filefn.is_open()){
                filefn << computedFn.format(CSVFormatfn);
                filefn.close();
              }

              // count_l += 1; 
            }
      }
    }


}



// // TEST SEMI-PARAMETRICO
// // test 3   
// //    domain:         c-shaped
// //    space sampling: locations != nodes
// //    time sampling:  locations != nodes
// //    missing data:   yes
// //    penalization:   simple laplacian
// //    covariates:     no
// //    BC:             no
// //    order FE:       1
// //    time penalization: separable (mass penalization)
// TEST(sqrpde_time_test_semiparam, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

//     // Marco 
//     // std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/space_time/Test_3"; 
//     // Ilenia 
//     std::string R_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/models/space_time/Test_3"; 

//     std::vector<double> alphas = {0.1};  // , 0.5, 0.9};   // , 0.9}; 
//     std::vector<std::string> data_types = {"all"}; 
//     std::string p_string = "50";   
//     std::string lambda_selection_type = "gcv_smooth_eps1e-1";   // gcv gcv_smooth manual 

//     // define temporal domain
//     unsigned int M = 7; 
//     std::string M_string = std::to_string(M);
//     std::string path_mesh = "/M_7";

//     double tf = fdapde::testing::pi;   // final time 
//     Triangulation<1, 1> time_mesh(0, tf, M-1);     // t0, tf, #subintervals 

//     // define spatial domain and regularizing PDE
//     MeshLoader<Triangulation<2,2>> domain("c_shaped_adj");

//     // import locs from files
//     DMatrix<double> space_locs = read_csv<double>(R_path + "/space_locs.csv");
//     DMatrix<double> time_locs = read_csv<double>(R_path + "/time_locs.csv");

//     // import covariate (they do not depend on the sim)
//     DMatrix<double> X = read_csv<double>(R_path + "/simulations/all_covariate/X.csv");
//     std::cout << "Dimension X = ( " << X.rows() << " , " << X.cols() << " )" << std::endl ; 

//     // define regularizing PDE in space 
//     auto Ld = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3 * time_mesh.n_nodes(), 1);
//     PDE<Triangulation<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

//     // define regularizing PDE in time
//     auto Lt = -bilaplacian<SPLINE>();
//     PDE<Triangulation<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

//     unsigned int n_sim = 1; 

//     // // lambdas_d_t for RMSE
//     // std::vector<SVector<2>> lambdas_try; 
//     // for(double xs = -4.0; xs <= -1.5; xs +=0.5)
//     //   for(double xt = -7.0; xt <= -6.0; xt +=1.0) 
//     //     lambdas_try.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));

//     for(auto data_type : data_types){
//       if(data_type == "all")
//         std::cout << "---------------------------------------------ALL DATA----------------------------" << std::endl; 
//       else 
//         std::cout << "---------------------------------------------MISSING DATA----------------------------" << std::endl; 
      

//       for(unsigned int sim = 1; sim <= n_sim; ++sim){
        
//           std::cout << "---------------------------Simulation #" << sim << "--------------------------" << std::endl; 
//           for(double alpha : alphas){

//             // std::size_t count_l = 1; 
//             // for(SVector<2> l : lambdas_try){
//             //   const std::string lambda_string = std::to_string(count_l); 

//               unsigned int alpha_int = alpha*100; 
//               std::string alpha_string = std::to_string(alpha_int); 

//               // load data from .csv files
//               DMatrix<double> y; 
//               if(data_type == "all")
//                 y = read_csv<double>(R_path + "/simulations/all_covariate/sim_" + std::to_string(sim) + "/y.csv");
//               else
//                 y = read_csv<double>(R_path + "/simulations/covariate_miss_strategy_" + data_type + "/p_" + p_string + "/sim_" + std::to_string(sim) + "/y.csv");

//               std::string solutions_path; 
//               if(data_type == "all")
//                 solutions_path = R_path + "/simulations/all_covariate/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + path_mesh + "/" + lambda_selection_type; 
//               else
//                 solutions_path = R_path + "/simulations/covariate_miss_strategy_" + data_type + "/p_" + p_string + "/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + 
//                                 path_mesh + "/" + lambda_selection_type ;

//               std::cout << "Solution path: " << solutions_path << std::endl ; 


//               // optima lambdas 
//               double lambda_D;  
//               double lambda_T;  
//               std::ifstream fileLambdaS_gcv(solutions_path + "/lambda_s_opt.csv");
//               if(fileLambdaS_gcv.is_open()){
//                 fileLambdaS_gcv >> lambda_D; 
//                 fileLambdaS_gcv.close();
//               }
//               std::ifstream fileLambdaT(solutions_path + "/lambda_t_opt.csv");
//               if(fileLambdaT.is_open()){
//                 fileLambdaT >> lambda_T; 
//                 fileLambdaT.close();
//               }
             
//              std::cout << "lambdaS = " << lambda_D << std::endl ;  
//              std::cout << "lambdaT = " << lambda_T << std::endl ;  

//               QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);
//               model.set_lambda_D(lambda_D);
//               model.set_lambda_T(lambda_T);
//               model.set_spatial_locations(space_locs);
//               model.set_temporal_locations(time_locs);

//               // set model's data
//               BlockFrame<double, int> df;
//               df.stack(OBSERVATIONS_BLK, y);
//               df.insert(DESIGN_MATRIX_BLK, X);
//               model.set_data(df);
//               // solve smoothing problem
//               model.init();
//               model.solve();

//               // Save C++ solution 
//               DMatrix<double> computedF = model.f();
//               const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//               std::ofstream filef(solutions_path + "/f.csv");
//               if (filef.is_open()){
//                 filef << computedF.format(CSVFormatf);
//                 filef.close();
//               }

//               DMatrix<double> computedFn = model.Psi()*model.f();
//               const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//               std::ofstream filefn(solutions_path + "/fn.csv");
//               if (filefn.is_open()){
//                 filefn << computedFn.format(CSVFormatfn);
//                 filefn.close();
//               }

//               DMatrix<double> computedBeta = model.beta();
//               const static Eigen::IOFormat CSVFormatbeta(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//               std::ofstream filebeta(solutions_path + "/beta.csv");
//               if (filebeta.is_open()){
//                 filebeta << computedBeta.format(CSVFormatbeta);
//                 filebeta.close();
//               }

//               //count_l += 1; 
//             }
//       }
//     }


// }






// ------------ PER PALU

// // test 1 
// //    domain:         c-shaped
// //    space sampling: locations != nodes
// //    time sampling:  locations != nodes
// //    missing data:   no
// //    penalization:   simple laplacian
// //    covariates:     no
// //    BC:             no
// //    order FE:       1
// //    time penalization: separable (mass penalization)
// TEST(qstrpde_test, laplacian_nonparametric_samplingatlocations_separable_monolithic) {
  
//   // define temporal domain
//   unsigned int M = 7; 
//   double tf = fdapde::testing::pi;   
//   Triangulation<1, 1> time_mesh(0, tf, M-1);     // t0, tf, #subintervals 

//   // define spatial domain and regularizing PDE
//   MeshLoader<Triangulation<2,2>> domain("c_shaped_adj");

//   // define regularizing PDE in space 
//   auto Ld = -laplacian<FEM>();
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3 * time_mesh.n_nodes(), 1);
//   PDE<Triangulation<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

//   // define regularizing PDE in time
//   auto Lt = -bilaplacian<SPLINE>();
//   PDE<Triangulation<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

//   // import locs from files
//   DMatrix<double> space_locs = read_csv<double>("../data/models/qstrpde/2D_test1/locs.csv");
//   DMatrix<double> time_locs = read_csv<double>("../data/models/qstrpde/2D_test1/time_locations.csv");

//   // load data from .csv files
//   DMatrix<double> y; 
//   y = read_csv<double>("../data/models/qstrpde/2D_test1/y.csv");

//   // alpha
//   double alpha = 0.5;

//   // lambdas 
//   double lambda_D = 1e-3; 
//   double lambda_T = 1e-6; 

//   // define model
//   QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);
//   model.set_lambda_D(lambda_D);
//   model.set_lambda_T(lambda_T);
//   model.set_spatial_locations(space_locs);
//   model.set_temporal_locations(time_locs);

//   // set model's data
//   BlockFrame<double, int> df;
//   df.stack(OBSERVATIONS_BLK, y);
//   model.set_data(df);

//   // solve smoothing problem
//   model.init();
//   model.solve();

//   // test correctness
//   EXPECT_TRUE(almost_equal(model.f(), "../data/models/qstrpde/2D_test1/sol.mtx"));

//     // Save C++ solution 
//   DMatrix<double> computedF = model.f();
//   Eigen::saveMarket(computedF, "../data/models/qstrpde/2D_test1/sol_modifica_palu.mtx");

// }



  // // Save C++ solution 
  // DMatrix<double> computedF = model.f();
  // Eigen::saveMarket(computedF, "../data/models/qstrpde/2D_test1/sol.mtx");



// ---------------------------------------------------OLD-----------------------------------------------------------------------


// /* test 1 
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    time penalization: separable (mass penalization)
//  */
// TEST(sqrpde_time_test, laplacian_nonparametric_samplingatnodes_separable_monolithic) {

//   // Parameters 
//   const std::string TestNumber = "1"; 

//   std::vector<double> alphas = {0.1};   // , 0.5, 0.9}; 

//   // number of simulations
//   unsigned int n_sim = 11;
  
//   // Marco
//   // std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/PACS_project_shared"; 
//   // Ilenia 
//   std::string path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared"; 

//   std::string path_test = path + "/models/space_time/Test_" + TestNumber ;


//   // define time domain
//   DVector<double> time_mesh;
//   time_mesh.resize(11);
//   std::size_t i = 0;
//   for(double x = 0; x <= 2; x+=0.2, ++i) time_mesh[i] = x;
  
//   // define spatial domain 
//   MeshLoader<Triangulation<2,2>> domain("unit_square_coarse");
//   // define regularizing PDE    
//   auto L = -laplacian<FEM>();
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3 * time_mesh.rows(), 1);
//   PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
  
//   // define statistical model
//   double lambda_D ; // smoothing in space
//   double lambda_T ; // smoothing in time

//   for(double alpha : alphas){

//     unsigned int alpha_int = alpha*100; 

//     std::cout << "------------------------------------------alpha=" << std::to_string(alpha_int) << "%-------------------------------------------------" << std::endl; 

//     SQRPDE<decltype(problem), SpaceTimeSeparable, GeoStatMeshNodes, MonolithicSolver> model(problem, time_mesh, alpha);

//     for(unsigned int sim = 11; sim <= n_sim; ++sim){
    
//     // // Read lambda
//     // std::ifstream fileLambdaS(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/GCV/Exact/LambdaS.csv");
//     // if (fileLambdaS.is_open()){
//     //   fileLambdaS >> lambdaS; 
//     //   fileLambdaS.close();
//     // }

//     // std::ifstream fileLambdaT(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/GCV/Exact/LambdaT.csv");
//     // if (fileLambdaT.is_open()){
//     //   fileLambdaT >> lambdaT; 
//     //   fileLambdaT.close();
//     // }

//     lambda_D = std::pow(10,-4.4);
//     lambda_T = std::pow(10,-4.4);

//     model.set_lambda_D(lambda_D);
//     model.set_lambda_T(lambda_T);

//     // import data from files
//     DMatrix<double> y = read_csv<double>(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/y.csv");

//     std::cout << " 1" << std::endl ; 
//     // set model data
//     BlockFrame<double, int> df;
//     std::cout << " 2" << std::endl ;
//     df.stack(OBSERVATIONS_BLK, y);
//     std::cout << " 3" << std::endl ;
//     model.set_data(df);
//     std::cout << " 4" << std::endl ;
    
//     // solve smoothing problem
//     model.init();
//     std::cout << " 5" << std::endl ;
//     model.solve();
//     std::cout << " 6" << std::endl ;

//     // Save C++ solution 
//     DMatrix<double> computedF = model.f();
//     std::cout << "Dim f = " << computedF.rows() << std::endl ; 
//     const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//     // std::ofstream filef(path_test + "/alpha_" + alpha_string + "/fCpp.csv");
//     std::ofstream filef(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/fCpp.csv");
//     if (filef.is_open()){
//           filef << computedF.format(CSVFormatf);
//           filef.close();
//     }

//     }

//   }

// }




// /* test 2 
//    domain:       c-shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    time penalization: separable (mass penalization)
//  */
// TEST(sqrpde_time_test, laplacian_semiparametric_samplingatlocations_separable_monolithic) {

//   // Parameters 
//   const std::string TestNumber = "2"; 

//   std::vector<double> alphas = {0.1};   // , 0.5, 0.9};

//   // number of simulations
//   unsigned int n_sim = 11;
  
//   // Marco
//   // std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/PACS_project_shared"; 
//   // Ilenia 
//   std::string path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared"; 

//   std::string path_test = path + "/models/space_time/Test_" + TestNumber ;


//   // define time domain
//   DVector<double> time_mesh;
//   time_mesh.resize(5);
//   for(std::size_t i = 0; i < 5; ++i)
//     time_mesh[i] = (fdapde::testing::pi/4)*i;

//   // define domain and regularizing PDE
//   MeshLoader<Triangulation<2,2>> domain("c_shaped");
//   // define regularizing PDE
//   auto L = -laplacian<FEM>();
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
//   PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
//   // define model
//   double lambda_D ;
//   double lambda_T ;

//   // import data from files
//   DMatrix<double> locs = read_csv<double>(path_test + "/locs.csv");
//   DMatrix<double> X    = read_csv<double>(path_test + "/X.csv");

//   for(double alpha : alphas){

//     unsigned int alpha_int = alpha*100; 

//     std::cout << "------------------------------------------alpha=" << std::to_string(alpha_int) << "%-------------------------------------------------" << std::endl; 

//     SQRPDE<decltype(problem), SpaceTimeSeparable, GeoStatLocations, MonolithicSolver> model(problem, time_mesh, alpha);

//     model.set_spatial_locations(locs);

//     for(unsigned int sim = 11; sim <= n_sim; ++sim){
    
//       // Read lambda
//       std::ifstream fileLambdaS(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/GCV/Exact/LambdaS.csv");
//       if (fileLambdaS.is_open()){
//         fileLambdaS >> lambda_D; 
//         fileLambdaS.close();
//       }

//       std::ifstream fileLambdaT(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/GCV/Exact/LambdaT.csv");
//       if (fileLambdaT.is_open()){
//         fileLambdaT >> lambda_T; 
//         fileLambdaT.close();
//       }

//       model.set_lambda_D(lambda_D);
//       model.set_lambda_T(lambda_T);

//       // load data from .csv files
//       DMatrix<double> y  = read_csv<double>(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/y.csv");

//       // set model's data
//       BlockFrame<double, int> df;
//       df.stack(OBSERVATIONS_BLK, y);
//       df.stack(DESIGN_MATRIX_BLK, X);
//       model.set_data(df);
//       // solve smoothing problem
//       model.init();
//       model.solve();
      
//       // Save C++ solution 
//       DMatrix<double> computedF = model.f();
//       const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//       std::ofstream filef(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/fCpp.csv");
//       if (filef.is_open()){
//             filef << computedF.format(CSVFormatf);
//             filef.close();
//       }

      
//       // Save C++ solution 
//       DMatrix<double> computedFn = model.Psi()*model.f();
//       const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//       std::ofstream filefn(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/fnCpp.csv");
//       if (filefn.is_open()){
//             filefn << computedFn.format(CSVFormatfn);
//             filefn.close();
//       }


//       // Save beta
//       DMatrix<double> computedbeta = model.beta();
//       const static Eigen::IOFormat CSVFormatbeta(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//       std::ofstream filebeta(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/betaCpp.csv");
//       if (filebeta.is_open()){
//             filebeta << computedbeta.format(CSVFormatbeta);
//             filebeta.close();
//       } 

//     }

//   }
  
// }






// /* test 3
//    domain:       quasicircular domain
//    sampling:     areal
//    penalization: non-costant coefficients PDE
//    covariates:   no
//    BC:           no
//    order FE:     1
//    time penalization: parabolic (monolithic solution)
//  */
// TEST(SQTRPDE, Test3_NonCostantCoefficientsPDE_NonParametric_Areal_Parabolic_Monolithic_EstimatedIC) {
//   // define time domain, we skip the first time instant because we are going to use the first block of data
//   // for the estimation of the initial condition
//   DVector<double> time_mesh;
//   time_mesh.resize(11);
//   for(std::size_t i = 0; i < 10; ++i) time_mesh[i] = 0.4*i;

//   // define domain and regularizing PDE
//   MeshLoader<Triangulation<2,2>> domain("quasi_circle");
//   // load PDE coefficients data
//   CSVReader<double> reader{};
//   CSVFile<double> diffFile; // diffusion tensor
//   diffFile = reader.parseFile("data/models/STQRPDE/2D_test3/K.csv");
//   DMatrix<double> diffData = diffFile.toEigen();
//   CSVFile<double> adveFile; // transport vector
//   adveFile = reader.parseFile("data/models/STQRPDE/2D_test3/b.csv");
//   DMatrix<double> adveData = adveFile.toEigen();

//   // define non-constant coefficients
//   SpaceVaryingDiffusion<2> diffCoeff;
//   diffCoeff.setData(diffData);
//   SpaceVaryingAdvection<2> adveCoeff;
//   adveCoeff.setData(adveData);
//   // parabolic PDE
//   auto L = dT() + Laplacian(diffCoeff.asParameter()) + Gradient(adveCoeff.asParameter());
  
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, time_mesh.rows()); 
//   PDE problem(domain.mesh, L, u); // definition of regularizing PDE
//   // define statistical model
//   double lambdaS = std::pow(0.1, 8); // smoothing in space
//   double lambdaT = std::pow(0.1, 8); // smoothing in time

//   CSVReader<int> int_reader{};
//   CSVFile<int> arealFile; // incidence matrix for specification of subdomains
//   arealFile = int_reader.parseFile("data/models/STQRPDE/2D_test3/incidence_matrix.csv");
//   DMatrix<int> areal = arealFile.toEigen();
  
//   double alpha = 0.5; 
//   SQRPDE<decltype(problem), fdaPDE::models::SpaceTimeParabolic, fdaPDE::models::Areal,
// 	 fdaPDE::models::MonolithicSolver> model(problem, time_mesh, alpha);
//   model.setLambdaS(lambdaS);
//   model.setLambdaT(lambdaT);
//   model.set_spatial_locations(areal);

//   // load data from .csv files
//   CSVFile<double> yFile; // observation file
//   yFile = reader.parseFile  ("data/models/STQRPDE/2D_test3/y.csv");
//   DMatrix<double> y = yFile.toEigen();

//   // set model data
//   BlockFrame<double, int> df;
//   df.stack (OBSERVATIONS_BLK, y);
//   model.setData(df);
  
//   // define initial condition estimator over grid of lambdas
//   InitialConditionEstimator ICestimator(model);
//   std::vector<SVector<1>> lambdas;
//   for(double x = -9; x <= 3; x += 0.1) lambdas.push_back(SVector<1>(std::pow(10,x))); 
//   // compute estimate
//   ICestimator.apply(lambdas);
//   DMatrix<double> ICestimate = ICestimator.get();
//   // // test computation initial condition
//   // CSVFile<double> ICfile;
//   // ICfile = reader.parseFile("data/models/STQRPDE/2D_test3/IC.csv");  
//   // DMatrix<double> expectedIC = ICfile.toEigen();
//   // EXPECT_TRUE( almost_equal(expectedIC, ICestimate) );

//   // set estimated initial condition
//   model.setInitialCondition(ICestimate);
//   model.shift_time(1); // shift data one time instant forward
  
//   model.init();
//   model.solve();
  
//   // //   **  test correctness of computed results  **   

//   // DMatrix<double> computedF;
//   // computedF.resize((model.n_temporal_locs()+1)*model.n_basis(), 1);
//   // computedF << model.s(), model.f();
  
//   // // estimate of spatial field \hat f (with estimatate of initial condition)
//   // SpMatrix<double> expectedSolution;
//   // Eigen::loadMarket(expectedSolution, "data/models/STQRPDE/2D_test3/sol.mtx");
//   // std::size_t N = computedF.rows();
//   // EXPECT_TRUE( almost_equal(DMatrix<double>(expectedSolution).topRows(N), computedF) );

//   const std::string path_test = "/mnt/c/Users/marco/PACS/Project/Code/Cpp/fdaPDE-fork/test/data/models/STQRPDE/2D_test3"; 
//   // Save C++ solution 
//   DMatrix<double> computedF = model.f();
//   const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//   std::ofstream filef(path_test + "/f.csv");
//   if(filef.is_open()){
//     filef << computedF.format(CSVFormatf);
//     filef.close();
//   }

//   DMatrix<double> computedFn = model.Psi()*computedF;
//   const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//   std::ofstream filefn(path_test + "/fn.csv");
//   if(filefn.is_open()){
//     filefn << computedFn.format(CSVFormatfn);
//     filefn.close();
//   }


// }

/* test 4
   domain:       unit square 
   sampling:     locations != nodes 
   penalization: laplacian 
   covariates:   yes
   BC:           no
   order FE:     1
   time penalization: parabolic
 */
// TEST(STQRPDE, Test4_Laplacian_SemiParametric_Locations_Parabolic_Monolithic_EstimatedIC){

//   // Solver type
//   std::string solver = "Iterative";   // Iterative 

//   // define time domain, we skip the first time instant because we are going to use the first block of data
//   // for the estimation of the initial condition
//   DVector<double> time_mesh;
//   time_mesh.resize(6);
//   for(std::size_t i = 0; i < 5; ++i) time_mesh[i] = 0.4*i;

// //   // define domain and regularizing PDE
// //   MeshLoader<Triangulation<2,2>> domain("unit_square_coarse");
  
// //   // parabolic PDE
// //   auto L = dT() + Laplacian();
  
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, time_mesh.rows()); 
//   PDE problem(domain.mesh, L, u); // definition of regularizing PDE
//   // define statistical model
//   double lambdaS = std::pow(0.1, 6); // smoothing in space
//   double lambdaT = std::pow(0.1, 6); // smoothing in time

//   std::string lambda_string; 
//   if(almost_equal(lambdaS, 1e-06) && almost_equal(lambdaT, 1e-06))
//     lambda_string = "1e-06"; 

//   std::vector<double> alphas = {0.1, 0.5, 0.9};

//   // load sample position
//   CSVReader<double> reader{};
//   CSVFile<double> locFile; // locations file
//   locFile = reader.parseFile("/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/space_time/Test_4/locs.csv");
//   DMatrix<double> loc = locFile.toEigen();

//   CSVFile<double> XFile; // design matrix
//   XFile = reader.parseFile("/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/space_time/Test_4/X.csv");
//   DMatrix<double> X = XFile.toEigen();

//   // number of simulations
//   unsigned int n_sim = 10; 

//   for(double alpha : alphas){

//     unsigned int alpha_int = alpha*100; 

//     std::cout << "------------------------------------------alpha=" << std::to_string(alpha_int) << "%-------------------------------------------------" << std::endl; 

//     for(unsigned int sim = 1; sim <= n_sim; ++sim){

//       const std::string path_test = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/space_time/Test_4/results_lambdas_" + lambda_string + "/alpha_" + std::to_string(alpha_int) + "/" + solver + "/sim_" + std::to_string(sim);

//       // load data from .csv files
//       CSVFile<double> yFile; // observation file
//       yFile = reader.parseFile(path_test + "/y.csv");
//       DMatrix<double> y = yFile.toEigen();

//       BlockFrame<double, int> df;
//       df.stack(OBSERVATIONS_BLK, y);
//       df.stack(DESIGN_MATRIX_BLK, X);

//       SQRPDE<decltype(problem), fdaPDE::models::SpaceTimeParabolic, fdaPDE::models::GeoStatLocations,
//             fdaPDE::models::IterativeSolver> model(problem, time_mesh, alpha);
//       model.setLambdaS(lambdaS);
//       model.setLambdaT(lambdaT);
//       model.set_spatial_locations(loc);
//       model.setData(df);

//       // define initial condition estimator over grid of lambdas
//       InitialConditionEstimator ICestimator(model);
//       std::vector<SVector<1>> lambdas_IC;
//       for(double x = -9; x <= 3; x += 0.1) lambdas_IC.push_back(SVector<1>(std::pow(10,x))); 
//       // compute estimate
//       ICestimator.apply(lambdas_IC);
//       DMatrix<double> ICestimate = ICestimator.get();

//       // set estimated initial condition
//       model.setInitialCondition(ICestimate);
//       model.shift_time(1); // shift data one time instant forward
      
//       model.init();
//       model.solve();
  
//       // Save C++ solution 
//       DMatrix<double> computedF = model.f();
//       const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//       std::ofstream filef(path_test + "/f.csv");
//       if(filef.is_open()){
//         filef << computedF.format(CSVFormatf);
//         filef.close();
//       }

//       DMatrix<double> computedFn = model.Psi()*computedF;
//       const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//       std::ofstream filefn(path_test + "/fn.csv");
//       if(filefn.is_open()){
//         filefn << computedFn.format(CSVFormatfn);
//         filefn.close();
//       }

//       DMatrix<double> computedBeta = model.beta();
//       const static Eigen::IOFormat CSVFormatbeta(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//       std::ofstream filebeta(path_test + "/beta.csv");
//       if(filebeta.is_open()){
//         filebeta << computedBeta.format(CSVFormatbeta);
//         filebeta.close();
//       }

//       unsigned int computedNiter = model.n_iter();
//       std::ofstream ofile_niter;
//       ofile_niter.open(path_test + "/n_iter.txt");
//       ofile_niter << computedNiter; 
//       ofile_niter.close(); 

//     }

//   }
// }
// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

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
using fdapde::core::laplacian;
using fdapde::core::PDE;
using fdapde::core::Triangulation;

#include "../../fdaPDE/models/regression/mqsrpde.h"
#include "../../fdaPDE/models/sampling_design.h"
using fdapde::models::MQSRPDE;
using fdapde::models::SpaceOnly;

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

// double RMSE_metric(DVector<double> v1, DVector<double> v2){
//     double res = 0.; 
//     if(v1.size() != v2.size())
//         std::cout << std::endl << "----------ERROR IN RMSE COMPUTATION---------" << std::endl; 
//     for(auto i = 0; i < v1.size(); ++i){
//         res += (v1[i]-v2[i])*(v1[i]-v2[i]); 
//     }
//     return std::sqrt(1./(v1.size())*res); 
// }


// // TEST TESI

// // test 7
// //    domain:       unit square
// //    sampling:     locations != nodes
// //    penalization: constant PDE coefficients 
// //    covariates:   no
// //    BC:           no
// //    order FE:     1
// //    GCV optimization: grid exact
// TEST(gcv_msqrpde_test7, pde_nonparametric_samplingatlocations_spaceonly_gridexact) {

//     // path test  
//     std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/multiple_quantiles/Tests/Test_7"; 
//     //std::string R_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/models/multiple_quantiles/Tests/Test_7"; 

//     const std::string pde_type = "_lap";    // "_lap" "_Ktrue" "_casc"

//     // define domain
//     MeshLoader<Triangulation<2, 2>> domain("unit_square_test7");

//     // rhs 
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);

//     // define regularizing PDE

//     // lap 
//     if(pde_type != "_lap")
//         std::cout << "ERROR: YOU WANT TO USE K = I BUT YOU ARE USING SOMETHING ELSE" << std::endl; 
//     auto L = -laplacian<FEM>();   
//     PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

//     // // K = K_true
//     // if(pde_type != "_Ktrue")
//     //     std::cout << "ERROR: YOU WANT TO USE K = K_true BUT YOU ARE USING SOMETHING ELSE" << std::endl;
//     // SMatrix<2> K = read_csv<double>(R_path + "/data/true/K_true.csv"); 
//     // auto L = -diffusion<FEM>(K);  
//     // PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

//     // define statistical model
//     std::vector<double> alphas = {0.01, 0.02, 0.05, 0.10, 0.25, 
//                                   0.50, 0.75, 0.90, 0.91, 0.92, 
//                                   0.93, 0.94, 0.95, 0.96, 0.97, 
//                                   0.98, 0.99};  

//     // define grid of lambda values
//     std::vector<std::string> lambda_selection_types = {"gcv_smooth_eps1e-1"}; // {"gcv", "gcv_smooth_eps1e-3", "gcv_smooth_eps1e-2", "gcv_smooth_eps1e-1.5", "gcv_smooth_eps1e-1"}; // {"gcv", "gcv_smooth_eps1e-3", "gcv_smooth_eps1e-2", "gcv_smooth_eps1e-1.5", "gcv_smooth_eps1e-1"};     
//     std::vector<double> lambdas_1; std::vector<double> lambdas_2; std::vector<double> lambdas_5;
//     std::vector<double> lambdas_10; std::vector<double> lambdas_25; std::vector<double> lambdas_50;
//     std::vector<double> lambdas_75; std::vector<double> lambdas_90; std::vector<double> lambdas_91; 
//     std::vector<double> lambdas_92; std::vector<double> lambdas_93;  std::vector<double> lambdas_94; 
//     std::vector<double> lambdas_95; std::vector<double> lambdas_96;  std::vector<double> lambdas_97; 
//     std::vector<double> lambdas_98; std::vector<double> lambdas_99; 
//     for(double x = -7.5; x <= -3.5; x += 0.1) lambdas_1.push_back(std::pow(10, x));
//     for(double x = -7.5; x <= -3.5; x += 0.1) lambdas_2.push_back(std::pow(10, x));
//     for(double x = -8.0; x <= -2.0; x += 0.1) lambdas_5.push_back(std::pow(10, x)); 
//     for(double x = -7.5; x <= -3.0; x += 0.1) lambdas_10.push_back(std::pow(10, x)); 
//     for(double x = -7.0; x <= -3.5; x += 0.1) lambdas_25.push_back(std::pow(10, x));
//     for(double x = -6.0; x <= -3.0; x += 0.1) lambdas_50.push_back(std::pow(10, x)); 
//     for(double x = -6.0; x <= -3.0; x += 0.1) lambdas_75.push_back(std::pow(10, x)); 
//     for(double x = -7.5; x <= -3.5; x += 0.1) lambdas_90.push_back(std::pow(10, x));
//     for(double x = -6.5; x <= -3.5; x += 0.1) lambdas_91.push_back(std::pow(10, x));
//     for(double x = -6.5; x <= -3.5; x += 0.1) lambdas_92.push_back(std::pow(10, x));
//     for(double x = -6.5; x <= -3.5; x += 0.1) lambdas_93.push_back(std::pow(10, x));
//     for(double x = -6.5; x <= -4.0; x += 0.1) lambdas_94.push_back(std::pow(10, x));
//     for(double x = -7.5; x <= -3.5; x += 0.1) lambdas_95.push_back(std::pow(10, x));
//     for(double x = -6.5; x <= -4.0; x += 0.1) lambdas_96.push_back(std::pow(10, x));
//     for(double x = -6.5; x <= -4.0; x += 0.1) lambdas_97.push_back(std::pow(10, x));
//     for(double x = -7.5; x <= -3.5; x += 0.1) lambdas_98.push_back(std::pow(10, x)); 
//     for(double x = -8.0; x <= -3.5; x += 0.1) lambdas_99.push_back(std::pow(10, x));
//     double best_lambda; 

//     // Read covariates and locations
//     DMatrix<double> loc = read_csv<double>(R_path + "/data" + "/locs.csv"); 

//     // Simulations 
//     std::vector<unsigned int> simulations = {25}; // {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}; 
//     for(auto sim : simulations){
//         std::cout << "--------------------Simulation #" << std::to_string(sim) << "-------------" << std::endl; 

//         std::string solutions_path_rmse = R_path + "/data/simulations/sim_" + std::to_string(sim) + "/single_est" + pde_type + "/RMSE"; 

//         // // K = K_est
//         // if(pde_type != "_casc")
//         //     std::cout << "ERROR: YOU WANT TO USE K = K_est BUT YOU ARE USING SOMETHING ELSE" << std::endl; 
//         // SMatrix<2> K = read_csv<double>(R_path + "/data/simulations/sim_" + std::to_string(sim) + "/K.csv"); 
//         // auto L = -diffusion<FEM>(K);   // anisotropic diffusion  
//         // PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

//         for(auto alpha : alphas){

//             unsigned int alpha_int = alpha*100; 
//             std::string alpha_string = std::to_string(alpha_int); 

//             std::cout << "------------------alpha=" << alpha_string << "-----------------" << std::endl; 

//             // load data from .csv files
//             DMatrix<double> y = read_csv<double>(R_path + "/data/simulations/sim_" + std::to_string(sim) + "/y.csv");
//             BlockFrame<double, int> df;
//             df.insert(OBSERVATIONS_BLK, y);

//             // GCV:
//             for(auto lambda_selection_type : lambda_selection_types){
                
//                 std::string solutions_path_gcv = R_path + "/data/simulations/sim_" + std::to_string(sim) + "/single_est" + pde_type + "/" + lambda_selection_type; 
                
//                 QSRPDE<SpaceOnly> model_gcv(problem, Sampling::pointwise, alpha);
//                 model_gcv.set_spatial_locations(loc);

//                 std::vector<double> lambdas;
//                 if(almost_equal(alpha, 0.01)){
//                     lambdas = lambdas_1; 
//                 }  
//                 if(almost_equal(alpha, 0.02)){
//                     lambdas = lambdas_2; 
//                 }  
//                 if(almost_equal(alpha, 0.05)){
//                     lambdas = lambdas_5; 
//                 }  
//                 if(almost_equal(alpha, 0.10)){
//                     lambdas = lambdas_10; 
//                 }  
//                 if(almost_equal(alpha, 0.25)){
//                     lambdas = lambdas_25; 
//                 }  
//                 if(almost_equal(alpha, 0.50)){
//                     lambdas = lambdas_50; 
//                 }  
//                 if(almost_equal(alpha, 0.75)){
//                     lambdas = lambdas_75; 
//                 }  
//                 if(almost_equal(alpha, 0.90)){
//                     lambdas = lambdas_90; 
//                 } 
//                 if(almost_equal(alpha, 0.91)){
//                     lambdas = lambdas_91; 
//                 }   
//                 if(almost_equal(alpha, 0.92)){
//                     lambdas = lambdas_92; 
//                 }  
//                 if(almost_equal(alpha, 0.93)){
//                     lambdas = lambdas_93; 
//                 } 
//                 if(almost_equal(alpha, 0.94)){
//                     lambdas = lambdas_94; 
//                 }   
//                 if(almost_equal(alpha, 0.95)){
//                     lambdas = lambdas_95; 
//                 } 
//                 if(almost_equal(alpha, 0.96)){
//                     lambdas = lambdas_96; 
//                 }    
//                 if(almost_equal(alpha, 0.97)){
//                     lambdas = lambdas_97; 
//                 }    
//                 if(almost_equal(alpha, 0.98)){
//                     lambdas = lambdas_98; 
//                 }  
//                 if(almost_equal(alpha, 0.99)){
//                     lambdas = lambdas_99; 
//                 }  
//                 // define lambda sequence as matrix 
//                 DMatrix<double> lambdas_mat;
//                 lambdas_mat.resize(lambdas.size(), 1); 
//                 for(auto i = 0; i < lambdas_mat.rows(); ++i){
//                    lambdas_mat(i,0) = lambdas[i]; 
//                  }
//                 // set model's data
//                 model_gcv.set_exact_gcv(lambda_selection_type == "gcv"); 

//                 if(lambda_selection_type == "gcv_smooth_eps1e-3"){
//                     model_gcv.set_eps_power(-3.0); 
//                 }
//                 if(lambda_selection_type == "gcv_smooth_eps1e-2"){
//                     model_gcv.set_eps_power(-2.0); 
//                 }
//                 if(lambda_selection_type == "gcv_smooth_eps1e-1.5"){
//                     model_gcv.set_eps_power(-1.5); 
//                 }
//                 if(lambda_selection_type == "gcv_smooth_eps1e-1"){
//                     model_gcv.set_eps_power(-1.0); 
//                 }
                
//                 model_gcv.set_data(df);
//                 model_gcv.init();

//                 // define GCV function and grid of \lambda_D values
//                 auto GCV = model_gcv.gcv<ExactEDF>();
//                 // optimize GCV
//                 Grid<fdapde::Dynamic> opt;
//                 opt.optimize(GCV, lambdas_mat);
                
//                 best_lambda = opt.optimum()(0,0);
        
//                 std::cout << "Best lambda is: " << std::setprecision(16) << best_lambda << std::endl; 

//                 // Save lambda sequence 
//                 std::ofstream fileLambdaS(solutions_path_gcv + "/lambdas_seq_alpha_" + alpha_string + ".csv");
//                 for(std::size_t i = 0; i < lambdas.size(); ++i) 
//                     fileLambdaS << std::setprecision(16) << lambdas[i] << "\n"; 
//                 fileLambdaS.close();

//                 // Save lambda GCVopt for all alphas
//                 std::ofstream fileLambdaoptS(solutions_path_gcv + "/lambdas_opt_alpha_" + alpha_string + ".csv");
//                 if(fileLambdaoptS.is_open()){
//                     fileLambdaoptS << std::setprecision(16) << best_lambda;
//                     fileLambdaoptS.close();
//                 }

//                 // Save GCV 
//                 std::ofstream fileGCV_scores(solutions_path_gcv + "/score_alpha_" + alpha_string + ".csv");
//                 for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
//                     fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.gcvs()[i]) << "\n"; 
//                 fileGCV_scores.close();
//             }

//         }


//     }
// }



// // test 8
// //    domain:       unit square
// //    sampling:     locations != nodes
// //    penalization: constant PDE coefficients 
// //    covariates:   yes
// //    BC:           no
// //    order FE:     1
// //    GCV optimization: grid exact
// TEST(gcv_msqrpde_test8, pde_nonparametric_samplingatlocations_spaceonly_gridexact) {

//     // path test  
//     std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/multiple_quantiles/Tests/Test_8"; 
//     //std::string R_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/models/multiple_quantiles/Tests/Test_8"; 

//     const std::string pde_type = "_casc";    // "_lap" "_Ktrue" "_casc"
//     const bool save = false; 

//     // define domain
//     MeshLoader<Triangulation<2, 2>> domain("unit_square_test8");

//     // rhs 
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);

//     // define regularizing PDE

//     // // lap 
//     // if(pde_type != "_lap")
//     //     std::cout << "ERROR: YOU WANT TO USE K = I BUT YOU ARE USING SOMETHING ELSE" << std::endl; 
//     // auto L = -laplacian<FEM>();   
//     // PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

//     // // K = K_true
//     // if(pde_type != "_Ktrue")
//     //     std::cout << "ERROR: YOU WANT TO USE K = K_true BUT YOU ARE USING SOMETHING ELSE" << std::endl;
//     // SMatrix<2> K = read_csv<double>(R_path + "/data/true/K_true.csv"); 
//     // auto L = -diffusion<FEM>(K);  
//     // PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

//     // define statistical model
//     std::vector<double> alphas = {0.50, 0.75}; 
//     // {0.01, 0.02, 0.05, 0.10, 
//     //                               0.25, 0.50, 0.75, 
//     //                               0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99};  

//     // define grid of lambda values
//     std::vector<std::string> lambda_selection_types = {"gcv_smooth_eps1e-1"}; // {"gcv", "gcv_smooth_eps1e-3", "gcv_smooth_eps1e-2", "gcv_smooth_eps1e-1.5", "gcv_smooth_eps1e-1"}; // {"gcv", "gcv_smooth_eps1e-3", "gcv_smooth_eps1e-2", "gcv_smooth_eps1e-1.5", "gcv_smooth_eps1e-1"};     
//     std::vector<double> lambdas_1; std::vector<double> lambdas_2; std::vector<double> lambdas_5;
//     std::vector<double> lambdas_10; std::vector<double> lambdas_25; std::vector<double> lambdas_50;
//     std::vector<double> lambdas_75; std::vector<double> lambdas_90; std::vector<double> lambdas_91; 
//     std::vector<double> lambdas_92; std::vector<double> lambdas_93;  std::vector<double> lambdas_94; 
//     std::vector<double> lambdas_95; std::vector<double> lambdas_96;  std::vector<double> lambdas_97; 
//     std::vector<double> lambdas_98; std::vector<double> lambdas_99; 
//     for(double x = -7.5; x <= -3.5; x += 0.1) lambdas_1.push_back(std::pow(10, x));
//     for(double x = -7.5; x <= -3.5; x += 0.1) lambdas_2.push_back(std::pow(10, x));
//     for(double x = -8.0; x <= -2.0; x += 0.1) lambdas_5.push_back(std::pow(10, x)); 
//     for(double x = -7.5; x <= -3.0; x += 0.1) lambdas_10.push_back(std::pow(10, x)); 
//     for(double x = -7.0; x <= -3.5; x += 0.1) lambdas_25.push_back(std::pow(10, x));
//     for(double x = -6.0; x <= -3.0; x += 1.0) lambdas_50.push_back(std::pow(10, x)); 
//     for(double x = -6.0; x <= -3.0; x += 1.0) lambdas_75.push_back(std::pow(10, x)); 
//     for(double x = -7.5; x <= -3.5; x += 0.1) lambdas_90.push_back(std::pow(10, x));
//     for(double x = -6.5; x <= -3.5; x += 0.1) lambdas_91.push_back(std::pow(10, x));
//     for(double x = -6.5; x <= -3.5; x += 0.1) lambdas_92.push_back(std::pow(10, x));
//     for(double x = -6.5; x <= -3.5; x += 0.1) lambdas_93.push_back(std::pow(10, x));
//     for(double x = -6.5; x <= -4.0; x += 0.1) lambdas_94.push_back(std::pow(10, x));
//     for(double x = -7.5; x <= -3.5; x += 0.1) lambdas_95.push_back(std::pow(10, x));
//     for(double x = -6.5; x <= -4.0; x += 0.1) lambdas_96.push_back(std::pow(10, x));
//     for(double x = -6.5; x <= -4.0; x += 0.1) lambdas_97.push_back(std::pow(10, x));
//     for(double x = -7.5; x <= -3.5; x += 0.1) lambdas_98.push_back(std::pow(10, x)); 
//     for(double x = -8.0; x <= -1.5; x += 0.1) lambdas_99.push_back(std::pow(10, x));
//     double best_lambda; 

//     // Read covariates and locations
//     DMatrix<double> loc = read_csv<double>(R_path + "/data" + "/locs.csv"); 
//     DMatrix<double> X = read_csv<double>(R_path + "/data" + "/X.csv"); 

//     // Simulations 
//     std::vector<unsigned int> simulations = {1}; // {1,2,3,4,5,6,7,8,9,10}; 
//     for(auto sim : simulations){

//         std::cout << "--------------------Simulation #" << std::to_string(sim) << "-------------" << std::endl; 

//         // K = K_est
//         if(pde_type != "_casc")
//             std::cout << "ERROR: YOU WANT TO USE K = K_est BUT YOU ARE USING SOMETHING ELSE" << std::endl; 
//         SMatrix<2> K = read_csv<double>(R_path + "/data/simulations/sim_" + std::to_string(sim) + "/K.csv"); 
//         auto L = -diffusion<FEM>(K);   // anisotropic diffusion  
//         PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

//         for(auto alpha : alphas){

//             unsigned int alpha_int = alpha*100; 
//             std::string alpha_string = std::to_string(alpha_int); 

//             std::cout << "------------------alpha=" << alpha_string << "-----------------" << std::endl; 

//             // load data from .csv files
//             DMatrix<double> y = read_csv<double>(R_path + "/data/simulations/sim_" + std::to_string(sim) + "/y.csv");
//             BlockFrame<double, int> df;
//             df.insert(OBSERVATIONS_BLK, y);
//             df.insert(DESIGN_MATRIX_BLK, X);

//             // GCV:
//             for(auto lambda_selection_type : lambda_selection_types){
                
//                 std::string solutions_path_gcv = R_path + "/data/simulations/sim_" + std::to_string(sim) + "/single_est" + pde_type + "/" + lambda_selection_type; 
                
//                 QSRPDE<SpaceOnly> model_gcv(problem, Sampling::pointwise, alpha);
//                 model_gcv.set_spatial_locations(loc);

//                 std::vector<double> lambdas;
//                 if(almost_equal(alpha, 0.01)){
//                     lambdas = lambdas_1; 
//                 }  
//                 if(almost_equal(alpha, 0.02)){
//                     lambdas = lambdas_2; 
//                 }  
//                 if(almost_equal(alpha, 0.05)){
//                     lambdas = lambdas_5; 
//                 }  
//                 if(almost_equal(alpha, 0.10)){
//                     lambdas = lambdas_10; 
//                 }  
//                 if(almost_equal(alpha, 0.25)){
//                     lambdas = lambdas_25; 
//                 }  
//                 if(almost_equal(alpha, 0.50)){
//                     lambdas = lambdas_50; 
//                 }  
//                 if(almost_equal(alpha, 0.75)){
//                     lambdas = lambdas_75; 
//                 }  
//                 if(almost_equal(alpha, 0.90)){
//                     lambdas = lambdas_90; 
//                 } 
//                 if(almost_equal(alpha, 0.91)){
//                     lambdas = lambdas_91; 
//                 }   
//                 if(almost_equal(alpha, 0.92)){
//                     lambdas = lambdas_92; 
//                 }  
//                 if(almost_equal(alpha, 0.93)){
//                     lambdas = lambdas_93; 
//                 } 
//                 if(almost_equal(alpha, 0.94)){
//                     lambdas = lambdas_94; 
//                 }   
//                 if(almost_equal(alpha, 0.95)){
//                     lambdas = lambdas_95; 
//                 } 
//                 if(almost_equal(alpha, 0.96)){
//                     lambdas = lambdas_96; 
//                 }    
//                 if(almost_equal(alpha, 0.97)){
//                     lambdas = lambdas_97; 
//                 }    
//                 if(almost_equal(alpha, 0.98)){
//                     lambdas = lambdas_98; 
//                 }  
//                 if(almost_equal(alpha, 0.99)){
//                     lambdas = lambdas_99; 
//                 }  
//                 // define lambda sequence as matrix 
//                 DMatrix<double> lambdas_mat;
//                 lambdas_mat.resize(lambdas.size(), 1); 
//                 for(auto i = 0; i < lambdas_mat.rows(); ++i){
//                    lambdas_mat(i,0) = lambdas[i]; 
//                 }

//                 if(lambda_selection_type == "gcv_smooth_eps1e-3"){
//                     model_gcv.set_eps_power(-3.0); 
//                 }
//                 if(lambda_selection_type == "gcv_smooth_eps1e-2"){
//                     model_gcv.set_eps_power(-2.0); 
//                 }
//                 if(lambda_selection_type == "gcv_smooth_eps1e-1.5"){
//                     model_gcv.set_eps_power(-1.5); 
//                 }
//                 if(lambda_selection_type == "gcv_smooth_eps1e-1"){
//                     model_gcv.set_eps_power(-1.0); 
//                 }
                
//                 model_gcv.set_data(df);
//                 model_gcv.init();

//                 // define GCV function and grid of \lambda_D values
//                 auto GCV = model_gcv.gcv<ExactEDF>();
//                 // optimize GCV
//                 Grid<fdapde::Dynamic> opt;
//                 opt.optimize(GCV, lambdas_mat);
                
//                 best_lambda = opt.optimum()(0,0);
        
//                 std::cout << "Best lambda is: " << std::setprecision(16) << best_lambda << std::endl; 

//                 if(save){
//                     // Save lambda sequence 
//                     std::ofstream fileLambdaS(solutions_path_gcv + "/lambdas_seq_alpha_" + alpha_string + ".csv");
//                     for(std::size_t i = 0; i < lambdas.size(); ++i) 
//                         fileLambdaS << std::setprecision(16) << lambdas[i] << "\n"; 
//                     fileLambdaS.close();

//                     // Save lambda GCVopt for all alphas
//                     std::ofstream fileLambdaoptS(solutions_path_gcv + "/lambdas_opt_alpha_" + alpha_string + ".csv");
//                     if(fileLambdaoptS.is_open()){
//                         fileLambdaoptS << std::setprecision(16) << best_lambda;
//                         fileLambdaoptS.close();
//                     }

//                     // Save GCV 
//                     std::ofstream fileGCV_scores(solutions_path_gcv + "/score_alpha_" + alpha_string + ".csv");
//                     for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
//                         fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.gcvs()[i]) << "\n"; 
//                     fileGCV_scores.close();
//                 }

//             }

//         }


//     }
// }



// //// 

// // TEST OBS RIPETUTE

// test 7 OBS RIPETUTE
//    domain:       unit square
//    sampling:     locations != nodes
//    penalization: constant PDE coefficients 
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
TEST(gcv_msqrpde_test7_obs_rip, pde_nonparametric_samplingatlocations_spaceonly_gridexact) {

    // path test  
    // std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/multiple_quantiles/Tests/Test_7"; 
    std::string R_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/models/multiple_quantiles/Tests/Test_7_obs_ripetute"; 

    const std::string pde_type = "_lap";    // "_lap" "_Ktrue" "_casc"
    
    // usare le obs ripetute?
    bool bool_obs_rip = true;

    // define domain
    MeshLoader<Triangulation<2, 2>> domain("unit_square_test7");

    // rhs 
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);

    // define regularizing PDE

    // lap 
    if(pde_type != "_lap")
        std::cout << "ERROR: YOU WANT TO USE K = I BUT YOU ARE USING SOMETHING ELSE" << std::endl; 
    auto L = -laplacian<FEM>();   
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

    // // K = K_true
    // if(pde_type != "_Ktrue")
    //     std::cout << "ERROR: YOU WANT TO USE K = K_true BUT YOU ARE USING SOMETHING ELSE" << std::endl;
    // SMatrix<2> K = read_csv<double>(R_path + "/data/true/K_true.csv"); 
    // auto L = -diffusion<FEM>(K);  
    // PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

    // // define statistical model
    // std::vector<double> alphas = {0.01, 0.02, 0.05, 0.10, 0.25, 
    //                               0.50, 0.75, 0.90, 0.91, 0.92, 
    //                               0.93, 0.94, 0.95, 0.96, 0.97, 
    //                               0.98, 0.99};  

    std::vector<double> alphas = {0.01, 0.02, 0.05, 0.10, 0.25, 
                                  0.75, 0.91, 0.92, 
                                  0.93, 0.94, 0.96, 0.97, 
                                  0.98};  

    // std::vector<double> alphas = {0.50, 0.90, 0.95, 0.99};  

    // define grid of lambda values
    std::vector<std::string> lambda_selection_types = {"gcv_smooth_eps1e-1"}; // {"gcv", "gcv_smooth_eps1e-3", "gcv_smooth_eps1e-2", "gcv_smooth_eps1e-1.5", "gcv_smooth_eps1e-1"}; // {"gcv", "gcv_smooth_eps1e-3", "gcv_smooth_eps1e-2", "gcv_smooth_eps1e-1.5", "gcv_smooth_eps1e-1"};     
    const std::string eps_string = "1e-1"; 
    

    std::vector<double> lambdas_1; std::vector<double> lambdas_2; std::vector<double> lambdas_5;
    std::vector<double> lambdas_10; std::vector<double> lambdas_25; std::vector<double> lambdas_50;
    std::vector<double> lambdas_75; std::vector<double> lambdas_90; std::vector<double> lambdas_91; 
    std::vector<double> lambdas_92; std::vector<double> lambdas_93;  std::vector<double> lambdas_94; 
    std::vector<double> lambdas_95; std::vector<double> lambdas_96;  std::vector<double> lambdas_97; 
    std::vector<double> lambdas_98; std::vector<double> lambdas_99; 

    // for(double x = -10.0; x <= -3.0; x += 0.2) lambdas_50.push_back(std::pow(10, x));
    // for(double x = -10.0; x <= -3.0; x += 0.2) lambdas_90.push_back(std::pow(10, x)); 
    // for(double x = -10.0; x <= -3.0; x += 0.2) lambdas_95.push_back(std::pow(10, x));
    // for(double x = -10.0; x <= -3.0; x += 0.2) lambdas_99.push_back(std::pow(10, x));

    for(double x = -10.0; x <= -3.0; x += 0.2) lambdas_1.push_back(std::pow(10, x));
    for(double x = -10.0; x <= -3.0; x += 0.2) lambdas_2.push_back(std::pow(10, x));
    for(double x = -10.0; x <= -3.0; x += 0.2) lambdas_5.push_back(std::pow(10, x)); 
    for(double x = -10.0; x <= -3.0; x += 0.2) lambdas_10.push_back(std::pow(10, x)); 
    for(double x = -10.0; x <= -3.0; x += 0.2) lambdas_25.push_back(std::pow(10, x));
    for(double x = -10.0; x <= -3.0; x += 0.2) lambdas_75.push_back(std::pow(10, x)); 
    for(double x = -10.0; x <= -3.0; x += 0.2) lambdas_91.push_back(std::pow(10, x));
    for(double x = -10.0; x <= -3.0; x += 0.2) lambdas_92.push_back(std::pow(10, x));
    for(double x = -10.0; x <= -3.0; x += 0.2) lambdas_93.push_back(std::pow(10, x));
    for(double x = -10.0; x <= -3.0; x += 0.2) lambdas_94.push_back(std::pow(10, x));
    for(double x = -10.0; x <= -3.0; x += 0.2) lambdas_96.push_back(std::pow(10, x));
    for(double x = -10.0; x <= -3.0; x += 0.2) lambdas_97.push_back(std::pow(10, x));
    for(double x = -10.0; x <= -3.0; x += 0.2) lambdas_98.push_back(std::pow(10, x)); 

    double best_lambda; 

    // Read covariates and locations
    DMatrix<double> loc ; 
    if(bool_obs_rip)
        loc = read_csv<double>(R_path + "/data" + "/locs_ripetute.csv"); 
    else
        loc = read_csv<double>(R_path + "/data" + "/locs.csv"); 

    // Simulations 
    const unsigned int n_sim = 10; 
    for(auto sim = 1; sim <= n_sim; ++sim){
        std::cout << "--------------------Simulation #" << std::to_string(sim) << "-------------" << std::endl; 

        // std::string solutions_path_rmse = R_path + "/data/simulations/sim_" + std::to_string(sim) + "/single_est" + pde_type + "/RMSE"; 

        // // K = K_est
        // if(pde_type != "_casc")
        //     std::cout << "ERROR: YOU WANT TO USE K = K_est BUT YOU ARE USING SOMETHING ELSE" << std::endl; 
        // SMatrix<2> K = read_csv<double>(R_path + "/data/simulations/sim_" + std::to_string(sim) + "/K.csv"); 
        // auto L = -diffusion<FEM>(K);   // anisotropic diffusion  
        // PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

        for(auto alpha : alphas){

            unsigned int alpha_int = alpha*100; 
            std::string alpha_string = std::to_string(alpha_int); 

            std::cout << "------------------alpha=" << alpha_string << "-----------------" << std::endl; 

            // load data from .csv files
            DMatrix<double> y; 
            if(bool_obs_rip)
                y = read_csv<double>(R_path + "/data/simulations/sim_" + std::to_string(sim) + "/y_ripetute.csv");
            else
                y = read_csv<double>(R_path + "/data/simulations/sim_" + std::to_string(sim) + "/y.csv");

            // check dimensions
            std::cout << "dim loc " << loc.rows() << " , " << loc.cols() << std::endl;
            std::cout << "dim y " << y.rows() << std::endl;


            BlockFrame<double, int> df;
            df.insert(OBSERVATIONS_BLK, y);

            // GCV:
            for(auto lambda_selection_type : lambda_selection_types){
                
                std::string solutions_path_gcv = R_path + "/data/simulations/sim_" + std::to_string(sim) + "/single_est" + pde_type + "/" + lambda_selection_type; 

                if(bool_obs_rip)
                    solutions_path_gcv = solutions_path_gcv + "/obs_ripetute"; 
                else
                    solutions_path_gcv = solutions_path_gcv + "/obs_singola"; 


                
                QSRPDE<SpaceOnly> model_gcv(problem, Sampling::pointwise, alpha);
                model_gcv.set_spatial_locations(loc);

                std::vector<double> lambdas;
                if(almost_equal(alpha, 0.01)){
                    lambdas = lambdas_1; 
                }  
                if(almost_equal(alpha, 0.02)){
                    lambdas = lambdas_2; 
                }  
                if(almost_equal(alpha, 0.05)){
                    lambdas = lambdas_5; 
                }  
                if(almost_equal(alpha, 0.10)){
                    lambdas = lambdas_10; 
                }  
                if(almost_equal(alpha, 0.25)){
                    lambdas = lambdas_25; 
                }  
                if(almost_equal(alpha, 0.50)){
                    lambdas = lambdas_50; 
                }  
                if(almost_equal(alpha, 0.75)){
                    lambdas = lambdas_75; 
                }  
                if(almost_equal(alpha, 0.90)){
                    lambdas = lambdas_90; 
                } 
                if(almost_equal(alpha, 0.91)){
                    lambdas = lambdas_91; 
                }   
                if(almost_equal(alpha, 0.92)){
                    lambdas = lambdas_92; 
                }  
                if(almost_equal(alpha, 0.93)){
                    lambdas = lambdas_93; 
                } 
                if(almost_equal(alpha, 0.94)){
                    lambdas = lambdas_94; 
                }   
                if(almost_equal(alpha, 0.95)){
                    lambdas = lambdas_95; 
                } 
                if(almost_equal(alpha, 0.96)){
                    lambdas = lambdas_96; 
                }    
                if(almost_equal(alpha, 0.97)){
                    lambdas = lambdas_97; 
                }    
                if(almost_equal(alpha, 0.98)){
                    lambdas = lambdas_98; 
                }  
                if(almost_equal(alpha, 0.99)){
                    lambdas = lambdas_99; 
                }  

                // define lambda sequence as matrix 
                DMatrix<double> lambdas_mat;
                lambdas_mat.resize(lambdas.size(), 1); 
                for(auto i = 0; i < lambdas_mat.rows(); ++i){
                    lambdas_mat(i,0) = lambdas[i]; 
                }


                // set model's data
                if(eps_string == "1e-1"){
                    model_gcv.set_eps_power(-1.0); 
                }
                
                model_gcv.set_data(df);
                model_gcv.init();

                // define GCV function and grid of \lambda_D values
                auto GCV = model_gcv.gcv<ExactEDF>();
                // optimize GCV
                Grid<fdapde::Dynamic> opt;
                opt.optimize(GCV, lambdas_mat);
                
                best_lambda = opt.optimum()(0,0);
        
                std::cout << "Best lambda is: " << std::setprecision(16) << best_lambda << std::endl; 

                // Save lambda sequence 
                std::ofstream fileLambdaS(solutions_path_gcv + "/lambdas_seq_alpha_" + alpha_string + ".csv");
                for(std::size_t i = 0; i < lambdas.size(); ++i) 
                    fileLambdaS << std::setprecision(16) << lambdas[i] << "\n"; 
                fileLambdaS.close();

                // Save lambda GCVopt for all alphas
                std::ofstream fileLambdaoptS(solutions_path_gcv + "/lambdas_opt_alpha_" + alpha_string + ".csv");
                if(fileLambdaoptS.is_open()){
                    fileLambdaoptS << std::setprecision(16) << best_lambda;
                    fileLambdaoptS.close();
                }

                // Save GCV scores
                std::ofstream fileGCV_scores(solutions_path_gcv + "/score_alpha_" + alpha_string + ".csv");
                for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
                    fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.gcvs()[i]) << "\n"; 
                fileGCV_scores.close();

                // Save GCV edf
                std::ofstream fileGCV_edfs(solutions_path_gcv + "/edf_alpha_" + alpha_string + ".csv");
                for(std::size_t i = 0; i < GCV.edfs().size(); ++i) 
                    fileGCV_edfs << std::setprecision(16) << GCV.edfs()[i] << "\n"; 
                fileGCV_edfs.close();
            }

        }


    }
}



//// 

// NUOVI TEST POST-TESI

// test 1
//    domain:       unit square
//    sampling:     locations != nodes
//    penalization: constant PDE coefficients 
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
TEST(gcv_msqrpde_test1, pde_nonparametric_samplingatlocations_spaceonly_gridexact) {

    // path test  
    std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/PhD/Codice/models/MQSRPDE/Tests/Test_1"; 

    // const std::string pde_type = "_lap";    // "_lap" "_Ktrue" "_casc"

    const unsigned int n_sim = 20;
    const std::string gcv_refinement = "fine";    // "lasco" "fine"
    double lambdas_step; 
    if(gcv_refinement == "lasco"){
        lambdas_step = 0.5;
    } 
    if(gcv_refinement == "fine"){
        lambdas_step = 0.1; 
    } 

    // define domain
    MeshLoader<Triangulation<2, 2>> domain("unit_square_test7");  // il dominio è lo stesso del vecchio test 7 

    // rhs 
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);

    // define regularizing PDE

    // // lap 
    // if(pde_type != "_lap")
    //     std::cout << "ERROR: YOU WANT TO USE K = I BUT YOU ARE USING SOMETHING ELSE" << std::endl; 
    // auto L = -laplacian<FEM>();   
    // PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

    // // K = K_true
    // if(pde_type != "_Ktrue")
    //     std::cout << "ERROR: YOU WANT TO USE K = K_true BUT YOU ARE USING SOMETHING ELSE" << std::endl;
    // SMatrix<2> K = read_csv<double>(R_path + "/data/true/K_true.csv"); 
    // auto L = -diffusion<FEM>(K);  
    // PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

    // define statistical model
    std::vector<double> alphas = {0.01, 0.02, 0.03,   // ATT: aggiunto 3%
                                  0.05, 0.10, 0.25, 
                                  0.50, 0.75, 0.90, 0.91, 0.92, 
                                  0.93, 0.94, 0.95, 0.96, 0.97, 
                                  0.98, 0.99};  

    // define grid of lambda values
    std::vector<std::string> lambda_selection_types = {"gcv_smooth_eps1e-1"}; // {"gcv", "gcv_smooth_eps1e-3", "gcv_smooth_eps1e-2", "gcv_smooth_eps1e-1.5", "gcv_smooth_eps1e-1"}; // {"gcv", "gcv_smooth_eps1e-3", "gcv_smooth_eps1e-2", "gcv_smooth_eps1e-1.5", "gcv_smooth_eps1e-1"};     
    std::vector<double> lambdas_1; std::vector<double> lambdas_2; std::vector<double> lambdas_3; 
    std::vector<double> lambdas_5;
    std::vector<double> lambdas_10; std::vector<double> lambdas_25; std::vector<double> lambdas_50;
    std::vector<double> lambdas_75; std::vector<double> lambdas_90; std::vector<double> lambdas_91; 
    std::vector<double> lambdas_92; std::vector<double> lambdas_93;  std::vector<double> lambdas_94; 
    std::vector<double> lambdas_95; std::vector<double> lambdas_96;  std::vector<double> lambdas_97; 
    std::vector<double> lambdas_98; std::vector<double> lambdas_99; 
    for(double x = -7.0; x <= -4.0; x += lambdas_step) lambdas_1.push_back(std::pow(10, x));
    for(double x = -7.0; x <= -4.0; x += lambdas_step) lambdas_2.push_back(std::pow(10, x));
    for(double x = -7.5; x <= -3.5; x += lambdas_step) lambdas_3.push_back(std::pow(10, x));
    for(double x = -7.5; x <= -3.5; x += lambdas_step) lambdas_5.push_back(std::pow(10, x)); 
    for(double x = -7.0; x <= -3.0; x += lambdas_step) lambdas_10.push_back(std::pow(10, x)); 
    for(double x = -6.0; x <= -3.0; x += lambdas_step) lambdas_25.push_back(std::pow(10, x));
    for(double x = -6.0; x <= -3.0; x += lambdas_step) lambdas_50.push_back(std::pow(10, x)); 
    for(double x = -6.0; x <= -3.0; x += lambdas_step) lambdas_75.push_back(std::pow(10, x)); 
    for(double x = -7.5; x <= -3.0; x += lambdas_step) lambdas_90.push_back(std::pow(10, x));
    for(double x = -7.5; x <= -3.0; x += lambdas_step) lambdas_91.push_back(std::pow(10, x));
    for(double x = -7.5; x <= -3.0; x += lambdas_step) lambdas_92.push_back(std::pow(10, x));
    for(double x = -7.5; x <= -3.0; x += lambdas_step) lambdas_93.push_back(std::pow(10, x));
    for(double x = -7.5; x <= -3.5; x += lambdas_step) lambdas_94.push_back(std::pow(10, x));
    for(double x = -7.5; x <= -3.0; x += lambdas_step) lambdas_95.push_back(std::pow(10, x));
    for(double x = -7.5; x <= -4.0; x += lambdas_step) lambdas_96.push_back(std::pow(10, x));
    for(double x = -6.5; x <= -3.0; x += lambdas_step) lambdas_97.push_back(std::pow(10, x));
    for(double x = -6.5; x <= -2.5; x += lambdas_step) lambdas_98.push_back(std::pow(10, x)); 
    for(double x = -6.5; x <= -2.0; x += lambdas_step) lambdas_99.push_back(std::pow(10, x));
    double best_lambda; 

    // Read covariates and locations
    DMatrix<double> loc = read_csv<double>(R_path + "/locs.csv"); 

    // Simulations  
    // std::vector<unsigned int> simulations = {25}; // {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}; 
    for(auto sim = 1; sim <= n_sim; ++sim){
        std::cout << "--------------------Simulation #" << std::to_string(sim) << "-------------" << std::endl; 

        // std::string solutions_path_rmse = R_path + "/data/simulations/sim_" + std::to_string(sim) + "/single_est" + pde_type + "/RMSE"; 
        std::string solutions_path_rmse = R_path + "/simulations/sim_" + std::to_string(sim) + "/single" + "/true_lambda"; 

        // // K = K_est
        // if(pde_type != "_casc")
        //     std::cout << "ERROR: YOU WANT TO USE K = K_est BUT YOU ARE USING SOMETHING ELSE" << std::endl; 
        SMatrix<2> K = read_csv<double>(R_path + "/simulations/sim_" + std::to_string(sim) + "/K.csv"); 

        auto L = -diffusion<FEM>(K);   // anisotropic diffusion  
        PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

        for(auto alpha : alphas){

            unsigned int alpha_int = alpha*100; 
            std::string alpha_string = std::to_string(alpha_int); 

            std::cout << "------------------alpha=" << alpha_string << "-----------------" << std::endl; 

            // load data from .csv files
            DMatrix<double> y = read_csv<double>(R_path + "/simulations/sim_" + std::to_string(sim) + "/y.csv");
            BlockFrame<double, int> df;
            df.insert(OBSERVATIONS_BLK, y);

            // GCV:
            for(auto lambda_selection_type : lambda_selection_types){
                
                // std::string solutions_path_gcv = R_path + "/data/simulations/sim_" + std::to_string(sim) + "/single_est" + pde_type + "/" + lambda_selection_type; 
                std::string solutions_path_gcv = R_path + "/simulations/sim_" + std::to_string(sim) + "/single" + "/" + lambda_selection_type + "/" + gcv_refinement; 
                       
                QSRPDE<SpaceOnly> model_gcv(problem, Sampling::pointwise, alpha);
                model_gcv.set_spatial_locations(loc);

                std::vector<double> lambdas;
                if(almost_equal(alpha, 0.01)){
                    lambdas = lambdas_1; 
                }  
                if(almost_equal(alpha, 0.02)){
                    lambdas = lambdas_2; 
                }
                if(almost_equal(alpha, 0.03)){
                    lambdas = lambdas_3; 
                }  
                if(almost_equal(alpha, 0.05)){
                    lambdas = lambdas_5; 
                }  
                if(almost_equal(alpha, 0.10)){
                    lambdas = lambdas_10; 
                }  
                if(almost_equal(alpha, 0.25)){
                    lambdas = lambdas_25; 
                }  
                if(almost_equal(alpha, 0.50)){
                    lambdas = lambdas_50; 
                }  
                if(almost_equal(alpha, 0.75)){
                    lambdas = lambdas_75; 
                }  
                if(almost_equal(alpha, 0.90)){
                    lambdas = lambdas_90; 
                } 
                if(almost_equal(alpha, 0.91)){
                    lambdas = lambdas_91; 
                }   
                if(almost_equal(alpha, 0.92)){
                    lambdas = lambdas_92; 
                }  
                if(almost_equal(alpha, 0.93)){
                    lambdas = lambdas_93; 
                } 
                if(almost_equal(alpha, 0.94)){
                    lambdas = lambdas_94; 
                }   
                if(almost_equal(alpha, 0.95)){
                    lambdas = lambdas_95; 
                } 
                if(almost_equal(alpha, 0.96)){
                    lambdas = lambdas_96; 
                }    
                if(almost_equal(alpha, 0.97)){
                    lambdas = lambdas_97; 
                }    
                if(almost_equal(alpha, 0.98)){
                    lambdas = lambdas_98; 
                }  
                if(almost_equal(alpha, 0.99)){
                    lambdas = lambdas_99; 
                }  

                // define lambda sequence as matrix 
                DMatrix<double> lambdas_mat;
                lambdas_mat.resize(lambdas.size(), 1); 
                for(auto i = 0; i < lambdas_mat.rows(); ++i){
                    lambdas_mat(i,0) = lambdas[i]; 
                }

                // set model's data
                // model_gcv.set_exact_gcv(lambda_selection_type == "gcv"); 

                if(lambda_selection_type == "gcv_smooth_eps1e-3"){
                    model_gcv.set_eps_power(-3.0); 
                }
                if(lambda_selection_type == "gcv_smooth_eps1e-2"){
                    model_gcv.set_eps_power(-2.0); 
                }
                if(lambda_selection_type == "gcv_smooth_eps1e-1.5"){
                    model_gcv.set_eps_power(-1.5); 
                }
                if(lambda_selection_type == "gcv_smooth_eps1e-1"){
                    model_gcv.set_eps_power(-1.0); 
                }
                
                model_gcv.set_data(df);
                model_gcv.init();

                // define GCV function and grid of \lambda_D values
                auto GCV = model_gcv.gcv<ExactEDF>();
                // optimize GCV
                Grid<fdapde::Dynamic> opt;
                opt.optimize(GCV, lambdas_mat);
                
                best_lambda = opt.optimum()(0,0);
        
                std::cout << "Best lambda is: " << std::setprecision(16) << best_lambda << std::endl; 

                // Save lambda sequence 
                std::ofstream fileLambdaS(solutions_path_gcv + "/lambdas_seq_alpha_" + alpha_string + ".csv");
                for(std::size_t i = 0; i < lambdas.size(); ++i) 
                    fileLambdaS << std::setprecision(16) << lambdas[i] << "\n"; 
                fileLambdaS.close();

                // Save lambda GCVopt for all alphas
                std::ofstream fileLambdaoptS(solutions_path_gcv + "/lambdas_opt_alpha_" + alpha_string + ".csv");
                if(fileLambdaoptS.is_open()){
                    fileLambdaoptS << std::setprecision(16) << best_lambda;
                    fileLambdaoptS.close();
                }

                // Save GCV 
                std::ofstream fileGCV_scores(solutions_path_gcv + "/score_alpha_" + alpha_string + ".csv");
                for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
                    fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.gcvs()[i]) << "\n"; 
                fileGCV_scores.close();
            }

        }


    }
}










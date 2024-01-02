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
using fdapde::core::FEM;
using fdapde::core::fem_order;
using fdapde::core::laplacian;
using fdapde::core::MatrixDataWrapper;
using fdapde::core::PDE;
using fdapde::core::VectorDataWrapper;

#include "../../fdaPDE/models/regression/msqrpde.h"
#include "../../fdaPDE/models/sampling_design.h"
using fdapde::models::Areal;
using fdapde::models::GeoStatLocations;
using fdapde::models::GeoStatMeshNodes;
using fdapde::models::MonolithicSolver;
using fdapde::models::SpaceOnly;
using fdapde::models::MSQRPDE;

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_csv;


// // test 1
// //    domain:       unit square [1,1] x [1,1]
// //    sampling:     locations = nodes
// //    penalization: simple laplacian
// //    covariates:   no
// //    BC:           no
// //    order FE:     1
// TEST(msqrpde_test1, laplacian_nonparametric_samplingatnodes) {

//     // path test  
//     std::string C_path = "/mnt/c/Users/marco/PACS/Project/Code/Cpp/fdaPDE-fork/test/data/models/msqrpde/2D_test1"; 
//     std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/multiple_quantiles/Tests/Test_1"; 

//     // define domain
//     MeshLoader<Mesh2D> domain("unit_square_44");

//     // define regularizing PDE
//     auto L = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
//     PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
//     // define statistical model
//     std::vector<double> alphas = {0.01, 0.05, 0.10, 0.90, 0.95, 0.99}; 

//     const std::string data_type = "hetero"; 

//     // Simulations 
//     const unsigned int M = 1; 
//     for(auto m = 1; m <= M; ++m){

//         MSQRPDE<decltype(problem), SpaceOnly, GeoStatMeshNodes, MonolithicSolver> model(problem, alphas);
//         // use optimal lambda to avoid possible numerical issues
//         DMatrix<double> lambdas = read_csv<double>(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/single_est/lambdas_opt.csv");   // the vector should be saved with the "R format"
//         model.setLambdas_D(lambdas);
//         // load data from .csv files
//         DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/z.csv"); 

//         // set model data
//         BlockFrame<double, int> df;
//         df.insert(OBSERVATIONS_BLK, y);
//         model.set_data(df);

//         // solve smoothing problem
//         model.init();
//         model.solve();

//         // Save solution
//         DMatrix<double> computedF = model.f();
//         const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filef(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/f_all.csv");
//         if (filef.is_open()){
//             filef << computedF.format(CSVFormatf);
//             filef.close();
//         }

//         // // debug 
//         // DMatrix<double> computedA = model.A_mult();
//         // const static Eigen::IOFormat CSVFormatA(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         // std::ofstream fileA(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/A.csv");
//         // if (fileA.is_open()){
//         //   fileA << computedA.format(CSVFormatA);
//         //   fileA.close();
//         // }

//         // DMatrix<double> computedW = model.W_mult();
//         // const static Eigen::IOFormat CSVFormatW(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         // std::ofstream fileW(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/W.csv");
//         // if (fileW.is_open()){
//         //   fileW << computedW.format(CSVFormatW);
//         //   fileW.close();
//         // }

//         // DMatrix<double> computedWbar = model.Wbar_mult();
//         // const static Eigen::IOFormat CSVFormatWbar(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         // std::ofstream fileWbar(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/Wbar.csv");
//         // if (fileWbar.is_open()){
//         //   fileWbar << computedWbar.format(CSVFormatWbar);
//         //   fileWbar.close();
//         // }

//         // DMatrix<double> computedDelta = model.Delta_mult();
//         // const static Eigen::IOFormat CSVFormatDelta(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         // std::ofstream fileDelta(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/Delta.csv");
//         // if (fileDelta.is_open()){
//         //   fileDelta << computedDelta.format(CSVFormatDelta);
//         //   fileDelta.close();
//         // }

//         // DMatrix<double> computedD_script = model.D_script();
//         // const static Eigen::IOFormat CSVFormatD_script(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         // std::ofstream fileD_script(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/D_script.csv");
//         // if (fileD_script.is_open()){
//         //   fileD_script << computedD_script.format(CSVFormatD_script);
//         //   fileD_script.close();
//         // }

//         // DMatrix<double> computedDscriptj = model.Dscriptj();
//         // const static Eigen::IOFormat CSVFormatDscriptj(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         // std::ofstream fileDscriptj(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/D_script_j.csv");
//         // if (fileDscriptj.is_open()){
//         //   fileDscriptj << computedDscriptj.format(CSVFormatDscriptj);
//         //   fileDscriptj.close();
//         // }

//         // DMatrix<double> computedPsi = model.Psi_mult();
//         // const static Eigen::IOFormat CSVFormatPsi(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         // std::ofstream filePsi(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/Psi.csv");
//         // if (filePsi.is_open()){
//         //   filePsi << computedPsi.format(CSVFormatPsi);
//         //   filePsi.close();
//         // }

//     }

// }



// // test 2
// //    domain:       c-shaped
// //    sampling:     locations != nodes
// //    penalization: simple laplacian
// //    covariates:   no
// //    BC:           no
// //    order FE:     1
// TEST(msqrpde_test2, laplacian_nonparametric_samplingatlocations) {

//     // path test  
//     std::string C_path = "/mnt/c/Users/marco/PACS/Project/Code/Cpp/fdaPDE-fork/test/data/models/msqrpde/2D_test2"; 
//     std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/multiple_quantiles/Tests/Test_2"; 

//     // define domain
//     MeshLoader<Mesh2D> domain("c_shaped_adj");

//     // define regularizing PDE
//     auto L = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
//     PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
//     // define statistical model
//     std::vector<double> alphas = {0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99}; 

//     const std::string data_type = "hetero"; 

//     // Read locs
//     DMatrix<double> loc = read_csv<double>(R_path + "/locs.csv"); 

//     // Simulations 
//     const unsigned int M = 1; 

//     // Simultaneous estimation
//     for(auto m = 1; m <= M; ++m){

//         std::cout << "--------------------Simulation #" << std::to_string(m) << "-------------" << std::endl; 

//         MSQRPDE<decltype(problem), SpaceOnly, GeoStatLocations, MonolithicSolver> model(problem, alphas);
//         model.set_spatial_locations(loc);

//         // use optimal lambda to avoid possible numerical issues
//         //std::ifstream fileLambdas(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/single_est/lambdas_opt.csv");
//         // if(fileLambdas.is_open()){
//         //     fileLambdas >> lambdas; 
//         //     fileLambdas.close();
//         // }

//         DMatrix<double> lambdas = read_csv<double>(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/single_est/lambdas_opt.csv"); 
//         model.setLambdas_D(lambdas);

//         // load data from .csv files
//         DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/y.csv");

//         // set model data
//         BlockFrame<double, int> df;
//         df.insert(OBSERVATIONS_BLK, y);
//         model.set_data(df);

//         // solve smoothing problem
//         model.init();
//         model.solve();

//         // Save solution
//         DMatrix<double> computedF = model.f();
//         const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filef(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/mult_est/f_all.csv");
//         if (filef.is_open()){
//             filef << computedF.format(CSVFormatf);
//             filef.close();
//         }

//         DMatrix<double> computedFn = model.Psi_mult()*model.f();
//         const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filefn(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/mult_est/fn_all.csv");
//         if (filefn.is_open()){
//             filefn << computedFn.format(CSVFormatfn);
//             filefn.close();
//         }

//     }

    
//     // Single estimations
//     for(auto m = 1; m <= M; ++m){
//         std::cout << "--------------------Simulation #" << std::to_string(m) << "-------------" << std::endl; 

//         // load data from .csv files
//         DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/y.csv");

//         // model data
//         BlockFrame<double, int> df;
//         df.insert(OBSERVATIONS_BLK, y);

//         unsigned int ind = 0; 
//         for(auto alpha : alphas){

//             unsigned int alpha_int = alpha*100; 
//             std::string alpha_string = std::to_string(alpha_int); 
//             std::cout << "------------------alpha=" << alpha_string << "-----------------" << std::endl; 

//             SQRPDE<decltype(problem), SpaceOnly, GeoStatLocations, MonolithicSolver> model(problem, alpha);
//             model.set_spatial_locations(loc);

//             // use optimal lambda to avoid possible numerical issues
//             // std::ifstream fileLambdas_single(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/single_est/lambdas_opt.csv");
//             // if(fileLambdas_single.is_open()){
//             //     fileLambdas_single >> lambdas_single; 
//             //     fileLambdas_single.close();
//             // } 

//             DMatrix<double> lambdas = read_csv<double>(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/single_est/lambdas_opt.csv"); 
//             model.set_lambda_D(lambdas(ind,0));
//             model.set_data(df);

//             model.init(); // init model
//             model.solve();

//             // Save solution
//             DMatrix<double> computedF = model.f();
//             const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//             std::ofstream filef(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/single_est/f_" +  alpha_string + ".csv");
//             if (filef.is_open()){
//                 filef << computedF.format(CSVFormatf);
//                 filef.close();
//             }

//             DMatrix<double> computedFn = model.Psi()*model.f();
//             const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//             std::ofstream filefn(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/single_est/fn_" +  alpha_string + ".csv");
//             if (filefn.is_open()){
//                 filefn << computedFn.format(CSVFormatfn);
//                 filefn.close();
//             }

//             ind += 1;
//         }


//     }


// }



// // test 3
// //    domain:       c-shaped
// //    sampling:     locations != nodes
// //    penalization: simple laplacian
// //    covariates:   yes
// //    BC:           no
// //    order FE:     1
// TEST(msqrpde_test3, laplacian_semiparametric_samplingatlocations) {

//     // path test  
//     std::string C_path = "/mnt/c/Users/marco/PACS/Project/Code/Cpp/fdaPDE-fork/test/data/models/msqrpde/2D_test3"; 
//     std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/multiple_quantiles/Tests/Test_3"; 

//     // define domain
//     MeshLoader<Mesh2D> domain("c_shaped_adj");

//     // define regularizing PDE
//     auto L = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
//     PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
//     // define statistical model
//     std::vector<double> alphas = {0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99}; 

//     const std::string data_type = "hetero"; 

//     // Read covariates
//     DMatrix<double> X = read_csv<double>(R_path + "/X.csv"); 
//     DMatrix<double> loc = read_csv<double>(R_path + "/locs.csv"); 

//     // Simulations 
//     const unsigned int M = 10; 

//     // Simultaneous estimation
//     for(auto m = 1; m <= M; ++m){

//         std::cout << "--------------------Simulation #" << std::to_string(m) << "-------------" << std::endl; 

//         MSQRPDE<decltype(problem), SpaceOnly, GeoStatLocations, MonolithicSolver> model(problem, alphas);
//         model.set_spatial_locations(loc);

//         // use optimal lambda to avoid possible numerical issues
//         std::ifstream fileLambdas(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/single_est/lambdas_opt.csv");
//         if(fileLambdas.is_open()){
//             fileLambdas >> lambdas; 
//             fileLambdas.close();
//         }
//         model.setLambdas_D(lambdas);

//         // load data from .csv files
//         DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/y.csv");

//         // set model data
//         BlockFrame<double, int> df;
//         df.insert(OBSERVATIONS_BLK, y);
//         df.insert(DESIGN_MATRIX_BLK, X);
//         model.set_data(df);

//         // solve smoothing problem
//         model.init();
//         model.solve();

//         // Save solution
//         DMatrix<double> computedF = model.f();
//         const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filef(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/mult_est/f_all.csv");
//         if (filef.is_open()){
//             filef << computedF.format(CSVFormatf);
//             filef.close();
//         }

//         DMatrix<double> computedFn = model.Psi_mult()*model.f();
//         const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filefn(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/mult_est/fn_all.csv");
//         if (filefn.is_open()){
//             filefn << computedFn.format(CSVFormatfn);
//             filefn.close();
//         }

//         DMatrix<double> computedBeta = model.beta();
//         const static Eigen::IOFormat CSVFormatBeta(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream fileBeta(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/mult_est/beta_all.csv");
//         if (fileBeta.is_open()){
//             fileBeta << computedBeta.format(CSVFormatBeta);
//             fileBeta.close();
//         }
//     }

//     // Single estimations
//     for(auto m = 1; m <= M; ++m){
//         std::cout << "--------------------Simulation #" << std::to_string(m) << "-------------" << std::endl; 

//         // load data from .csv files
//         DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/y.csv");

//         // model data
//         BlockFrame<double, int> df;
//         df.insert(OBSERVATIONS_BLK, y);
//         df.insert(DESIGN_MATRIX_BLK, X);

//         unsigned int ind = 0; 
//         for(auto alpha : alphas){

//             unsigned int alpha_int = alpha*100; 
//             std::string alpha_string = std::to_string(alpha_int); 
//             std::cout << "------------------alpha=" << alpha_string << "-----------------" << std::endl; 

//             SQRPDE<decltype(problem), SpaceOnly, GeoStatLocations, MonolithicSolver> model(problem, alpha);
//             model.set_spatial_locations(loc);

//             // use optimal lambda to avoid possible numerical issues
//             DMatrix<double> lambdas = read_csv<double>(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/single_est/lambdas_opt.csv");   // the vector should be saved with the "R format"
//             model.set_lambda_D(lambdas(ind,0));
 
//             model.set_data(df);

//             model.init(); // init model
//             model.solve();

//             // Save solution
//             DMatrix<double> computedF = model.f();
//             const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//             std::ofstream filef(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/single_est/f_" +  alpha_string + ".csv");
//             if (filef.is_open()){
//                 filef << computedF.format(CSVFormatf);
//                 filef.close();
//             }

//             DMatrix<double> computedFn = model.Psi()*model.f();
//             const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//             std::ofstream filefn(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/single_est/fn_" +  alpha_string + ".csv");
//             if (filefn.is_open()){
//                 filefn << computedFn.format(CSVFormatfn);
//                 filefn.close();
//             }

//             DMatrix<double> computedBeta = model.beta();
//             const static Eigen::IOFormat CSVFormatBeta(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//             std::ofstream fileBeta(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/single_est/beta_" +  alpha_string + ".csv");
//             if (fileBeta.is_open()){
//                 fileBeta << computedBeta.format(CSVFormatBeta);
//                 fileBeta.close();
//             }



//             ind += 1;
//         }


//     }


// }


// // test 6 (for GCV comparison)
// //    domain:         unit square 
// //    sampling: locations != nodes
// //    penalization:   simple laplacian
// //    covariates:     no
// //    BC:             no
// //    order FE:       1

// TEST(sqrpde_test, laplacian_nonparametric_samplingatlocations) {

//     // path test  
//     std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/multiple_quantiles/Tests/Test_6"; 

//     // define domain
//     MeshLoader<Mesh2D> domain("unit_square_test6");
//     const std::string data_type = "hetero_3"; 
//     std::vector<std::string> lambda_selection_types = {"gcv", "gcv_smooth_eps1e-3", "gcv_smooth_eps1e-2", "gcv_smooth_eps1e-1.5", "gcv_smooth_eps1e-1"};     
//     const std::string pde_type = "_lap";    // "_lap" ""

//     // rhs
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
//     // define regularizing PDE
//     // lap 
//     auto L = -laplacian<FEM>();   
//     PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

//     // // K = K_true
//     // SMatrix<2> K;
//     // K << 6, 4, 4, 6;
//     // auto L = -diffusion<FEM>(K);   // anisotropic diffusion  
//     // PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

//     // define statistical model
//     std::vector<double> alphas = {0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99};  

//     // Read locs
//     DMatrix<double> loc = read_csv<double>(R_path + "/locs.csv"); 

//     // Simulations 
//     unsigned int n_sim = 20; 

//     // simulations 
//     for(unsigned int sim = 1; sim <= n_sim; ++sim){
//         std::cout << std::endl << "---------------Simulation #" << sim << "--------------" << std::endl; 
//         // load data from .csv files
//         DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/y.csv");

//         for(double alpha : alphas){

//             unsigned int alpha_int = alpha*100; 
//             std::string alpha_string = std::to_string(alpha_int);
//             std::cout << "--------alpha=" << alpha_string << "%" << std::endl;  

//             BlockFrame<double, int> df;

//             for(auto lambda_selection_type : lambda_selection_types){
                
//                 std::string true_path = R_path + "/data_" + data_type + "/true/f_true_" + alpha_string + ".csv"; 
//                 std::string solutions_path = R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/single_est" + pde_type + "/" + lambda_selection_type; 
//                 double lambda;
//                 std::ifstream fileLambdaS(solutions_path + "/lambdas_opt_alpha_" + alpha_string + ".csv");
//                 if(fileLambdaS.is_open()){
//                     fileLambdaS >> lambda; 
//                     fileLambdaS.close();
//                 }
                

//                 SQRPDE<decltype(problem), SpaceOnly, GeoStatLocations, MonolithicSolver> model(problem, alpha);
//                 model.set_lambda_D(lambda);
//                 model.set_spatial_locations(loc);
//                 // set model data
//                 df.insert(OBSERVATIONS_BLK, y);
//                 model.set_data(df);
//                 // solve smoothing problem
//                 model.init();
//                 model.solve();

//                 // Save C++ solution 
//                 DMatrix<double> computedF = model.f();
//                 const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//                 std::ofstream filef(solutions_path + "/f_" + alpha_string + ".csv");
//                 if(filef.is_open()){
//                     filef << computedF.format(CSVFormatf);
//                     filef.close();
//                 }

//                 DMatrix<double> computedFn =  model.Psi()*model.f();
//                 const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//                 std::ofstream filefn(solutions_path + "/fn_" + alpha_string + ".csv");
//                 if(filefn.is_open()){
//                     filefn << computedFn.format(CSVFormatfn);
//                     filefn.close();
//                 }


//             }

//         }
//       }

// }


// test 6 (run multiple & PP)
//    domain:       unit square
//    sampling:     locations != nodes
//    penalization: constant coefficients PDE
//    covariates:   no
//    BC:           no
//    order FE:     1
TEST(msqrpde_test6, laplacian_nonparametric_samplingatlocations) {

    // path test   
    std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/multiple_quantiles/Tests/Test_6"; 

    // define domain
    MeshLoader<Mesh2D> domain("unit_square_test6");
    const std::string data_type = "hetero_3";
    const std::string lambda_selection = "gcv_smooth_eps1e-1"; 
    const std::string pde_type = "";    // "_lap" ""
    const bool single_est = false;
    const bool mult_est = true; 

    // rhs 
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);

    // define regularizing PDE
    // // lap 
    // auto L = -laplacian<FEM>();   
    // PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

    // K = K_true
    SMatrix<2> K;
    K << 6, 4, 4, 6;
    auto L = -diffusion<FEM>(K);   // anisotropic diffusion  
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

    // define statistical model
    std::vector<double> alphas = {0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99};  

    // Read locs
    DMatrix<double> loc = read_csv<double>(R_path + "/locs.csv"); 

    // Simulations 
    const unsigned int n_sim = 20; 

    // Single estimations
    if(single_est){
        std::cout << "-----------------------SINGLE running---------------" << std::endl;
        for(auto sim = 1; sim <= n_sim; ++sim){

                std::cout << "--------------------Simulation #" << std::to_string(sim) << "-------------" << std::endl; 

                // load data from .csv files
                DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/y.csv");
                unsigned int idx = 0; 

                for(double alpha : alphas){
                    SQRPDE<decltype(problem), SpaceOnly, GeoStatLocations, MonolithicSolver> model(problem, alpha);
                    model.set_spatial_locations(loc);
                    unsigned int alpha_int = alphas[idx]*100;  
                    double lambda; 
                    std::ifstream fileLambda(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/single_est" + pde_type +  "/" + lambda_selection + "/lambdas_opt_alpha_" + std::to_string(alpha_int) + ".csv");
                    if(fileLambda.is_open()){
                        fileLambda >> lambda; 
                        fileLambda.close();
                    }
                    model.set_lambda_D(lambda);

                    // set model data
                    BlockFrame<double, int> df;
                    df.insert(OBSERVATIONS_BLK, y);
                    model.set_data(df);

                    // solve smoothing problem
                    model.init();
                    model.solve();

                    // Save solution
                    DMatrix<double> computedF = model.f();
                    const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
                    std::ofstream filef(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/single_est" + pde_type + "/" + lambda_selection + "/f_" + std::to_string(alpha_int) + ".csv");
                    if(filef.is_open()){
                        filef << computedF.format(CSVFormatf);
                        filef.close();
                    }

                    DMatrix<double> computedFn = model.Psi()*model.f();
                    const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
                    std::ofstream filefn(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/single_est" + pde_type + "/" + lambda_selection + "/fn_" + std::to_string(alpha_int) + ".csv");
                    if(filefn.is_open()){
                        filefn << computedFn.format(CSVFormatfn);
                        filefn.close();
                    }

                    idx++;
                }

        }


    }


    // Simultaneous estimations
    if(mult_est){
        for(std::string method : {"PP_new"}){   // "mult", "PP", "PP_new"
            bool force; 
            bool processing; 
            if(method == "mult"){
                processing = false; 
                force = false; 
                std::cout << "-------------------------MULTIPLE running-----------------" << std::endl;
            }
            if(method == "PP"){
                processing = true; 
                force = false; 
                std::cout << "-------------------------PP running-----------------" << std::endl;
            }
            if(method == "PP_new"){
                processing = true; 
                force = true; 
                std::cout << "-------------------------PP new running-----------------" << std::endl;
            }

            for(auto sim = 1; sim <= n_sim; ++sim){

                std::cout << "--------------------Simulation #" << std::to_string(sim) << "-------------" << std::endl; 

                MSQRPDE<decltype(problem), SpaceOnly, GeoStatLocations, MonolithicSolver> model(problem, alphas);
                model.set_spatial_locations(loc);
                model.set_preprocess_option(processing); 
                model.set_forcing_option(force);

                // use optimal lambda to avoid possible numerical issues
                DMatrix<double> lambdas;
                DVector<double> lambdas_temp; 
                lambdas_temp.resize(alphas.size());
                for(std::size_t idx = 0; idx < alphas.size(); ++idx){
                    unsigned int alpha_int = alphas[idx]*100;  
                    std::ifstream fileLambdas(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/single_est" + pde_type + "/" + lambda_selection + "/lambdas_opt_alpha_" + std::to_string(alpha_int) + ".csv");
                    if(fileLambdas.is_open()){
                        fileLambdas >> lambdas_temp(idx); 
                        fileLambdas.close();
                    }
                }
                lambdas = lambdas_temp;
                // std::cout << "lam dim = " << lambdas.rows() << " " << lambdas.cols() << std::endl; 
                // for(auto i = 0; i < lambdas.rows(); ++i)
                //     std::cout << lambdas(i,0) << std::endl; 

                //DMatrix<double> lambdas = read_csv<double>(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/single_est/lambdas_opt.csv"); 
                
                model.setLambdas_D(lambdas);

                // load data from .csv files
                DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/y.csv");

                // set model data
                BlockFrame<double, int> df;
                df.insert(OBSERVATIONS_BLK, y);
                model.set_data(df);

                // solve smoothing problem
                model.init();
                model.solve();

                // Save solution
                if(method == "mult"){
                    DMatrix<double> computedF = model.f();
                    const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
                    std::ofstream filef(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/mult_est" + pde_type + "/f_all.csv");
                    if(filef.is_open()){
                        filef << computedF.format(CSVFormatf);
                        filef.close();
                    }

                    DMatrix<double> computedFn = model.Psi_mult()*model.f();
                    const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
                    std::ofstream filefn(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/mult_est" + pde_type + "/fn_all.csv");
                    if(filefn.is_open()){
                        filefn << computedFn.format(CSVFormatfn);
                        filefn.close();
                    }
                } 
                if(method == "PP"){
                    DMatrix<double> computedF = model.f();
                    const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
                    std::ofstream filef(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/competitors/f_all" + pde_type + "_postproc.csv");
                    if(filef.is_open()){
                        filef << computedF.format(CSVFormatf);
                        filef.close();
                    }

                    DMatrix<double> computedFn = model.Psi_mult()*model.f();
                    const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
                    std::ofstream filefn(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/competitors/fn_all" + pde_type + "_postproc.csv");
                    if(filefn.is_open()){
                        filefn << computedFn.format(CSVFormatfn);
                        filefn.close();
                    }

                }
                if(method == "PP_new"){
                    DMatrix<double> computedF = model.f();
                    const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
                    std::ofstream filef(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/competitors/f_all" + pde_type + "_postproc_new.csv");
                    if(filef.is_open()){
                        filef << computedF.format(CSVFormatf);
                        filef.close();
                    }

                    DMatrix<double> computedFn = model.Psi_mult()*model.f();
                    const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
                    std::ofstream filefn(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/competitors/fn_all" + pde_type + "_postproc_new.csv");
                    if(filefn.is_open()){
                        filefn << computedFn.format(CSVFormatfn);
                        filefn.close();
                    }

                }


            }

        }

    }


}

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

double RMSE_metric(DVector<double> v1, DVector<double> v2){
    double res = 0.; 
    if(v1.size() != v2.size())
        std::cout << std::endl << "----------ERROR IN RMSE COMPUTATION---------" << std::endl; 
    for(auto i = 0; i < v1.size(); ++i){
        res += (v1[i]-v2[i])*(v1[i]-v2[i]); 
    }
    return std::sqrt(1./(v1.size())*res); 
}

// // test 1
// //    domain:       unit square [0,1] x [0,1] 
// //    sampling:     locations = nodes
// //    penalization: simple laplacian
// //    covariates:   no
// //    BC:           no
// //    order FE:     1
// //    GCV optimization: grid exact
// TEST(gcv_msqrpde_test1, laplacian_nonparametric_samplingatnodes_spaceonly_gridexact) {
//     // path test  
//     std::string C_path = "/mnt/c/Users/marco/OneDrive Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/fdaPDE_official-fork/fdaPDE-cpp-sqrpde/test/data/models/msqrpde/2D_test1"; 
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

//     // define grid of lambda values
//     std::vector<SVector<1>> lambdas;
//     for(double x = -7.2; x <= -6.3; x +=0.1) lambdas.push_back(SVector<1>(std::pow(10,x)));
//     DVector<double> best_lambda;
//     best_lambda.resize(alphas.size());  

//     // Simulations 
//     const unsigned int M = 1; 

//     for(auto m = 1; m <= M; ++m){
//         unsigned int ind = 0; 
//         std::ofstream fileGCV(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/single_est/gcv_scores.csv");
//         for(auto alpha : alphas){

//             std::cout << "------------------alpha=" << std::to_string(alpha) << "-----------------" << std::endl; 

//             SQRPDE<decltype(problem), SpaceOnly, GeoStatMeshNodes, MonolithicSolver> model(problem, alpha);

//             // load data from .csv files
//             DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/z.csv"); 

//             // set model data
//             BlockFrame<double, int> df;
//             df.insert(OBSERVATIONS_BLK, y);
//             model.set_data(df);

//             model.init(); // init model

//             // define GCV function and optimize
//             GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
//             Grid<1> opt;

//             ScalarField<1, decltype(GCV)> obj(GCV);
//             opt.optimize(obj, lambdas); // optimize gcv field
//             std::cout << "opt: " << opt.optimum()[0] << std::endl; 
//             best_lambda[ind] = opt.optimum()[0];
            
//             std::cout << "Best lambda is: " << std::setprecision(16) << best_lambda[ind] << std::endl; 
//             ind++;

//             // gcv scores
//             for(std::size_t i = 0; i < GCV.values().size(); ++i) 
//                 fileGCV << std::setprecision(16) << std::sqrt(GCV.values()[i]) << "\n"; 

//         }

//         const static Eigen::IOFormat CSVFormatL(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream fileL(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/single_est/lambdas_opt.csv");
//         if (fileL.is_open()){
//             fileL << best_lambda.format(CSVFormatL);
//             fileL.close();
//         }

//         fileGCV.close(); 
//     }
// }



// // test 2
// //    domain:       c-shaped
// //    sampling:     locations != nodes
// //    penalization: simple laplacian 
// //    covariates:   no
// //    BC:           no
// //    order FE:     1
// //    GCV optimization: grid exact
// TEST(gcv_msqrpde_test2, spacevaryingpde_nonparametric_samplingatlocations_spaceonly_gridexact) {

//     // path test  
//     std::string C_path = "/mnt/c/Users/marco/OneDrive Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/fdaPDE_official-fork/fdaPDE-cpp-sqrpde/test/data/models/msqrpde/2D_test2"; 
//     std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/multiple_quantiles/Tests/Test_2"; 
//     // define domain
//     MeshLoader<Mesh2D> domain("c_shaped_adj");
//     // define regularizing PDE
//     auto L = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
//     PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
//     // define statistical model
//     std::vector<double> alphas =  {0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99}; 

//     const std::string data_type = "hetero"; 

//     // define grid of lambda values
//     std::vector<SVector<1>> lambdas;
//     for(double x = -4.0; x <= -0.1; x +=0.2) lambdas.push_back(SVector<1>(std::pow(10,x)));
//     DVector<double> best_lambda;
//     best_lambda.resize(alphas.size());  

//     // Read covariates and locations
//     DMatrix<double> loc = read_csv<double>(R_path + "/locs.csv"); 

//     // Simulations 
//     const unsigned int M = 1; 
//     for(auto m = 1; m <= M; ++m){
//         std::cout << "--------------------Simulation #" << std::to_string(m) << "-------------" << std::endl; 
//         unsigned int ind = 0; 
//         std::ofstream fileGCV(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/single_est/gcv_scores.csv");
//         for(auto alpha : alphas){

//             unsigned int alpha_int = alpha*100; 
//             std::string alpha_string = std::to_string(alpha_int); 

//             std::cout << "------------------alpha=" << alpha_string << "-----------------" << std::endl; 

//             SQRPDE<decltype(problem), SpaceOnly, GeoStatLocations, MonolithicSolver> model(problem, alpha);
//             model.set_spatial_locations(loc);

//             // load data from .csv files
//             DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/y.csv");

//             // set model data
//             BlockFrame<double, int> df;
//             df.insert(OBSERVATIONS_BLK, y);
//             model.set_data(df);

//             model.init(); // init model

//             // define GCV function and optimize
//             GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
//             Grid<1> opt;

//             ScalarField<1, decltype(GCV)> obj(GCV);
//             opt.optimize(obj, lambdas); // optimize gcv field
//             best_lambda[ind] = opt.optimum()[0];
            
//             std::cout << "Best lambda is: " << std::setprecision(16) << best_lambda[ind] << std::endl; 
//             ind++;

//             // gcv scores
//             for(std::size_t i = 0; i < GCV.values().size(); ++i){
//                 fileGCV << std::setprecision(16) << std::sqrt(GCV.values()[i]) << "\n"; 
//             }

//         }

//         const static Eigen::IOFormat CSVFormatL(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream fileL(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/single_est/lambdas_opt.csv");
//         if (fileL.is_open()){
//             fileL << best_lambda.format(CSVFormatL);
//             fileL.close();
//         }

//         fileGCV.close(); 
//     }

// }



// // test 3
// //    domain:       c-shaped
// //    sampling:     locations != nodes
// //    penalization: simple laplacian 
// //    covariates:   yes
// //    BC:           no
// //    order FE:     1
// //    GCV optimization: grid exact
// TEST(gcv_msqrpde_test3, spacevaryingpde_semiparametric_samplingatlocations_spaceonly_gridexact) {

//     // path test  
//     std::string C_path = "/mnt/c/Users/marco/OneDrive Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/fdaPDE_official-fork/fdaPDE-cpp-sqrpde/test/data/models/msqrpde/2D_test3"; 
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

//     // define grid of lambda values
//     std::vector<SVector<1>> lambdas;
//     for(double x = -5.0; x <= -0.1; x +=0.2) lambdas.push_back(SVector<1>(std::pow(10,x)));
//     DVector<double> best_lambda;
//     best_lambda.resize(alphas.size());  

//     // Read covariates and locations
//     DMatrix<double> X = read_csv<double>(R_path + "/X.csv"); 
//     DMatrix<double> loc = read_csv<double>(R_path + "/locs.csv"); 

//     // Simulations 
//     const unsigned int M = 10; 
//     for(auto m = 1; m <= M; ++m){
//         std::cout << "--------------------Simulation #" << std::to_string(m) << "-------------" << std::endl; 
//         unsigned int ind = 0; 
//         std::ofstream fileGCV(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/single_est/gcv_scores.csv");
//         for(auto alpha : alphas){

//             unsigned int alpha_int = alpha*100; 
//             std::string alpha_string = std::to_string(alpha_int); 

//             std::cout << "------------------alpha=" << alpha_string << "-----------------" << std::endl; 

//             SQRPDE<decltype(problem), SpaceOnly, GeoStatLocations, MonolithicSolver> model(problem, alpha);
//             model.set_spatial_locations(loc);

//             // load data from .csv files
//             DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/y.csv");

//             // set model data
//             BlockFrame<double, int> df;
//             df.insert(OBSERVATIONS_BLK, y);
//             df.insert(DESIGN_MATRIX_BLK, X);
//             model.set_data(df);

//             model.init(); // init model

//             // define GCV function and optimize
//             GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
//             Grid<1> opt;

//             ScalarField<1, decltype(GCV)> obj(GCV);
//             opt.optimize(obj, lambdas); // optimize gcv field
//             best_lambda[ind] = opt.optimum()[0];
            
//             std::cout << "Best lambda is: " << std::setprecision(16) << best_lambda[ind] << std::endl; 
//             ind++;

//             // gcv scores
//             for(std::size_t i = 0; i < GCV.values().size(); ++i){
//                 fileGCV << std::setprecision(16) << std::sqrt(GCV.values()[i]) << "\n"; 
//             }

//         }

//         const static Eigen::IOFormat CSVFormatL(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream fileL(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/single_est/lambdas_opt.csv");
//         if (fileL.is_open()){
//             fileL << best_lambda.format(CSVFormatL);
//             fileL.close();
//         }

//         fileGCV.close(); 
//     }

// }



// // test 6
// //    domain:       unit square
// //    sampling:     locations != nodes
// //    penalization: constant PDE coefficients 
// //    covariates:   no
// //    BC:           no
// //    order FE:     1
// //    GCV optimization: grid exact
// TEST(gcv_msqrpde_test6, pde_nonparametric_samplingatlocations_spaceonly_gridexact) {

//     // path test  
//     // std::string C_path = "/mnt/c/Users/marco/OneDrive Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/fdaPDE_official-fork/fdaPDE-cpp-sqrpde/test/data/models/msqrpde/2D_test6"; 
//     // std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/multiple_quantiles/Tests/Test_6"; 

//     std::string R_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/models/multiple_quantiles/Tests/Test_6"; 
    
//     const std::string data_type = "hetero_3";  

//     // define domain
//     MeshLoader<Mesh2D> domain("unit_square_test6");
//     // define regularizing PDE
//     //auto L = -laplacian<FEM>();   // hetero 

//     // hetero_2
//     SMatrix<2> K;
//     K << 6, 4, 4, 6;
//     auto L = -diffusion<FEM>(K);   // anisotropic diffusion

          
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
//     PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
//     // define statistical model
//     std::vector<double> alphas = {0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99};   // quantiles_1
//     //std::vector<double> alphas = {0.90, 0.91, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99};   // quantiles_2

//     // define grid of lambda values
//     std::vector<SVector<1>> lambdas;
//     for(double x = -8.0; x <= -5.2; x +=0.05) lambdas.push_back(SVector<1>(std::pow(10,x)));
//     DVector<double> best_lambda;
//     best_lambda.resize(alphas.size());  

//     // Read covariates and locations
//     DMatrix<double> loc = read_csv<double>(R_path + "/locs.csv"); 

//     // Simulations 
//     const unsigned int M = 20; 
//     for(auto m = 1; m <= M; ++m){
//         std::cout << "--------------------Simulation #" << std::to_string(m) << "-------------" << std::endl; 
//         unsigned int ind = 0; 
//         std::ofstream fileGCV(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/single_est/gcv_scores.csv");
//         for(auto alpha : alphas){

//             unsigned int alpha_int = alpha*100; 
//             std::string alpha_string = std::to_string(alpha_int); 

//             std::cout << "------------------alpha=" << alpha_string << "-----------------" << std::endl; 

//             SQRPDE<decltype(problem), SpaceOnly, GeoStatLocations, MonolithicSolver> model(problem, alpha);
//             model.set_spatial_locations(loc);

//             // load data from .csv files
//             DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/y.csv");

//             // set model data
//             BlockFrame<double, int> df;
//             df.insert(OBSERVATIONS_BLK, y);
//             model.set_data(df);

//             model.init(); // init model

//             // define GCV function and optimize
//             GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
//             Grid<1> opt;

//             ScalarField<1, decltype(GCV)> obj(GCV);
//             opt.optimize(obj, lambdas); // optimize gcv field
//             best_lambda[ind] = opt.optimum()[0];
            
//             std::cout << "Best lambda is: " << std::setprecision(16) << best_lambda[ind] << std::endl; 
//             ind++;

//             // gcv scores
//             for(std::size_t i = 0; i < GCV.values().size(); ++i){
//                 fileGCV << std::setprecision(16) << std::sqrt(GCV.values()[i]) << "\n"; 
//             }

//         }

//         const static Eigen::IOFormat CSVFormatL(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream fileL(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(m) + "/single_est/lambdas_opt.csv");
//         if (fileL.is_open()){
//             fileL << best_lambda.format(CSVFormatL);
//             fileL.close();
//         }

//         fileGCV.close(); 
//     }

// }



// New test 6
// test 6
//    domain:       unit square
//    sampling:     locations != nodes
//    penalization: constant PDE coefficients 
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
TEST(gcv_msqrpde_test6, pde_nonparametric_samplingatlocations_spaceonly_gridexact) {

    // path test  
    std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/multiple_quantiles/Tests/Test_6"; 
    //std::string R_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/models/multiple_quantiles/Tests/Test_6"; 

    const std::string data_type = "hetero_3";  
    const std::string pde_type = "_lap";    // "_lap" ""

    // define domain
    MeshLoader<Mesh2D> domain("unit_square_test6");

    // rhs 
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);

    // define regularizing PDE

    // lap 
    auto L = -laplacian<FEM>();   
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

    // // K = K_true
    // SMatrix<2> K;
    // K << 6, 4, 4, 6;
    // auto L = -diffusion<FEM>(K);   // anisotropic diffusion  
    // PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);


    // define statistical model
    std::vector<double> alphas = {0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99};  

    // define grid of lambda values
    std::vector<std::string> lambda_selection_types = {"gcv", "gcv_smooth_eps1e-3", "gcv_smooth_eps1e-2", "gcv_smooth_eps1e-1.5", "gcv_smooth_eps1e-1"}; // {"gcv", "gcv_smooth_eps1e-3", "gcv_smooth_eps1e-2", "gcv_smooth_eps1e-1.5", "gcv_smooth_eps1e-1"};     
    bool compute_rmse = false; 
    bool compute_gcv = true; 
    std::vector<DVector<double>> lambdas_1;
    std::vector<DVector<double>> lambdas_5;
    std::vector<DVector<double>> lambdas_10;
    std::vector<DVector<double>> lambdas_25;
    std::vector<DVector<double>> lambdas_50;
    std::vector<DVector<double>> lambdas_75;
    std::vector<DVector<double>> lambdas_90;
    std::vector<DVector<double>> lambdas_95;
    std::vector<DVector<double>> lambdas_99;
    std::vector<double> lambdas_rmse;
    for(double x = -7.5; x <= -3.5; x += 0.05) lambdas_1.push_back(SVector<1>(std::pow(10, x)));
    for(double x = -8.0; x <= -4.0; x += 0.05) lambdas_5.push_back(SVector<1>(std::pow(10, x))); 
    for(double x = -7.5; x <= -3.0; x += 0.05) lambdas_10.push_back(SVector<1>(std::pow(10, x))); 
    for(double x = -7.0; x <= -3.5; x += 0.05) lambdas_25.push_back(SVector<1>(std::pow(10, x)));
    for(double x = -6.0; x <= -3.0; x += 0.05) lambdas_50.push_back(SVector<1>(std::pow(10, x))); 
    for(double x = -6.0; x <= -3.0; x += 0.05) lambdas_75.push_back(SVector<1>(std::pow(10, x))); 
    for(double x = -7.5; x <= -3.5; x += 0.05) lambdas_90.push_back(SVector<1>(std::pow(10, x)));
    for(double x = -7.5; x <= -3.5; x += 0.05) lambdas_95.push_back(SVector<1>(std::pow(10, x))); 
    for(double x = -8.0; x <= -3.5; x += 0.05) lambdas_99.push_back(SVector<1>(std::pow(10, x)));

    for(double x = -9.0; x <= -4.0; x += 0.05) lambdas_rmse.push_back(std::pow(10, x));   // for all alphas 

    double best_lambda;

    // Read covariates and locations
    DMatrix<double> loc = read_csv<double>(R_path + "/locs.csv"); 

    // Simulations 
    const unsigned int n_sim = 20; 
    for(auto sim = 1; sim <= n_sim; ++sim){
        std::cout << "--------------------Simulation #" << std::to_string(sim) << "-------------" << std::endl; 

        std::string solutions_path_rmse = R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/single_est" + pde_type + "/RMSE"; 
        for(auto alpha : alphas){

            unsigned int alpha_int = alpha*100; 
            std::string alpha_string = std::to_string(alpha_int); 

            std::cout << "------------------alpha=" << alpha_string << "-----------------" << std::endl; 

            // load data from .csv files
            DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/y.csv");
            BlockFrame<double, int> df;
            df.insert(OBSERVATIONS_BLK, y);

            if(compute_gcv){

                // GCV:
                for(auto lambda_selection_type : lambda_selection_types){
                    
                    std::string solutions_path_gcv = R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/single_est" + pde_type + "/" + lambda_selection_type; 
                    
                    QSRPDE<SpaceOnly> model_gcv(problem, Sampling::pointwise, alpha);
                    model_gcv.set_spatial_locations(loc);

                    std::vector<DVector<double>> lambdas;
                    if(almost_equal(alpha, 0.01)){
                        lambdas = lambdas_1; 
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
                    if(almost_equal(alpha, 0.95)){
                        lambdas = lambdas_95; 
                    }  
                    if(almost_equal(alpha, 0.99)){
                        lambdas = lambdas_99; 
                    }  

                    // set model's data
                    model_gcv.set_exact_gcv(lambda_selection_type == "gcv"); 

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
                    opt.optimize(GCV, lambdas);
                    
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

            if(compute_rmse){
                std::cout << "-----RMSE computation-----" << std::endl; 
                // RMSE
                DMatrix<double> f_true = read_csv<double>(R_path + "/data_" + data_type + "/true/f_true_" + alpha_string + ".csv");

                std::vector<double> rmse_score; 
                rmse_score.resize(lambdas_rmse.size()); 
                double count_l = 0; 
                for(auto lambda : lambdas_rmse){
                    QSRPDE<SpaceOnly> model(problem, Sampling::pointwise, alpha);
                    // set model's data
                    model.set_spatial_locations(loc);
                    model.set_lambda_D(lambda);           
                    
                    model.set_data(df);
                    model.init();
                    model.solve();

                    rmse_score[count_l] = RMSE_metric(model.f(), f_true); 

                    count_l = count_l+1; 
                }

                auto min_idx = std::distance(std::begin(rmse_score), std::min_element(std::begin(rmse_score), std::end(rmse_score))); 
                
                // Save lambda sequence 
                std::ofstream fileLambdaS_rmse(solutions_path_rmse + "/lambdas_seq_alpha_" + alpha_string + ".csv");
                for(std::size_t i = 0; i < lambdas_rmse.size(); ++i) 
                    fileLambdaS_rmse << std::setprecision(16) << lambdas_rmse[i] << "\n"; 
                fileLambdaS_rmse.close();

                // Save lambda RMSEopt for all alphas
                std::ofstream fileLambdaoptS_rmse(solutions_path_rmse + "/lambdas_opt_alpha_" + alpha_string + ".csv");
                if(fileLambdaoptS_rmse.is_open()){
                    fileLambdaoptS_rmse << std::setprecision(16) << lambdas_rmse[min_idx]; ;
                    fileLambdaoptS_rmse.close();
                }

                // Save score 
                std::ofstream fileRMSE_scores(solutions_path_rmse + "/score_alpha_" + alpha_string + ".csv");
                for(std::size_t i = 0; i < rmse_score.size(); ++i) 
                    fileRMSE_scores << std::setprecision(16) << rmse_score[i] << "\n"; 
                fileRMSE_scores.close();
            
            }

        }


    }
}

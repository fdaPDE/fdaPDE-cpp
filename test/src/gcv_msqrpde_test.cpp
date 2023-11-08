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
using fdapde::core::MatrixDataWrapper;
using fdapde::core::PDE;
using fdapde::core::VectorDataWrapper;

#include "../../fdaPDE/models/regression/msqrpde.h"
#include "../../fdaPDE/models/sampling_design.h"
using fdapde::models::Areal;
using fdapde::models::GeoStatLocations;
using fdapde::models::GeoStatMeshNodes;
using fdapde::models::MSQRPDE;
using fdapde::models::SpaceOnly;
using fdapde::models::MonolithicSolver;

#include "../../fdaPDE/calibration/gcv.h"
using fdapde::calibration::ExactEDF;
using fdapde::calibration::ExactGCV;
using fdapde::calibration::GCV;
using fdapde::calibration::StochasticEDF;

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_mtx;
using fdapde::testing::read_csv;

// test 1
//    domain:       unit square [0,1] x [0,1] 
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
TEST(gcv_msqrpde_test1, laplacian_nonparametric_samplingatnodes_spaceonly_gridexact) {
    // path test  
    std::string C_path = "/mnt/c/Users/marco/OneDrive Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/fdaPDE_official-fork/fdaPDE-cpp-sqrpde/test/data/models/msqrpde/2D_test1"; 
    std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/multiple_quantiles/Tests/Test_1"; 
    // define domain
    MeshLoader<Mesh2D> domain("unit_square_44");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define statistical model
    std::vector<double> alphas = {0.01, 0.05, 0.10, 0.90, 0.95, 0.99}; 

    const std::string data_type = "hetero"; 

    // define grid of lambda values
    std::vector<SVector<1>> lambdas;
    for(double x = -7.2; x <= -6.3; x +=0.1) lambdas.push_back(SVector<1>(std::pow(10,x)));
    DVector<double> best_lambda;
    best_lambda.resize(alphas.size());  

    // Simulations 
    const unsigned int M = 10; 

    for(auto m = 1; m <= M; ++m){
        unsigned int ind = 0; 
        std::ofstream fileGCV(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/single_est/gcv_scores.csv");
        for(auto alpha : alphas){

            std::cout << "------------------alpha=" << std::to_string(alpha) << "-----------------" << std::endl; 

            SQRPDE<decltype(problem), SpaceOnly, GeoStatMeshNodes, MonolithicSolver> model(problem, alpha);

            // load data from .csv files
            DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/y.csv"); 

            // set model data
            BlockFrame<double, int> df;
            df.insert(OBSERVATIONS_BLK, y);
            model.set_data(df);

            model.init(); // init model

            // define GCV function and optimize
            GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
            Grid<1> opt;

            ScalarField<1, decltype(GCV)> obj(GCV);
            opt.optimize(obj, lambdas); // optimize gcv field
            std::cout << "opt: " << opt.optimum()[0] << std::endl; 
            best_lambda[ind] = opt.optimum()[0];
            
            std::cout << "Best lambda is: " << std::setprecision(16) << best_lambda[ind] << std::endl; 
            ind++;

            // gcv scores
            for(std::size_t i = 0; i < GCV.values().size(); ++i) 
                fileGCV << std::setprecision(16) << std::sqrt(GCV.values()[i]) << "\n"; 

        }

        const static Eigen::IOFormat CSVFormatL(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        std::ofstream fileL(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/single_est/lambdas_opt.csv");
        if (fileL.is_open()){
            fileL << best_lambda.format(CSVFormatL);
            fileL.close();
        }

        fileGCV.close(); 
    }
}



// test 2
//    domain:       unit square [0,1] x [0,1] 
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
TEST(gcv_msqrpde_test2, laplacian_semiparametric_samplingatnodes_spaceonly_gridexact) {

    // path test  
    std::string C_path = "/mnt/c/Users/marco/OneDrive Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/fdaPDE_official-fork/fdaPDE-cpp-sqrpde/test/data/models/msqrpde/2D_test2"; 
    std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/multiple_quantiles/Tests/Test_2"; 
    // define domain
    MeshLoader<Mesh2D> domain("unit_square_32");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define statistical model
    std::vector<double> alphas = {0.01, 0.05, 0.10, 0.90, 0.95, 0.99}; 

    const std::string data_type = "hetero"; 

    // define grid of lambda values
    std::vector<SVector<1>> lambdas;
    for(double x = -6.9; x <= -5.8; x +=0.1) lambdas.push_back(SVector<1>(std::pow(10,x)));
    DVector<double> best_lambda;
    best_lambda.resize(alphas.size());  

    // Read covariates
    DMatrix<double> X = read_csv<double>(R_path + "/data_" + data_type + "/X.csv"); 

    // Simulations 
    const unsigned int M = 10; 

    for(auto m = 1; m <= M; ++m){
        std::cout << "--------------------Simulation #" << std::to_string(m) << "-------------" << std::endl; 
        unsigned int ind = 0; 
        std::ofstream fileGCV(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/single_est/gcv_scores.csv");
        for(auto alpha : alphas){

            std::cout << "------------------alpha=" << std::to_string(alpha) << "-----------------" << std::endl; 

            SQRPDE<decltype(problem), SpaceOnly, GeoStatMeshNodes, MonolithicSolver> model(problem, alpha);

            // load data from .csv files
            DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/y.csv");

            // set model data
            BlockFrame<double, int> df;
            df.insert(OBSERVATIONS_BLK, y);
            df.insert(DESIGN_MATRIX_BLK, X);
            model.set_data(df);

            model.init(); // init model

            // define GCV function and optimize
            GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
            Grid<1> opt;

            ScalarField<1, decltype(GCV)> obj(GCV);
            opt.optimize(obj, lambdas); // optimize gcv field
            std::cout << "opt: " << opt.optimum()[0] << std::endl; 
            best_lambda[ind] = opt.optimum()[0];
            
            std::cout << "Best lambda is: " << std::setprecision(16) << best_lambda[ind] << std::endl; 
            ind++;

            // gcv scores
            for(std::size_t i = 0; i < GCV.values().size(); ++i) 
                fileGCV << std::setprecision(16) << std::sqrt(GCV.values()[i]) << "\n"; 

        }

        const static Eigen::IOFormat CSVFormatL(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        std::ofstream fileL(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/single_est/lambdas_opt.csv");
        if (fileL.is_open()){
            fileL << best_lambda.format(CSVFormatL);
            fileL.close();
        }

        fileGCV.close(); 
    }

}



// test 3
//    domain:       c-shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian 
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
TEST(gcv_msqrpde_test3, spacevaryingpde_semiparametric_samplingatlocations_spaceonly_gridexact) {

    // path test  
    std::string C_path = "/mnt/c/Users/marco/OneDrive Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/fdaPDE_official-fork/fdaPDE-cpp-sqrpde/test/data/models/msqrpde/2D_test3"; 
    std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/multiple_quantiles/Tests/Test_2"; 
    // define domain
    MeshLoader<Mesh2D> domain("c_shaped_631");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define statistical model
    std::vector<double> alphas = {0.01, 0.05, 0.10, 0.90, 0.95, 0.99}; 

    const std::string data_type = "hetero"; 

    // define grid of lambda values
    std::vector<SVector<1>> lambdas;
    for(double x = -6.9; x <= -5.8; x +=0.1) lambdas.push_back(SVector<1>(std::pow(10,x)));
    DVector<double> best_lambda;
    best_lambda.resize(alphas.size());  

    // Read covariates and locations
    DMatrix<double> X = read_csv<double>(R_path + "/data_" + data_type + "/X.csv"); 
    DMatrix<double> loc = read_csv<double>(R_path + "/data_" + data_type + "/locs.csv"); 

    // Simulations 
    const unsigned int M = 10; 
    for(auto m = 1; m <= M; ++m){
        std::cout << "--------------------Simulation #" << std::to_string(m) << "-------------" << std::endl; 
        unsigned int ind = 0; 
        std::ofstream fileGCV(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/single_est/gcv_scores.csv");
        for(auto alpha : alphas){

            std::cout << "------------------alpha=" << std::to_string(alpha) << "-----------------" << std::endl; 

            SQRPDE<decltype(problem), SpaceOnly, GeoStatMeshNodes, MonolithicSolver> model(problem, alpha);
            model.set_spatial_locations(loc);

            // load data from .csv files
            DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/y.csv");

            // set model data
            BlockFrame<double, int> df;
            df.insert(OBSERVATIONS_BLK, y);
            df.insert(DESIGN_MATRIX_BLK, X);
            model.set_data(df);

            model.init(); // init model

            // define GCV function and optimize
            GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
            Grid<1> opt;

            ScalarField<1, decltype(GCV)> obj(GCV);
            opt.optimize(obj, lambdas); // optimize gcv field
            std::cout << "opt: " << opt.optimum()[0] << std::endl; 
            best_lambda[ind] = opt.optimum()[0];
            
            std::cout << "Best lambda is: " << std::setprecision(16) << best_lambda[ind] << std::endl; 
            ind++;

            // gcv scores
            for(std::size_t i = 0; i < GCV.values().size(); ++i){
                fileGCV << std::setprecision(16) << std::sqrt(GCV.values()[i]) << "\n"; 
            }

        }

        const static Eigen::IOFormat CSVFormatL(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        std::ofstream fileL(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/single_est/lambdas_opt.csv");
        if (fileL.is_open()){
            fileL << best_lambda.format(CSVFormatL);
            fileL.close();
        }

        fileGCV.close(); 
    }

}




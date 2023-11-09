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


// test 1
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
TEST(msqrpde_test1, laplacian_nonparametric_samplingatnodes) {

    // path test  
    std::string C_path = "/mnt/c/Users/marco/PACS/Project/Code/Cpp/fdaPDE-fork/test/data/models/msqrpde/2D_test1"; 
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

    // Simulations 
    const unsigned int M = 1; 
    for(auto m = 1; m <= M; ++m){

        MSQRPDE<decltype(problem), SpaceOnly, GeoStatMeshNodes, MonolithicSolver> model(problem, alphas);
        // use optimal lambda to avoid possible numerical issues
        DMatrix<double> lambdas = read_csv<double>(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/single_est/lambdas_opt.csv");   // the vector should be saved with the "R format"
        model.setLambdas_D(lambdas);
        // load data from .csv files
        DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/z.csv"); 

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
        std::ofstream filef(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/f_all.csv");
        if (filef.is_open()){
            filef << computedF.format(CSVFormatf);
            filef.close();
        }

        // // debug 
        // DMatrix<double> computedA = model.A_mult();
        // const static Eigen::IOFormat CSVFormatA(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        // std::ofstream fileA(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/A.csv");
        // if (fileA.is_open()){
        //   fileA << computedA.format(CSVFormatA);
        //   fileA.close();
        // }

        // DMatrix<double> computedW = model.W_mult();
        // const static Eigen::IOFormat CSVFormatW(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        // std::ofstream fileW(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/W.csv");
        // if (fileW.is_open()){
        //   fileW << computedW.format(CSVFormatW);
        //   fileW.close();
        // }

        // DMatrix<double> computedWbar = model.Wbar_mult();
        // const static Eigen::IOFormat CSVFormatWbar(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        // std::ofstream fileWbar(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/Wbar.csv");
        // if (fileWbar.is_open()){
        //   fileWbar << computedWbar.format(CSVFormatWbar);
        //   fileWbar.close();
        // }

        // DMatrix<double> computedDelta = model.Delta_mult();
        // const static Eigen::IOFormat CSVFormatDelta(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        // std::ofstream fileDelta(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/Delta.csv");
        // if (fileDelta.is_open()){
        //   fileDelta << computedDelta.format(CSVFormatDelta);
        //   fileDelta.close();
        // }

        // DMatrix<double> computedD_script = model.D_script();
        // const static Eigen::IOFormat CSVFormatD_script(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        // std::ofstream fileD_script(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/D_script.csv");
        // if (fileD_script.is_open()){
        //   fileD_script << computedD_script.format(CSVFormatD_script);
        //   fileD_script.close();
        // }

        // DMatrix<double> computedDscriptj = model.Dscriptj();
        // const static Eigen::IOFormat CSVFormatDscriptj(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        // std::ofstream fileDscriptj(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/D_script_j.csv");
        // if (fileDscriptj.is_open()){
        //   fileDscriptj << computedDscriptj.format(CSVFormatDscriptj);
        //   fileDscriptj.close();
        // }

        // DMatrix<double> computedPsi = model.Psi_mult();
        // const static Eigen::IOFormat CSVFormatPsi(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        // std::ofstream filePsi(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/Psi.csv");
        // if (filePsi.is_open()){
        //   filePsi << computedPsi.format(CSVFormatPsi);
        //   filePsi.close();
        // }

    }

}



// test 2
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
TEST(msqrpde_test2, laplacian_semiparametric_samplingatnodes) {

    // path test  
    std::string C_path = "/mnt/c/Users/marco/PACS/Project/Code/Cpp/fdaPDE-fork/test/data/models/msqrpde/2D_test2"; 
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

    // Read covariates
    DMatrix<double> X = read_csv<double>(R_path + "/data_" + data_type + "/X.csv"); 

    // Simulations 
    const unsigned int M = 1; 
    for(auto m = 1; m <= M; ++m){

        std::cout << "--------------------Simulation #" << std::to_string(m) << "-------------" << std::endl; 

        MSQRPDE<decltype(problem), SpaceOnly, GeoStatMeshNodes, MonolithicSolver> model(problem, alphas);
        // use optimal lambda to avoid possible numerical issues
        DMatrix<double> lambdas = read_csv<double>(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/single_est/lambdas_opt.csv");   // the vector should be saved with the "R format"
        model.setLambdas_D(lambdas);
        // load data from .csv files
        DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/z.csv"); 

        // set model data
        BlockFrame<double, int> df;
        df.insert(OBSERVATIONS_BLK, y);
        df.insert(DESIGN_MATRIX_BLK, X);
        model.set_data(df);

        // solve smoothing problem
        model.init();
        model.solve();

        // Save solution
        DMatrix<double> computedF = model.f();
        const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        std::ofstream filef(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/f_all.csv");
        if (filef.is_open()){
            filef << computedF.format(CSVFormatf);
            filef.close();
        }

        DMatrix<double> computedBeta = model.beta();
        const static Eigen::IOFormat CSVFormatBeta(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        std::ofstream fileBeta(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/beta_all.csv");
        if (fileBeta.is_open()){
            fileBeta << computedBeta.format(CSVFormatBeta);
            fileBeta.close();
        }
    }
}

// test 3
//    domain:       c-shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
TEST(msqrpde_test3, laplacian_semiparametric_samplingatlocations) {

    // path test  
    std::string C_path = "/mnt/c/Users/marco/PACS/Project/Code/Cpp/fdaPDE-fork/test/data/models/msqrpde/2D_test3"; 
    std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/multiple_quantiles/Tests/Test_3"; 

    // define domain
    MeshLoader<Mesh2D> domain("c_shaped_631");

    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define statistical model
    std::vector<double> alphas = {0.01, 0.05, 0.10, 0.90, 0.95, 0.99}; 

    const std::string data_type = "hetero"; 

    // Read covariates
    DMatrix<double> X = read_csv<double>(R_path + "/data_" + data_type + "/X.csv"); 
    DMatrix<double> loc = read_csv<double>(R_path + "/data_" + data_type + "/locs.csv"); 

    // Simulations 
    const unsigned int M = 1; 
    for(auto m = 1; m <= M; ++m){

        std::cout << "--------------------Simulation #" << std::to_string(m) << "-------------" << std::endl; 

        MSQRPDE<decltype(problem), SpaceOnly, GeoStatLocations, MonolithicSolver> model(problem, alphas);
        model.set_spatial_locations(loc);

        // use optimal lambda to avoid possible numerical issues
        DMatrix<double> lambdas = read_csv<double>(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/single_est/lambdas_opt.csv");   // the vector should be saved with the "R format"
        model.setLambdas_D(lambdas);

        // load data from .csv files
        DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/z.csv");

        // set model data
        BlockFrame<double, int> df;
        df.insert(OBSERVATIONS_BLK, y);
        df.insert(DESIGN_MATRIX_BLK, X);
        model.set_data(df);

        // solve smoothing problem
        model.init();
        model.solve();

        // Save solution
        DMatrix<double> computedF = model.f();
        const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        std::ofstream filef(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/f_all.csv");
        if (filef.is_open()){
            filef << computedF.format(CSVFormatf);
            filef.close();
        }

        DMatrix<double> computedFn = model.Psi_mult()*model.f();
        const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        std::ofstream filefn(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/fn_all.csv");
        if (filefn.is_open()){
            filefn << computedFn.format(CSVFormatfn);
            filefn.close();
        }

        DMatrix<double> computedBeta = model.beta();
        const static Eigen::IOFormat CSVFormatBeta(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        std::ofstream fileBeta(R_path + "/data_" + data_type + "/sim_" + std::to_string(m) + "/mult_est/beta_all.csv");
        if (fileBeta.is_open()){
            fileBeta << computedBeta.format(CSVFormatBeta);
            fileBeta.close();
        }
  }
}



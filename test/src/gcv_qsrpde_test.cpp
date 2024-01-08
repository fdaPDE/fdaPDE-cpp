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
using fdapde::core::laplacian;
using fdapde::core::MatrixDataWrapper;
using fdapde::core::PDE;
using fdapde::core::VectorDataWrapper;

#include "../../fdaPDE/models/sampling_design.h"
#include "../../fdaPDE/models/regression/qsrpde.h"
#include "../../fdaPDE/models/regression/gcv.h"
using fdapde::models::QSRPDE;
using fdapde::models::SpaceOnly;
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

double RMSE(DVector<double> v1, DVector<double> v2){
    double res = 0.; 
    if(v1.size() != v2.size())
        std::cout << std::endl << "----------ERROR IN RMSE COMPUTATION---------" << std::endl; 
    for(auto i = 0; i < v1.size(); ++i){
        res += (v1[i]-v2[i])*(v1[i]-v2[i]); 
    }
    return std::sqrt(1./(v1.size())*res); 
}

// test 1
//    domain:       unit square [1,1] x [1,1] (coarse)
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
TEST(gcv_qsrpde_test, laplacian_nonparametric_samplingatnodes_spaceonly_gridexact) {
    // define domain
    MeshLoader<Mesh2D> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test1/y.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double alpha = 0.1;
    QSRPDE<SpaceOnly> model(problem, Sampling::mesh_nodes, alpha);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    auto GCV = model.gcv<ExactEDF>();
    std::vector<DVector<double>> lambdas;
    for (double x = -8.0; x <= -5.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas);
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/gcv/qsrpde/2D_test1/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/gcv/qsrpde/2D_test1/gcvs.mtx"));
}

// test 2
//    domain:       unit square [1,1] x [1,1] (coarse)
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid stochastic
TEST(gcv_qsrpde_test, laplacian_nonparametric_samplingatnodes_spaceonly_gridstochastic) {
    // define domain
    MeshLoader<Mesh2D> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test2/y.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double alpha = 0.1;
    QSRPDE<SpaceOnly> model(problem, Sampling::mesh_nodes, alpha);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    std::size_t seed = 438172;
    auto GCV = model.gcv<StochasticEDF>(1000, seed);
    std::vector<DVector<double>> lambdas;
    for (double x = -8.0; x <= -5.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas);
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/gcv/qsrpde/2D_test2/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/gcv/qsrpde/2D_test2/gcvs.mtx"));
}

// test 3
//    domain:       c-shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
TEST(gcv_qsrpde_test, laplacian_semiparametric_samplingatlocations_gridexact) {
    // define domain
    MeshLoader<Mesh2D> domain("c_shaped");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/gcv/qsrpde/2D_test3/locs.csv");
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test3/y.csv");
    DMatrix<double> X = read_csv<double>("../data/gcv/qsrpde/2D_test3/X.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double alpha = 0.9;
    QSRPDE<SpaceOnly> model(problem, Sampling::pointwise, alpha);
    model.set_spatial_locations(locs);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    auto GCV = model.gcv<ExactEDF>();
    std::vector<DVector<double>> lambdas;
    for (double x = -5.0; x <= -3.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas);
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/gcv/qsrpde/2D_test3/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/gcv/qsrpde/2D_test3/gcvs.mtx"));
}

// test 4
//    domain:       c-shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    GCV optimization: grid stochastic
TEST(gcv_qsrpde_test, laplacian_semiparametric_samplingatlocations_gridstochastic) {
    // define domain
    MeshLoader<Mesh2D> domain("c_shaped");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/gcv/qsrpde/2D_test4/locs.csv");
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test4/y.csv");
    DMatrix<double> X = read_csv<double>("../data/gcv/qsrpde/2D_test4/X.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double alpha = 0.9;
    QSRPDE<SpaceOnly> model(problem, Sampling::pointwise, alpha);
    model.set_spatial_locations(locs);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D value
    std::size_t seed = 66546513;
    auto GCV = model.gcv<StochasticEDF>(1000, seed);
    std::vector<DVector<double>> lambdas;
    for (double x = -5.0; x <= -3.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas);
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/gcv/qsrpde/2D_test4/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/gcv/qsrpde/2D_test4/gcvs.mtx"));
}

// test 5
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes
//    penalization: costant coefficients PDE
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
TEST(gcv_qsrpde_test, costantcoefficientspde_nonparametric_samplingatnodes_gridexact) {
    // define domain
    MeshLoader<Mesh2D> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test5/y.csv");
    // define regularizing PDE
    SMatrix<2> K;
    K << 1, 0, 0, 4;
    auto L = -diffusion<FEM>(K);   // anisotropic diffusion
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double alpha = 0.1;
    QSRPDE<SpaceOnly> model(problem, Sampling::mesh_nodes, alpha);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    auto GCV = model.gcv<ExactEDF>();
    std::vector<DVector<double>> lambdas;
    for (double x = -7.0; x <= -5.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas);   // optimize gcv field
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/gcv/qsrpde/2D_test5/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/gcv/qsrpde/2D_test5/gcvs.mtx"));
}

// test 6
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes
//    penalization: costant coefficients PDE
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid stochastic
TEST(gcv_qsrpde_test, costantcoefficientspde_nonparametric_samplingatnodes_gridstochastic) {
    // define domain
    MeshLoader<Mesh2D> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test6/y.csv");
    // define regularizing PDE
    SMatrix<2> K;
    K << 1, 0, 0, 4;
    auto L = -diffusion<FEM>(K);   // anisotropic diffusion
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double alpha = 0.1;
    QSRPDE<SpaceOnly> model(problem, Sampling::mesh_nodes, alpha);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    std::size_t seed = 438172;
    auto GCV = model.gcv<StochasticEDF>(1000, seed);
    std::vector<DVector<double>> lambdas;
    for (double x = -7.0; x <= -5.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas);   // optimize gcv field
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/gcv/qsrpde/2D_test6/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/gcv/qsrpde/2D_test6/gcvs.mtx"));
}

// test 7
//    domain:       c-shaped
//    sampling:     areal
//    penalization: simple laplacian
//    covariates:   no
//    BC:           yes
//    order FE:     1
//    GCV optimization: grid exact
TEST(gcv_qsrpde_test, laplacian_semiparametric_samplingareal_gridexact) {
    // define domain and regularizing PDE
    MeshLoader<Mesh2D> domain("c_shaped_areal");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test7/y.csv");
    DMatrix<double> X = read_csv<double>("../data/gcv/qsrpde/2D_test7/X.csv");
    DMatrix<double> subdomains = read_csv<double>("../data/gcv/qsrpde/2D_test7/incidence_matrix.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double alpha = 0.5;
    QSRPDE<SpaceOnly> model(problem, Sampling::areal, alpha);
    model.set_spatial_locations(subdomains);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    auto GCV = model.gcv<ExactEDF>();
    std::vector<DVector<double>> lambdas;
    for (double x = -4.0; x <= -1.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas);   // optimize gcv field
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/gcv/qsrpde/2D_test7/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/gcv/qsrpde/2D_test7/gcvs.mtx"));
}

// test 8
//    domain:       c-shaped
//    sampling:     areal
//    penalization: simple laplacian
//    covariates:   no
//    BC:           yes
//    order FE:     1
//    GCV optimization: grid stochastic
TEST(gcv_qsrpde_test, laplacian_semiparametric_samplingareal_gridstochastic) {
    // define domain and regularizing PDE
    // define domain and regularizing PDE
    MeshLoader<Mesh2D> domain("c_shaped_areal");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/gcv/qsrpde/2D_test8/y.csv");
    DMatrix<double> X = read_csv<double>("../data/gcv/qsrpde/2D_test8/X.csv");
    DMatrix<double> subdomains = read_csv<double>("../data/gcv/qsrpde/2D_test8/incidence_matrix.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double alpha = 0.5;
    QSRPDE<SpaceOnly> model(problem, Sampling::areal, alpha);
    model.set_spatial_locations(subdomains);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    std::size_t seed = 438172;
    auto GCV = model.gcv<StochasticEDF>(100, seed);
    std::vector<DVector<double>> lambdas;
    for (double x = -4.0; x <= -1.0; x += 0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas);
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/gcv/qsrpde/2D_test8/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/gcv/qsrpde/2D_test8/gcvs.mtx"));
}



// test 5  (gcv test)
//    domain:         unit square 
//    sampling: locations != nodes
//    penalization:   simple laplacian
//    covariates:     no
//    BC:             no
//    order FE:       1

TEST(gcv_qsrpde_test, laplacian_nonparametric_samplingatlocations){

    // Marco 
    std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/space_only/Test_5"; 

    std::vector<double> alphas = {0.1, 0.5, 0.9};  
    unsigned int n_sim = 20;
    std::vector<std::string> lambda_selection_types ={"gcv", "gcv_smooth_eps1e-3", "gcv_smooth_eps1e-2", "gcv_smooth_eps1e-1.5", "gcv_smooth_eps1e-1"};     
    bool compute_rmse = false; 
    bool compute_gcv = true; 

    std::vector<std::string> data_types = {"hetero_4"};  // {"hetero_1", "hetero_2", "hetero_3", "hetero_4"};
    for(auto data_type : data_types){
            
        std::vector<DVector<double>> lambdas_10;
        std::vector<DVector<double>> lambdas_50;
        std::vector<DVector<double>> lambdas_90;
        std::vector<double> lambdas_rmse;

        if(data_type == "hetero_1"){
            for(double x = -6.5; x <= -3.0; x += 0.05) lambdas_10.push_back(SVector<1>(std::pow(10, x)));
            for(double x = -6.5; x <= -3.0; x += 0.05) lambdas_50.push_back(SVector<1>(std::pow(10, x))); 
            for(double x = -6.5; x <= -3.0; x += 0.05) lambdas_90.push_back(SVector<1>(std::pow(10, x)));
            for(double x = -7.0; x <= -5.0; x += 0.05) lambdas_rmse.push_back(std::pow(10, x));  
        }

        if(data_type == "hetero_2"){
            for(double x = -6.5; x <= -2.5; x += 0.05) lambdas_10.push_back(SVector<1>(std::pow(10, x)));
            for(double x = -6.0; x <= -2.0; x += 0.05) lambdas_50.push_back(SVector<1>(std::pow(10, x))); 
            for(double x = -6.5; x <= -2.0; x += 0.05) lambdas_90.push_back(SVector<1>(std::pow(10, x))); 
            for(double x = -7.5; x <= -4.0; x += 0.05) lambdas_rmse.push_back(std::pow(10, x));
        }

        if(data_type == "hetero_3"){
            for(double x = -6.0; x <= -1.0; x += 0.05) lambdas_10.push_back(SVector<1>(std::pow(10, x)));
            for(double x = -5.4; x <= -1.0; x += 0.05) lambdas_50.push_back(SVector<1>(std::pow(10, x))); 
            for(double x = -5.5; x <= -0.2; x += 0.05) lambdas_90.push_back(SVector<1>(std::pow(10, x))); 
            for(double x = -9.0; x <= -4.5; x += 0.05) lambdas_rmse.push_back(std::pow(10, x));
        }

        if(data_type == "hetero_4"){
            for(double x = -4.8; x <= 2.5; x += 0.05) lambdas_10.push_back(SVector<1>(std::pow(10, x)));
            for(double x = -0.5; x <= 2.5; x += 0.05) lambdas_50.push_back(SVector<1>(std::pow(10, x))); 
            for(double x = -5.2; x <= 2.5; x += 0.05) lambdas_90.push_back(SVector<1>(std::pow(10, x))); 
            for(double x = -5.7; x <= -2.0; x += 0.05) lambdas_rmse.push_back(std::pow(10, x)); 
        }

        // define spatial domain and regularizing PDE 
        MeshLoader<Mesh2D> domain("unit_square_17"); 

        // import locs from files
        DMatrix<double> locs = read_csv<double>(R_path + "/data_" + data_type + "/locs.csv");

        // define regularizing PDE
        auto L = -laplacian<FEM>();
        DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements()*3, 1);
        PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

        // simulations 
        for(unsigned int sim = 1; sim <= n_sim; ++sim){
            std::cout << std::endl << "---------------Simulation #" << sim << "--------------" << std::endl; 
            // load data from .csv files
            DMatrix<double> y = read_csv<double>(R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/y.csv");

            for(double alpha : alphas){

                unsigned int alpha_int = alpha*100; 
                std::string alpha_string = std::to_string(alpha_int);
                std::cout << "--------alpha=" << alpha_string << "%" << std::endl;  

                BlockFrame<double, int> df;

                if(compute_gcv){

                    // GCV:
                    for(auto lambda_selection_type : lambda_selection_types){
                        QSRPDE<SpaceOnly> model_gcv(problem, Sampling::pointwise, alpha);
                
                        std::vector<DVector<double>> lambdas;
                        if(almost_equal(alpha, 0.1)){
                            lambdas = lambdas_10; 
                        }  
                        if(almost_equal(alpha, 0.5)){
                            lambdas = lambdas_50; 
                        }  
                        if(almost_equal(alpha, 0.9)){
                            lambdas = lambdas_90; 
                        }  

                        // set model's data
                        model_gcv.set_spatial_locations(locs);
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
                        
                        df.stack(OBSERVATIONS_BLK, y);
                        model_gcv.set_data(df);
                        model_gcv.init();

                        // define GCV function and grid of \lambda_D values
                        auto GCV = model_gcv.gcv<ExactEDF>();
                        // optimize GCV
                        Grid<fdapde::Dynamic> opt;
                        opt.optimize(GCV, lambdas);
                        SVector<1> best_lambda = opt.optimum();

                        std::string solutions_path_gcv = R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + "/" + lambda_selection_type; 

                        // Save lambda sequence 
                        std::ofstream fileLambdaS(solutions_path_gcv + "/lambdas_seq.csv");
                        for(std::size_t i = 0; i < lambdas.size(); ++i) 
                            fileLambdaS << std::setprecision(16) << lambdas[i] << "\n"; 
                        fileLambdaS.close();
                        // Save lambda GCVopt 
                        std::ofstream fileLambdaoptS(solutions_path_gcv + "/lambda_s_opt.csv");
                        if(fileLambdaoptS.is_open()){
                            fileLambdaoptS << std::setprecision(16) << best_lambda[0];
                            fileLambdaoptS.close();
                        }
                        // Save lambda GCV 
                        std::ofstream fileGCV_scores(solutions_path_gcv + "/score.csv");
                        for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
                            fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.gcvs()[i]) << "\n"; 
                        fileGCV_scores.close();
                    }

                }

                if(compute_rmse){
                    std::cout << "-----RMSE computation-----" << std::endl; 
                    // RMSE
                    DMatrix<double> f_true = read_csv<double>(R_path + "/data_" + data_type + "/true/f_true_" + alpha_string + ".csv");
                    std::string solutions_path_rmse = R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + "/RMSE"; 
                    std::vector<double> rmse_score; 
                    rmse_score.resize(lambdas_rmse.size()); 
                    double count_l = 0; 
                    for(auto lambda : lambdas_rmse){
                        QSRPDE<SpaceOnly> model(problem, Sampling::pointwise, alpha);
                        // set model's data
                        model.set_spatial_locations(locs);
                        model.set_lambda_D(lambda);            
                        df.stack(OBSERVATIONS_BLK, y);
                        model.set_data(df);
                        model.init();
                        model.solve();

                        rmse_score[count_l] = RMSE(model.f(), f_true); 

                        count_l = count_l+1; 
                    }

                    auto min_idx = std::distance(std::begin(rmse_score), std::min_element(std::begin(rmse_score), std::end(rmse_score))); 
                    
                    // Save lambda sequence 
                    std::ofstream fileLambdaS_rmse(solutions_path_rmse + "/lambdas_seq.csv");
                    for(std::size_t i = 0; i < lambdas_rmse.size(); ++i) 
                        fileLambdaS_rmse << std::setprecision(16) << lambdas_rmse[i] << "\n"; 
                    fileLambdaS_rmse.close();

                    // Save score 
                    std::ofstream fileRMSE_scores(solutions_path_rmse + "/score.csv");
                    for(std::size_t i = 0; i < rmse_score.size(); ++i) 
                        fileRMSE_scores << std::setprecision(16) << rmse_score[i] << "\n"; 
                    fileRMSE_scores.close();

                    // Save lambda RMSEopt 
                    std::ofstream fileLambdaoptS_rmse(solutions_path_rmse + "/lambda_s_opt.csv");
                    if(fileLambdaoptS_rmse.is_open()){
                        fileLambdaoptS_rmse << std::setprecision(16) << lambdas_rmse[min_idx]; ;
                        fileLambdaoptS_rmse.close();
                    }
                
                }


            }
        }

    }
}

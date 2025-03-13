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
using fdapde::core::FEM;
using fdapde::core::fem_order;
using fdapde::core::laplacian;
using fdapde::core::diffusion;
using fdapde::core::PDE;
using fdapde::core::Triangulation;
using fdapde::core::bilaplacian;
using fdapde::core::SPLINE;
using fdapde::core::spline_order;

#include "../../fdaPDE/models/regression/qsrpde.h"
#include "../../fdaPDE/models/sampling_design.h"
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
using fdapde::testing::read_csv;

// test 1
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
TEST(qsrpde_test, laplacian_nonparametric_samplingatnodes) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/models/qsrpde/2D_test1/y.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double lambda = 1.778279 * std::pow(0.1, 4);
    double alpha = 0.1;
    QSRPDE<SpaceOnly> model(problem, Sampling::mesh_nodes, alpha);
    model.set_lambda_D(lambda);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f(), "../data/models/qsrpde/2D_test1/sol.mtx"));
}

// test 2
//    domain:       c-shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
TEST(qsrpde_test, laplacian_semiparametric_samplingatlocations) {
    // define domain and regularizing PDE
    MeshLoader<Triangulation<2, 2>> domain("c_shaped");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/models/qsrpde/2D_test2/locs.csv");
    DMatrix<double> y    = read_csv<double>("../data/models/qsrpde/2D_test2/y.csv");
    DMatrix<double> X    = read_csv<double>("../data/models/qsrpde/2D_test2/X.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define statistical model
    double alpha = 0.9;
    double lambda = 3.162277660168379 * std::pow(0.1, 4);   // use optimal lambda to avoid possible numerical issues
    QSRPDE<SpaceOnly> model(problem, Sampling::pointwise, alpha);
    model.set_lambda_D(lambda);
    model.set_spatial_locations(locs);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f()   , "../data/models/qsrpde/2D_test2/sol.mtx" ));
    EXPECT_TRUE(almost_equal(model.beta(), "../data/models/qsrpde/2D_test2/beta.mtx"));
}

// test 3
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes
//    penalization: costant coefficients PDE
//    covariates:   no
//    BC:           no
//    order FE:     1
TEST(qsrpde_test, costantcoefficientspde_nonparametric_samplingatnodes) {
    // define domain and regularizing PDE
    MeshLoader<Triangulation<2, 2>> domain("unit_square_coarse");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/models/qsrpde/2D_test3/y.csv");
    // define regularizing PDE
    SMatrix<2> K;
    K << 1, 0, 0, 4;
    auto L = -diffusion<FEM>(K);   // anisotropic diffusion
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define statistical model
    double alpha = 0.1;
    double lambda = 5.623413251903491 * pow(0.1, 4);
    QSRPDE<SpaceOnly> model(problem, Sampling::mesh_nodes, alpha);
    model.set_lambda_D(lambda);
    // set model data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f(), "../data/models/qsrpde/2D_test3/sol.mtx"));
}

// test 4
//    domain:       c-shaped
//    sampling:     areal
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
TEST(qsrpde_test, laplacian_semiparametric_samplingareal) {
    // define domain and regularizing PDE
    MeshLoader<Triangulation<2, 2>> domain("c_shaped_areal");
    // import data from files
    DMatrix<double> y = read_csv<double>("../data/models/qsrpde/2D_test4/y.csv");
    DMatrix<double> X = read_csv<double>("../data/models/qsrpde/2D_test4/X.csv");
    DMatrix<double> subdomains = read_csv<double>("../data/models/qsrpde/2D_test4/incidence_matrix.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define statistical model
    double alpha = 0.5;
    double lambda = 5.623413251903491 * std::pow(0.1, 3);   // use optimal lambda to avoid possible numerical issues
    QSRPDE<SpaceOnly> model(problem, Sampling::areal, alpha);
    model.set_lambda_D(lambda);
    model.set_spatial_locations(subdomains);
    // set model data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.f()   , "../data/models/qsrpde/2D_test4/sol.mtx" ));
    EXPECT_TRUE(almost_equal(model.beta(), "../data/models/qsrpde/2D_test4/beta.mtx"));
}

// test 5
//    domain:         c-shaped
//    space sampling: locations != nodes
//    time sampling:  locations != nodes
//    missing data:   no
//    penalization:   simple laplacian
//    covariates:     no
//    BC:             no
//    order FE:       1
//    time penalization: separable (mass penalization)
TEST(qsrpde_test, laplacian_nonparametric_samplingatlocations_separable_monolithic) {
    // define temporal and spatial domain
    Triangulation<1, 1> time_mesh(0, fdapde::testing::pi, 6);   // interval [0, \pi] with 7 knots
    MeshLoader<Triangulation<2, 2>> domain("c_shaped_adj");
    // import data from files
    DMatrix<double> space_locs = read_csv<double>("../data/models/qsrpde/2D_test5/locs.csv");
    DMatrix<double> time_locs  = read_csv<double>("../data/models/qsrpde/2D_test5/time_locations.csv");
    DMatrix<double> y          = read_csv<double>("../data/models/qsrpde/2D_test5/y.csv");
    // define regularizing PDE in space
    auto Ld = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3 * time_mesh.n_nodes(), 1);
    PDE<Triangulation<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);
    // define regularizing PDE in time
    auto Lt = -bilaplacian<SPLINE>();
    PDE<Triangulation<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);
    // define model
    double alpha = 0.5;
    double lambda_D = 1e-3;
    double lambda_T = 1e-6;
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
    // test correctness
    EXPECT_TRUE(almost_equal(model.f(), "../data/models/qsrpde/2D_test5/sol.mtx"));
}



/////

// CONFRONTO METODI CV

// test cv
//    domain:       unit square
//    sampling:     locations != nodes
//    penalization: laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
TEST(sqrpde_test_cv, pde_nonparametric_samplingatlocations_spaceonly_gridexact) {

    std::string test_number = "2";    // "1" "2"

    // path test  
    std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/PhD/Codice/models/QSRPDE/Tests_cv/Test_" + test_number; 

    // define statistical model
    std::vector<double> alphas = {0.50, 0.95}; // {0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99};   

    // define grid of lambda values
    std::vector<std::string> lambda_selection_types = {"gcv_smooth_eps1e-1"}; // {"gcv_smooth_eps1e-3", "gcv_smooth_eps1e-1"};   
    const int eps_power = -1.0;  
    
    // methods 
    std::vector<std::string> score_types = {"gcv", "k-fold", "10-fold"};   // "gcv" "gacv" "gacv_star" "k-fold" "10-fold"
    
    bool compute_rmse = true;
    bool compute_gcv = false;     // if you want to compute either gcv, gacv or gacv*

    bool force_lambdas_longer = false;
    bool lambda_long = false; 

    const unsigned int n_sim = 10;

    // define domain
    std::string domain_str; 
    if(test_number == "1"){
        domain_str = "unit_square_15"; 
    }
    if(test_number == "2"){
        domain_str = "unit_square_25"; 
    }
    MeshLoader<Triangulation<2, 2>> domain(domain_str);

    // rhs 
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);

    // define regularizing PDE
    auto L = -laplacian<FEM>();   
    PDE<Triangulation<2, 2>, decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
 
    double best_lambda; 

    // Read covariates and locations
    DMatrix<double> loc = read_csv<double>(R_path + "/locs.csv"); 

    // Simulations 
    for(auto sim = 1; sim <= n_sim; ++sim){
        std::cout << "--------------------Simulation #" << std::to_string(sim) << "-------------" << std::endl; 
        
        // load data from .csv files
        DMatrix<double> y = read_csv<double>(R_path + "/simulations/sim_" + std::to_string(sim) + "/y.csv");
        BlockFrame<double, int> df;
        df.insert(OBSERVATIONS_BLK, y);

        for(auto alpha : alphas){

            unsigned int alpha_int = alpha*100; 
            std::string alpha_string = std::to_string(alpha_int); 

            std::cout << "------------------alpha=" << alpha_string << "-----------------" << std::endl; 

            std::string solutions_path_rmse = R_path + "/simulations/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + "/RMSE"; 

            if(compute_gcv){
                for(auto lambda_selection_type : lambda_selection_types){
                    for(auto score_type : score_types){

                        std::cout << "------------------score=" << score_type << "-----------------" << std::endl;

                        std::string solutions_path_gcv = R_path + "/simulations/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + "/" + lambda_selection_type + "/" + score_type; 
                            
                        QSRPDE<SpaceOnly> model(problem, Sampling::pointwise, alpha);
                        
                        // Read lambda
                        double lambda; 
                        if(lambda_long){
                            if(score_type == "k-fold"){
                                std::ifstream fileLambda(solutions_path_gcv + "/lambda_s_opt_long.csv");
                                if(fileLambda.is_open()){
                                    fileLambda >> lambda; 
                                    fileLambda.close();
                                }
                            } else{
                                std::ifstream fileLambda(solutions_path_gcv + "/lambda_s_opt.csv");
                                if(fileLambda.is_open()){
                                    fileLambda >> lambda; 
                                    fileLambda.close();
                                }
                            }
                        } else{
                            std::ifstream fileLambda(solutions_path_gcv + "/lambda_s_opt.csv");
                            if(fileLambda.is_open()){
                                fileLambda >> lambda; 
                                fileLambda.close();
                            }
                        }


                        std::cout << "lambda=" << std::setprecision(16) << lambda << std::endl; 


                        model.set_spatial_locations(loc);         
                        model.set_data(df);
                        model.set_lambda_D(lambda);

                        model.init();

                        // solve smoothing problem
                        model.init();; 
                        model.solve();

                        // Save solution
                        DMatrix<double> computedF = model.f();
                        const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
                        std::ofstream filef(solutions_path_gcv + "/f.csv");
                        if(filef.is_open()){
                            filef << computedF.format(CSVFormatf);
                            filef.close();
                        }
                        DMatrix<double> computedFn = model.Psi()*model.f();
                        const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
                        std::ofstream filefn(solutions_path_gcv + "/fn.csv");
                        if(filefn.is_open()){
                            filefn << computedFn.format(CSVFormatfn);
                            filefn.close();
                        }

                    }
                        
                }
            }

            if(compute_rmse){
                std::cout << "------------------score=RMSE selection-----------------" << std::endl;

                std::string solutions_path_gcv = R_path + "/simulations/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + "/RMSE"; 
                    
                QSRPDE<SpaceOnly> model(problem, Sampling::pointwise, alpha);
                
                // Read lambda
                double lambda; 
                std::ifstream fileLambda(solutions_path_gcv + "/lambda_s_opt.csv");
                if(fileLambda.is_open()){
                    fileLambda >> lambda; 
                    fileLambda.close();
                }
            

                std::cout << "lambda=" << std::setprecision(16) << lambda << std::endl; 


                model.set_spatial_locations(loc);         
                model.set_data(df);
                model.set_lambda_D(lambda);

                model.init();

                // solve smoothing problem
                model.init();; 
                model.solve();

                // Save solution
                DMatrix<double> computedF = model.f();
                const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
                std::ofstream filef(solutions_path_gcv + "/f.csv");
                if(filef.is_open()){
                    filef << computedF.format(CSVFormatf);
                    filef.close();
                }
                DMatrix<double> computedFn = model.Psi()*model.f();
                const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
                std::ofstream filefn(solutions_path_gcv + "/fn.csv");
                if(filefn.is_open()){
                    filefn << computedFn.format(CSVFormatfn);
                    filefn.close();
                }
            }

                    

        }


    }
}

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

#include "../../fdaPDE/models/regression/sqrpde.h"
#include "../../fdaPDE/models/sampling_design.h"
using fdapde::models::Areal;
using fdapde::models::GeoStatLocations;
using fdapde::models::GeoStatMeshNodes;
using fdapde::models::MonolithicSolver;
using fdapde::models::SpaceOnly;
using fdapde::models::SQRPDE;

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
// TEST(sqrpde_test, laplacian_nonparametric_samplingatnodes) {
//     // define domain
//     MeshLoader<Mesh2D> domain("unit_square_coarse");
//     // import data from files
//     DMatrix<double> y = read_csv<double>("../data/models/sqrpde/2D_test1/y.csv");
//     // define regularizing PDE
//     auto L = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
//     PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
//     // define model
//     double lambda = 1.778279 * std::pow(0.1, 4);
//     double alpha = 0.1;
//     SQRPDE<decltype(problem), SpaceOnly, GeoStatMeshNodes, MonolithicSolver> model(problem, alpha);
//     model.set_lambda_D(lambda);
//     // set model's data
//     BlockFrame<double, int> df;
//     df.insert(OBSERVATIONS_BLK, y);
//     model.set_data(df);
//     // solve smoothing problem
//     model.init();
//     model.solve();
//     // test correctness
//     EXPECT_TRUE(almost_equal(model.f(), "../data/models/sqrpde/2D_test1/sol.mtx"));
// }

// // test 2
// //    domain:       c-shaped
// //    sampling:     locations != nodes
// //    penalization: simple laplacian
// //    covariates:   yes
// //    BC:           no
// //    order FE:     1
// TEST(sqrpde_test, laplacian_semiparametric_samplingatlocations) {
//     // define domain and regularizing PDE
//     MeshLoader<Mesh2D> domain("c_shaped");
//     // import data from files
//     DMatrix<double> locs = read_csv<double>("../data/models/sqrpde/2D_test2/locs.csv");
//     DMatrix<double> y    = read_csv<double>("../data/models/sqrpde/2D_test2/y.csv");
//     DMatrix<double> X    = read_csv<double>("../data/models/sqrpde/2D_test2/X.csv");
//     // define regularizing PDE
//     auto L = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
//     PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
//     // define statistical model
//     double alpha = 0.9;
//     double lambda = 3.162277660168379 * std::pow(0.1, 4);   // use optimal lambda to avoid possible numerical issues
//     SQRPDE<decltype(problem), SpaceOnly, GeoStatLocations, MonolithicSolver> model(problem, alpha);
//     model.set_lambda_D(lambda);
//     model.set_spatial_locations(locs);
//     // set model's data
//     BlockFrame<double, int> df;
//     df.insert(OBSERVATIONS_BLK, y);
//     df.insert(DESIGN_MATRIX_BLK, X);
//     model.set_data(df);
//     // solve smoothing problem
//     model.init();
//     model.solve();
//     // test correctness
//     EXPECT_TRUE(almost_equal(model.f()   , "../data/models/sqrpde/2D_test2/sol.mtx" ));
//     EXPECT_TRUE(almost_equal(model.beta(), "../data/models/sqrpde/2D_test2/beta.mtx"));
// }

// // test 3
// //    domain:       unit square [1,1] x [1,1]
// //    sampling:     locations = nodes
// //    penalization: costant coefficients PDE
// //    covariates:   no
// //    BC:           no
// //    order FE:     1
// TEST(sqrpde_test, costantcoefficientspde_nonparametric_samplingatnodes) {
//     // define domain and regularizing PDE
//     MeshLoader<Mesh2D> domain("unit_square_coarse");
//     // import data from files
//     DMatrix<double> y = read_csv<double>("../data/models/sqrpde/2D_test3/y.csv");
//     // define regularizing PDE
//     SMatrix<2> K;
//     K << 1, 0, 0, 4;
//     auto L = -diffusion<FEM>(K);   // anisotropic diffusion
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
//     PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
//     // define statistical model
//     double alpha = 0.1;
//     double lambda = 5.623413251903491 * pow(0.1, 4);
//     SQRPDE<decltype(problem), SpaceOnly, GeoStatMeshNodes, MonolithicSolver> model(problem, alpha);
//     model.set_lambda_D(lambda);
//     // set model data
//     BlockFrame<double, int> df;
//     df.insert(OBSERVATIONS_BLK, y);
//     model.set_data(df);
//     // solve smoothing problem
//     model.init();
//     model.solve();
//     // test correctness
//     EXPECT_TRUE(almost_equal(model.f(), "../data/models/sqrpde/2D_test3/sol.mtx"));
// }

// // test 4
// //    domain:       c-shaped
// //    sampling:     areal
// //    penalization: simple laplacian
// //    covariates:   yes
// //    BC:           no
// //    order FE:     1
// TEST(sqrpde_test, laplacian_semiparametric_samplingareal) {
//     // define domain and regularizing PDE
//     MeshLoader<Mesh2D> domain("c_shaped_areal");
//     // import data from files
//     DMatrix<double> y = read_csv<double>("../data/models/sqrpde/2D_test4/y.csv");
//     DMatrix<double> X = read_csv<double>("../data/models/sqrpde/2D_test4/X.csv");
//     DMatrix<int> subdomains = read_csv<int>("../data/models/sqrpde/2D_test4/incidence_matrix.csv");
//     // define regularizing PDE
//     auto L = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
//     PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
//     // define statistical model
//     double alpha = 0.5;
//     double lambda = 5.623413251903491 * std::pow(0.1, 3);   // use optimal lambda to avoid possible numerical issues
//     SQRPDE<decltype(problem), SpaceOnly, Areal, MonolithicSolver> model(problem, alpha);
//     model.set_lambda_D(lambda);
//     model.set_spatial_locations(subdomains);
//     // set model data
//     BlockFrame<double, int> df;
//     df.insert(OBSERVATIONS_BLK, y);
//     df.insert(DESIGN_MATRIX_BLK, X);
//     model.set_data(df);
//     // solve smoothing problem
//     model.init();
//     model.solve();
//     // test correctness
//     EXPECT_TRUE(almost_equal(model.f()   , "../data/models/sqrpde/2D_test4/sol.mtx" ));
//     EXPECT_TRUE(almost_equal(model.beta(), "../data/models/sqrpde/2D_test4/beta.mtx"));
// }




// test 5  
//    domain:         unit square 
//    sampling: locations != nodes
//    penalization:   simple laplacian
//    covariates:     no
//    BC:             no
//    order FE:       1

TEST(sqrpde_test, laplacian_nonparametric_samplingatlocations) {

    // Marco 
    std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/space_only/Test_5"; 

    std::vector<double> alphas = {0.1, 0.5, 0.9}; 
    unsigned int n_sim = 20; 
    std::vector<std::string> lambda_selection_types = {"gcv", "gcv_smooth_eps1e-3", "gcv_smooth_eps1e-2",  "gcv_smooth_eps1e-1.5", "gcv_smooth_eps1e-1"}; 

    std::vector<std::string> data_types = {"hetero_4"};
    for(auto data_type : data_types){
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

                for(auto lambda_selection_type : lambda_selection_types){
                    
                    std::string true_path = R_path + "/data_" + data_type + "/true/f_true_" + alpha_string + ".csv"; 
                    std::string solutions_path = R_path + "/data_" + data_type + "/simulations/sim_" + std::to_string(sim) + "/alpha_" + alpha_string + "/" + lambda_selection_type; 
                    double lambda_s;
                    std::ifstream fileLambdaS(solutions_path + "/lambda_s_opt.csv");
                    if(fileLambdaS.is_open()){
                        fileLambdaS >> lambda_s; 
                        fileLambdaS.close();
                    }

                    SQRPDE<decltype(problem), SpaceOnly, GeoStatLocations, MonolithicSolver> model(problem, alpha);
                    model.set_lambda_D(lambda_s);
                    model.set_spatial_locations(locs);
                    // set model data
                    df.insert(OBSERVATIONS_BLK, y);
                    model.set_data(df);
                    // solve smoothing problem
                    model.init();
                    model.solve();

                    // Save C++ solution 
                    DMatrix<double> computedF = model.f();
                    const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
                    std::ofstream filef(solutions_path + "/f.csv");
                    if(filef.is_open()){
                        filef << computedF.format(CSVFormatf);
                        filef.close();
                    }

                }

            }
        }
    }

}




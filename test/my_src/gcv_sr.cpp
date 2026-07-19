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

#include <fdaPDE/fdapde.h>

#include <iostream>
#include <string>
#include <vector>
#include "my_fun.h"   // include le funzioni per COSP 


using namespace fdapde;
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using sparse_matrix_t = Eigen::SparseMatrix<double>;
using fdapde::GridSearch; 


int test_1(); 


int main(){
    test_1();
    return 0; 
}



// Test 1 (COSP)
//    mesh:         unit square
//    sampling:     locations != nodes + areal 
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid stochastic 


int test_1() {

    const bool verbose = false; 

    const unsigned int sim_start = 1; 
    const unsigned int n_sim = 50; 

    const std::string R_path = "../../OneDrive - Politecnico di Milano/Corsi/PhD/Codice/models/COSP/space-only/Test_1/";

    const bool hetero_fit = true;
    const bool naive_fit = false;    // caso sigma_p = sigma_A (=> W=I)
    // NOTA: il naive fit corrisponde all'inizializzazione del caso etero
    const bool only_point_fit = false; 
    const bool only_area_fit = false; 


    // define lambda sequences 
    std::vector<double> lambdas_hetero;
    std::vector<double> lambdas_naive;
    std::vector<double> lambdas_only_point;
    std::vector<double> lambdas_only_area;


    // // sequenze lasche 
    // for(double xs = -9.0; xs <= -2.0; xs += 0.50)   
    //     lambdas_hetero.push_back(std::pow(10,xs));

    // for(double xs = -9.0; xs <= -2.0; xs += 0.50)   
    //     lambdas_naive.push_back(std::pow(10,xs));

    // for(double xs = -9.0; xs <= -2.0; xs += 0.50)   
    //     lambdas_only_point.push_back(std::pow(10,xs));

    // for(double xs = -9.0; xs <= -2.0; xs += 0.50)   
    //     lambdas_only_area.push_back(std::pow(10,xs));


    // sequenze fini
    for(double xs = -5.0; xs <= -2.0; xs += 0.10)   
        lambdas_hetero.push_back(std::pow(10,xs));

    for(double xs = -5.0; xs <= -2.0; xs += 0.10)   
        lambdas_naive.push_back(std::pow(10,xs));

    for(double xs = -7.0; xs <= -4.0; xs += 0.10)   
        lambdas_only_point.push_back(std::pow(10,xs));

    for(double xs = -7.0; xs <= -4.0; xs += 0.10)   
        lambdas_only_area.push_back(std::pow(10,xs));


    // geometry
    std::string mesh_path = "../my_data/mesh/unit_square_21/";

    Eigen::Matrix<double, Dynamic, Dynamic> points = read_csv<double>("my_data/mesh/unit_square_21/points.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> elements = read_csv<int>("my_data/mesh/unit_square_21/elements.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> boundary = read_csv<int>("my_data/mesh/unit_square_21/boundary.csv").as_matrix();

    elements.array() -= 1; // non necessario
    
    Triangulation<2, 2> D(points, elements, boundary);

    // seeds 
    unsigned int random_seed = 37;
    unsigned int random_seed_p = 17; 
    unsigned int random_seed_A = 23; 
        
    // physics
    std::cout << "Define the physics..." << std::endl;
    FeSpace Vh(D, P1<1>);
    TrialFunction f_trial(Vh);
    TestFunction  v_test(Vh);
    auto a = integral(D)(dot(grad(f_trial), grad(v_test)));
    ZeroField<2> u_fun;
    auto F = integral(D)(u_fun * v_test);

 
    // --- Draft implementation (outside fdaPDE architecture) --- 

    // stiffness and mass matrices and forcing vector
    sparse_matrix_t R1 = a.assemble();
    sparse_matrix_t R0 = integral(D)(f_trial * v_test).assemble();
    vector_t u = F.assemble();
    unsigned int n_dofs = Vh.n_dofs();

    for(auto sim = sim_start; sim <= n_sim; ++sim){

        std::cout << std::endl; 
        std::cout << "SIMULATION " << sim << " / " << n_sim << std::endl; 

        std::string solution_path = R_path + "results/sim_" + std::to_string(sim) + "/hetero/";
        std::string solution_path_naive = R_path + "results/sim_" + std::to_string(sim) + "/naive/";
        std::string solution_path_only_point = R_path + "results/sim_" + std::to_string(sim) + "/only_point/";
        std::string solution_path_only_area = R_path + "results/sim_" + std::to_string(sim) + "/only_area/";

        // data
        std::cout << "Define the data structure..." << std::endl;
        GeoFrame geo_data(D);
        // pointwise layer
        auto& l_point = geo_data.insert_scalar_layer<POINT>("l1", R_path + "locs.csv");    
        // areal layer
        auto& l_areal = geo_data.insert_scalar_layer<POLYGON>("l2", R_path + "incidence_mat.csv");

        // // print layers 
        // std::cout << "Printing layers:" << std::endl;
        // std::cout << l_point << std::endl;
        // std::cout << l_areal << std::endl;

        l_point.load_csv<double>(R_path + "results/sim_" + std::to_string(sim) + "/response_point.csv");
        l_areal.load_csv<double>(R_path + "results/sim_" + std::to_string(sim) + "/response_areal.csv");

        // sample sizes 
        unsigned int n_points = l_point.rows();
        std::cout << "Number of pointwise locations: " << n_points << std::endl;
        unsigned int n_areal = l_areal.rows();
        std::cout << "Number of areal locations: " << n_areal << std::endl;

        // total sample size 
        const unsigned int n = n_points + n_areal;

        // extract data
        vector_t y_point = l_point.col<double>(0).as_matrix();
        vector_t y_areal = l_areal.col<double>(0).as_matrix();
        vector_t y(n_points + n_areal);
        y << y_point, y_areal;

        // Psi matrices and D_areal (containes the areal measures)
        sparse_matrix_t Psi_point = internals::point_basis_eval(Vh, l_point.geometry<0>());
        auto [Psi_areal, D_areal] = internals::areal_basis_eval(Vh, l_areal.geometry<0>());
        
        // construct Psi matrix by stacking point over areal Psi  
        std::vector<Eigen::Triplet<double>> triplets_Psi;
        for (int k = 0; k < Psi_point.outerSize(); ++k) {
            for (sparse_matrix_t::InnerIterator it(Psi_point, k); it; ++it) {
                triplets_Psi.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
            }
        }
        for (int k = 0; k < Psi_areal.outerSize(); ++k) {
            for (sparse_matrix_t::InnerIterator it(Psi_areal, k); it; ++it) {
                triplets_Psi.push_back(Eigen::Triplet<double>(it.row() + n_points, it.col(), it.value()));
            }
        }
        sparse_matrix_t Psi(n_points + n_areal, n_dofs);
        Psi.setFromTriplets(triplets_Psi.begin(), triplets_Psi.end());

        // construct diagonal matrix D by stacking point over areal D
        std::vector<Eigen::Triplet<double>> triplets_D;
        for (int i = 0; i < n_points; ++i) {
            triplets_D.push_back(Eigen::Triplet<double>(i, i, 1.0));
        }
        for (int i = 0; i < n_areal; ++i) {
            triplets_D.push_back(Eigen::Triplet<double>(i + n_points, i + n_points, D_areal[i]));
        }
        sparse_matrix_t D_matrix(n_points + n_areal, n_points + n_areal);
        D_matrix.setFromTriplets(triplets_D.begin(), triplets_D.end());


        // run GCV hetero
        if(hetero_fit){

            std::cout << "Running hetero fit..." << std::endl;

            // loop over lambdas: store GCV values 
            vector_t gcv_values(lambdas_hetero.size());
            for (int i = 0; i < lambdas_hetero.size(); ++i) {

                gcv_values[i] = gcv_score_fun(true, n_points, n_areal, n_dofs, D_matrix, Psi, R1, R0, y, u, lambdas_hetero[i], 
                                            random_seed_p, random_seed_A, random_seed, verbose);

                if(verbose){
                    std::cout << "Lambda: " << lambdas_hetero[i] << ", GCV score: " << gcv_values[i] << std::endl;
                }
                    
            }
            Eigen::Index min_index;
            gcv_values.minCoeff(&min_index);
            double lambda_optimum = lambdas_hetero[min_index];
            
            // save
            write_csv(solution_path + "lambdas_seq.csv", lambdas_hetero);
            std::ofstream file_lambda_opt(solution_path + "lambda_opt.csv");
            if(file_lambda_opt.is_open()){
                file_lambda_opt << lambda_optimum << "\n"; 
                file_lambda_opt.close();
            }
            write_csv(solution_path + "gcv_scores.csv", gcv_values);
        }

        // run GCV naive
        if(naive_fit){

            std::cout << "Running naive fit..." << std::endl;

            // loop over lambdas: store GCV values 
            vector_t gcv_values(lambdas_naive.size());
            for (int i = 0; i < lambdas_naive.size(); ++i) {

                gcv_values[i] = gcv_score_fun(false, n_points, n_areal, n_dofs, D_matrix, Psi, R1, R0, y, u, lambdas_naive[i], 
                                            random_seed_p, random_seed_A, random_seed, verbose);

                if(verbose){
                    std::cout << "Lambda: " << lambdas_naive[i] << ", GCV score: " << gcv_values[i] << std::endl;
                }
            }
            Eigen::Index min_index;
            gcv_values.minCoeff(&min_index);
            double lambda_optimum = lambdas_naive[min_index];

            // save
            write_csv(solution_path_naive + "lambdas_seq.csv", lambdas_naive);
            std::ofstream file_lambda_opt(solution_path_naive + "lambda_opt.csv");
            if(file_lambda_opt.is_open()){
                file_lambda_opt << lambda_optimum << "\n"; 
                file_lambda_opt.close();
            }
            write_csv(solution_path_naive + "gcv_scores.csv", gcv_values);
        }

        // run GCV only-point
        if(only_point_fit){

            std::cout << "Running only-point fit..." << std::endl;

            // data
            GeoFrame geo_data_only_point(D);
            // pointwise layer
            auto& l_only_point = geo_data_only_point.insert_scalar_layer<POINT>("l1", R_path + "locs.csv");
            l_only_point.load_csv<double>(R_path + "results/sim_" + std::to_string(sim) + "/response_point.csv");


            GridSearch<1> optimizer;
            SRPDE model_only_point("y ~ f", geo_data_only_point, fe_ls_elliptic(a, F));
            optimizer.optimize(model_only_point.gcv(100, random_seed_p), lambdas_only_point);

            Eigen::Matrix<double, Dynamic, 1> best_lambda = optimizer.optimum();
            std::cout << "Best lambdas is: " << std::setprecision(16) << best_lambda << std::endl; 

            // save 
            write_csv(solution_path_only_point + "lambdas_seq.csv", lambdas_only_point);
            std::ofstream fileLambdaoptS(solution_path_only_point + "lambda_opt.csv");
            if(fileLambdaoptS.is_open()){
                fileLambdaoptS << std::setprecision(16) << best_lambda(0,0);
                fileLambdaoptS.close();
            }
            write_csv(solution_path_only_point + "gcv_scores.csv", optimizer.values());
        }

        // run GCV only-area
        if(only_area_fit){

            std::cout << "Running only-area fit..." << std::endl;

            // data
            GeoFrame geo_data_only_area(D);
            // areal layer
            auto& l_only_area = geo_data_only_area.insert_scalar_layer<POLYGON>("l2", R_path + "incidence_mat.csv");
            l_only_area.load_csv<double>(R_path + "results/sim_" + std::to_string(sim) + "/response_areal.csv");

            GridSearch<1> optimizer;
            SRPDE model_only_area("y ~ f", geo_data_only_area, fe_ls_elliptic(a, F));
            optimizer.optimize(model_only_area.gcv(100, random_seed_A), lambdas_only_area);

            Eigen::Matrix<double, Dynamic, 1> best_lambda = optimizer.optimum();
            std::cout << "Best lambdas is: " << std::setprecision(16) << best_lambda << std::endl; 

            // save 
            write_csv(solution_path_only_area + "lambdas_seq.csv", lambdas_only_area);
            std::ofstream fileLambdaoptS(solution_path_only_area + "lambda_opt.csv");
            if(fileLambdaoptS.is_open()){
                fileLambdaoptS << std::setprecision(16) << best_lambda(0,0);
                fileLambdaoptS.close();
            }
            write_csv(solution_path_only_area + "gcv_scores.csv", optimizer.values());

        }



    }
    


    return 0;
}

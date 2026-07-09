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
#include "my_fun.h"  // include le funzioni per COSP 

using namespace fdapde;
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using sparse_matrix_t = Eigen::SparseMatrix<double>;


int test_20(); 


int main(){
    test_20();
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

int test_20() {

    const unsigned int sim_start = 1; 
    const unsigned int n_sim = 50; 

    const std::string R_path = "../../OneDrive - Politecnico di Milano/Corsi/PhD/Codice/models/COSP/space-only/Test_1/";

    const bool naive_fit = true;    // caso sigma_p = sigma_A (=> W=I)
    const bool hetero_fit = true;
    // NOTA: il naive fit corrisponde all'inizializzazione del caso etero
    const bool only_point_fit = true; 
    const bool only_area_fit = true; 
    
    // geometry
    std::string mesh_path = "../my_data/mesh/unit_square_21/";
    // Triangulation<2, 2> D(mesh_path + "points.csv", mesh_path + "elements.csv", mesh_path + "boundary.csv", true, true);

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

   
    // modeling
    // std::cout << "Fitting model..." << std::endl;
    //SRPDE m("y ~ f", data, fe_ls_elliptic(a, F));
    //m.fit(0.0001428571428571429);
    //std::cout << "Model fitted." << std::endl;


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
        for(int k = 0; k < Psi_point.outerSize(); ++k) {
            for(sparse_matrix_t::InnerIterator it(Psi_point, k); it; ++it) {
                triplets_Psi.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
            }
        }
        for(int k = 0; k < Psi_areal.outerSize(); ++k) {
            for(sparse_matrix_t::InnerIterator it(Psi_areal, k); it; ++it) {
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


        // run hetero
        if(hetero_fit){

            // read optimal lambda 
            double best_lambda; 
            std::ifstream fileLambda(solution_path + "lambda_opt.csv");
            if(fileLambda.is_open()){
                fileLambda >> best_lambda; 
                fileLambda.close();
            }

            auto results_hetero = run_hetero_srpde(n_points, n_areal, n_dofs, D_matrix, Psi, R1, R0, y, u, best_lambda, 
                random_seed_p, random_seed_A, random_seed);

            // save results at convergence 
            write_csv(solution_path + "f.csv", results_hetero.f);
            write_csv(solution_path + "fn.csv", results_hetero.fitted);
            write_csv(solution_path + "g.csv", results_hetero.g);

            std::ofstream file_sigma_sq_p(solution_path + "sigma_sq_p.csv");
            if(file_sigma_sq_p.is_open()){
                file_sigma_sq_p << results_hetero.sigma_sq_p << "\n"; 
                file_sigma_sq_p.close();
            }

            std::ofstream file_sigma_sq_a(solution_path + "sigma_sq_a.csv");
            if(file_sigma_sq_a.is_open()){
                file_sigma_sq_a << results_hetero.sigma_sq_A << "\n"; 
                file_sigma_sq_a.close();    
            }

            std::ofstream file_niter(solution_path + "n_iter.csv");
            if(file_niter.is_open()){
                file_niter << results_hetero.n_iter << "\n"; 
                file_niter.close();
            }

            write_csv(solution_path + "obj_history.csv", results_hetero.obj_history);
            write_csv(solution_path + "loss_p_history.csv", results_hetero.loss_p_history);
            write_csv(solution_path + "loss_a_history.csv", results_hetero.loss_A_history);
            write_csv(solution_path + "loss_global_history.csv", results_hetero.loss_global_history);
            write_csv(solution_path + "sigma_sq_p_history.csv", results_hetero.sigma_sq_p_history);
            write_csv(solution_path + "sigma_sq_a_history.csv", results_hetero.sigma_sq_A_history);
        }

        // run naive
        if(naive_fit){

            // read optimal lambda 
            double best_lambda; 
            std::ifstream fileLambda(solution_path_naive + "lambda_opt.csv");
            if(fileLambda.is_open()){
                fileLambda >> best_lambda; 
                fileLambda.close();
            }


            auto results_naive = run_naive_srpde(n_points, n_areal, n_dofs, D_matrix, Psi, R1, R0, y, u, best_lambda, random_seed);

            // save results 
            write_csv(solution_path_naive + "f.csv", results_naive.f);
            write_csv(solution_path_naive + "fn.csv", results_naive.fitted);  // nonparametric case 
            write_csv(solution_path_naive + "g.csv", results_naive.g);

            std::ofstream file_sigma_sq(solution_path_naive + "sigma_sq.csv");
            if(file_sigma_sq.is_open()){
                file_sigma_sq << results_naive.sigma_sq << "\n"; 
                file_sigma_sq.close();
            }


        }

        // run GCV only-point
        if(only_point_fit){
            std::cout << "Running only-point fit..." << std::endl;

            // data
            std::cout << "Define the only-point data structure..." << std::endl;
            GeoFrame geo_data_only_point(D);
            // pointwise layer
            auto& l_only_point = geo_data_only_point.insert_scalar_layer<POINT>("l1", R_path + "locs.csv");
            l_only_point.load_csv<double>(R_path + "results/sim_" + std::to_string(sim) + "/response_point.csv");

            // std::cout << "Printing layers:" << std::endl;
            // std::cout << l_only_point << std::endl;

            // read optimal lambda 
            double best_lambda; 
            std::ifstream fileLambda(solution_path_only_point + "lambda_opt.csv");
            if(fileLambda.is_open()){
                fileLambda >> best_lambda; 
                fileLambda.close();
            }
    
            std::cout << "Defining model only-point" << std::endl;
            SRPDE model_only_point("y ~ f", geo_data_only_point, fe_ls_elliptic(a, F));
            model_only_point.fit(best_lambda);

            write_csv(solution_path_only_point + "f.csv", model_only_point.f());
        
            // compute sigma_sq and save 
            vector_t y_point = l_only_point.col<double>(0).as_matrix();
            double computedsigmahat = (y_point - model_only_point.fitted()).squaredNorm() / (n_points - model_only_point.edf());
            std::ofstream filesigmahat(solution_path_only_point + "sigma_sq.csv");
            if(filesigmahat.is_open()){
                filesigmahat << std::setprecision(16) << computedsigmahat << "\n"; 
                filesigmahat.close();
            }

        
        }

        // run GCV only-area
        if(only_area_fit){
            std::cout << "Running only-area fit..." << std::endl;

            // data
            std::cout << "Define the only-area data structure..." << std::endl;
            GeoFrame geo_data_only_area(D);
            // areal layer
            auto& l_only_area = geo_data_only_area.insert_scalar_layer<POLYGON>("l2", R_path + "incidence_mat.csv");
            l_only_area.load_csv<double>(R_path + "results/sim_" + std::to_string(sim) + "/response_areal.csv");

            // read optimal lambda 
            double best_lambda; 
            std::ifstream fileLambda(solution_path_only_area + "lambda_opt.csv");
            if(fileLambda.is_open()){
                fileLambda >> best_lambda; 
                fileLambda.close();
            }

            std::cout << "Defining model only-area" << std::endl;
            SRPDE model_only_area("y ~ f", geo_data_only_area, fe_ls_elliptic(a, F));
            model_only_area.fit(best_lambda);

            write_csv(solution_path_only_area + "f.csv", model_only_area.f());
        
            // compute sigma_sq and save 
            vector_t y_area = l_only_area.col<double>(0).as_matrix();
            double computedsigmahat = (y_area - model_only_area.fitted()).squaredNorm() / (n_areal - model_only_area.edf());
            std::ofstream filesigmahat(solution_path_only_area + "sigma_sq.csv");
            if(filesigmahat.is_open()){
                filesigmahat << std::setprecision(16) << computedsigmahat << "\n"; 
                filesigmahat.close();
            }


        }


    }


    return 0;
}

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

using namespace fdapde;
using namespace std::chrono;  // to measure computational times 

int test_06(); 

int main(){
    test_06();
    return 0; 
}




// test 6
//    mesh:         unit square
//    sampling:     locations != nodes
//    penalization: anisotropic diffusion
//    time penalization: separable
//    covariates:   yes
//    BC:           no
//    order FE:     1
int test_06() {

    const unsigned int sim_start = 1; 
    const unsigned int n_sim = 30; 

    // run SRPDE and/or MSRPDE ? 
    const bool run_srpde = true;     // stprde
    const bool run_msrpde = true;     // mixed-effects anisotropic
    const bool run_msr_iso = true;   // mixed-effects isotropic
    const bool run_srpde_d = true;   // strpde con dummies

    bool likelihood_dataloss_type; // false = fpirls data loss, true = likelihood
    bool sigma_edf_type;           // false = sigma senza edf nelle iterazioni, true = sigma con edf

    const std::string trial_number = "14"; 
    const std::string R_path = "../../../OneDrive - Politecnico di Milano/Corsi/PhD/Codice/models/MSRPDE/Tests/space-time/Test_6/trial_" + trial_number + "/";

    if(trial_number == "13" || trial_number == "14"){
        likelihood_dataloss_type = false;
    } else{
        likelihood_dataloss_type = true;
    }

    if(trial_number == "14"){
        sigma_edf_type = true;
    } else{
        sigma_edf_type = false;
    }

    unsigned int M; 
    if(trial_number == "12" || trial_number == "13" || trial_number == "14"){
        M = 8; 
    }
    
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 1, M);  // ATT qui non bisogna fare più M-1 come nella vecchia lib!! Vuole direttamente il numero di nodi, cioè M!!! 
    
    // geometry 
    using PointT = Eigen::Matrix<double, 2, 1>;

    
    std::string N_string; 
    if(trial_number == "12" || trial_number == "13" || trial_number == "14"){
        N_string = "476"; 
    }
    Eigen::Matrix<double, Dynamic, Dynamic> points = read_csv<double>("my_data/mesh/unit_square_reduced_censoring_" + N_string + "/points.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> elements = read_csv<int>("my_data/mesh/unit_square_reduced_censoring_" + N_string + "/elements.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> boundary = read_csv<int>("my_data/mesh/unit_square_reduced_censoring_" + N_string + "/boundary.csv").as_matrix();

    elements.array() -= 1; // non necessario
    
    Triangulation<2, 2> D(points, elements, boundary);

    const unsigned int max_fpirls_iter = 15;

    // time penalty 
    BsSpace Bh(T, 3);   // cubic B-splines in time
    TrialFunction g(Bh);
    TestFunction  w(Bh);
    auto a_T = integral(T)(dxx(g) * dxx(w));
    ZeroField<1> u_T;
    auto F_T = integral(T)(u_T * w);    

    std::vector<double> lambdas_d; std::vector<double> lambdas_t; 
    std::vector<Eigen::Matrix<double, Dynamic, 1>> lambdas_d_t;
    if(trial_number == "12" || trial_number == "13" || trial_number == "14"){
        for(double xs = -4.0-3.0; xs <= -2.0-3.0; xs += 0.25)   // traslato di 3 ordini (n*m = 1100) rispetto alla lib vecchia  -> inoltre, accorciata sequenza
        lambdas_d.push_back(std::pow(10,xs));

        for(double xt = -4.0-3.0; xt <= -4.0-3.0; xt += 2.0)    // traslato di 3 ordini (n*m = 1100) rispetto alla lib vecchia
            lambdas_t.push_back(std::pow(10,xt));
    } 

    for(auto i = 0; i < lambdas_d.size(); ++i)
        for(auto j = 0; j < lambdas_t.size(); ++j) 
            lambdas_d_t.push_back(Eigen::Matrix<double, 2, 1>(lambdas_d[i], lambdas_t[j]));

    Eigen::Matrix<double, Dynamic, 2> lambdas_mat(lambdas_d.size()*lambdas_t.size(), 2);
    for(int i = 0; i < lambdas_d.size(); ++i) { 
        for (int j = 0; j < lambdas_t.size(); ++j) {
            lambdas_mat(i * lambdas_t.size() + j, 0) = lambdas_d[i];
            lambdas_mat(i * lambdas_t.size() + j, 1) = lambdas_t[j];
        }
    }

    // Simulations MSRPDE  
    if(run_msrpde){

        for(auto sim = sim_start; sim <= n_sim; ++sim){


            std::cout << "--------------------Simulation GCV MSRPDE #" << std::to_string(sim) << "-------------" << std::endl; 

            // data 
            GeoFrame data_msrpde(D, T);
            auto& l_msrpde = data_msrpde.insert_scalar_layer<POINT, POINT>("layer", std::pair{R_path + "space_locs.csv", R_path + "time_locs.csv"});
            l_msrpde.load_csv<double>(R_path + "X.csv");
            l_msrpde.load_csv<double>(R_path + "ids_groups.csv"); 


            // load data from .csv files
            l_msrpde.load_csv<double>(R_path + "simulations/sim_" + std::to_string(sim) + "/y_cpp.csv");
                    
            std::string solutions_path_gcv = R_path + "simulations/sim_" + std::to_string(sim) + "/fit_newlib/"; 
            std::string solution_path = R_path + "simulations/sim_" + std::to_string(sim) + "/fit_newlib/"; 
    
            // physics 
            FeSpace Vh(D, P1<1>);   // functional space definition
            Eigen::Matrix<double, 2, 2> K = read_csv<double>(R_path + "simulations/sim_" + std::to_string(sim) + "/K.csv").as_matrix(); 
            TrialFunction f(Vh);
            TestFunction v(Vh);
            auto a_D = integral(D)(dot(K * grad(f), grad(v)));
            // homogeneous forcing linear form
            ZeroField<2> u_D;
            auto F_D = integral(D)(u_D * v);
    
            // Start measuring time
            auto start_time_gcv = high_resolution_clock::now();
    
            // modeling
            MSRPDE m("y ~ x1 + x2 + 1|g + f", data_msrpde, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

            m.set_fpirls_max_iter(max_fpirls_iter);
            m.set_likelihood_dataloss_type(likelihood_dataloss_type);
            m.set_compute_sigma_with_edf(sigma_edf_type); 

            // calibration
            GridSearch<2> opt;   // dimension 2 for space-time problems 
            opt.optimize(m.gcv(100, 1234), lambdas_mat);  // stochastic GCV 

            // Stop measuring time
            auto stop_time_gcv = high_resolution_clock::now();
            auto duration_gcv = duration_cast<milliseconds>(stop_time_gcv - start_time_gcv).count();
            std::cout << "Execution time GCV: " << duration_gcv << " ms" << std::endl;

            Eigen::Matrix<double, Dynamic, 1> best_lambda = opt.optimum();
            std::cout << "Best lambdas are: " << std::setprecision(16) << best_lambda << std::endl; 
    
            // Save lambda sequence 
            write_csv(solutions_path_gcv + "lambdas_seq_S.csv", lambdas_d);
            write_csv(solutions_path_gcv + "lambdas_seq_T.csv", lambdas_t);

            std::ofstream fileLambdaoptS(solutions_path_gcv + "lambda_s_opt.csv");
            if(fileLambdaoptS.is_open()){
              fileLambdaoptS << std::setprecision(16) << best_lambda(0,0);
              fileLambdaoptS.close();
            }
            std::ofstream fileLambdaoptT(solutions_path_gcv + "lambda_t_opt.csv");
            if(fileLambdaoptT.is_open()){
              fileLambdaoptT << std::setprecision(16) << best_lambda(1,0);
              fileLambdaoptT.close();
            }
    
            write_csv(solutions_path_gcv + "score.csv", opt.values());

            std::ofstream file_time_gcv(solutions_path_gcv + "time_gcv.csv"); 
            if(file_time_gcv.is_open()){
                file_time_gcv << duration_gcv << "\n";
                file_time_gcv.close();
            }


            std::cout << "End GCV MSRPDE" << std::endl; 
            
    
        }
    

    }

    // Simulations STRPDE  
    if(run_srpde){

        for(auto sim = sim_start; sim <= n_sim; ++sim){

            std::cout << "--------------------Simulation GCV SRPDE #" << std::to_string(sim) << "-------------" << std::endl; 

            // data 
            GeoFrame data_srpde(D, T);
            auto& l_srpde = data_srpde.insert_scalar_layer<POINT, POINT>("layer", std::pair{R_path + "space_locs.csv", R_path + "time_locs.csv"});
            l_srpde.load_csv<double>(R_path + "X.csv");

            // load data from .csv files
            l_srpde.load_csv<double>(R_path + "simulations/sim_" + std::to_string(sim) + "/y_cpp.csv");
                    
            std::string solutions_path_gcv = R_path + "simulations/sim_" + std::to_string(sim) + "/fit_newlib_srpde/"; 
            std::string solution_path = R_path + "simulations/sim_" + std::to_string(sim) + "/fit_newlib_srpde/"; 
    
            // physics 
            FeSpace Vh(D, P1<1>);   // functional space definition
            Eigen::Matrix<double, 2, 2> K = read_csv<double>(R_path + "simulations/sim_" + std::to_string(sim) + "/K.csv").as_matrix(); 
            TrialFunction f(Vh);
            TestFunction v(Vh);
            auto a_D = integral(D)(dot(K * grad(f), grad(v)));
            // homogeneous forcing linear form
            ZeroField<2> u_D;
            auto F_D = integral(D)(u_D * v);
    
            // Start measuring time
            auto start_time_gcv = high_resolution_clock::now();
    
            // modeling
            SRPDE m("y ~ x1 + x2 + f", data_srpde, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

            // calibration
            GridSearch<2> opt;   // dimension 2 for space-time problems 
            opt.optimize(m.gcv(100, 1234), lambdas_mat);  // stochastic GCV 

            // Stop measuring time
            auto stop_time_gcv = high_resolution_clock::now();
            auto duration_gcv = duration_cast<milliseconds>(stop_time_gcv - start_time_gcv).count();
            std::cout << "Execution time GCV: " << duration_gcv << " ms" << std::endl;

            Eigen::Matrix<double, Dynamic, 1> best_lambda = opt.optimum();
            std::cout << "Best lambdas are: " << std::setprecision(16) << best_lambda << std::endl; 
    
            // Save lambda sequence 
            write_csv(solutions_path_gcv + "lambdas_seq_S.csv", lambdas_d);
            write_csv(solutions_path_gcv + "lambdas_seq_T.csv", lambdas_t);

            std::ofstream fileLambdaoptS(solutions_path_gcv + "lambda_s_opt.csv");
            if(fileLambdaoptS.is_open()){
              fileLambdaoptS << std::setprecision(16) << best_lambda(0,0);
              fileLambdaoptS.close();
            }
            std::ofstream fileLambdaoptT(solutions_path_gcv + "lambda_t_opt.csv");
            if(fileLambdaoptT.is_open()){
              fileLambdaoptT << std::setprecision(16) << best_lambda(1,0);
              fileLambdaoptT.close();
            }
    
            write_csv(solutions_path_gcv + "score.csv", opt.values());

            std::ofstream file_time_gcv(solutions_path_gcv + "time_gcv.csv"); 
            if(file_time_gcv.is_open()){
                file_time_gcv << duration_gcv << "\n";
                file_time_gcv.close();
            }
            
    
        }
    

    }

    // Simulations MSR-ISO  
    if(run_msr_iso){

        for(auto sim = sim_start; sim <= n_sim; ++sim){

            std::cout << "--------------------Simulation GCV MSR-ISO #" << std::to_string(sim) << "-------------" << std::endl; 
   
            // data 
            GeoFrame data_msr_iso(D, T);
            auto& l_msr_iso = data_msr_iso.insert_scalar_layer<POINT, POINT>("layer", std::pair{R_path + "space_locs.csv", R_path + "time_locs.csv"});
            l_msr_iso.load_csv<double>(R_path + "X.csv");
            l_msr_iso.load_csv<double>(R_path + "ids_groups.csv");              

            // load data from .csv files
            l_msr_iso.load_csv<double>(R_path + "simulations/sim_" + std::to_string(sim) + "/y_cpp.csv");
                    
            std::string solutions_path_gcv = R_path + "simulations/sim_" + std::to_string(sim) + "/fit_newlib_iso/"; 
            std::string solution_path = R_path + "simulations/sim_" + std::to_string(sim) + "/fit_newlib_iso/"; 
    
            // physics 
            FeSpace Vh(D, P1<1>);   // functional space definition
            Eigen::Matrix<double, 2, 2> K; 
            K << 1, 0, 0, 1; 
            TrialFunction f(Vh);
            TestFunction v(Vh);
            auto a_D = integral(D)(dot(K * grad(f), grad(v)));
            // homogeneous forcing linear form
            ZeroField<2> u_D;
            auto F_D = integral(D)(u_D * v);
    
            // Start measuring time
            auto start_time_gcv = high_resolution_clock::now();
    
            // modeling
            MSRPDE m("y ~ x1 + x2 + 1|g + f", data_msr_iso, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

            m.set_fpirls_max_iter(max_fpirls_iter);
            m.set_likelihood_dataloss_type(likelihood_dataloss_type);
            m.set_compute_sigma_with_edf(sigma_edf_type);

            // calibration
            GridSearch<2> opt;   // dimension 2 for space-time problems 
            opt.optimize(m.gcv(100, 1234), lambdas_mat);  // stochastic GCV 

            // Stop measuring time
            auto stop_time_gcv = high_resolution_clock::now();
            auto duration_gcv = duration_cast<milliseconds>(stop_time_gcv - start_time_gcv).count();
            std::cout << "Execution time GCV: " << duration_gcv << " ms" << std::endl;

            Eigen::Matrix<double, Dynamic, 1> best_lambda = opt.optimum();
            std::cout << "Best lambdas are: " << std::setprecision(16) << best_lambda << std::endl; 
    
            // Save lambda sequence 
            write_csv(solutions_path_gcv + "lambdas_seq_S.csv", lambdas_d);
            write_csv(solutions_path_gcv + "lambdas_seq_T.csv", lambdas_t);

            std::ofstream fileLambdaoptS(solutions_path_gcv + "lambda_s_opt.csv");
            if(fileLambdaoptS.is_open()){
              fileLambdaoptS << std::setprecision(16) << best_lambda(0,0);
              fileLambdaoptS.close();
            }
            std::ofstream fileLambdaoptT(solutions_path_gcv + "lambda_t_opt.csv");
            if(fileLambdaoptT.is_open()){
              fileLambdaoptT << std::setprecision(16) << best_lambda(1,0);
              fileLambdaoptT.close();
            }
    
            write_csv(solutions_path_gcv + "score.csv", opt.values());

            std::ofstream file_time_gcv(solutions_path_gcv + "time_gcv.csv"); 
            if(file_time_gcv.is_open()){
                file_time_gcv << duration_gcv << "\n";
                file_time_gcv.close();
            }
            
    
        }
    

    }
   
    // Simulations STRPDE DUMMIES 
    if(run_srpde_d){

        for(auto sim = sim_start; sim <= n_sim; ++sim){

            std::cout << "--------------------Simulation GCV SRPDE-DUMMIES #" << std::to_string(sim) << "-------------" << std::endl; 

            // data 
            GeoFrame data_srpde_d(D, T);
            auto& l_srpde_d = data_srpde_d.insert_scalar_layer<POINT, POINT>("layer", std::pair{R_path + "space_locs.csv", R_path + "time_locs.csv"});
            l_srpde_d.load_csv<double>(R_path + "X_dummies.csv");

            // load data from .csv files
            l_srpde_d.load_csv<double>(R_path + "simulations/sim_" + std::to_string(sim) + "/y_cpp.csv");
            std::string solutions_path_gcv = R_path + "simulations/sim_" + std::to_string(sim) + "/fit_newlib_srpde_d/"; 
            std::string solution_path = R_path + "simulations/sim_" + std::to_string(sim) + "/fit_newlib_srpde_d/"; 
    
            // physics 
            FeSpace Vh(D, P1<1>);   // functional space definition
            Eigen::Matrix<double, 2, 2> K = read_csv<double>(R_path + "simulations/sim_" + std::to_string(sim) + "/K_dummies.csv").as_matrix(); 
            TrialFunction f(Vh);
            TestFunction v(Vh);
            auto a_D = integral(D)(dot(K * grad(f), grad(v)));
            // homogeneous forcing linear form
            ZeroField<2> u_D;
            auto F_D = integral(D)(u_D * v);
    
            // Start measuring time
            auto start_time_gcv = high_resolution_clock::now();
    
            // modeling
            SRPDE m("y ~ x1 + x2 + dummy1 + dummy2 + dummy3 + dummy4 + dummy5 + f", data_srpde_d, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

            // calibration
            GridSearch<2> opt;   // dimension 2 for space-time problems 
            opt.optimize(m.gcv(100, 1234), lambdas_mat);  // stochastic GCV 

            // Stop measuring time
            auto stop_time_gcv = high_resolution_clock::now();
            auto duration_gcv = duration_cast<milliseconds>(stop_time_gcv - start_time_gcv).count();
            std::cout << "Execution time GCV: " << duration_gcv << " ms" << std::endl;

            Eigen::Matrix<double, Dynamic, 1> best_lambda = opt.optimum();
            std::cout << "Best lambdas are: " << std::setprecision(16) << best_lambda << std::endl; 
    
            // Save lambda sequence 
            write_csv(solutions_path_gcv + "lambdas_seq_S.csv", lambdas_d);
            write_csv(solutions_path_gcv + "lambdas_seq_T.csv", lambdas_t);

            std::ofstream fileLambdaoptS(solutions_path_gcv + "lambda_s_opt.csv");
            if(fileLambdaoptS.is_open()){
              fileLambdaoptS << std::setprecision(16) << best_lambda(0,0);
              fileLambdaoptS.close();
            }
            std::ofstream fileLambdaoptT(solutions_path_gcv + "lambda_t_opt.csv");
            if(fileLambdaoptT.is_open()){
              fileLambdaoptT << std::setprecision(16) << best_lambda(1,0);
              fileLambdaoptT.close();
            }
    
            write_csv(solutions_path_gcv + "score.csv", opt.values());

            std::ofstream file_time_gcv(solutions_path_gcv + "time_gcv.csv"); 
            if(file_time_gcv.is_open()){
                file_time_gcv << duration_gcv << "\n";
                file_time_gcv.close();
            }
            
    
        }
    

    }
    

    return 0;
}

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


// int test_06(); 
int test_06_scalability(); 
// int test_07();
// int test_08();  

int main(){
    // test_06();
    test_06_scalability();
    // test_07();
    // test_08();
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
//    missing:      no
int test_06() {

    const unsigned int sim_start = 1; 
    const unsigned int n_sim = 30; 

    // run SRPDE and/or MSRPDE ? 
    const bool run_srpde = true;    // stprde
    const bool run_msrpde = true;   // mixed-effects anisotropic
    const bool run_msr_iso = true;  // mixed-effects isotropic
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

    // Simulations MSRPDE  
    if(run_msrpde){

        for(auto sim = sim_start; sim <= n_sim; ++sim){

            std::cout << "--------------------Simulation RUN MSRPDE #" << std::to_string(sim) << "-------------" << std::endl; 
    
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
            std::cout << "K = " << K << std::endl;
            
            TrialFunction f(Vh);
            TestFunction v(Vh);
            auto a_D = integral(D)(dot(K * grad(f), grad(v)));
            // homogeneous forcing linear form
            ZeroField<2> u_D;
            auto F_D = integral(D)(u_D * v);

            // read lambdas
            double lambda_D;  
            double lambda_T;  
    
            std::ifstream fileLambdaS_gcv(solution_path + "/lambda_s_opt.csv");
            if(fileLambdaS_gcv.is_open()){
                fileLambdaS_gcv >> lambda_D; 
                fileLambdaS_gcv.close();
            }
            std::ifstream fileLambdaT(solution_path + "/lambda_t_opt.csv");
            if(fileLambdaT.is_open()){
                fileLambdaT >> lambda_T; 
                fileLambdaT.close();
            }

            std::cout << "Optimal lambda_D: " << std::setprecision(16) << lambda_D << std::endl;
            std::cout << "Optimal lambda_T: " << std::setprecision(16) << lambda_T << std::endl;

            // Start measuring time
            auto start_time_run = high_resolution_clock::now();
    
            // modeling
            MSRPDE m("y ~ x1 + x2 + 1|g + f", data_msrpde, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

            m.set_fpirls_max_iter(max_fpirls_iter);
            m.set_likelihood_dataloss_type(likelihood_dataloss_type);
            m.set_compute_sigma_with_edf(sigma_edf_type);

            // fit at optimal smoothing level
            m.fit(lambda_D, lambda_T);

            // Stop measuring time
            auto stop_time_run = high_resolution_clock::now();
            auto duration_run = duration_cast<milliseconds>(stop_time_run - start_time_run).count();
            std::cout << "Execution time RUN: " << duration_run << " ms" << std::endl;

            // Save results 
            write_csv(solution_path + "f.csv", m.f());
            write_csv(solution_path + "fn.csv", m.fn());
            write_csv(solution_path + "beta.csv", m.beta());

            Eigen::Matrix<double, Dynamic, Dynamic> computed_b;
            computed_b.resize(m.b_hat().size(), m.n_random_covs());  
            for(int i=0; i<m.b_hat().size(); ++i){
                computed_b.row(i) = m.b_hat()[i].transpose();   // NOTE: .transpose() is important to have the correct shape and save all the values in the case with >1 RE; 
            }
            write_csv(solution_path + "b_random.csv", computed_b);
    
            double computedsigmahat = std::sqrt(m.sigma_sq_hat());
            std::ofstream filesigmahat(solution_path + "/sigma_hat.csv");
            if(filesigmahat.is_open()){
                filesigmahat << std::setprecision(16) << computedsigmahat << "\n"; 
                filesigmahat.close();
            }
    
            write_csv(solution_path + "Sigma_b_hat.csv", m.Sigma_b());

    
            std::ofstream file_time_run(solution_path + "/time_run.csv"); 
            if(file_time_run.is_open()){
                file_time_run << duration_run << "\n"; 
                file_time_run.close();
            }

            std::ofstream filen_iter(solution_path + "/n_iter.csv");
            if(filen_iter.is_open()){
                filen_iter << m.n_iter() << "\n"; 
                filen_iter.close();
            }



        }
    

    }

    // Simulations STRPDE  
    if(run_srpde){

        for(auto sim = sim_start; sim <= n_sim; ++sim){

            std::cout << "--------------------Simulation RUN STRPDE #" << std::to_string(sim) << "-------------" << std::endl; 
    
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
            std::cout << "K = " << K << std::endl;
            
            TrialFunction f(Vh);
            TestFunction v(Vh);
            auto a_D = integral(D)(dot(K * grad(f), grad(v)));
            // homogeneous forcing linear form
            ZeroField<2> u_D;
            auto F_D = integral(D)(u_D * v);

            // read lambdas
            double lambda_D;  
            double lambda_T;  
    
            std::ifstream fileLambdaS_gcv(solution_path + "/lambda_s_opt.csv");
            if(fileLambdaS_gcv.is_open()){
                fileLambdaS_gcv >> lambda_D; 
                fileLambdaS_gcv.close();
            }
            std::ifstream fileLambdaT(solution_path + "/lambda_t_opt.csv");
            if(fileLambdaT.is_open()){
                fileLambdaT >> lambda_T; 
                fileLambdaT.close();
            }

            std::cout << "Optimal lambda_D: " << std::setprecision(16) << lambda_D << std::endl;
            std::cout << "Optimal lambda_T: " << std::setprecision(16) << lambda_T << std::endl;

            // Start measuring time
            auto start_time_run = high_resolution_clock::now();
    
            // modeling
            SRPDE m("y ~ x1 + x2 + f", data_srpde, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

            // fit at optimal smoothing level
            m.fit(lambda_D, lambda_T);

            // Stop measuring time
            auto stop_time_run = high_resolution_clock::now();
            auto duration_run = duration_cast<milliseconds>(stop_time_run - start_time_run).count();
            std::cout << "Execution time RUN: " << duration_run << " ms" << std::endl;

            // Save results 
            write_csv(solution_path + "f.csv", m.f());
            write_csv(solution_path + "fn.csv", m.fn());
            write_csv(solution_path + "beta.csv", m.beta());


            std::ofstream file_time_run(solution_path + "/time_run.csv"); 
            if(file_time_run.is_open()){
                file_time_run << duration_run << "\n"; 
                file_time_run.close();
            }



        }
    

    }

    // Simulations MSR-ISO  
    if(run_msr_iso){

        for(auto sim = sim_start; sim <= n_sim; ++sim){

            std::cout << "--------------------Simulation RUN MSR-ISO #" << std::to_string(sim) << "-------------" << std::endl;

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
            std::cout << "K = " << K << std::endl;
            
            TrialFunction f(Vh);
            TestFunction v(Vh);
            auto a_D = integral(D)(dot(K * grad(f), grad(v)));
            // homogeneous forcing linear form
            ZeroField<2> u_D;
            auto F_D = integral(D)(u_D * v);

            // read lambdas
            double lambda_D;  
            double lambda_T;  
    
            std::ifstream fileLambdaS_gcv(solution_path + "/lambda_s_opt.csv");
            if(fileLambdaS_gcv.is_open()){
                fileLambdaS_gcv >> lambda_D; 
                fileLambdaS_gcv.close();
            }
            std::ifstream fileLambdaT(solution_path + "/lambda_t_opt.csv");
            if(fileLambdaT.is_open()){
                fileLambdaT >> lambda_T; 
                fileLambdaT.close();
            }

            std::cout << "Optimal lambda_D: " << std::setprecision(16) << lambda_D << std::endl;
            std::cout << "Optimal lambda_T: " << std::setprecision(16) << lambda_T << std::endl;

            // Start measuring time
            auto start_time_run = high_resolution_clock::now();
    
            // modeling
            MSRPDE m("y ~ x1 + x2 + 1|g + f", data_msr_iso, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

            m.set_fpirls_max_iter(max_fpirls_iter);
            m.set_likelihood_dataloss_type(likelihood_dataloss_type);
            m.set_compute_sigma_with_edf(sigma_edf_type);

            // fit at optimal smoothing level
            m.fit(lambda_D, lambda_T);

            // Stop measuring time
            auto stop_time_run = high_resolution_clock::now();
            auto duration_run = duration_cast<milliseconds>(stop_time_run - start_time_run).count();
            std::cout << "Execution time RUN: " << duration_run << " ms" << std::endl;

            // Save results 
            write_csv(solution_path + "f.csv", m.f());
            write_csv(solution_path + "fn.csv", m.fn());
            write_csv(solution_path + "beta.csv", m.beta());

            Eigen::Matrix<double, Dynamic, Dynamic> computed_b;
            computed_b.resize(m.b_hat().size(), m.n_random_covs());  
            for(int i=0; i<m.b_hat().size(); ++i){
                computed_b.row(i) = m.b_hat()[i].transpose();   // NOTE: .transpose() is important to have the correct shape and save all the values in the case with >1 RE; 
            }
            write_csv(solution_path + "b_random.csv", computed_b);
    

            double computedsigmahat = std::sqrt(m.sigma_sq_hat());
            std::ofstream filesigmahat(solution_path + "/sigma_hat.csv");
            if(filesigmahat.is_open()){
                filesigmahat << std::setprecision(16) << computedsigmahat << "\n"; 
                filesigmahat.close();
            }
    
            write_csv(solution_path + "Sigma_b_hat.csv", m.Sigma_b());

    
            std::ofstream file_time_run(solution_path + "/time_run.csv"); 
            if(file_time_run.is_open()){
                file_time_run << duration_run << "\n"; 
                file_time_run.close();
            }

            std::ofstream filen_iter(solution_path + "/n_iter.csv");
            if(filen_iter.is_open()){
                filen_iter << m.n_iter() << "\n"; 
                filen_iter.close();
            }




        }
    

    }

    // Simulations STRPDE DUMMIES  
    if(run_srpde_d){

        for(auto sim = sim_start; sim <= n_sim; ++sim){

            std::cout << "--------------------Simulation RUN STRPDE DUMMIES  #" << std::to_string(sim) << "-------------" << std::endl; 
    
            // data 
            GeoFrame data_srpde(D, T);
            auto& l_srpde = data_srpde.insert_scalar_layer<POINT, POINT>("layer", std::pair{R_path + "space_locs.csv", R_path + "time_locs.csv"});
            l_srpde.load_csv<double>(R_path + "X_dummies.csv");

            // load data from .csv files
            l_srpde.load_csv<double>(R_path + "simulations/sim_" + std::to_string(sim) + "/y_cpp.csv");
                    
            std::string solutions_path_gcv = R_path + "simulations/sim_" + std::to_string(sim) + "/fit_newlib_srpde_d/"; 
            std::string solution_path = R_path + "simulations/sim_" + std::to_string(sim) + "/fit_newlib_srpde_d/"; 
    
            // physics 
            FeSpace Vh(D, P1<1>);   // functional space definition
            Eigen::Matrix<double, 2, 2> K = read_csv<double>(R_path + "simulations/sim_" + std::to_string(sim) + "/K_dummies.csv").as_matrix(); 
            std::cout << "K = " << K << std::endl;
            
            TrialFunction f(Vh);
            TestFunction v(Vh);
            auto a_D = integral(D)(dot(K * grad(f), grad(v)));
            // homogeneous forcing linear form
            ZeroField<2> u_D;
            auto F_D = integral(D)(u_D * v);

            // read lambdas
            double lambda_D;  
            double lambda_T;  
    
            std::ifstream fileLambdaS_gcv(solution_path + "/lambda_s_opt.csv");
            if(fileLambdaS_gcv.is_open()){
                fileLambdaS_gcv >> lambda_D; 
                fileLambdaS_gcv.close();
            }
            std::ifstream fileLambdaT(solution_path + "/lambda_t_opt.csv");
            if(fileLambdaT.is_open()){
                fileLambdaT >> lambda_T; 
                fileLambdaT.close();
            }

            std::cout << "Optimal lambda_D: " << std::setprecision(16) << lambda_D << std::endl;
            std::cout << "Optimal lambda_T: " << std::setprecision(16) << lambda_T << std::endl;

            // Start measuring time
            auto start_time_run = high_resolution_clock::now();
    
            // modeling
            SRPDE m("y ~ x1 + x2 + dummy1 + dummy2 + dummy3 + dummy4 + dummy5 + f", data_srpde, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

            // fit at optimal smoothing level
            m.fit(lambda_D, lambda_T);

            // Stop measuring time
            auto stop_time_run = high_resolution_clock::now();
            auto duration_run = duration_cast<milliseconds>(stop_time_run - start_time_run).count();
            std::cout << "Execution time RUN: " << duration_run << " ms" << std::endl;

            // Save results 
            write_csv(solution_path + "f.csv", m.f());
            write_csv(solution_path + "fn.csv", m.fn());
            write_csv(solution_path + "beta.csv", m.beta());

            std::ofstream file_time_run(solution_path + "/time_run.csv"); 
            if(file_time_run.is_open()){
                file_time_run << duration_run << "\n"; 
                file_time_run.close();
            }



        }
    

    }



    return 0;
}


// test 6 scalability (for cluster)
//    mesh:         unit square
//    sampling:     locations != nodes
//    penalization: anisotropic diffusion
//    time penalization: separable
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    missing:      no
int test_06_scalability() {

    const bool scale_n = true; 
    std::vector<unsigned int> nn_vec; 
    if(scale_n){
        nn_vec = {100, 200, 400, 800, 1600, 3200}; 
    }

    const unsigned int sim_start = 1; 
    const unsigned int n_sim = 50; 

    const bool run_msrpde = true;   // mixed-effects anisotropic
    const bool run_msr_iso = true;  // mixed-effects isotropic

    bool likelihood_dataloss_type; // false = fpirls data loss, true = likelihood
    bool sigma_edf_type;           // false = sigma senza edf nelle iterazioni, true = sigma con edf

    const std::string trial_number = "14"; 
    std::string R_path = "/u/desanctis/R_scripts/Test_6_scalability/trial_" + trial_number + "/";
    if(scale_n){
        R_path += "scale_n/";
    }


    likelihood_dataloss_type = false;   // false: FPIRLS data loss; true: likelihood data loss 
    sigma_edf_type = true;
  

    unsigned int M; 
    M = 8; 
    
    
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 1, M);  // ATT qui non bisogna fare più M-1 come nella vecchia lib!! Vuole direttamente il numero di nodi, cioè M!!! 
    
    
    // geometry 
    using PointT = Eigen::Matrix<double, 2, 1>;

    
    std::string N_string; 
    N_string = "476"; 
    
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


    // // Read true solutions for later RMSE computation
    // std::string true_path = "/u/desanctis/R_scripts/Test_6_scalability/trial_" + trial_number + "/true"; 
    // std::string X_eval_path = "/u/desanctis/R_scripts/Test_6_scalability/trial_" + trial_number; 
    
    // vector_t betas_true;  
    // std::ifstream file_betas_true(true_path + "/betas.csv");
    // if(file_betas_true.is_open()){
    //     file_betas_true >> betas_true; 
    //     file_betas_true.close();
    // }

    // vector_t f_true_eval;  
    // std::ifstream file_f_true_eval(true_path + "/spate.sim/f_grf.csv");
    // if(file_f_true_eval.is_open()){
    //     file_f_true_eval >> f_true_eval; 
    //     file_f_true_eval.close();
    // }

    // matrix_t X_eval; 
    // std::ifstream file_X_eval(X_eval_path + "/X_eval.csv");
    // if(file_X_eval.is_open()){
    //     file_X_eval >> X_eval; 
    //     file_X_eval.close();
    // }

    // vector_t mu_true_eval = f_true_eval + X_eval * betas_true;


    // Simulations MSRPDE  
    if(run_msrpde){


        if(scale_n){

            for(unsigned int nn : nn_vec){

                std::string path_nn = R_path + "n_" + std::to_string(nn) + "/";
                std::cout << "======== Runnig n = " << nn << " =========" << std::endl;

                for(auto sim = sim_start; sim <= n_sim; ++sim){

                    std::cout << "Simulation RUN MSRPDE #" << std::to_string(sim) << std::endl; 
            
                    // data 
                    GeoFrame data_msrpde(D, T);
                    auto& l_msrpde = data_msrpde.insert_scalar_layer<POINT, POINT>("layer", std::pair{path_nn + "space_locs.csv", R_path + "time_locs.csv"});
                    // NOTA: nel caso scale_n = true, le space_locs.csv sono in /n_***, mentre le time_locs.csv sono sempre quelle in R_path
                    
                    l_msrpde.load_csv<double>(path_nn + "X.csv");
                    l_msrpde.load_csv<double>(path_nn + "ids_groups.csv");

                    // load data from .csv files
                    l_msrpde.load_csv<double>(path_nn + "simulations/sim_" + std::to_string(sim) + "/y_cpp.csv");
                            
                    std::string solutions_path_gcv = path_nn + "simulations/sim_" + std::to_string(sim) + "/fit_newlib/"; 
                    std::string solution_path = path_nn + "simulations/sim_" + std::to_string(sim) + "/fit_newlib/"; 
            
                    // physics 
                    FeSpace Vh(D, P1<1>);   // functional space definition
                    Eigen::Matrix<double, 2, 2> K = read_csv<double>(path_nn + "simulations/sim_" + std::to_string(sim) + "/K.csv").as_matrix(); 
                    std::cout << "K = " << K << std::endl;
                    
                    TrialFunction f(Vh);
                    TestFunction v(Vh);
                    auto a_D = integral(D)(dot(K * grad(f), grad(v)));
                    // homogeneous forcing linear form
                    ZeroField<2> u_D;
                    auto F_D = integral(D)(u_D * v);

                    // read lambdas
                    double lambda_D;  
                    double lambda_T;  
            
                    std::ifstream fileLambdaS_gcv(solution_path + "/lambda_s_opt.csv");
                    if(fileLambdaS_gcv.is_open()){
                        fileLambdaS_gcv >> lambda_D; 
                        fileLambdaS_gcv.close();
                    }
                    std::ifstream fileLambdaT(solution_path + "/lambda_t_opt.csv");
                    if(fileLambdaT.is_open()){
                        fileLambdaT >> lambda_T; 
                        fileLambdaT.close();
                    }

                    // std::cout << "Optimal lambda_D: " << std::setprecision(16) << lambda_D << std::endl;
                    // std::cout << "Optimal lambda_T: " << std::setprecision(16) << lambda_T << std::endl;

                    // Start measuring time
                    auto start_time_run = high_resolution_clock::now();
            
                    // modeling
                    MSRPDE m("y ~ x1 + x2 + 1|g + f", data_msrpde, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

                    m.set_fpirls_max_iter(max_fpirls_iter);
                    m.set_likelihood_dataloss_type(likelihood_dataloss_type);
                    m.set_compute_sigma_with_edf(sigma_edf_type);

                    // fit at optimal smoothing level
                    m.fit(lambda_D, lambda_T);

                    // Stop measuring time
                    auto stop_time_run = high_resolution_clock::now();
                    auto duration_run = duration_cast<milliseconds>(stop_time_run - start_time_run).count();
                    std::cout << "Execution time RUN: " << duration_run << " ms" << std::endl;

                    // Save results 
                    write_csv(solution_path + "f.csv", m.f());    
                    write_csv(solution_path + "beta.csv", m.beta());

                    // Eigen::Matrix<double, Dynamic, Dynamic> computed_b;
                    // computed_b.resize(m.b_hat().size(), m.n_random_covs());  
                    // for(int i=0; i<m.b_hat().size(); ++i){
                    //     computed_b.row(i) = m.b_hat()[i].transpose();   // NOTE: .transpose() is important to have the correct shape and save all the values in the case with >1 RE; 
                    // }
                    // write_csv(solution_path + "b_random.csv", computed_b);
            
                    double computedsigmahat = std::sqrt(m.sigma_sq_hat());
                    std::ofstream filesigmahat(solution_path + "/sigma_hat.csv");
                    if(filesigmahat.is_open()){
                        filesigmahat << std::setprecision(16) << computedsigmahat << "\n"; 
                        filesigmahat.close();
                    }
            
                    write_csv(solution_path + "Sigma_b_hat.csv", m.Sigma_b());

            
                    std::ofstream file_time_run(solution_path + "/time_run.csv"); 
                    if(file_time_run.is_open()){
                        file_time_run << duration_run << "\n"; 
                        file_time_run.close();
                    }

                    std::ofstream filen_iter(solution_path + "/n_iter.csv");
                    if(filen_iter.is_open()){
                        filen_iter << m.n_iter() << "\n"; 
                        filen_iter.close();
                    }


                }
    
            }

        }



    }

    // Simulations MSR-ISO  
    if(run_msr_iso){


        if(scale_n){

            for(unsigned int nn : nn_vec){

                std::string path_nn = R_path + "n_" + std::to_string(nn) + "/";
                std::cout << "======== Runnig n = " << nn << " =========" << std::endl;


                for(auto sim = sim_start; sim <= n_sim; ++sim){

                    std::cout << "--------------------Simulation RUN MSR-ISO #" << std::to_string(sim) << "-------------" << std::endl;

                    // data 
                    GeoFrame data_msr_iso(D, T);
                    auto& l_msr_iso = data_msr_iso.insert_scalar_layer<POINT, POINT>("layer", std::pair{path_nn + "space_locs.csv", R_path + "time_locs.csv"});
                    // NOTA: nel caso scale_n = true, le space_locs.csv sono in /n_***, mentre le time_locs.csv sono sempre quelle in R_path
                    
                    l_msr_iso.load_csv<double>(path_nn + "X.csv");
                    l_msr_iso.load_csv<double>(path_nn + "ids_groups.csv");

                    // load data from .csv files
                    l_msr_iso.load_csv<double>(path_nn + "simulations/sim_" + std::to_string(sim) + "/y_cpp.csv");
                            
                    std::string solutions_path_gcv = path_nn + "simulations/sim_" + std::to_string(sim) + "/fit_newlib_iso/"; 
                    std::string solution_path = path_nn + "simulations/sim_" + std::to_string(sim) + "/fit_newlib_iso/"; 
            
                    // physics 
                    FeSpace Vh(D, P1<1>);   // functional space definition
                    Eigen::Matrix<double, 2, 2> K; 
                    K << 1, 0, 0, 1;        // isotropic case 
                    // std::cout << "K = " << K << std::endl;
                    
                    TrialFunction f(Vh);
                    TestFunction v(Vh);
                    auto a_D = integral(D)(dot(K * grad(f), grad(v)));
                    // homogeneous forcing linear form
                    ZeroField<2> u_D;
                    auto F_D = integral(D)(u_D * v);

                    // read lambdas
                    double lambda_D;  
                    double lambda_T;  
            
                    std::ifstream fileLambdaS_gcv(solution_path + "/lambda_s_opt.csv");
                    if(fileLambdaS_gcv.is_open()){
                        fileLambdaS_gcv >> lambda_D; 
                        fileLambdaS_gcv.close();
                    }
                    std::ifstream fileLambdaT(solution_path + "/lambda_t_opt.csv");
                    if(fileLambdaT.is_open()){
                        fileLambdaT >> lambda_T; 
                        fileLambdaT.close();
                    }

                    // std::cout << "Optimal lambda_D: " << std::setprecision(16) << lambda_D << std::endl;
                    // std::cout << "Optimal lambda_T: " << std::setprecision(16) << lambda_T << std::endl;

                    // Start measuring time
                    auto start_time_run = high_resolution_clock::now();
            
                    // modeling
                    MSRPDE m("y ~ x1 + x2 + 1|g + f", data_msr_iso, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

                    m.set_fpirls_max_iter(max_fpirls_iter);
                    m.set_likelihood_dataloss_type(likelihood_dataloss_type);
                    m.set_compute_sigma_with_edf(sigma_edf_type);

                    // fit at optimal smoothing level
                    m.fit(lambda_D, lambda_T);

                    // Stop measuring time
                    auto stop_time_run = high_resolution_clock::now();
                    auto duration_run = duration_cast<milliseconds>(stop_time_run - start_time_run).count();
                    std::cout << "Execution time RUN: " << duration_run << " ms" << std::endl;

                    // Save results 
                    write_csv(solution_path + "f.csv", m.f());
                    write_csv(solution_path + "beta.csv", m.beta());

                    // Eigen::Matrix<double, Dynamic, Dynamic> computed_b;
                    // computed_b.resize(m.b_hat().size(), m.n_random_covs());  
                    // for(int i=0; i<m.b_hat().size(); ++i){
                    //     computed_b.row(i) = m.b_hat()[i].transpose();   // NOTE: .transpose() is important to have the correct shape and save all the values in the case with >1 RE; 
                    // }
                    // write_csv(solution_path + "b_random.csv", computed_b);
            

                    double computedsigmahat = std::sqrt(m.sigma_sq_hat());
                    std::ofstream filesigmahat(solution_path + "/sigma_hat.csv");
                    if(filesigmahat.is_open()){
                        filesigmahat << std::setprecision(16) << computedsigmahat << "\n"; 
                        filesigmahat.close();
                    }
            
                    write_csv(solution_path + "Sigma_b_hat.csv", m.Sigma_b());

            
                    std::ofstream file_time_run(solution_path + "/time_run.csv"); 
                    if(file_time_run.is_open()){
                        file_time_run << duration_run << "\n"; 
                        file_time_run.close();
                    }

                    std::ofstream filen_iter(solution_path + "/n_iter.csv");
                    if(filen_iter.is_open()){
                        filen_iter << m.n_iter() << "\n"; 
                        filen_iter.close();
                    }



                }
    

            }

        }


    }



    return 0;
}


// test 7
//    mesh:         unit square
//    sampling:     locations != nodes
//    penalization: anisotropic diffusion
//    time penalization: separable
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    missing:      yes
int test_07() {

    const std::string schema = "d"; 

    const unsigned int sim_start = 1; 
    const unsigned int n_sim = 30; 

    // run SRPDE and/or MSRPDE ? 
    const bool run_srpde = false;    // stprde
    const bool run_msrpde = true;   // mixed-effects anisotropic
    const bool run_msr_iso = false;  // mixed-effects isotropic
    const bool run_srpde_d = false;   // strpde con dummies

    bool likelihood_dataloss_type; // false = fpirls data loss, true = likelihood
    bool sigma_edf_type;           // false = sigma senza edf nelle iterazioni, true = sigma con edf

    const std::string trial_number = "2"; 
    const std::string R_path = "../../../OneDrive - Politecnico di Milano/Corsi/PhD/Codice/models/MSRPDE/Tests/space-time/Test_7/trial_" + trial_number + "/miss_" + schema + "/";

    if(trial_number == "1" || trial_number == "2" || trial_number == "3" || trial_number == "4" || trial_number == "5"){
        likelihood_dataloss_type = false;
    } 
    if(trial_number == "1" || trial_number == "2" || trial_number == "3" || trial_number == "4" || trial_number == "5"){
        sigma_edf_type = true;
    } 

    unsigned int M; 
    if(schema == "a"){
        if(trial_number == "1"){
            M = 8; 
        }
    }
    if(schema == "b"){
        if(trial_number == "1"|| trial_number == "3"){
            M = 8; 
        }
    }
    if(schema == "d"){
        if(trial_number == "1" || trial_number == "4" || trial_number == "5"){
            M = 8; 
        }
        if(trial_number == "2"){
            M = 16; 
        }
    }
    
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 1, M);  // ATT qui non bisogna fare più M-1 come nella vecchia lib!! Vuole direttamente il numero di nodi, cioè M!!! 
    
    
    // geometry 
    using PointT = Eigen::Matrix<double, 2, 1>;

    
    std::string N_string; 
    if(trial_number == "1" || trial_number == "2" || trial_number == "3" || trial_number == "4" || trial_number == "5"){
        N_string = "476"; 
    }
    Eigen::Matrix<double, Dynamic, Dynamic> points = read_csv<double>("my_data/mesh/unit_square_reduced_censoring_" + N_string + "/points.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> elements = read_csv<int>("my_data/mesh/unit_square_reduced_censoring_" + N_string + "/elements.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> boundary = read_csv<int>("my_data/mesh/unit_square_reduced_censoring_" + N_string + "/boundary.csv").as_matrix();

    elements.array() -= 1; // non necessario
    
    Triangulation<2, 2> D(points, elements, boundary);

    unsigned int max_fpirls_iter = 15; 

    // time penalty 
    BsSpace Bh(T, 3);   // cubic B-splines in time
    TrialFunction g(Bh);
    TestFunction  w(Bh);
    auto a_T = integral(T)(dxx(g) * dxx(w));
    ZeroField<1> u_T;
    auto F_T = integral(T)(u_T * w);  
    
    
    std::cout << "R_path: " << R_path << std::endl;

    // Simulations MSRPDE  
    if(run_msrpde){

        for(auto sim = sim_start; sim <= n_sim; ++sim){

            std::cout << "--------------------Simulation RUN MSRPDE #" << std::to_string(sim) << "-------------" << std::endl; 
    
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
            std::cout << "K = " << K << std::endl;
            
            TrialFunction f(Vh);
            TestFunction v(Vh);
            auto a_D = integral(D)(dot(K * grad(f), grad(v)));
            // homogeneous forcing linear form
            ZeroField<2> u_D;
            auto F_D = integral(D)(u_D * v);

            // read lambdas
            double lambda_D;  
            double lambda_T;  
    
            std::ifstream fileLambdaS_gcv(solution_path + "/lambda_s_opt.csv");
            if(fileLambdaS_gcv.is_open()){
                fileLambdaS_gcv >> lambda_D; 
                fileLambdaS_gcv.close();
            }
            std::ifstream fileLambdaT(solution_path + "/lambda_t_opt.csv");
            if(fileLambdaT.is_open()){
                fileLambdaT >> lambda_T; 
                fileLambdaT.close();
            }

            std::cout << "Optimal lambda_D: " << std::setprecision(16) << lambda_D << std::endl;
            std::cout << "Optimal lambda_T: " << std::setprecision(16) << lambda_T << std::endl;

            // Start measuring time
            auto start_time_run = high_resolution_clock::now();
    
            // modeling
            MSRPDE m("y ~ x1 + x2 + 1|g + f", data_msrpde, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

            m.set_fpirls_max_iter(max_fpirls_iter);
            m.set_likelihood_dataloss_type(likelihood_dataloss_type);
            m.set_compute_sigma_with_edf(sigma_edf_type);

            // fit at optimal smoothing level
            m.fit(lambda_D, lambda_T);

            // Stop measuring time
            auto stop_time_run = high_resolution_clock::now();
            auto duration_run = duration_cast<milliseconds>(stop_time_run - start_time_run).count();
            std::cout << "Execution time RUN: " << duration_run << " ms" << std::endl;

            // Save results 
            write_csv(solution_path + "f.csv", m.f());
            write_csv(solution_path + "fn.csv", m.fn());
            write_csv(solution_path + "beta.csv", m.beta());

            Eigen::Matrix<double, Dynamic, Dynamic> computed_b;
            computed_b.resize(m.b_hat().size(), m.n_random_covs());  
            for(int i=0; i<m.b_hat().size(); ++i){
                computed_b.row(i) = m.b_hat()[i].transpose();   // NOTE: .transpose() is important to have the correct shape and save all the values in the case with >1 RE; 
            }
            write_csv(solution_path + "b_random.csv", computed_b);
    
            double computedsigmahat = std::sqrt(m.sigma_sq_hat());
            std::ofstream filesigmahat(solution_path + "/sigma_hat.csv");
            if(filesigmahat.is_open()){
                filesigmahat << std::setprecision(16) << computedsigmahat << "\n"; 
                filesigmahat.close();
            }
    
            write_csv(solution_path + "Sigma_b_hat.csv", m.Sigma_b());

    
            std::ofstream file_time_run(solution_path + "/time_run.csv"); 
            if(file_time_run.is_open()){
                file_time_run << duration_run << "\n"; 
                file_time_run.close();
            }

            std::ofstream filen_iter(solution_path + "/n_iter.csv");
            if(filen_iter.is_open()){
                filen_iter << m.n_iter() << "\n"; 
                filen_iter.close();
            }



        }
    

    }

    // Simulations STRPDE  
    if(run_srpde){

        for(auto sim = sim_start; sim <= n_sim; ++sim){

            std::cout << "--------------------Simulation RUN STRPDE #" << std::to_string(sim) << "-------------" << std::endl; 
    
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
            std::cout << "K = " << K << std::endl;
            
            TrialFunction f(Vh);
            TestFunction v(Vh);
            auto a_D = integral(D)(dot(K * grad(f), grad(v)));
            // homogeneous forcing linear form
            ZeroField<2> u_D;
            auto F_D = integral(D)(u_D * v);

            // read lambdas
            double lambda_D;  
            double lambda_T;  
    
            std::ifstream fileLambdaS_gcv(solution_path + "/lambda_s_opt.csv");
            if(fileLambdaS_gcv.is_open()){
                fileLambdaS_gcv >> lambda_D; 
                fileLambdaS_gcv.close();
            }
            std::ifstream fileLambdaT(solution_path + "/lambda_t_opt.csv");
            if(fileLambdaT.is_open()){
                fileLambdaT >> lambda_T; 
                fileLambdaT.close();
            }

            std::cout << "Optimal lambda_D: " << std::setprecision(16) << lambda_D << std::endl;
            std::cout << "Optimal lambda_T: " << std::setprecision(16) << lambda_T << std::endl;

            // Start measuring time
            auto start_time_run = high_resolution_clock::now();
    
            // modeling
            SRPDE m("y ~ x1 + x2 + f", data_srpde, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

            // fit at optimal smoothing level
            m.fit(lambda_D, lambda_T);

            // Stop measuring time
            auto stop_time_run = high_resolution_clock::now();
            auto duration_run = duration_cast<milliseconds>(stop_time_run - start_time_run).count();
            std::cout << "Execution time RUN: " << duration_run << " ms" << std::endl;

            // Save results 
            write_csv(solution_path + "f.csv", m.f());
            write_csv(solution_path + "fn.csv", m.fn());
            write_csv(solution_path + "beta.csv", m.beta());


            std::ofstream file_time_run(solution_path + "/time_run.csv"); 
            if(file_time_run.is_open()){
                file_time_run << duration_run << "\n"; 
                file_time_run.close();
            }



        }
    

    }

    // Simulations MSR-ISO  
    if(run_msr_iso){

        for(auto sim = sim_start; sim <= n_sim; ++sim){

            std::cout << "--------------------Simulation RUN MSR-ISO #" << std::to_string(sim) << "-------------" << std::endl;

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
            std::cout << "K = " << K << std::endl;
            
            TrialFunction f(Vh);
            TestFunction v(Vh);
            auto a_D = integral(D)(dot(K * grad(f), grad(v)));
            // homogeneous forcing linear form
            ZeroField<2> u_D;
            auto F_D = integral(D)(u_D * v);

            // read lambdas
            double lambda_D;  
            double lambda_T;  
    
            std::ifstream fileLambdaS_gcv(solution_path + "/lambda_s_opt.csv");
            if(fileLambdaS_gcv.is_open()){
                fileLambdaS_gcv >> lambda_D; 
                fileLambdaS_gcv.close();
            }
            std::ifstream fileLambdaT(solution_path + "/lambda_t_opt.csv");
            if(fileLambdaT.is_open()){
                fileLambdaT >> lambda_T; 
                fileLambdaT.close();
            }

            std::cout << "Optimal lambda_D: " << std::setprecision(16) << lambda_D << std::endl;
            std::cout << "Optimal lambda_T: " << std::setprecision(16) << lambda_T << std::endl;

            // Start measuring time
            auto start_time_run = high_resolution_clock::now();
    
            // modeling
            MSRPDE m("y ~ x1 + x2 + 1|g + f", data_msr_iso, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

            m.set_fpirls_max_iter(max_fpirls_iter);
            m.set_likelihood_dataloss_type(likelihood_dataloss_type);
            m.set_compute_sigma_with_edf(sigma_edf_type);

            // fit at optimal smoothing level
            m.fit(lambda_D, lambda_T);

            // Stop measuring time
            auto stop_time_run = high_resolution_clock::now();
            auto duration_run = duration_cast<milliseconds>(stop_time_run - start_time_run).count();
            std::cout << "Execution time RUN: " << duration_run << " ms" << std::endl;

            // Save results 
            write_csv(solution_path + "f.csv", m.f());
            write_csv(solution_path + "fn.csv", m.fn());
            write_csv(solution_path + "beta.csv", m.beta());

            Eigen::Matrix<double, Dynamic, Dynamic> computed_b;
            computed_b.resize(m.b_hat().size(), m.n_random_covs());  
            for(int i=0; i<m.b_hat().size(); ++i){
                computed_b.row(i) = m.b_hat()[i].transpose();   // NOTE: .transpose() is important to have the correct shape and save all the values in the case with >1 RE; 
            }
            write_csv(solution_path + "b_random.csv", computed_b);
    

            double computedsigmahat = std::sqrt(m.sigma_sq_hat());
            std::ofstream filesigmahat(solution_path + "/sigma_hat.csv");
            if(filesigmahat.is_open()){
                filesigmahat << std::setprecision(16) << computedsigmahat << "\n"; 
                filesigmahat.close();
            }
    
            write_csv(solution_path + "Sigma_b_hat.csv", m.Sigma_b());

    
            std::ofstream file_time_run(solution_path + "/time_run.csv"); 
            if(file_time_run.is_open()){
                file_time_run << duration_run << "\n"; 
                file_time_run.close();
            }

            std::ofstream filen_iter(solution_path + "/n_iter.csv");
            if(filen_iter.is_open()){
                filen_iter << m.n_iter() << "\n"; 
                filen_iter.close();
            }




        }
    

    }

    // Simulations STRPDE DUMMIES  
    if(run_srpde_d){

        for(auto sim = sim_start; sim <= n_sim; ++sim){

            std::cout << "--------------------Simulation RUN STRPDE DUMMIES  #" << std::to_string(sim) << "-------------" << std::endl; 
    
            // data 
            GeoFrame data_srpde(D, T);
            auto& l_srpde = data_srpde.insert_scalar_layer<POINT, POINT>("layer", std::pair{R_path + "space_locs.csv", R_path + "time_locs.csv"});
            l_srpde.load_csv<double>(R_path + "X_dummies.csv");

            // load data from .csv files
            l_srpde.load_csv<double>(R_path + "simulations/sim_" + std::to_string(sim) + "/y_cpp.csv");
                    
            std::string solutions_path_gcv = R_path + "simulations/sim_" + std::to_string(sim) + "/fit_newlib_srpde_d/"; 
            std::string solution_path = R_path + "simulations/sim_" + std::to_string(sim) + "/fit_newlib_srpde_d/"; 
    
            // physics 
            FeSpace Vh(D, P1<1>);   // functional space definition
            Eigen::Matrix<double, 2, 2> K = read_csv<double>(R_path + "simulations/sim_" + std::to_string(sim) + "/K_dummies.csv").as_matrix(); 
            std::cout << "K = " << K << std::endl;
            
            TrialFunction f(Vh);
            TestFunction v(Vh);
            auto a_D = integral(D)(dot(K * grad(f), grad(v)));
            // homogeneous forcing linear form
            ZeroField<2> u_D;
            auto F_D = integral(D)(u_D * v);

            // read lambdas
            double lambda_D;  
            double lambda_T;  
    
            std::ifstream fileLambdaS_gcv(solution_path + "/lambda_s_opt.csv");
            if(fileLambdaS_gcv.is_open()){
                fileLambdaS_gcv >> lambda_D; 
                fileLambdaS_gcv.close();
            }
            std::ifstream fileLambdaT(solution_path + "/lambda_t_opt.csv");
            if(fileLambdaT.is_open()){
                fileLambdaT >> lambda_T; 
                fileLambdaT.close();
            }

            std::cout << "Optimal lambda_D: " << std::setprecision(16) << lambda_D << std::endl;
            std::cout << "Optimal lambda_T: " << std::setprecision(16) << lambda_T << std::endl;

            // Start measuring time
            auto start_time_run = high_resolution_clock::now();
    
            // modeling
            SRPDE m("y ~ x1 + x2 + dummy1 + dummy2 + dummy3 + dummy4 + dummy5 + f", data_srpde, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

            // fit at optimal smoothing level
            m.fit(lambda_D, lambda_T);

            // Stop measuring time
            auto stop_time_run = high_resolution_clock::now();
            auto duration_run = duration_cast<milliseconds>(stop_time_run - start_time_run).count();
            std::cout << "Execution time RUN: " << duration_run << " ms" << std::endl;

            // Save results 
            write_csv(solution_path + "f.csv", m.f());
            write_csv(solution_path + "fn.csv", m.fn());
            write_csv(solution_path + "beta.csv", m.beta());

            std::ofstream file_time_run(solution_path + "/time_run.csv"); 
            if(file_time_run.is_open()){
                file_time_run << duration_run << "\n"; 
                file_time_run.close();
            }



        }
    

    }



    return 0;
}



// test 8
//    mesh:         unit square
//    sampling:     locations != nodes
//    penalization: anisotropic diffusion
//    time penalization: separable
//    covariates:   yes (q=2, p=3)
//    BC:           no
//    order FE:     1
//    missing:      no
int test_08() {

    const unsigned int sim_start = 1; 
    const unsigned int n_sim = 30; 

    // run SRPDE and/or MSRPDE ? 
    const bool run_srpde = false;     // stprde
    const bool run_msrpde = true;     // mixed-effects anisotropic
    const bool run_msr_iso = false;   // mixed-effects isotropic
    const bool run_srpde_d = false;   // strpde con dummies

    bool likelihood_dataloss_type; // false = fpirls data loss, true = likelihood
    bool sigma_edf_type;           // false = sigma senza edf nelle iterazioni, true = sigma con edf

    const std::string trial_number = "2"; 
    const std::string R_path = "../../../OneDrive - Politecnico di Milano/Corsi/PhD/Codice/models/MSRPDE/Tests/space-time/Test_8/trial_" + trial_number + "/";

    std::string formula_mixed; 
    if(trial_number == "1" || trial_number == "2" || trial_number == "3" || trial_number == "4"){
        formula_mixed = "y ~ x1 + x2 + 1|g + x1|g + x2|g + f";
    } 
    if(trial_number == "5"){
        formula_mixed = "y ~ x1 + x2 + 1|g + x1|g + f";
    }


    // likelihood data loss? 
    if(trial_number == "1" || trial_number == "2" || trial_number == "3" || trial_number == "4" || trial_number == "5"){
        likelihood_dataloss_type = false;
    }

    // sigma with edf?
    if(trial_number == "1" || trial_number == "2" || trial_number == "3" || trial_number == "4" || trial_number == "5"){
        sigma_edf_type = true;
    } 


    unsigned int M; 
    if(trial_number == "1" || trial_number == "2" || trial_number == "3" || trial_number == "4" || trial_number == "5"){
        M = 8; 
    }
    
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 1, M);  // ATT qui non bisogna fare più M-1 come nella vecchia lib!! Vuole direttamente il numero di nodi, cioè M!!! 
    
    
    // geometry 
    using PointT = Eigen::Matrix<double, 2, 1>;

    
    std::string N_string; 
    if(trial_number == "1" || trial_number == "2" || trial_number == "3" || trial_number == "4" || trial_number == "5"){ 
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

    // Simulations MSRPDE  
    if(run_msrpde){

        for(auto sim = sim_start; sim <= n_sim; ++sim){

            std::cout << "--------------------Simulation RUN MSRPDE #" << std::to_string(sim) << "-------------" << std::endl; 
    
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
            std::cout << "K = " << K << std::endl;
            
            TrialFunction f(Vh);
            TestFunction v(Vh);
            auto a_D = integral(D)(dot(K * grad(f), grad(v)));
            // homogeneous forcing linear form
            ZeroField<2> u_D;
            auto F_D = integral(D)(u_D * v);

            // read lambdas
            double lambda_D;  
            double lambda_T;  
    
            std::ifstream fileLambdaS_gcv(solution_path + "/lambda_s_opt.csv");
            if(fileLambdaS_gcv.is_open()){
                fileLambdaS_gcv >> lambda_D; 
                fileLambdaS_gcv.close();
            }
            std::ifstream fileLambdaT(solution_path + "/lambda_t_opt.csv");
            if(fileLambdaT.is_open()){
                fileLambdaT >> lambda_T; 
                fileLambdaT.close();
            }

            std::cout << "Optimal lambda_D: " << std::setprecision(16) << lambda_D << std::endl;
            std::cout << "Optimal lambda_T: " << std::setprecision(16) << lambda_T << std::endl;

            // Start measuring time
            auto start_time_run = high_resolution_clock::now();
    
            // modeling
            MSRPDE m(formula_mixed, data_msrpde, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

            m.set_fpirls_max_iter(max_fpirls_iter);
            m.set_likelihood_dataloss_type(likelihood_dataloss_type);
            m.set_compute_sigma_with_edf(sigma_edf_type);

            // fit at optimal smoothing level
            m.fit(lambda_D, lambda_T);

            // Stop measuring time
            auto stop_time_run = high_resolution_clock::now();
            auto duration_run = duration_cast<milliseconds>(stop_time_run - start_time_run).count();
            std::cout << "Execution time RUN: " << duration_run << " ms" << std::endl;

            // Save results 
            write_csv(solution_path + "f.csv", m.f());
            write_csv(solution_path + "fn.csv", m.fn());
            write_csv(solution_path + "beta.csv", m.beta());

            Eigen::Matrix<double, Dynamic, Dynamic> computed_b;
            computed_b.resize(m.b_hat().size(), m.n_random_covs());  
            for(int i=0; i<m.b_hat().size(); ++i){
                computed_b.row(i) = m.b_hat()[i].transpose();   // NOTE: .transpose() is important to have the correct shape and save all the values in the case with >1 RE 
            }
            write_csv(solution_path + "b_random.csv", computed_b);
    
            double computedsigmahat = std::sqrt(m.sigma_sq_hat());
            std::ofstream filesigmahat(solution_path + "/sigma_hat.csv");
            if(filesigmahat.is_open()){
                filesigmahat << std::setprecision(16) << computedsigmahat << "\n"; 
                filesigmahat.close();
            }
    
            write_csv(solution_path + "Sigma_b_hat.csv", m.Sigma_b());
    
            std::ofstream file_time_run(solution_path + "/time_run.csv"); 
            if(file_time_run.is_open()){
                file_time_run << duration_run << "\n"; 
                file_time_run.close();
            }

            std::ofstream filen_iter(solution_path + "/n_iter.csv");
            if(filen_iter.is_open()){
                filen_iter << m.n_iter() << "\n"; 
                filen_iter.close();
            }

            // Save Delta J at convergence 
            double DeltaJ = std::abs(m.Jnew_debug()-m.Jold_debug()); 
            std::ofstream fileDeltaJ(solution_path + "/DeltaJ.csv");
            if(fileDeltaJ.is_open()){
                fileDeltaJ << std::setprecision(16) << DeltaJ << "\n"; 
                fileDeltaJ.close();
            }

            // Debug: Delta_debug: p x n_iter matrix which at each column has the Delta at each iteration
            write_csv(solution_path + "Delta_debug.csv", m.Delta_debug());


        }
    

    }

    // Simulations STRPDE  
    if(run_srpde){

        for(auto sim = sim_start; sim <= n_sim; ++sim){

            std::cout << "--------------------Simulation RUN STRPDE #" << std::to_string(sim) << "-------------" << std::endl; 
    
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
            std::cout << "K = " << K << std::endl;
            
            TrialFunction f(Vh);
            TestFunction v(Vh);
            auto a_D = integral(D)(dot(K * grad(f), grad(v)));
            // homogeneous forcing linear form
            ZeroField<2> u_D;
            auto F_D = integral(D)(u_D * v);

            // read lambdas
            double lambda_D;  
            double lambda_T;  
    
            std::ifstream fileLambdaS_gcv(solution_path + "/lambda_s_opt.csv");
            if(fileLambdaS_gcv.is_open()){
                fileLambdaS_gcv >> lambda_D; 
                fileLambdaS_gcv.close();
            }
            std::ifstream fileLambdaT(solution_path + "/lambda_t_opt.csv");
            if(fileLambdaT.is_open()){
                fileLambdaT >> lambda_T; 
                fileLambdaT.close();
            }

            std::cout << "Optimal lambda_D: " << std::setprecision(16) << lambda_D << std::endl;
            std::cout << "Optimal lambda_T: " << std::setprecision(16) << lambda_T << std::endl;

            // Start measuring time
            auto start_time_run = high_resolution_clock::now();
    
            // modeling
            SRPDE m("y ~ x1 + x2 + f", data_srpde, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

            // fit at optimal smoothing level
            m.fit(lambda_D, lambda_T);

            // Stop measuring time
            auto stop_time_run = high_resolution_clock::now();
            auto duration_run = duration_cast<milliseconds>(stop_time_run - start_time_run).count();
            std::cout << "Execution time RUN: " << duration_run << " ms" << std::endl;

            // Save results 
            write_csv(solution_path + "f.csv", m.f());
            write_csv(solution_path + "fn.csv", m.fn());
            write_csv(solution_path + "beta.csv", m.beta());


            std::ofstream file_time_run(solution_path + "/time_run.csv"); 
            if(file_time_run.is_open()){
                file_time_run << duration_run << "\n"; 
                file_time_run.close();
            }



        }
    

    }

    // Simulations MSR-ISO  
    if(run_msr_iso){

        for(auto sim = sim_start; sim <= n_sim; ++sim){

            std::cout << "--------------------Simulation RUN MSR-ISO #" << std::to_string(sim) << "-------------" << std::endl;

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
            std::cout << "K = " << K << std::endl;
            
            TrialFunction f(Vh);
            TestFunction v(Vh);
            auto a_D = integral(D)(dot(K * grad(f), grad(v)));
            // homogeneous forcing linear form
            ZeroField<2> u_D;
            auto F_D = integral(D)(u_D * v);

            // read lambdas
            double lambda_D;  
            double lambda_T;  
    
            std::ifstream fileLambdaS_gcv(solution_path + "/lambda_s_opt.csv");
            if(fileLambdaS_gcv.is_open()){
                fileLambdaS_gcv >> lambda_D; 
                fileLambdaS_gcv.close();
            }
            std::ifstream fileLambdaT(solution_path + "/lambda_t_opt.csv");
            if(fileLambdaT.is_open()){
                fileLambdaT >> lambda_T; 
                fileLambdaT.close();
            }

            std::cout << "Optimal lambda_D: " << std::setprecision(16) << lambda_D << std::endl;
            std::cout << "Optimal lambda_T: " << std::setprecision(16) << lambda_T << std::endl;

            // Start measuring time
            auto start_time_run = high_resolution_clock::now();
    
            // modeling
            MSRPDE m(formula_mixed, data_msr_iso, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

            m.set_fpirls_max_iter(max_fpirls_iter);
            m.set_likelihood_dataloss_type(likelihood_dataloss_type);
            m.set_compute_sigma_with_edf(sigma_edf_type);

            // fit at optimal smoothing level
            m.fit(lambda_D, lambda_T);

            // Stop measuring time
            auto stop_time_run = high_resolution_clock::now();
            auto duration_run = duration_cast<milliseconds>(stop_time_run - start_time_run).count();
            std::cout << "Execution time RUN: " << duration_run << " ms" << std::endl;

            // Save results 
            write_csv(solution_path + "f.csv", m.f());
            write_csv(solution_path + "fn.csv", m.fn());
            write_csv(solution_path + "beta.csv", m.beta());

            Eigen::Matrix<double, Dynamic, Dynamic> computed_b;
            computed_b.resize(m.b_hat().size(), m.n_random_covs());  
            for(int i=0; i<m.b_hat().size(); ++i){
                computed_b.row(i) = m.b_hat()[i].transpose();   // NOTE: .transpose() is important to have the correct shape and save all the values in the case with >1 RE; 
            }
            write_csv(solution_path + "b_random.csv", computed_b);
    

            double computedsigmahat = std::sqrt(m.sigma_sq_hat());
            std::ofstream filesigmahat(solution_path + "/sigma_hat.csv");
            if(filesigmahat.is_open()){
                filesigmahat << std::setprecision(16) << computedsigmahat << "\n"; 
                filesigmahat.close();
            }
    
            write_csv(solution_path + "Sigma_b_hat.csv", m.Sigma_b());

    
            std::ofstream file_time_run(solution_path + "/time_run.csv"); 
            if(file_time_run.is_open()){
                file_time_run << duration_run << "\n"; 
                file_time_run.close();
            }

            std::ofstream filen_iter(solution_path + "/n_iter.csv");
            if(filen_iter.is_open()){
                filen_iter << m.n_iter() << "\n"; 
                filen_iter.close();
            }

            // Save Delta J at convergence 
            double DeltaJ = std::abs(m.Jnew_debug()-m.Jold_debug()); 
            std::ofstream fileDeltaJ(solution_path + "/DeltaJ.csv");
            if(fileDeltaJ.is_open()){
                fileDeltaJ << std::setprecision(16) << DeltaJ << "\n"; 
                fileDeltaJ.close();
            }





        }
    

    }

    // Simulations STRPDE DUMMIES  
    if(run_srpde_d){

        for(auto sim = sim_start; sim <= n_sim; ++sim){

            std::cout << "--------------------Simulation RUN STRPDE DUMMIES  #" << std::to_string(sim) << "-------------" << std::endl; 
    
            // data 
            GeoFrame data_srpde(D, T);
            auto& l_srpde = data_srpde.insert_scalar_layer<POINT, POINT>("layer", std::pair{R_path + "space_locs.csv", R_path + "time_locs.csv"});
            l_srpde.load_csv<double>(R_path + "X_dummies.csv");

            // load data from .csv files
            l_srpde.load_csv<double>(R_path + "simulations/sim_" + std::to_string(sim) + "/y_cpp.csv");
                    
            std::string solutions_path_gcv = R_path + "simulations/sim_" + std::to_string(sim) + "/fit_newlib_srpde_d/"; 
            std::string solution_path = R_path + "simulations/sim_" + std::to_string(sim) + "/fit_newlib_srpde_d/"; 
    
            // physics 
            FeSpace Vh(D, P1<1>);   // functional space definition
            Eigen::Matrix<double, 2, 2> K = read_csv<double>(R_path + "simulations/sim_" + std::to_string(sim) + "/K_dummies.csv").as_matrix(); 
            std::cout << "K = " << K << std::endl;
            
            TrialFunction f(Vh);
            TestFunction v(Vh);
            auto a_D = integral(D)(dot(K * grad(f), grad(v)));
            // homogeneous forcing linear form
            ZeroField<2> u_D;
            auto F_D = integral(D)(u_D * v);

            // read lambdas
            double lambda_D;  
            double lambda_T;  
    
            std::ifstream fileLambdaS_gcv(solution_path + "/lambda_s_opt.csv");
            if(fileLambdaS_gcv.is_open()){
                fileLambdaS_gcv >> lambda_D; 
                fileLambdaS_gcv.close();
            }
            std::ifstream fileLambdaT(solution_path + "/lambda_t_opt.csv");
            if(fileLambdaT.is_open()){
                fileLambdaT >> lambda_T; 
                fileLambdaT.close();
            }

            std::cout << "Optimal lambda_D: " << std::setprecision(16) << lambda_D << std::endl;
            std::cout << "Optimal lambda_T: " << std::setprecision(16) << lambda_T << std::endl;

            // Start measuring time
            auto start_time_run = high_resolution_clock::now();
    
            // modeling
            SRPDE m("y ~ x1 + x2 + dummy1 + dummy2 + dummy3 + dummy4 + dummy5 + f", data_srpde, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

            // fit at optimal smoothing level
            m.fit(lambda_D, lambda_T);

            // Stop measuring time
            auto stop_time_run = high_resolution_clock::now();
            auto duration_run = duration_cast<milliseconds>(stop_time_run - start_time_run).count();
            std::cout << "Execution time RUN: " << duration_run << " ms" << std::endl;

            // Save results 
            write_csv(solution_path + "f.csv", m.f());
            write_csv(solution_path + "fn.csv", m.fn());
            write_csv(solution_path + "beta.csv", m.beta());

            std::ofstream file_time_run(solution_path + "/time_run.csv"); 
            if(file_time_run.is_open()){
                file_time_run << duration_run << "\n"; 
                file_time_run.close();
            }



        }
    

    }



    return 0;
}

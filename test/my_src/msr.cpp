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
int test_07(); 

int main(){
    // test_06();
    test_07();
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
                computed_b.row(i) = m.b_hat()[i]; 
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
                computed_b.row(i) = m.b_hat()[i]; 
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

    const std::string schema = "a"; 

    const unsigned int sim_start = 1; 
    const unsigned int n_sim = 30; 

    // run SRPDE and/or MSRPDE ? 
    const bool run_srpde = true;    // stprde
    const bool run_msrpde = false;   // mixed-effects anisotropic
    const bool run_msr_iso = false;  // mixed-effects isotropic
    const bool run_srpde_d = true;   // strpde con dummies

    bool likelihood_dataloss_type; // false = fpirls data loss, true = likelihood
    bool sigma_edf_type;           // false = sigma senza edf nelle iterazioni, true = sigma con edf

    const std::string trial_number = "1"; 
    const std::string R_path = "../../../OneDrive - Politecnico di Milano/Corsi/PhD/Codice/models/MSRPDE/Tests/space-time/Test_7/trial_" + trial_number + "/miss_" + schema + "/";

    if(trial_number == "1"){
        likelihood_dataloss_type = false;
    } 
    if(trial_number == "1"){
        sigma_edf_type = true;
    } 

    unsigned int M; 
    if(trial_number == "1"){
        M = 8; 
    }
    
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 1, M);  // ATT qui non bisogna fare più M-1 come nella vecchia lib!! Vuole direttamente il numero di nodi, cioè M!!! 
    
    
    // geometry 
    using PointT = Eigen::Matrix<double, 2, 1>;

    
    std::string N_string; 
    if(trial_number == "1"){
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
                computed_b.row(i) = m.b_hat()[i]; 
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
                computed_b.row(i) = m.b_hat()[i]; 
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

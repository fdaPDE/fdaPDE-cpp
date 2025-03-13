#include <fdaPDE/core.h>
#include <gtest/gtest.h>   // testing framework

#include <cstddef>
using fdapde::core::advection;
using fdapde::core::diffusion;
using fdapde::core::laplacian;
using fdapde::core::bilaplacian;
using fdapde::core::fem_order;
using fdapde::core::FEM;
using fdapde::core::Grid; 
using fdapde::core::Mesh; 
using fdapde::core::SPLINE;
using fdapde::core::spline_order;
using fdapde::core::PDE;

#include "../../fdaPDE/models/sampling_design.h"
#include "../../fdaPDE/models/regression/strpde.h"
#include "../../fdaPDE/models/regression/srpde.h"
#include "../../fdaPDE/models/regression/qsrpde.h"
using fdapde::models::STRPDE;
using fdapde::models::SRPDE;
using fdapde::models::QSRPDE;
using fdapde::models::SpaceTimeSeparable;
using fdapde::models::SpaceTimeParabolic;

#include "../../fdaPDE/models/regression/gcv.h"
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
using fdapde::models::not_nan; 


// gcv  

TEST(case_study_gcv, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

    std::string pollutant = "NO2";   // PM10 NO2

    const std::string type_locs = "all";   // "reduced" "all"

    const std::string type_space_mesh = "coarse";   // coarse ---> CAMBIA MESH SOTTO !!! 
    
    const std::string type_data = "/data_original";  // "data_original"  ---> MODIFICA ANCHE GIU'lettura dati !!!

    const std::string time_step = "by_week";   // "" "by_day" "by_week"

    std::string gcv_type = "stochastic";   // "exact" "stochastic"  ---> MODIFICA ANCHE GIU'!
    std::size_t seed = 438172;
    unsigned int MC_run = 100; 

    std::string est_type = "quantile";    // mean quantile
    std::vector<double> alphas = {0.99};     //  {0.5, 0.9, 0.95};     

    std::string lambda_selection_type = "gcv_smooth_eps1e-0.5";

    // // Marco 
    // std::string path_data = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant; 
    
    // std::string solutions_path; 
    // if(est_type == "mean")
    //     solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/STRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 
    // if(est_type == "quantile")
    //     solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/QSTRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 

    // Ilenia 
    std::string path_data = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant + "/" + time_step ; 
    std::string solutions_path; 
    if(est_type == "mean")
        solutions_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/STRPDE/" + pollutant + "/" + time_step + "/gcv_Cpp_" + gcv_type; 
    if(est_type == "quantile")
        solutions_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/QSTRPDE/" + pollutant +"/" + time_step + "/gcv_Cpp_" + gcv_type; 


    // lambdas sequence 
    std::vector<DVector<double>> lambdas_d_t;
    std::vector<double> lambdas_d;
    std::vector<double> lambdas_t;

    std::vector<double> lambdas50_d; std::vector<double> lambdas50_t; std::vector<DVector<double>> lambdas50_d_t; 
    std::vector<double> lambdas90_d; std::vector<double> lambdas90_t; std::vector<DVector<double>> lambdas90_d_t; 
    std::vector<double> lambdas95_d; std::vector<double> lambdas95_t; std::vector<DVector<double>> lambdas95_d_t;

    if(est_type == "mean"){
        for(double xs = -4.0; xs <= -1.0; xs += 0.05)
            lambdas_d.push_back(std::pow(10,xs));

        for(double xt = +3.0; xt <= +3.1; xt += 2.0)
            lambdas_t.push_back(std::pow(10,xt));    

        for(auto i = 0; i < lambdas_d.size(); ++i)
            for(auto j = 0; j < lambdas_t.size(); ++j) 
                lambdas_d_t.push_back(SVector<2>(lambdas_d[i], lambdas_t[j]));
    }
    if(est_type == "quantile"){

        // 50% 
        for(double xs = -10.0; xs <= -3.0; xs += 1.0)
            lambdas50_d.push_back(std::pow(10,xs));

        for(double xt = -5.0; xt <= -2.9; xt += 2.0)
            lambdas50_t.push_back(std::pow(10,xt));    

        for(auto i = 0; i < lambdas50_d.size(); ++i)
            for(auto j = 0; j < lambdas50_t.size(); ++j) 
                lambdas50_d_t.push_back(SVector<2>(lambdas50_d[i], lambdas50_t[j]));
        

        for(double xs = -11.0; xs <= +2.0; xs += 1.0)
            lambdas_d.push_back(std::pow(10,xs));

        for(double xt = 0.0; xt <= +1.1; xt += 2.0)
            lambdas_t.push_back(std::pow(10,xt)); 
        // 90%
        for(double xs = -8.0; xs <= -1.0; xs +=1.0)
            for(double xt = -3.0; xt <= +1.1; xt +=2.0)  
            lambdas90_d_t.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));
        // 90%
        for(double xs = -8.0; xs <= -1.0; xs +=1.0)
            for(double xt = -3.0; xt <= +1.1; xt +=2.0)  
            lambdas95_d_t.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));

    }


    // define temporal domain
    unsigned int M = 22;     // DA CAMBIARE !  
    std::string M_string = std::to_string(M);
    double tf;
    if(time_step == "by_week"){
        tf = 357.0;    // final time 
    }else{
        tf = 364.0;    // final time 
    }
        
    
    Mesh<1, 1> time_mesh(0, tf, M-1);
    // std::cout << " # Time mesh elements = " << time_mesh.n_elements()  << std::endl ; 
    // std::cout << "Time mesh elements = " << time_mesh.elements()  << std::endl ; 
    // std::cout << "# Time mesh nodes = " << time_mesh.n_nodes()  << std::endl ; 
    // std::cout << "Time mesh nodes = " << time_mesh.nodes()  << std::endl ; 

    // define spatial domain and regularizing PDE

    MeshLoader<Mesh2D> domain("mesh_lombardia_coarse");        // ATT: da cambiare in base a type_space_mesh ma con if non funziona

    // if(type_space_mesh == "fine"){
    //     MeshLoader<Mesh2D> domain("mesh_lombardia");
    // }else{
    //     MeshLoader<Mesh2D> domain("mesh_lombardia_coarse");
    // }
    
    // import data and locs from files
    DMatrix<double> y; DMatrix<double> space_locs; DMatrix<double> time_locs; 
    if(type_locs == "reduced"){
        if(pollutant == "NO2"){
            y = read_csv<double>(path_data + "/NO2_max_2022_Cpp_reduced.csv");
        } else{
            y = read_csv<double>(path_data + "/PM10_2022_Cpp_reduced.csv");
        }
        time_locs = read_csv<double>(path_data + "/time_locations_reduced.csv");
    } else{
        if(pollutant == "NO2"){
            y = read_csv<double>(path_data + "/NO2_max_2022_Cpp.csv");   // ATT: CAMBIA TIPO DATI NORMALIZED/RESCALED/NO
            time_locs = read_csv<double>(path_data + "/time_locations.csv");
        } else{
            y = read_csv<double>(path_data + "/PM10_2022.csv");  
            time_locs = read_csv<double>(path_data + "/time_locs.csv"); 
        }     
    }
    space_locs = read_csv<double>(path_data + "/locs.csv");
    
    // check dimensions
    std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
    std::cout << "dim time loc " << time_locs.rows() << " " << time_locs.cols() << std::endl;
    std::cout << "dim y " << y.rows() << " " << y.cols() << std::endl;

    BlockFrame<double, int> df;
    df.stack(OBSERVATIONS_BLK, y);
   
    // define regularizing PDE in space 
    auto Ld = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
    PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

    // define regularizing PDE in time
    auto Lt = -bilaplacian<SPLINE>();
    PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

    if(est_type == "mean"){

        STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::pointwise);
    
        // set model's data
        model.set_spatial_locations(space_locs);
        model.set_temporal_locations(time_locs);
        
        model.set_data(df);
        model.init();

        // define GCV function and grid of \lambda_D values

        // stochastic
        auto GCV = model.gcv<StochasticEDF>(MC_run, seed);
        // // exact
        // auto GCV = model.gcv<ExactEDF>();

        // optimize GCV
        Grid<fdapde::Dynamic> opt;
        opt.optimize(GCV, lambdas_d_t);
        SVector<2> best_lambda = opt.optimum();

        // Save lambda sequence 
        std::ofstream fileLambda_S_Seq(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + type_data +"/lambdas_S_seq.csv");
        for(std::size_t i = 0; i < lambdas_d.size(); ++i) 
            fileLambda_S_Seq << std::setprecision(16) << lambdas_d[i] << "\n"; 
        fileLambda_S_Seq.close();


        std::ofstream fileLambda_T_Seq(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + type_data +"/lambdas_T_seq.csv");
        for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
            fileLambda_T_Seq << std::setprecision(16) << lambdas_t[i] << "\n"; 
        fileLambda_T_Seq.close();

        // Save Lambda opt
        std::ofstream fileLambdaoptS(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + type_data +"/lambda_s_opt.csv");
        if(fileLambdaoptS.is_open()){
            fileLambdaoptS << std::setprecision(16) << best_lambda[0];
            fileLambdaoptS.close();
        }
        std::ofstream fileLambdaoptT(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + type_data + "/lambda_t_opt.csv");
        if (fileLambdaoptT.is_open()){
            fileLambdaoptT << std::setprecision(16) << best_lambda[1];
            fileLambdaoptT.close();
        }

        // Save GCV scores
        std::ofstream fileGCV_scores(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + type_data +"/gcv_scores.csv");
        for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
            fileGCV_scores << std::setprecision(16) << GCV.gcvs()[i] << "\n"; 
        fileGCV_scores.close();

        // Save edfs
        std::ofstream fileEDF(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + type_data +"/edfs.csv");
        for(std::size_t i = 0; i < GCV.edfs().size(); ++i) 
            fileEDF << std::setprecision(16) << GCV.edfs()[i] << "\n"; 
        fileEDF.close();

        // Save matrix A
        // SpMatrix<double> computedA = model.A();

        // Eigen::saveMarket(computedA, solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + type_data +"/A");


    }

    if(est_type == "quantile"){
        
        for(double alpha : alphas){

            unsigned int alpha_int = alpha*100; 
            std::string alpha_string = std::to_string(alpha_int); 

            std::cout << "----------------ALPHA = " << alpha_string << "----------------" << std::endl; 

            QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);

            if(alpha_string == "50"){
                lambdas_d = lambdas50_d; 
                lambdas_t = lambdas50_t;
                lambdas_d_t = lambdas50_d_t;
            }
                 
            if(alpha_string == "90"){
                lambdas_d = lambdas90_d; 
                lambdas_t = lambdas90_t;
                lambdas_d_t = lambdas90_d_t;
            }
                 
            if(alpha_string == "95"){
                lambdas_d = lambdas95_d; 
                lambdas_t = lambdas95_t;
                lambdas_d_t = lambdas95_d_t;
            }
    
            // lambdas_d = lambdas_d_t


            // set model's data
            model.set_spatial_locations(space_locs);
            model.set_temporal_locations(time_locs);

            model.set_exact_gcv(lambda_selection_type == "gcv");
            if(lambda_selection_type == "gcv_smooth_eps1e-3")
                model.set_eps_power(-3.0);
            if(lambda_selection_type == "gcv_smooth_eps1e-2")
                model.set_eps_power(-2.0);
            if(lambda_selection_type == "gcv_smooth_eps1e-1.5")
                model.set_eps_power(-1.5);
            if(lambda_selection_type == "gcv_smooth_eps1e-1")
                model.set_eps_power(-1.0);
            if(lambda_selection_type == "gcv_smooth_eps1e+0")
                model.set_eps_power(0.0);
            
            model.set_data(df);
            model.init();

            // define GCV function and grid of \lambda_D values

            // stochastic
            auto GCV = model.gcv<StochasticEDF>(MC_run, seed);
            // // exact
            // auto GCV = model.gcv<ExactEDF>();

            // optimize GCV
            Grid<fdapde::Dynamic> opt;
            opt.optimize(GCV, lambdas_d_t);
            SVector<2> best_lambda = opt.optimum();

            // Save lambda sequences 
            std::ofstream fileLambda_S_Seq(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + "/" + 
                                            lambda_selection_type + "/alpha_" + alpha_string + type_data + "/lambdas_S_seq.csv");
            for(std::size_t i = 0; i < lambdas_d.size(); ++i) 
                fileLambda_S_Seq << std::setprecision(16) << lambdas_d[i] << "\n"; 
            fileLambda_S_Seq.close();
            for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
                std::cout << lambdas_t[i] << "\n"; 
            std::ofstream fileLambda_T_Seq(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + "/" + 
                                            lambda_selection_type + "/alpha_" + alpha_string + type_data + "/lambdas_T_seq.csv");
            for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
                fileLambda_T_Seq << std::setprecision(16) << lambdas_t[i] << "\n"; 
            fileLambda_T_Seq.close();

            // Save Lambda opt
            std::ofstream fileLambdaoptS(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + "/" + 
                                            lambda_selection_type + "/alpha_" + alpha_string + type_data + "/lambda_s_opt.csv");
            if(fileLambdaoptS.is_open()){
                fileLambdaoptS << std::setprecision(16) << best_lambda[0];
                fileLambdaoptS.close();
            }
            std::ofstream fileLambdaoptT(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + "/" + 
                                            lambda_selection_type + "/alpha_" + alpha_string + type_data + "/lambda_t_opt.csv");
            if (fileLambdaoptT.is_open()){
                fileLambdaoptT << std::setprecision(16) << best_lambda[1];
                fileLambdaoptT.close();
            }
            // Save GCV scores
            std::ofstream fileGCV_scores(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + "/" + 
                                            lambda_selection_type + "/alpha_" + alpha_string + type_data + "/gcv_scores.csv");
            for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
                fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.gcvs()[i]) << "\n"; 

            fileGCV_scores.close();

            // Save edfs
            std::ofstream fileEDF(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + "/" + 
                                            lambda_selection_type + "/alpha_" + alpha_string + type_data + "/edfs.csv");
            for(std::size_t i = 0; i < GCV.edfs().size(); ++i) 
                fileEDF << std::setprecision(16) << GCV.edfs()[i] << "\n"; 
            fileEDF.close();

        }

    }

      
}



// // gcv exact

// TEST(case_study_gcv, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

//     std::string pollutant = "NO2";   // PM10 NO2

//     const std::string type_locs = "all";   // "reduced" "all"

//     const std::string type_space_mesh = "coarse";   // coarse ---> CAMBIA MESH SOTTO !!! 
    
//     const std::string type_data = "/data_rescaled";  // "data_original"  ---> MODIFICA ANCHE GIU'lettura dati !!!

//     const std::string time_step = "by_week";   // "" "by_day" "by_week"

//     std::string gcv_type = "exact";   // "exact" "stochastic"  ---> MODIFICA ANCHE GIU'!

//     std::string est_type = "quantile";    // mean quantile
//     std::vector<double> alphas = {0.5};     //  {0.5, 0.9, 0.95};     

//     std::string lambda_selection_type = "gcv_smooth_eps1e+0";

//     // // Marco 
//     // std::string path_data = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant; 
    
//     // std::string solutions_path; 
//     // if(est_type == "mean")
//     //     solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/STRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 
//     // if(est_type == "quantile")
//     //     solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/QSTRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 

//     // Ilenia 
//     std::string path_data = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant + "/" + time_step ; 
//     std::string solutions_path; 
//     if(est_type == "mean")
//         solutions_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/STRPDE/" + pollutant + "/" + time_step + "/gcv_Cpp_" + gcv_type; 
//     if(est_type == "quantile")
//         solutions_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/QSTRPDE/" + pollutant +"/" + time_step + "/gcv_Cpp_" + gcv_type; 


//     // lambdas sequence 
//     std::vector<DVector<double>> lambdas_d_t;
//     std::vector<double> lambdas_d;
//     std::vector<double> lambdas_t;

//     std::vector<DVector<double>> lambdas50_d_t;
//     std::vector<DVector<double>> lambdas90_d_t;
//     std::vector<DVector<double>> lambdas95_d_t;

//     if(est_type == "mean"){
//         for(double xs = -4.0; xs <= -1.0; xs += 0.05)
//             lambdas_d.push_back(std::pow(10,xs));

//         for(double xt = +3.0; xt <= +3.1; xt += 2.0)
//             lambdas_t.push_back(std::pow(10,xt));    

//         for(auto i = 0; i < lambdas_d.size(); ++i)
//             for(auto j = 0; j < lambdas_t.size(); ++j) 
//                 lambdas_d_t.push_back(SVector<2>(lambdas_d[i], lambdas_t[j]));
//     }
//     if(est_type == "quantile"){

//         // 50% 
//         for(double xs = -8.0; xs <= -0.9; xs +=1.0)
//             for(double xt = 0.0; xt <= +1.1; xt += 2.0) 
//             lambdas50_d_t.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));

//         for(double xs = -11.0; xs <= +2.0; xs += 1.0)
//             lambdas_d.push_back(std::pow(10,xs));

//         for(double xt = 0.0; xt <= +1.1; xt += 2.0)
//             lambdas_t.push_back(std::pow(10,xt)); 
//         // 90%
//         for(double xs = -8.0; xs <= -1.0; xs +=1.0)
//             for(double xt = -3.0; xt <= +1.1; xt +=2.0)  
//             lambdas90_d_t.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));
//         // 90%
//         for(double xs = -8.0; xs <= -1.0; xs +=1.0)
//             for(double xt = -3.0; xt <= +1.1; xt +=2.0)  
//             lambdas95_d_t.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));

//     }


//     // define temporal domain
//     unsigned int M = 22;     // DA CAMBIARE !  
//     std::string M_string = std::to_string(M);
//     double tf;
//     if(time_step == "by_week"){
//         tf = 357.0;    // final time 
//     }else{
//         tf = 364.0;    // final time 
//     }
        
    
//     Mesh<1, 1> time_mesh(0, tf, M-1);
//     // std::cout << " # Time mesh elements = " << time_mesh.n_elements()  << std::endl ; 
//     // std::cout << "Time mesh elements = " << time_mesh.elements()  << std::endl ; 
//     // std::cout << "# Time mesh nodes = " << time_mesh.n_nodes()  << std::endl ; 
//     // std::cout << "Time mesh nodes = " << time_mesh.nodes()  << std::endl ; 

//     // define spatial domain and regularizing PDE

//     MeshLoader<Mesh2D> domain("mesh_lombardia_coarse");        // ATT: da cambiare in base a type_space_mesh ma con if non funziona

//     // if(type_space_mesh == "fine"){
//     //     MeshLoader<Mesh2D> domain("mesh_lombardia");
//     // }else{
//     //     MeshLoader<Mesh2D> domain("mesh_lombardia_coarse");
//     // }
    
//     // import data and locs from files
//     DMatrix<double> y; DMatrix<double> space_locs; DMatrix<double> time_locs; 
//     if(type_locs == "reduced"){
//         if(pollutant == "NO2"){
//             y = read_csv<double>(path_data + "/NO2_max_2022_Cpp_reduced.csv");
//         } else{
//             y = read_csv<double>(path_data + "/PM10_2022_Cpp_reduced.csv");
//         }
//         time_locs = read_csv<double>(path_data + "/time_locations_reduced.csv");
//     } else{
//         if(pollutant == "NO2"){
//             y = read_csv<double>(path_data + "/NO2_max_2022_Cpp_rescaled.csv");   // ATT: CAMBIA TIPO DATI NORMALIZED/NO
//             time_locs = read_csv<double>(path_data + "/time_locations.csv");
//         } else{
//             y = read_csv<double>(path_data + "/PM10_2022.csv");  
//             time_locs = read_csv<double>(path_data + "/time_locs.csv"); 
//         }     
//     }
//     space_locs = read_csv<double>(path_data + "/locs.csv");
    
//     // check dimensions
//     std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
//     std::cout << "dim time loc " << time_locs.rows() << " " << time_locs.cols() << std::endl;
//     std::cout << "dim y " << y.rows() << " " << y.cols() << std::endl;

//     BlockFrame<double, int> df;
//     df.stack(OBSERVATIONS_BLK, y);
   
//     // define regularizing PDE in space 
//     auto Ld = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
//     PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

//     // define regularizing PDE in time
//     auto Lt = -bilaplacian<SPLINE>();
//     PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

//     if(est_type == "mean"){

//         STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::pointwise);
    
//         // set model's data
//         model.set_spatial_locations(space_locs);
//         model.set_temporal_locations(time_locs);
        
//         model.set_data(df);
//         model.init();

//         // define GCV function and grid of \lambda_D values

//         // stochastic
//         // auto GCV = model.gcv<StochasticEDF>(MC_run, seed);
//         // // exact
//         auto GCV = model.gcv<ExactEDF>();

//         // optimize GCV
//         Grid<fdapde::Dynamic> opt;
//         opt.optimize(GCV, lambdas_d_t);
//         SVector<2> best_lambda = opt.optimum();

//         // Save lambda sequence 
//         std::ofstream fileLambda_S_Seq(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + type_data +"/lambdas_S_seq.csv");
//         for(std::size_t i = 0; i < lambdas_d.size(); ++i) 
//             fileLambda_S_Seq << std::setprecision(16) << lambdas_d[i] << "\n"; 
//         fileLambda_S_Seq.close();


//         std::ofstream fileLambda_T_Seq(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + type_data +"/lambdas_T_seq.csv");
//         for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
//             fileLambda_T_Seq << std::setprecision(16) << lambdas_t[i] << "\n"; 
//         fileLambda_T_Seq.close();

//         // Save Lambda opt
//         std::ofstream fileLambdaoptS(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + type_data +"/lambda_s_opt.csv");
//         if(fileLambdaoptS.is_open()){
//             fileLambdaoptS << std::setprecision(16) << best_lambda[0];
//             fileLambdaoptS.close();
//         }
//         std::ofstream fileLambdaoptT(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + type_data + "/lambda_t_opt.csv");
//         if (fileLambdaoptT.is_open()){
//             fileLambdaoptT << std::setprecision(16) << best_lambda[1];
//             fileLambdaoptT.close();
//         }

//         // Save GCV scores
//         std::ofstream fileGCV_scores(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + type_data +"/gcv_scores.csv");
//         for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
//             fileGCV_scores << std::setprecision(16) << GCV.gcvs()[i] << "\n"; 
//         fileGCV_scores.close();

//         // Save edfs
//         std::ofstream fileEDF(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + type_data +"/edfs.csv");
//         for(std::size_t i = 0; i < GCV.edfs().size(); ++i) 
//             fileEDF << std::setprecision(16) << GCV.edfs()[i] << "\n"; 
//         fileEDF.close();

//         // Save matrix A
//         // SpMatrix<double> computedA = model.A();

//         // Eigen::saveMarket(computedA, solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + type_data +"/A");


//     }

//     if(est_type == "quantile"){
        
//         for(double alpha : alphas){

//             unsigned int alpha_int = alpha*100; 
//             std::string alpha_string = std::to_string(alpha_int); 

//             std::cout << "----------------ALPHA = " << alpha_string << "----------------" << std::endl; 

//             QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);

//             if(almost_equal(alpha, 0.5)){
//             lambdas_d_t = lambdas50_d_t; 
//             }  
//             if(almost_equal(alpha, 0.9)){
//                 lambdas_d_t = lambdas95_d_t; 
//             }  
//             if(almost_equal(alpha, 0.95)){
//                 lambdas_d_t = lambdas95_d_t; 
//             }
    
//             // lambdas_d = lambdas_d_t


//             // set model's data
//             model.set_spatial_locations(space_locs);
//             model.set_temporal_locations(time_locs);

//             model.set_exact_gcv(lambda_selection_type == "gcv");
//             if(lambda_selection_type == "gcv_smooth_eps1e-3")
//                 model.set_eps_power(-3.0);
//             if(lambda_selection_type == "gcv_smooth_eps1e-2")
//                 model.set_eps_power(-2.0);
//             if(lambda_selection_type == "gcv_smooth_eps1e-1.5")
//                 model.set_eps_power(-1.5);
//             if(lambda_selection_type == "gcv_smooth_eps1e-1")
//                 model.set_eps_power(-1.0);
//             if(lambda_selection_type == "gcv_smooth_eps1e+0")
//                 model.set_eps_power(0.0);
            
//             model.set_data(df);
//             model.init();

//             // define GCV function and grid of \lambda_D values

//             // stochastic
//             auto GCV = model.gcv<ExactEDF>();
//             // // exact
//             // auto GCV = model.gcv<ExactEDF>();

//             std::cout << solutions_path << std::endl ; 

//             // optimize GCV
//             Grid<fdapde::Dynamic> opt;
//             opt.optimize(GCV, lambdas_d_t);
//             SVector<2> best_lambda = opt.optimum();

//             // Save lambda sequence 
//             std::ofstream fileLambda_S_Seq(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + "/" + 
//                                             lambda_selection_type + "/alpha_" + alpha_string + type_data + "/lambdas_S_seq.csv");
//             for(std::size_t i = 0; i < lambdas_d.size(); ++i) 
//                 fileLambda_S_Seq << std::setprecision(16) << lambdas_d[i] << "\n"; 
//             fileLambda_S_Seq.close();

//             for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
//                 std::cout << lambdas_t[i] << "\n"; 

//             std::ofstream fileLambda_T_Seq(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + "/" + 
//                                             lambda_selection_type + "/alpha_" + alpha_string + type_data + "/lambdas_T_seq.csv");
//             for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
//                 fileLambda_T_Seq << std::setprecision(16) << lambdas_t[i] << "\n"; 
//             fileLambda_T_Seq.close();

//             // Save Lambda opt
//             std::ofstream fileLambdaoptS(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + "/" + 
//                                             lambda_selection_type + "/alpha_" + alpha_string + type_data + "/lambda_s_opt.csv");
//             if(fileLambdaoptS.is_open()){
//                 fileLambdaoptS << std::setprecision(16) << best_lambda[0];
//                 fileLambdaoptS.close();
//             }
//             std::ofstream fileLambdaoptT(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + "/" + 
//                                             lambda_selection_type + "/alpha_" + alpha_string + type_data + "/lambda_t_opt.csv");
//             if (fileLambdaoptT.is_open()){
//                 fileLambdaoptT << std::setprecision(16) << best_lambda[1];
//                 fileLambdaoptT.close();
//             }
//             // Save GCV scores
//             std::ofstream fileGCV_scores(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + "/" + 
//                                             lambda_selection_type + "/alpha_" + alpha_string + type_data + "/gcv_scores.csv");
//             for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
//                 fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.gcvs()[i]) << "\n"; 

//             fileGCV_scores.close();

//             // Save edfs
//             std::ofstream fileEDF(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + "/" + 
//                                             lambda_selection_type + "/alpha_" + alpha_string + type_data + "/edfs.csv");
//             for(std::size_t i = 0; i < GCV.edfs().size(); ++i) 
//                 fileEDF << std::setprecision(16) << GCV.edfs()[i] << "\n"; 
//             fileEDF.close();

//         }

//     }

      
// }







// // run  

// TEST(case_study_run, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

//     std::string pollutant = "NO2";   // PM10 NO2

//     std::string time_step = "by_week";

//     std::string gcv_type = "stochastic";   // "exact" "stochastic"  

//     std::string est_type = "quantile";    // mean quantile
//     std::vector<double> alphas = {0.5};     // , 0.9, 0.95};     

//     std::string lambda_selection_type = "gcv_smooth_eps1e+0";  // _smooth_eps1e-1

//     // // Marco 
//     // std::string path_data = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant; 
//     // std::string solutions_path; 
//     // if(est_type == "mean")
//     //     solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/STRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 
//     // if(est_type == "quantile")
//     //     solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/QSTRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 

//     // Ilenia 
//     std::string path_data = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant + "/" + time_step ; 
//     std::string solutions_path; 
//     if(est_type == "mean")
//         solutions_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/STRPDE/" + pollutant + "/" + time_step + "/gcv_Cpp_" + gcv_type; 
//     if(est_type == "quantile")
//         solutions_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/QSTRPDE/" + pollutant + "/" + time_step + "/gcv_Cpp_" + gcv_type; 


//     const std::string type_locs = "all";   // "reduced" "all"
//     const std::string type_data = "/data_original";     // CAMBIA SOTTO I DATI
//     const std::string type_space_mesh = "coarse";

//     // define temporal domain
//     unsigned int M = 22;         // DA CAMBIARE 
//     std::string M_string = std::to_string(M);
//     double tf;
//     if(time_step == "by_week"){
//         tf = 357.0;    // final time 
//     }else{
//         tf = 364.0;    // final time 
//     }
        
//     Mesh<1, 1> time_mesh(0, tf, M-1);
//     // define spatial domain and regularizing PDE
//     MeshLoader<Mesh2D> domain("mesh_lombardia_coarse");

//     // import data and locs from files
//     DMatrix<double> y; DMatrix<double> space_locs; DMatrix<double> time_locs; 
//     if(type_locs == "reduced"){
//         if(pollutant == "NO2"){
//             y = read_csv<double>(path_data + "/NO2_max_2022_Cpp_reduced.csv");
//         } else{
//             y = read_csv<double>(path_data + "/PM10_2022_Cpp_reduced.csv");
//         }
//         time_locs = read_csv<double>(path_data + "/time_locations_reduced.csv");
//     } else{
//         if(pollutant == "NO2"){
//             y = read_csv<double>(path_data + "/NO2_max_2022_Cp.csv");       // ATT: CAMBIA TIPO DATI NORMALIZED/NO
//             time_locs = read_csv<double>(path_data + "/time_locations.csv");
//         } else{
//             y = read_csv<double>(path_data + "/PM10_2022.csv");  
//             time_locs = read_csv<double>(path_data + "/time_locs.csv"); 
//         }     
//     }
//     space_locs = read_csv<double>(path_data + "/locs.csv");

//     // check dimensions
//     std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
//     std::cout << "dim time loc " << time_locs.rows() << " " << time_locs.cols() << std::endl;
//     std::cout << "dim y " << y.rows() << " " << y.cols() << std::endl;

//     BlockFrame<double, int> df;
//     df.stack(OBSERVATIONS_BLK, y);
   
//     // define regularizing PDE in space 
//     auto Ld = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
//     PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

//     // define regularizing PDE in time
//     auto Lt = -bilaplacian<SPLINE>();
//     PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);


//     if(est_type == "mean"){

//         STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::pointwise);
    
//         // set model's data
//         model.set_spatial_locations(space_locs);
//         model.set_temporal_locations(time_locs);

//         // // Impose lambdas 
//         // double lambda_T = 1e-3; 
//         // double lambda_S = 1e-3; 

//         // Read optima lambdas 
//         double lambda_T ; 
//         double lambda_S ; 
        
//         std::ifstream fileLambdaT_opt(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + type_data + "/lambda_t_opt.csv");
//         if(fileLambdaT_opt.is_open()){
//             fileLambdaT_opt >> lambda_T; 
//             fileLambdaT_opt.close();
//         }
//         std::ifstream fileLambdaS_opt(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + type_data + "/lambda_s_opt.csv");
//         if(fileLambdaS_opt.is_open()){
//             fileLambdaS_opt >> lambda_S; 
//             fileLambdaS_opt.close();
//         }

//         std::cout << "lambda S " << lambda_S << std::endl;
//         std::cout << "lambda T " << lambda_T << std::endl;

//         model.set_lambda_D(lambda_S);
//         model.set_lambda_T(lambda_T);
        
//         model.set_data(df);

//         model.init();
//         model.solve();

//         // Save C++ solution 
//         DMatrix<double> computedF = model.f();
//         const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filef(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + type_data + "/f.csv");
//         if (filef.is_open()){
//             filef << computedF.format(CSVFormatf);
//             filef.close();
//         }

//         DMatrix<double> computedFn = model.Psi()*model.f();
//         const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filefn(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + type_data + "/fn.csv");
//         if (filefn.is_open()){
//             filefn << computedFn.format(CSVFormatfn);
//             filefn.close();
//         }

 

//     }

//     if(est_type == "quantile"){
        
//         for(double alpha : alphas){

//             unsigned int alpha_int = alpha*100; 
//             std::string alpha_string = std::to_string(alpha_int); 

            

//             // Read optima lambdas 
//             double lambda_T; 
//             double lambda_S; 

//             std::ifstream fileLambdaT_opt(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + "/" + 
//                                             lambda_selection_type + "/alpha_" + alpha_string + type_data + "/lambda_t_opt.csv");
//             if(fileLambdaT_opt.is_open()){
//                 fileLambdaT_opt >> lambda_T; 
//                 fileLambdaT_opt.close();
//             }
//             std::ifstream fileLambdaS_opt(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + "/" + 
//                                             lambda_selection_type + "/alpha_" + alpha_string + type_data + "/lambda_s_opt.csv");
//             if(fileLambdaS_opt.is_open()){
//                 fileLambdaS_opt >> lambda_S; 
//                 fileLambdaS_opt.close();
//             }

//             // // Impose lambdas 
//             // double lambda_T = 1e+0; 
//             // double lambda_S = 1e-6;

//             // std::vector<double> seq_lambda_S;
//             // for(double xs = -5.0; xs <= +2.1; xs += 7.0)
//             //     seq_lambda_S.push_back(std::pow(10,xs));

//             // std::size_t lambda_count = 1;

//             // for(double lambda_S_choice : seq_lambda_S ){

//             //     lambda_S = lambda_S_choice;
//                 std::cout << "lambda S " << lambda_S << std::endl;
//                 std::cout << "lambda T " << lambda_T << std::endl;

//                 QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);
    
//                 // set model's data
//                 model.set_spatial_locations(space_locs);
//                 model.set_temporal_locations(time_locs);

//                 model.set_lambda_D(lambda_S);
//                 model.set_lambda_T(lambda_T);
                
//                 model.set_data(df);

//                 model.init();
//                 model.solve();

//                 // Save C++ solution 
//                 DMatrix<double> computedF = model.f();
//                 const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//                 std::ofstream filef(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + "/" + 
//                                             lambda_selection_type + "/alpha_" + alpha_string + type_data + "/f.csv");  //  + std::to_string(lambda_count) + ".csv");

//                 if (filef.is_open()){
//                     filef << computedF.format(CSVFormatf);
//                     filef.close();
//                 }

//                 DMatrix<double> computedFn = model.Psi()*model.f();
//                 const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//                 std::ofstream filefn(solutions_path + "/" + type_locs + "/N_" + type_space_mesh + "/M_" + M_string + "/" + 
//                                             lambda_selection_type + "/alpha_" + alpha_string + type_data + "/fn.csv");  // _" + std::to_string(lambda_count) + ".csv");
//                 if (filefn.is_open()){
//                     filefn << computedFn.format(CSVFormatfn);
//                     filefn.close();
//                 }

//                 // lambda_count = lambda_count +1 ; 
//             // }



//         }

//     }

        
// }



// // run  for condition number of A

// TEST(case_study_run, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

//     std::string pollutant = "NO2";   // PM10 NO2

//     std::string time_step = "by_week";

//     std::string gcv_type = "stochastic";   // "exact" "stochastic"  

//     std::string est_type = "mean";    // mean quantile
//     std::vector<double> alphas = {0.5};     // {0.5, 0.9, 0.95};     

//     std::string lambda_selection_type = "gcv_smooth_eps1e-1";

//     // // Marco 
//     // std::string path_data = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant; 
//     // std::string solutions_path; 
//     // if(est_type == "mean")
//     //     solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/STRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 
//     // if(est_type == "quantile")
//     //     solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/QSTRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 

//     // Ilenia 
//     std::string path_data = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant + "/" + time_step ; 
//     std::string solutions_path; 
//     if(est_type == "mean")
//         solutions_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/STRPDE/" + pollutant + "/" + time_step + "/gcv_Cpp_" + gcv_type; 
//     if(est_type == "quantile")
//         solutions_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/QSTRPDE/" + pollutant + "/" + time_step + "/gcv_Cpp_" + gcv_type; 


//     const std::string type_locs = "all";   // "reduced" "all"
//     const std::string type_data = "/data_normalized";
//     const std::string type_space_mesh = "coarse";

//     // define temporal domain
//     unsigned int M = 22;         // DA CAMBIARE 
//     std::string M_string = std::to_string(M);
//     double tf;
//     if(time_step == "by_week"){
//         tf = 357.0;    // final time 
//     }else{
//         tf = 364.0;    // final time 
//     }
        
//     Mesh<1, 1> time_mesh(0, tf, M-1);
//     // define spatial domain and regularizing PDE
//     MeshLoader<Mesh2D> domain("mesh_lombardia_coarse");

//     // import data and locs from files
//     DMatrix<double> y; DMatrix<double> space_locs; DMatrix<double> time_locs; 
//     if(type_locs == "reduced"){
//         if(pollutant == "NO2"){
//             y = read_csv<double>(path_data + "/NO2_max_2022_Cpp_reduced.csv");
//         } else{
//             y = read_csv<double>(path_data + "/PM10_2022_Cpp_reduced.csv");
//         }
//         time_locs = read_csv<double>(path_data + "/time_locations_reduced.csv");
//     } else{
//         if(pollutant == "NO2"){
//             y = read_csv<double>(path_data + "/NO2_max_2022_Cpp_normalized.csv");       // ATT: CAMBIA TIPO DATI NORMALIZED/NO
//             time_locs = read_csv<double>(path_data + "/time_locations.csv");
//         } else{
//             y = read_csv<double>(path_data + "/PM10_2022.csv");  
//             time_locs = read_csv<double>(path_data + "/time_locs.csv"); 
//         }     
//     }
//     space_locs = read_csv<double>(path_data + "/locs.csv");

//     // check dimensions
//     std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
//     std::cout << "dim time loc " << time_locs.rows() << " " << time_locs.cols() << std::endl;
//     std::cout << "dim y " << y.rows() << " " << y.cols() << std::endl;

//     BlockFrame<double, int> df;
//     df.stack(OBSERVATIONS_BLK, y);
   
//     // define regularizing PDE in space 
//     auto Ld = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
//     PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

//     // define regularizing PDE in time
//     auto Lt = -bilaplacian<SPLINE>();
//     PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);


//     if(est_type == "mean"){

//         STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::pointwise);
    
//         // set model's data
//         model.set_spatial_locations(space_locs);
//         model.set_temporal_locations(time_locs);

//         // Impose lambdas 
//         double lambda_T = 1e+0; 
//         double lambda_S ;

//         std::vector<double> seq_lambda_S;
//         for(double xs = -11.0; xs <= +2.0; xs += 1.0)
//             seq_lambda_S.push_back(std::pow(10,xs));

//         std::size_t count_lambda = 1;

//         for(double lambda_S_choice : seq_lambda_S ){

//             lambda_S = lambda_S_choice;

//             std::cout << "lambda S " << lambda_S << std::endl;
//             std::cout << "lambda T " << lambda_T << std::endl;

//             model.set_lambda_D(lambda_S);
//             model.set_lambda_T(lambda_T);
            
//             model.set_data(df);

//             model.init();
//             model.solve();

//             // Save matrix A 
//             SpMatrix<double> computedA = model.A();

//             Eigen::saveMarket(computedA, solutions_path + "/" + type_locs + "/N_" + type_space_mesh + 
//                                 "/M_" + M_string + type_data + "/A_" + std::to_string(count_lambda) + ".mtx");

//             count_lambda = count_lambda +1;
        
//         }

        

 

//     }
     
// }


// ---------- PM10 restricted 
// // gcv time

// TEST(case_study_gcv, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

//     const std::string eps_string = "1e+0";   //  "0" "1e-1" "1e-0.5" "1e+0" "1e+0.5" "1e+1" "1e+2"

//     std::string gcv_type = "stochastic";   // "exact" "stochastic"  ---> MODIFICA ANCHE GIU'!
//     std::size_t seed = 438172;
//     unsigned int MC_run = 100; 
//     const std::string num_months = "one_month";   // "one_month" "two_months"
//     const std::string scaling = "_rescale";   // "" "_rescale"
//     const std::string model_type = "nonparametric";  // "nonparametric" "parametric"
//     const std::string mesh_type = "square";  // "square" "esagoni"

//     std::string est_type = "quantile";    // mean quantile
//     std::vector<double> alphas = {0.9}; 

//     // Marco 
//     std::string path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia_restricted"; 
//     std::string path_data = path + "/data/" + num_months; 
//     std::string solutions_path; 

//     // Ilenia 
//     // std::string path_data = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia_restricted/data"; 
//     // std::string solutions_path; 

//     if(est_type == "mean")
//         solutions_path = path + "/results/STRPDE/" + num_months + "/" + model_type + "/" + mesh_type; 
//     if(est_type == "quantile")
//         solutions_path = path + "/results/QSTRPDE/" + num_months + "/" + model_type + "/" + mesh_type; 


//     // lambdas sequence 
//     std::vector<DVector<double>> lambdas_d_t; std::vector<double> lambdas_d; std::vector<double> lambdas_t;
//     std::vector<double> lambdas50_d; std::vector<double> lambdas50_t; std::vector<DVector<double>> lambdas50_d_t; 
//     std::vector<double> lambdas90_d; std::vector<double> lambdas90_t; std::vector<DVector<double>> lambdas90_d_t; 
//     std::vector<double> lambdas95_d; std::vector<double> lambdas95_t; std::vector<DVector<double>> lambdas95_d_t; 

//     if(est_type == "mean"){
//         for(double xs = -3.5; xs <= 0.0; xs += 0.1)
//             lambdas_d.push_back(std::pow(10,xs));

//         for(double xt = -7.0; xt <= -1.0; xt += 2.0)
//             lambdas_t.push_back(std::pow(10,xt));    

//         for(auto i = 0; i < lambdas_d.size(); ++i)
//             for(auto j = 0; j < lambdas_t.size(); ++j) 
//                 lambdas_d_t.push_back(SVector<2>(lambdas_d[i], lambdas_t[j]));
//     }
//     if(est_type == "quantile"){

//         // 50% 
//         {
//             for(double xs = -4.4; xs <= -3.2; xs += 0.2)
//                 lambdas50_d.push_back(std::pow(10,xs));

//             for(double xt = -3.0; xt <= -2.0; xt += 2.0)
//                 lambdas50_t.push_back(std::pow(10,xt));    

//             for(auto i = 0; i < lambdas50_d.size(); ++i)
//                 for(auto j = 0; j < lambdas50_t.size(); ++j) 
//                     lambdas50_d_t.push_back(SVector<2>(lambdas50_d[i], lambdas50_t[j]));
//         }

//         // 90% 
//         {
//             for(double xs = -8.0; xs <= 0.0; xs += 1.0)
//                 lambdas90_d.push_back(std::pow(10,xs));

//             for(double xt = -3.0; xt <= -2.0; xt += 2.0)
//                 lambdas90_t.push_back(std::pow(10,xt));    
  
//             for(auto i = 0; i < lambdas90_d.size(); ++i)
//                 for(auto j = 0; j < lambdas90_t.size(); ++j) 
//                     lambdas90_d_t.push_back(SVector<2>(lambdas90_d[i], lambdas90_t[j]));
//         }

//         // 95% 
//         {
//             for(double xs = -8.0; xs <= -1.0; xs += 1.0)
//                 lambdas95_d.push_back(std::pow(10,xs));

//             for(double xt = -5.0; xt <= 1.0; xt += 2.0)
//                 lambdas95_t.push_back(std::pow(10,xt));    

//             for(auto i = 0; i < lambdas95_d.size(); ++i)
//                 for(auto j = 0; j < lambdas95_t.size(); ++j) 
//                     lambdas95_d_t.push_back(SVector<2>(lambdas95_d[i], lambdas95_t[j]));
//         }

//     }

//     // define temporal domain
//     unsigned int M; double tf; 
//     if(num_months == "two_months"){
//         M = 30; 
//         tf = 60.0; 
//     }
//     if(num_months == "one_month"){
//         M = 11;    // old: 22
//         tf = 21.0; 
//     }
//     std::string M_string = std::to_string(M);
//     Mesh<1, 1> time_mesh(0, tf, M-1);
//     // define spatial domain and regularizing PDE
//     MeshLoader<Mesh2D> domain("mesh_lombardia_" + mesh_type);

//     // import data and locs from files
//     DMatrix<double> y; DMatrix<double> space_locs; DMatrix<double> time_locs; DMatrix<double> X;  

//     y = read_csv<double>(path_data + "/y" + scaling + ".csv");  
//     time_locs = read_csv<double>(path_data + "/time_locs.csv"); 
//     space_locs = read_csv<double>(path_data + "/locs.csv");
//     if(model_type == "parametric")
//         X = read_csv<double>(path_data + "/X.csv");
    
//     // check dimensions
//     std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
//     std::cout << "dim time loc " << time_locs.rows() << " " << time_locs.cols() << std::endl;
//     std::cout << "dim y " << y.rows() << " " << y.cols() << std::endl;
//     std::cout << "dim X " << X.rows() << " " << X.cols() << std::endl;

//     BlockFrame<double, int> df;
//     df.stack(OBSERVATIONS_BLK, y);
//     if(model_type == "parametric")
//         df.stack(DESIGN_MATRIX_BLK, X);
   
//     // define regularizing PDE in space 
//     auto Ld = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
//     PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

//     // define regularizing PDE in time
//     auto Lt = -bilaplacian<SPLINE>();
//     PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

//     std::cout << "-----------------------------GCV STARTS------------------------" << std::endl; 

//     if(est_type == "mean"){

//         STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::pointwise);
        
//         // set model's data
//         model.set_spatial_locations(space_locs);
//         model.set_temporal_locations(time_locs);
        
//         model.set_data(df);
//         model.init();

//         // define GCV function and grid of \lambda_D values

//         // stochastic
//         auto GCV = model.gcv<StochasticEDF>(MC_run, seed);
//         // // exact
//         // auto GCV = model.gcv<ExactEDF>();

//         // optimize GCV
//         Grid<fdapde::Dynamic> opt;
//         opt.optimize(GCV, lambdas_d_t);
//         SVector<2> best_lambda = opt.optimum();

//         // Save lambda sequence 
//         std::ofstream fileLambda_S_Seq(solutions_path + "/lambdas_S_seq.csv");
//         for(std::size_t i = 0; i < lambdas_d.size(); ++i) 
//             fileLambda_S_Seq << std::setprecision(16) << lambdas_d[i] << "\n"; 
//         fileLambda_S_Seq.close();

//         std::ofstream fileLambda_T_Seq(solutions_path + "/lambdas_T_seq.csv");
//         for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
//             fileLambda_T_Seq << std::setprecision(16) << lambdas_t[i] << "\n"; 
//         fileLambda_T_Seq.close();

//         // Save Lambda opt
//         std::ofstream fileLambdaoptS(solutions_path + "/lambda_s_opt.csv");
//         if(fileLambdaoptS.is_open()){
//             fileLambdaoptS << std::setprecision(16) << best_lambda[0];
//             fileLambdaoptS.close();
//         }
//         std::ofstream fileLambdaoptT(solutions_path + "/lambda_t_opt.csv");
//         if (fileLambdaoptT.is_open()){
//             fileLambdaoptT << std::setprecision(16) << best_lambda[1];
//             fileLambdaoptT.close();
//         }
//         // Save GCV scores
//         std::ofstream fileGCV_scores(solutions_path + "/gcv_scores.csv");
//         for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
//             fileGCV_scores << std::setprecision(16) << GCV.gcvs()[i] << "\n"; 
//         fileGCV_scores.close();

//         // Save edfs
//         std::ofstream fileEDF(solutions_path + "/edfs.csv");
//         for(std::size_t i = 0; i < GCV.edfs().size(); ++i) 
//             fileEDF << std::setprecision(16) << GCV.edfs()[i] << "\n"; 
//         fileEDF.close();

//     }

//     if(est_type == "quantile"){
        
//         for(double alpha : alphas){

//             unsigned int alpha_int = alpha*100; 
//             std::string alpha_string = std::to_string(alpha_int); 

//             std::cout << "----------------ALPHA = " << alpha_string << "----------------" << std::endl; 

//             QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);
    
//             // set model's data
//             model.set_spatial_locations(space_locs);
//             model.set_temporal_locations(time_locs);

//             model.set_exact_gcv(eps_string == "0");   // true: no lisciamento curva gcv
//             if(eps_string == "1e+2"){
//                 model.set_eps_power(2.0); 
//             }
//             if(eps_string == "1e+1"){
//                 model.set_eps_power(1.0); 
//             }
//             if(eps_string == "1e+0.75"){
//                 model.set_eps_power(0.75); 
//             }
//             if(eps_string == "1e+0.5"){
//                 model.set_eps_power(0.5); 
//             }
//             if(eps_string == "1e+0"){
//                 model.set_eps_power(0.0); 
//             }
//             if(eps_string == "1e-0.5"){
//                 model.set_eps_power(-0.5); 
//             }
//             if(eps_string == "1e-1"){
//                 model.set_eps_power(-1.0); 
//             }
//             if(eps_string == "1e-1.5"){
//                 model.set_eps_power(-1.5); 
//             }
//             if(eps_string == "1e-2"){
//                 model.set_eps_power(-2.0); 
//             }
//             if(eps_string == "1e-3"){
//                 model.set_eps_power(-3.0); 
//             }
            
//             model.set_data(df);
//             model.init();

//             // define GCV function and grid of \lambda_D values

//             // stochastic
//             auto GCV = model.gcv<StochasticEDF>(MC_run, seed);
//             // // exact
//             // auto GCV = model.gcv<ExactEDF>();

//             if(alpha_string == "50"){
//                 lambdas_d = lambdas50_d; 
//                 lambdas_t = lambdas50_t;
//                 lambdas_d_t = lambdas50_d_t;
//             }
                 
//             if(alpha_string == "90"){
//                 lambdas_d = lambdas90_d; 
//                 lambdas_t = lambdas90_t;
//                 lambdas_d_t = lambdas90_d_t;
//             }
                 
//             if(alpha_string == "95"){
//                 lambdas_d = lambdas95_d; 
//                 lambdas_t = lambdas95_t;
//                 lambdas_d_t = lambdas95_d_t;
//             }
//             // optimize GCV
//             Grid<fdapde::Dynamic> opt;
//             opt.optimize(GCV, lambdas_d_t);
//             SVector<2> best_lambda = opt.optimum();

//             // Save lambda sequences 
//             std::ofstream fileLambda_S_Seq(solutions_path + "/alpha_" + alpha_string  + "/eps_" + eps_string + "/lambdas_S_seq.csv");
//             for(std::size_t i = 0; i < lambdas_d.size(); ++i) 
//                 fileLambda_S_Seq << std::setprecision(16) << lambdas_d[i] << "\n"; 
//             fileLambda_S_Seq.close();
//             for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
//                 std::cout << lambdas_t[i] << "\n"; 
//             std::ofstream fileLambda_T_Seq(solutions_path + "/alpha_" + alpha_string  + "/eps_" + eps_string + "/lambdas_T_seq.csv");
//             for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
//                 fileLambda_T_Seq << std::setprecision(16) << lambdas_t[i] << "\n"; 
//             fileLambda_T_Seq.close();

//             // Save Lambda opt
//             std::ofstream fileLambdaoptS(solutions_path + "/alpha_" + alpha_string  + "/eps_" + eps_string + "/lambda_s_opt.csv");
//             if(fileLambdaoptS.is_open()){
//                 fileLambdaoptS << std::setprecision(16) << best_lambda[0];
//                 fileLambdaoptS.close();
//             }
//             std::ofstream fileLambdaoptT(solutions_path + "/alpha_" + alpha_string  + "/eps_" + eps_string + "/lambda_t_opt.csv");
//             if (fileLambdaoptT.is_open()){
//                 fileLambdaoptT << std::setprecision(16) << best_lambda[1];
//                 fileLambdaoptT.close();
//             }
//             // Save GCV scores
//             std::ofstream fileGCV_scores(solutions_path + "/alpha_" + alpha_string  + "/eps_" + eps_string + "/gcv_scores.csv");
//             for(std::size_t i = 0; i < GCV.gcvs().size(); ++i){
//                 fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.gcvs()[i]) << "\n"; 
//                 std::cout << "----- " << std::sqrt(GCV.gcvs()[i]) << std::endl ; 
//             }
                
//             fileGCV_scores.close();
//             // Save edfs
//             std::ofstream fileEDF(solutions_path + "/alpha_" + alpha_string  + "/eps_" + eps_string + "/edfs.csv");
//             for(std::size_t i = 0; i < GCV.edfs().size(); ++i) 
//                 fileEDF << std::setprecision(16) << GCV.edfs()[i] << "\n"; 
//             fileEDF.close();  
//             // Save gcv numerator
//             std::ofstream fileGCVnum(solutions_path + "/alpha_" + alpha_string  + "/eps_" + eps_string + "/gcv_num.csv");
//             for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
//                 fileGCVnum << std::setprecision(16) << std::sqrt(GCV.gcvs()[i]*model.n_obs()) << "\n"; 
//             fileGCVnum.close();

//         }

//     }

      
// }

// // run time

// TEST(case_study_run, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

//     const std::string eps_string = "1e+0";   // "1e-1" "0" "1e+0" "1e+1"

//     std::string est_type = "mean";    // mean quantile
//     std::vector<double> alphas = {0.5}; 
//     const std::string num_months = "one_month";   // "one_month" "two_months"
//     const std::string model_type = "nonparametric";  // "nonparametric" "parametric"
//     const std::string scaling = "_rescale";   // "" "_rescale"
//     const std::string mesh_type = "square";  // "square" "esagoni"

//     // Marco 
//     std::string path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia_restricted"; 
//     std::string path_data = path + "/data/" + num_months; 
//     std::string solutions_path; 

//     // Ilenia 
//     // std::string path_data = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia_restricted/data"; 
//     // std::string solutions_path; 

//     if(est_type == "mean")
//         solutions_path = path + "/results/STRPDE/" + num_months + "/" + model_type + "/" + mesh_type;
//     if(est_type == "quantile")
//         solutions_path = path + "/results/QSTRPDE/" + num_months + "/" + model_type + "/" + mesh_type; 

//     // define temporal domain
//     unsigned int M; double tf; 
//     if(num_months == "two_months"){
//         M = 30; 
//         tf = 60.0; 
//     }
//     if(num_months == "one_month"){
//         M = 11;   // old: 22
//         tf = 21.0; 
//     }
//     std::string M_string = std::to_string(M);
//     Mesh<1, 1> time_mesh(0, tf, M-1);
//     // define spatial domain and regularizing PDE
//     MeshLoader<Mesh2D> domain("mesh_lombardia_" + mesh_type);

//     // import data and locs from files
//     DMatrix<double> y; DMatrix<double> space_locs; DMatrix<double> time_locs; DMatrix<double> X; 

//     y = read_csv<double>(path_data + "/y" + scaling + ".csv"); 
//     time_locs = read_csv<double>(path_data + "/time_locs.csv"); 
//     space_locs = read_csv<double>(path_data + "/locs.csv");
//     if(model_type == "parametric")
//         X = read_csv<double>(path_data + "/X.csv");
    
//     // check dimensions
//     std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
//     std::cout << "dim time loc " << time_locs.rows() << " " << time_locs.cols() << std::endl;
//     std::cout << "dim y " << y.rows() << " " << y.cols() << std::endl;
//     std::cout << "dim X " << X.rows() << " " << X.cols() << std::endl;

//     BlockFrame<double, int> df;
//     df.stack(OBSERVATIONS_BLK, y);
//     if(model_type == "parametric")
//         df.stack(DESIGN_MATRIX_BLK, X);
   
//     // define regularizing PDE in space 
//     auto Ld = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
//     PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

//     // define regularizing PDE in time
//     auto Lt = -bilaplacian<SPLINE>();
//     PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

//     std::cout << "--------------------------------RUN STARTS--------------------------------" << std::endl; 
//     if(est_type == "mean"){

//         STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::pointwise);
    
//         // set model's data
//         model.set_spatial_locations(space_locs);
//         model.set_temporal_locations(time_locs);

//         // Read optima lambdas 
//         double lambda_T; double lambda_S; 
//         std::ifstream fileLambdaT_opt(solutions_path + "/lambda_t_opt.csv");
//         if(fileLambdaT_opt.is_open()){
//             fileLambdaT_opt >> lambda_T; 
//             fileLambdaT_opt.close();
//         }
//         std::ifstream fileLambdaS_opt(solutions_path + "/lambda_s_opt.csv");
//         if(fileLambdaS_opt.is_open()){
//             fileLambdaS_opt >> lambda_S; 
//             fileLambdaS_opt.close();
//         }

//         std::cout << "lambda S " << lambda_S << std::endl;
//         std::cout << "lambda T " << lambda_T << std::endl;

//         model.set_lambda_D(lambda_S);
//         model.set_lambda_T(lambda_T);
        
//         model.set_data(df);

//         model.init();
//         model.solve();

//         // Save C++ solution 
//         DMatrix<double> computedF = model.f();
//         const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filef(solutions_path + "/f.csv");
//         if (filef.is_open()){
//             filef << computedF.format(CSVFormatf);
//             filef.close();
//         }

//         DMatrix<double> computedFn = model.Psi()*model.f();
//         const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filefn(solutions_path + "/fn.csv");
//         if (filefn.is_open()){
//             filefn << computedFn.format(CSVFormatfn);
//             filefn.close();
//         }

//         if(model_type == "parametric"){
//             DMatrix<double> computedBeta = model.beta();
//             const static Eigen::IOFormat CSVFormatBeta(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//             std::ofstream filebeta(solutions_path + "/beta.csv");
//             if (filebeta.is_open()){
//                 filebeta << computedBeta.format(CSVFormatBeta);
//                 filebeta.close();
//             }
//         }

//     }

//     if(est_type == "quantile"){
        
//         for(double alpha : alphas){

//             unsigned int alpha_int = alpha*100; 
//             std::string alpha_string = std::to_string(alpha_int); 

//             std::cout << "----------------ALPHA = " << alpha_string << "----------------" << std::endl; 

//             QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);
    
//             // set model's data
//             model.set_spatial_locations(space_locs);
//             model.set_temporal_locations(time_locs);

//             // Read optima lambdas 
//             double lambda_T; double lambda_S; 
//             std::ifstream fileLambdaT_opt(solutions_path + "/alpha_" + alpha_string + "/eps_" + eps_string + "/lambda_t_opt.csv");
//             if(fileLambdaT_opt.is_open()){
//                 fileLambdaT_opt >> lambda_T; 
//                 fileLambdaT_opt.close();
//             }
//             std::ifstream fileLambdaS_opt(solutions_path + "/alpha_" + alpha_string + "/eps_" + eps_string + "/lambda_s_opt.csv");
//             if(fileLambdaS_opt.is_open()){
//                 fileLambdaS_opt >> lambda_S; 
//                 fileLambdaS_opt.close();
//             }

//             // double lambda_T=1e-5; 
//             // double lambda_S=1e-8; 

//             std::cout << "lambda S " << lambda_S << std::endl;
//             std::cout << "lambda T " << lambda_T << std::endl;

//             model.set_lambda_D(lambda_S);
//             model.set_lambda_T(lambda_T);
            
//             model.set_data(df);

//             model.init();
//             model.solve();

//             // Save C++ solution 
//             DMatrix<double> computedF = model.f();
//             const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//             std::ofstream filef(solutions_path + "/alpha_" + alpha_string + "/eps_" + eps_string + "/f.csv");
//             if (filef.is_open()){
//                 filef << computedF.format(CSVFormatf);
//                 filef.close();
//             }

//             DMatrix<double> computedFn = model.Psi()*model.f();
//             const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//             std::ofstream filefn(solutions_path + "/alpha_" + alpha_string + "/eps_" + eps_string + "/fn.csv");
//             if (filefn.is_open()){
//                 filefn << computedFn.format(CSVFormatfn);
//                 filefn.close();
//             }

//             if(model_type == "parametric"){
//                 DMatrix<double> computedBeta = model.beta();
//                 const static Eigen::IOFormat CSVFormatBeta(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//                 std::ofstream filebeta(solutions_path + "/alpha_" + alpha_string + "/eps_" + eps_string + "/beta.csv");
//                 if (filebeta.is_open()){
//                     filebeta << computedBeta.format(CSVFormatBeta);
//                     filebeta.close();
//                 }
//             }


//         }

//     }        
// }


// ---------- NO2 restricted 

// gcv time
TEST(case_study_gcv, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

    const std::string eps_string = "1e-0.5";   // "1e-0.25" "0"  "1e+0" "1e+0.5" "1e+1" "1e+2"

    std::string gcv_type = "stochastic";   // "exact" "stochastic"  ---> MODIFICA ANCHE GIU'!
    std::size_t seed = 438172;
    unsigned int MC_run = 100; 
    const std::string num_months = "one_month";   // "one_month" "two_months"   
    const std::string model_type = "nonparametric";  // "nonparametric" "parametric"
    const std::string covariate_type = "radiation";
    const std::string mesh_type = "convex_hull";  // "square" "esagoni" "convex_hull"
    const std::string pollutant = "NO2";
    const bool return_smoothing = false;    // se true => Smetti exact gcv!! 

    std::string est_type = "quantile";    // mean quantile
    std::vector<double> alphas = {0.05}; //, 0.9, 0.95}; 

    // Marco 
    std::string path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia_restricted"; 
    std::string path_data = path + "/data/" + num_months + "/" + pollutant;  
    std::string solutions_path; 

    // // Ilenia 
    // std::string path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia_restricted"; 
    // std::string path_data = path + "/data/" + num_months + "/" + pollutant; 
    // std::string solutions_path;

    if(est_type == "mean"){
        if(model_type == "nonparametric"){
            solutions_path = path + "/results/STRPDE/" + num_months + "/" + model_type + "/" + mesh_type + "/" + pollutant;
        } else{
            solutions_path = path + "/results/STRPDE/" + num_months + "/" + model_type + "/" + mesh_type + "/" + pollutant + "/" + covariate_type;
        }
    }

    if(est_type == "quantile"){
        if(model_type == "nonparametric"){
            solutions_path = path + "/results/QSTRPDE/" + num_months + "/" + model_type + "/" + mesh_type + "/" + pollutant;
        } else{
            solutions_path = path + "/results/QSTRPDE/" + num_months + "/" + model_type + "/" + mesh_type + "/" + pollutant + "/" + covariate_type;
        }
    }

    // lambdas sequence 
    std::vector<DVector<double>> lambdas_d_t; std::vector<double> lambdas_d; std::vector<double> lambdas_t;
     std::vector<double> lambdas5_d; std::vector<double> lambdas5_t; std::vector<DVector<double>> lambdas5_d_t; 
    std::vector<double> lambdas50_d; std::vector<double> lambdas50_t; std::vector<DVector<double>> lambdas50_d_t; 
    std::vector<double> lambdas90_d; std::vector<double> lambdas90_t; std::vector<DVector<double>> lambdas90_d_t; 
    std::vector<double> lambdas95_d; std::vector<double> lambdas95_t; std::vector<DVector<double>> lambdas95_d_t; 

    if(est_type == "mean"){
        if(!return_smoothing){
            for(double xs = -3.0; xs <= -1.0; xs += 0.1)
                lambdas_d.push_back(std::pow(10,xs));
            for(double xt = -1.0; xt <= 0; xt += 2.0)
                lambdas_t.push_back(std::pow(10,xt));    

        } else{
            double lambda_S; double lambda_T;  
            std::ifstream fileLambdaS_opt(solutions_path + "/lambda_s_opt.csv");
            if(fileLambdaS_opt.is_open()){
                fileLambdaS_opt >> lambda_S; 
                fileLambdaS_opt.close();
            }
            lambdas_d.push_back(std::pow(10,lambda_S)); 

            std::ifstream fileLambdaT_opt(solutions_path + "/lambda_t_opt.csv");
            if(fileLambdaT_opt.is_open()){
                fileLambdaT_opt >> lambda_T; 
                fileLambdaT_opt.close();
            }
            lambdas_t.push_back(std::pow(10,lambda_T)); 
        }


        for(auto i = 0; i < lambdas_d.size(); ++i)
            for(auto j = 0; j < lambdas_t.size(); ++j) 
                lambdas_d_t.push_back(SVector<2>(lambdas_d[i], lambdas_t[j]));
    }

    if(return_smoothing && lambdas_d.size() > 1){
        std::cout << "ERROR: you want S, but you are providing more lambdas" << std::endl; 
    } 

    if(est_type == "quantile"){

        // 5% 
        {
            for(double xs = -9.5; xs <= -6.0; xs += 0.1)
                lambdas5_d.push_back(std::pow(10,xs));

            for(double xt = -3.0; xt <= -2.0; xt += 2.0)
                lambdas5_t.push_back(std::pow(10,xt));    

            for(auto i = 0; i < lambdas5_d.size(); ++i)
                for(auto j = 0; j < lambdas5_t.size(); ++j) 
                    lambdas5_d_t.push_back(SVector<2>(lambdas5_d[i], lambdas5_t[j]));
        }

        // 50% 
        {
            for(double xs = -7.5; xs <= -5.0; xs += 0.1)
                lambdas50_d.push_back(std::pow(10,xs));

            for(double xt = -3.0; xt <= -0.9; xt += 2.0)
                lambdas50_t.push_back(std::pow(10,xt));    

            for(auto i = 0; i < lambdas50_d.size(); ++i)
                for(auto j = 0; j < lambdas50_t.size(); ++j) 
                    lambdas50_d_t.push_back(SVector<2>(lambdas50_d[i], lambdas50_t[j]));
        }

        // 90% 
        {
            for(double xs = -9.0; xs <= -5.0; xs += 0.05)
                lambdas90_d.push_back(std::pow(10,xs));

            for(double xt = -3.0; xt <= -2.0; xt += 2.0)
                lambdas90_t.push_back(std::pow(10,xt));    

            for(auto i = 0; i < lambdas90_d.size(); ++i)
                for(auto j = 0; j < lambdas90_t.size(); ++j) 
                    lambdas90_d_t.push_back(SVector<2>(lambdas90_d[i], lambdas90_t[j]));
        }

        // 95% 
        {
            for(double xs = -9.5; xs <= -6.0; xs += 0.1)
                lambdas95_d.push_back(std::pow(10,xs));

            for(double xt = -3.0; xt <= -2.0; xt += 2.0)
                lambdas95_t.push_back(std::pow(10,xt));    

            for(auto i = 0; i < lambdas95_d.size(); ++i)
                for(auto j = 0; j < lambdas95_t.size(); ++j) 
                    lambdas95_d_t.push_back(SVector<2>(lambdas95_d[i], lambdas95_t[j]));
        }

    }

    // define temporal domain
    unsigned int M; double tf; 
    if(num_months == "two_months"){
        M = 30; 
        tf = 61.0; 
    }
    if(num_months == "one_month"){
        M = 11; 
        tf = 22.0; 
    }
    std::string M_string = std::to_string(M);
    Mesh<1, 1> time_mesh(0, tf, M-1);
    // define spatial domain and regularizing PDE
    MeshLoader<Mesh2D> domain("mesh_lombardia_" + mesh_type);

    // import data and locs from files
    DMatrix<double> y; DMatrix<double> space_locs; DMatrix<double> time_locs; DMatrix<double> X;  

    y = read_csv<double>(path_data + "/y_rescale.csv");     // ATT
    time_locs = read_csv<double>(path_data + "/time_locs.csv"); 
    space_locs = read_csv<double>(path_data + "/locs.csv");
    if(model_type == "parametric")
        X = read_csv<double>(path_data + "/X_" + covariate_type + ".csv");
    
    // check dimensions
    std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
    std::cout << "dim time loc " << time_locs.rows() << " " << time_locs.cols() << std::endl;
    std::cout << "dim y " << y.rows() << " " << y.cols() << std::endl;
    std::cout << "dim X " << X.rows() << " " << X.cols() << std::endl;

    BlockFrame<double, int> df;
    df.stack(OBSERVATIONS_BLK, y);
    if(model_type == "parametric")
        df.insert(DESIGN_MATRIX_BLK, X);
   
    // define regularizing PDE in space 
    auto Ld = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
    PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

    // define regularizing PDE in time
    auto Lt = -bilaplacian<SPLINE>();
    PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

    std::cout << "-----------------------------GCV STARTS------------------------" << std::endl; 

    if(est_type == "mean"){

        STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::pointwise);
        
        // set model's data
        model.set_spatial_locations(space_locs);
        model.set_temporal_locations(time_locs);
        
        model.set_data(df);
        model.init();

        // define GCV function and grid of \lambda_D values

        // stochastic
        auto GCV = model.gcv<StochasticEDF>(MC_run, seed);
        // // exact
        // auto GCV = model.gcv<ExactEDF>();

        // optimize GCV
        Grid<fdapde::Dynamic> opt;
        opt.optimize(GCV, lambdas_d_t);
        SVector<2> best_lambda = opt.optimum();


        if(!return_smoothing){
            // Save lambda sequence 
            std::ofstream fileLambda_S_Seq(solutions_path + "/lambdas_S_seq.csv");
            for(std::size_t i = 0; i < lambdas_d.size(); ++i) 
                fileLambda_S_Seq << std::setprecision(16) << lambdas_d[i] << "\n"; 
            fileLambda_S_Seq.close();

            std::ofstream fileLambda_T_Seq(solutions_path + "/lambdas_T_seq.csv");
            for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
                fileLambda_T_Seq << std::setprecision(16) << lambdas_t[i] << "\n"; 
            fileLambda_T_Seq.close();

            // Save Lambda opt
            std::ofstream fileLambdaoptS(solutions_path + "/lambda_s_opt.csv");
            if(fileLambdaoptS.is_open()){
                fileLambdaoptS << std::setprecision(16) << best_lambda[0];
                fileLambdaoptS.close();
            }
            std::ofstream fileLambdaoptT(solutions_path + "/lambda_t_opt.csv");
            if (fileLambdaoptT.is_open()){
                fileLambdaoptT << std::setprecision(16) << best_lambda[1];
                fileLambdaoptT.close();
            }
            // Save GCV scores
            std::ofstream fileGCV_scores(solutions_path + "/gcv_scores.csv");
            for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
                fileGCV_scores << std::setprecision(16) << GCV.gcvs()[i] << "\n"; 
            fileGCV_scores.close();

            // Save edfs
            std::ofstream fileEDF(solutions_path + "/edfs.csv");
            for(std::size_t i = 0; i < GCV.edfs().size(); ++i) 
                fileEDF << std::setprecision(16) << GCV.edfs()[i] << "\n"; 
            fileEDF.close();
        }

        if(return_smoothing){
            // Save S
            DMatrix<double> computedS = GCV.S_get_gcv();
            Eigen::saveMarket(computedS, solutions_path + "/S.mtx");
        }

    }

    if(est_type == "quantile"){
        
        for(double alpha : alphas){

            unsigned int alpha_int = alpha*100; 
            std::string alpha_string = std::to_string(alpha_int); 

            std::cout << "----------------ALPHA = " << alpha_string << "----------------" << std::endl; 

            QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);
    
            // set model's data
            model.set_spatial_locations(space_locs);
            model.set_temporal_locations(time_locs);

            if(eps_string == "1e+2"){
                model.set_eps_power(2.0); 
            }
            if(eps_string == "1e+1"){
                model.set_eps_power(1.0); 
            }
            if(eps_string == "1e+0.75"){
                model.set_eps_power(0.75); 
            }
            if(eps_string == "1e+0.5"){
                model.set_eps_power(0.5); 
            }
            if(eps_string == "1e+0"){
                model.set_eps_power(0.0); 
            }
            if(eps_string == "1e-1"){
                model.set_eps_power(-1.0); 
            }
            if(eps_string == "1e-0.5"){
                model.set_eps_power(-0.5); 
            }
            if(eps_string == "1e-0.25"){
                model.set_eps_power(-0.25); 
            }
            if(eps_string == "1e-1.5"){
                model.set_eps_power(-1.5); 
            }
            if(eps_string == "1e-2"){
                model.set_eps_power(-2.0); 
            }
            if(eps_string == "1e-3"){
                model.set_eps_power(-3.0); 
            }
            
            model.set_data(df);
            model.init();

            // define GCV function and grid of \lambda_D values

            // stochastic
            auto GCV = model.gcv<StochasticEDF>(MC_run, seed);
            // // exact
            // auto GCV = model.gcv<ExactEDF>();

                             
            if(alpha_string == "5"){
                lambdas_d = lambdas5_d; 
                lambdas_t = lambdas5_t;
                lambdas_d_t = lambdas5_d_t;
            }

            if(alpha_string == "50"){
                lambdas_d = lambdas50_d; 
                lambdas_t = lambdas50_t;
                lambdas_d_t = lambdas50_d_t;
            }
                 
            if(alpha_string == "90"){
                lambdas_d = lambdas90_d; 
                lambdas_t = lambdas90_t;
                lambdas_d_t = lambdas90_d_t;
            }
                 
            if(alpha_string == "95"){
                lambdas_d = lambdas95_d; 
                lambdas_t = lambdas95_t;
                lambdas_d_t = lambdas95_d_t;
            }
            // optimize GCV
            Grid<fdapde::Dynamic> opt;
            opt.optimize(GCV, lambdas_d_t);
            SVector<2> best_lambda = opt.optimum();

            // Save lambda sequences 
            std::ofstream fileLambda_S_Seq(solutions_path + "/alpha_" + alpha_string  + "/eps_" + eps_string + "/lambdas_S_seq.csv");
            for(std::size_t i = 0; i < lambdas_d.size(); ++i) 
                fileLambda_S_Seq << std::setprecision(16) << lambdas_d[i] << "\n"; 
            fileLambda_S_Seq.close();
            for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
                std::cout << lambdas_t[i] << "\n"; 
            std::ofstream fileLambda_T_Seq(solutions_path + "/alpha_" + alpha_string  + "/eps_" + eps_string + "/lambdas_T_seq.csv");
            for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
                fileLambda_T_Seq << std::setprecision(16) << lambdas_t[i] << "\n"; 
            fileLambda_T_Seq.close();

            // Save Lambda opt
            std::ofstream fileLambdaoptS(solutions_path + "/alpha_" + alpha_string  + "/eps_" + eps_string + "/lambda_s_opt.csv");
            if(fileLambdaoptS.is_open()){
                fileLambdaoptS << std::setprecision(16) << best_lambda[0];
                fileLambdaoptS.close();
            }
            std::ofstream fileLambdaoptT(solutions_path + "/alpha_" + alpha_string  + "/eps_" + eps_string + "/lambda_t_opt.csv");
            if (fileLambdaoptT.is_open()){
                fileLambdaoptT << std::setprecision(16) << best_lambda[1];
                fileLambdaoptT.close();
            }

            // Save GCV scores
            std::ofstream fileGCV_scores(solutions_path + "/alpha_" + alpha_string  + "/eps_" + eps_string + "/gcv_scores.csv");
            for(std::size_t i = 0; i < GCV.gcvs().size(); ++i){
                fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.gcvs()[i]) << "\n"; 
                std::cout << "----- " << std::sqrt(GCV.gcvs()[i]) << std::endl ; 
            }
                
            fileGCV_scores.close();

            // Save edfs
            std::ofstream fileEDF(solutions_path + "/alpha_" + alpha_string  + "/eps_" + eps_string + "/edfs.csv");
            for(std::size_t i = 0; i < GCV.edfs().size(); ++i) 
                fileEDF << std::setprecision(16) << GCV.edfs()[i] << "\n"; 
            fileEDF.close();  
            // Save gcv numerator
            std::ofstream fileGCVnum(solutions_path + "/alpha_" + alpha_string  + "/eps_" + eps_string + "/gcv_num.csv");
            for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
                fileGCVnum << std::setprecision(16) << std::sqrt(GCV.gcvs()[i]*model.n_obs()) << "\n"; 
            fileGCVnum.close();

        }

    }

      
}


// run time
TEST(case_study_run, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

    const std::string eps_string = "1e-0.5";   // "0" "1e+0" "1e+1"

    std::string est_type = "quantile";    // mean quantile
    std::vector<double> alphas = {0.05}; // {0.5, 0.90, 0.95}; 
    const std::string num_months = "one_month";   // "one_month" "two_months"
    const std::string model_type = "nonparametric";  // "nonparametric" "parametric"
    const std::string covariate_type = "radiation";
    const std::string mesh_type = "convex_hull";  // "square" "esagoni" "convex_hull"
    const std::string pollutant = "NO2"; 

    // Marco 
    std::string path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia_restricted"; 
    std::string path_data = path + "/data/" + num_months + "/" + pollutant; 
    std::string solutions_path; 

    // // Ilenia 
    // std::string path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia_restricted"; 
    // std::string path_data = path + "/data/" + num_months + "/" + pollutant; 
    // std::string solutions_path;

    if(est_type == "mean"){
        if(model_type == "nonparametric"){
            solutions_path = path + "/results/STRPDE/" + num_months + "/" + model_type + "/" + mesh_type + "/" + pollutant;
        } else{
            solutions_path = path + "/results/STRPDE/" + num_months + "/" + model_type + "/" + mesh_type + "/" + pollutant + "/" + covariate_type;
        }
    }
        
    if(est_type == "quantile"){
        if(model_type == "nonparametric"){
            solutions_path = path + "/results/QSTRPDE/" + num_months + "/" + model_type + "/" + mesh_type + "/" + pollutant;
        } else{
            solutions_path = path + "/results/QSTRPDE/" + num_months + "/" + model_type + "/" + mesh_type + "/" + pollutant + "/" + covariate_type;
        }
    }

    // define temporal domain
    unsigned int M; double tf; 
    if(num_months == "two_months"){
        M = 30; 
        tf = 61.0; 
    }
    if(num_months == "one_month"){
        M = 11; 
        tf = 22.0; 
    }
    std::string M_string = std::to_string(M);
    Mesh<1, 1> time_mesh(0, tf, M-1);
    // define spatial domain and regularizing PDE
    MeshLoader<Mesh2D> domain("mesh_lombardia_" + mesh_type);

    // import data and locs from files
    DMatrix<double> y; DMatrix<double> space_locs; DMatrix<double> time_locs; DMatrix<double> X; 

    y = read_csv<double>(path_data + "/y_rescale.csv");     // ATT
    time_locs = read_csv<double>(path_data + "/time_locs.csv"); 
    space_locs = read_csv<double>(path_data + "/locs.csv");
    if(model_type == "parametric")
        X = read_csv<double>(path_data + "/X_" + covariate_type + ".csv");
    
    // check dimensions
    std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
    std::cout << "dim time loc " << time_locs.rows() << " " << time_locs.cols() << std::endl;
    std::cout << "dim y " << y.rows() << " " << y.cols() << std::endl;
    std::cout << "dim X " << X.rows() << " " << X.cols() << std::endl;

    BlockFrame<double, int> df;
    df.stack(OBSERVATIONS_BLK, y);
    if(model_type == "parametric")
        df.insert(DESIGN_MATRIX_BLK, X);
   
    // define regularizing PDE in space 
    auto Ld = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
    PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

    // define regularizing PDE in time
    auto Lt = -bilaplacian<SPLINE>();
    PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

    std::cout << "--------------------------------RUN STARTS--------------------------------" << std::endl; 
    if(est_type == "mean"){

        STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::pointwise);
    
        // set model's data
        model.set_spatial_locations(space_locs);
        model.set_temporal_locations(time_locs);

        // Read optima lambdas 
        double lambda_T; double lambda_S; 
        std::ifstream fileLambdaT_opt(solutions_path + "/lambda_t_opt.csv");
        if(fileLambdaT_opt.is_open()){
            fileLambdaT_opt >> lambda_T; 
            fileLambdaT_opt.close();
        }
        std::ifstream fileLambdaS_opt(solutions_path + "/lambda_s_opt.csv");
        if(fileLambdaS_opt.is_open()){
            fileLambdaS_opt >> lambda_S; 
            fileLambdaS_opt.close();
        }

        std::cout << "lambda S " << lambda_S << std::endl;
        std::cout << "lambda T " << lambda_T << std::endl;

        model.set_lambda_D(lambda_S);
        model.set_lambda_T(lambda_T);
        
        model.set_data(df);

        model.init();
        model.solve();

        // Save C++ solution 
        DMatrix<double> computedF = model.f();
        const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        std::ofstream filef(solutions_path + "/f.csv");
        if (filef.is_open()){
            filef << computedF.format(CSVFormatf);
            filef.close();
        }

        DMatrix<double> computedFn = model.Psi(not_nan())*model.f();
        const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        std::ofstream filefn(solutions_path + "/fn.csv");
        if (filefn.is_open()){
            filefn << computedFn.format(CSVFormatfn);
            filefn.close();
        }

        if(model_type == "parametric"){
            DMatrix<double> computedBeta = model.beta();
            const static Eigen::IOFormat CSVFormatBeta(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
            std::ofstream filebeta(solutions_path + "/beta.csv");
            if (filebeta.is_open()){
                filebeta << computedBeta.format(CSVFormatBeta);
                filebeta.close();
            }
        }

    }

    if(est_type == "quantile"){
        
        for(double alpha : alphas){

            unsigned int alpha_int = alpha*100; 
            std::string alpha_string = std::to_string(alpha_int); 

            std::cout << "----------------ALPHA = " << alpha_string << "----------------" << std::endl; 

            QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);
    
            // set model's data
            model.set_spatial_locations(space_locs);
            model.set_temporal_locations(time_locs);

            // Read optima lambdas 
            double lambda_T; double lambda_S; 
            std::ifstream fileLambdaT_opt(solutions_path + "/alpha_" + alpha_string + "/eps_" + eps_string + "/lambda_t_opt.csv");
            if(fileLambdaT_opt.is_open()){
                fileLambdaT_opt >> lambda_T; 
                fileLambdaT_opt.close();
            }
            std::ifstream fileLambdaS_opt(solutions_path + "/alpha_" + alpha_string + "/eps_" + eps_string + "/lambda_s_opt.csv");
            if(fileLambdaS_opt.is_open()){
                fileLambdaS_opt >> lambda_S; 
                fileLambdaS_opt.close();
            }

            std::cout << "lambda S " << lambda_S << std::endl;
            std::cout << "lambda T " << lambda_T << std::endl;

            model.set_lambda_D(lambda_S);
            model.set_lambda_T(lambda_T);
            
            model.set_data(df);

            model.init();
            model.solve();

            // Save C++ solution 
            DMatrix<double> computedF = model.f();
            const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
            std::ofstream filef(solutions_path + "/alpha_" + alpha_string + "/eps_" + eps_string + "/f.csv");
            if (filef.is_open()){
                filef << computedF.format(CSVFormatf);
                filef.close();
            }

            DMatrix<double> computedFn = model.Psi(not_nan())*model.f();
            const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
            std::ofstream filefn(solutions_path + "/alpha_" + alpha_string + "/eps_" + eps_string + "/fn.csv");
            if (filefn.is_open()){
                filefn << computedFn.format(CSVFormatfn);
                filefn.close();
            }

            if(model_type == "parametric"){
                DMatrix<double> computedBeta = model.beta();
                const static Eigen::IOFormat CSVFormatBeta(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
                std::ofstream filebeta(solutions_path + "/alpha_" + alpha_string + "/eps_" + eps_string + "/beta.csv");
                if (filebeta.is_open()){
                    filebeta << computedBeta.format(CSVFormatBeta);
                    filebeta.close();
                }
            }


        }

    }        
}


// ---------- NO2 mapping piogge 

// // gcv time

// TEST(case_study_gcv_piogge, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

//     std::size_t seed = 438172;
//     unsigned int MC_run = 100; 
//     const std::string mesh_type = "esagoni";  // "square" "esagoni"
//     const std::string pollutant = "NO2";

//     // Marco 
//     std::string path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia_restricted"; 
//     std::string path_data = path + "/raw_data/covariate/location_mapping";  
//     std::string solutions_path = path_data; 

//     // // Ilenia 
//     // std::string path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia_restricted"; 
//     // std::string path_data = path + "/data/" + num_months + "/" + pollutant; 
//     // std::string solutions_path;

//     // lambdas sequence 
//     std::vector<DVector<double>> lambdas_d_t; std::vector<double> lambdas_d; std::vector<double> lambdas_t;
//     for(double xs = -2.0; xs <= +0.5; xs += 0.05)
//         lambdas_d.push_back(std::pow(10,xs));

//     for(double xt = -1.0; xt <= 0.; xt += 2.0)
//         lambdas_t.push_back(std::pow(10,xt));    

//     for(auto i = 0; i < lambdas_d.size(); ++i)
//         for(auto j = 0; j < lambdas_t.size(); ++j) 
//             lambdas_d_t.push_back(SVector<2>(lambdas_d[i], lambdas_t[j]));
    

//     // define temporal domain
//     unsigned int M; double tf; 
//     M = 11; 
//     tf = 22.0; 

//     std::string M_string = std::to_string(M);
//     Mesh<1, 1> time_mesh(0, tf, M-1);
//     // define spatial domain and regularizing PDE
//     MeshLoader<Mesh2D> domain("mesh_lombardia_esagoni");

//     // import data and locs from files
//     DMatrix<double> y; DMatrix<double> space_locs; DMatrix<double> time_locs;  

//     y = read_csv<double>(path_data + "/original_y.csv");   
//     time_locs = read_csv<double>(path_data + "/time_locs.csv"); 
//     space_locs = read_csv<double>(path_data + "/original_locs.csv");
    
//     // check dimensions
//     std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
//     std::cout << "dim time loc " << time_locs.rows() << " " << time_locs.cols() << std::endl;
//     std::cout << "dim y " << y.rows() << " " << y.cols() << std::endl;

//     BlockFrame<double, int> df;
//     df.stack(OBSERVATIONS_BLK, y);
   
//     // define regularizing PDE in space 
//     auto Ld = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
//     PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

//     // define regularizing PDE in time
//     auto Lt = -bilaplacian<SPLINE>();
//     PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

//     std::cout << "-----------------------------GCV STARTS------------------------" << std::endl; 

//     STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::pointwise);
    
//     // set model's data
//     model.set_spatial_locations(space_locs);
//     model.set_temporal_locations(time_locs);
    
//     model.set_data(df);
//     model.init();

//     // define GCV function and grid of \lambda_D values

//     // stochastic
//     auto GCV = model.gcv<StochasticEDF>(MC_run, seed);
//     // // exact
//     // auto GCV = model.gcv<ExactEDF>();

//     // optimize GCV
//     Grid<fdapde::Dynamic> opt;
//     opt.optimize(GCV, lambdas_d_t);
//     SVector<2> best_lambda = opt.optimum();

//     // Save lambda sequence 
//     std::ofstream fileLambda_S_Seq(solutions_path + "/lambdas_S_seq.csv");
//     for(std::size_t i = 0; i < lambdas_d.size(); ++i) 
//         fileLambda_S_Seq << std::setprecision(16) << lambdas_d[i] << "\n"; 
//     fileLambda_S_Seq.close();

//     std::ofstream fileLambda_T_Seq(solutions_path + "/lambdas_T_seq.csv");
//     for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
//         fileLambda_T_Seq << std::setprecision(16) << lambdas_t[i] << "\n"; 
//     fileLambda_T_Seq.close();

//     // Save Lambda opt
//     std::ofstream fileLambdaoptS(solutions_path + "/lambda_s_opt.csv");
//     if(fileLambdaoptS.is_open()){
//         fileLambdaoptS << std::setprecision(16) << best_lambda[0];
//         fileLambdaoptS.close();
//     }
//     std::ofstream fileLambdaoptT(solutions_path + "/lambda_t_opt.csv");
//     if (fileLambdaoptT.is_open()){
//         fileLambdaoptT << std::setprecision(16) << best_lambda[1];
//         fileLambdaoptT.close();
//     }
//     // Save GCV scores
//     std::ofstream fileGCV_scores(solutions_path + "/gcv_scores.csv");
//     for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
//         fileGCV_scores << std::setprecision(16) << GCV.gcvs()[i] << "\n"; 
//     fileGCV_scores.close();

//     // Save edfs
//     std::ofstream fileEDF(solutions_path + "/edfs.csv");
//     for(std::size_t i = 0; i < GCV.edfs().size(); ++i) 
//         fileEDF << std::setprecision(16) << GCV.edfs()[i] << "\n"; 
//     fileEDF.close();


// }

// // run time

// TEST(case_study_run_piogge, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

//     const std::string mesh_type = "esagoni";  // "square" "esagoni"
//     const std::string pollutant = "NO2"; 

//     // Marco 
//     std::string path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia_restricted"; 
//     std::string path_data = path + "/raw_data/covariate/location_mapping";  
//     std::string solutions_path = path_data; 

//     // // Ilenia 
//     // std::string path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia_restricted"; 
//     // std::string path_data = path + "/data/" + num_months + "/" + pollutant; 
//     // std::string solutions_path;

//     // define temporal domain
//     unsigned int M; double tf; 
//     M = 11; 
//     tf = 22.0; 

//     std::string M_string = std::to_string(M);
//     Mesh<1, 1> time_mesh(0, tf, M-1);
//     // define spatial domain and regularizing PDE
//     MeshLoader<Mesh2D> domain("mesh_lombardia_" + mesh_type);

//     // import data and locs from files
//     DMatrix<double> y; DMatrix<double> space_locs; DMatrix<double> time_locs;

//     y = read_csv<double>(path_data + "/original_y.csv");   
//     time_locs = read_csv<double>(path_data + "/time_locs.csv"); 
//     space_locs = read_csv<double>(path_data + "/original_locs.csv");
    
//     // check dimensions
//     std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
//     std::cout << "dim time loc " << time_locs.rows() << " " << time_locs.cols() << std::endl;

//     BlockFrame<double, int> df;
//     df.stack(OBSERVATIONS_BLK, y);
   
//     // define regularizing PDE in space 
//     auto Ld = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
//     PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

//     // define regularizing PDE in time
//     auto Lt = -bilaplacian<SPLINE>();
//     PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

//     std::cout << "--------------------------------RUN STARTS--------------------------------" << std::endl; 

//     STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::pointwise);

//     // set model's data
//     model.set_spatial_locations(space_locs);
//     model.set_temporal_locations(time_locs);

//     // Read optima lambdas 
//     double lambda_T; double lambda_S; 
//     std::ifstream fileLambdaT_opt(solutions_path + "/lambda_t_opt.csv");
//     if(fileLambdaT_opt.is_open()){
//         fileLambdaT_opt >> lambda_T; 
//         fileLambdaT_opt.close();
//     }
//     std::ifstream fileLambdaS_opt(solutions_path + "/lambda_s_opt.csv");
//     if(fileLambdaS_opt.is_open()){
//         fileLambdaS_opt >> lambda_S; 
//         fileLambdaS_opt.close();
//     }

//     std::cout << "lambda S " << lambda_S << std::endl;
//     std::cout << "lambda T " << lambda_T << std::endl;

//     model.set_lambda_D(lambda_S);
//     model.set_lambda_T(lambda_T);
    
//     model.set_data(df);

//     model.init();
//     model.solve();

//     // Save C++ solution 
//     DMatrix<double> computedF = model.f();
//     const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//     std::ofstream filef(solutions_path + "/f.csv");
//     if (filef.is_open()){
//         filef << computedF.format(CSVFormatf);
//         filef.close();
//     }

//     DMatrix<double> computedFn = model.Psi()*model.f();
//     const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//     std::ofstream filefn(solutions_path + "/fn.csv");
//     if (filefn.is_open()){
//         filefn << computedFn.format(CSVFormatfn);
//         filefn.close();
//     }

// }








// -------- Lake Victoria
// // Lake victoria -> anche lui crasha per dimensione
// TEST(strpde_case_study, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

//     // Marco
//     const std::string path_data = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/models/space_time/Arnone_Vicini/dati_LakeVictoria"; 

//     // Ilenia 
//     // std::string path_data = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/models/space_time/Arnone_Vicini/dati_LakeVictoria";
 
//     // define temporal domain
//     std::cout << "read " << std::endl;
//     DMatrix<double> time_points  = read_csv<double>(path_data + "/time_mesh.csv");
//     std::cout << "ok " << std::endl;
//     std::cout << time_points.rows() << " , " << time_points.cols() <<  std::endl;
//     std::cout << " value=" << time_points.coeff(0,0) <<  std::endl; 

//     // define temporal domain
//     unsigned int M = 80; 
//     std::string M_string = std::to_string(M);
//     double tf = 6119;    // final time 
//     Mesh<1, 1> time_mesh(0, tf, M-1);

//     // unsigned int M = 2;
//     // time_mesh.resize(M);
//     // time_mesh[0] = 0.00;
//     // time_mesh[0] = 182.00;
//     // time_mesh[0] = 364.00;


//     // define spatial domain and regularizing PDE
//     MeshLoader<Mesh2D> domain("lake_victoria");
//     // import data from files
//     DMatrix<double> time_locs  = read_csv<double>(path_data + "/time_locations.csv");
//     DMatrix<double> space_locs = read_csv<double>(path_data + "/locs.csv");
//     DMatrix<double> y          = read_csv<double>(path_data + "/y.csv");
//     std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
//     std::cout << "dim time loc " << time_locs.rows() << " " << time_locs.cols() << std::endl;
//     std::cout << "dim y " << y.rows() << " " << y.cols() << std::endl;

//     // define regularizing PDE in space 
//     auto Ld = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
//     PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

//     // define regularizing PDE in time
//     auto Lt = -bilaplacian<SPLINE>();
//     PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);
    
//     // ----GCV---

//     std::vector<DVector<double>> lambdas_d_t;
//     for(double xs = -5.0; xs <= -2.8; xs +=0.05)
//         for(double xt = -7.0; xt <= -6.0; xt +=1.0) 
//             lambdas_d_t.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));
    
//     // define model
//     STRPDE<SpaceTimeSeparable, fdapde::monolithic> model_gcv(space_penalty, time_penalty, Sampling::pointwise);

//     model_gcv.set_spatial_locations(space_locs);
//     model_gcv.set_temporal_locations(time_locs);
//     std::cout << "here 4" << std::endl;

//     BlockFrame<double, int> df_gcv;
//     df_gcv.stack(OBSERVATIONS_BLK, y);
//     model_gcv.set_data(df_gcv);
//     model_gcv.init();

//     // define GCV function and grid of \lambda_D values
//     auto GCV = model_gcv.gcv<ExactEDF>();
//     // optimize GCV
//     Grid<fdapde::Dynamic> opt;
//     opt.optimize(GCV, lambdas_d_t);
//     SVector<2> best_lambda = opt.optimum();
//     std::cout << "here 10" << std::endl;

//     // std::string solutions_path = "C:/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/STRPDE/NO2/gcv_3";

//     std::string solutions_path = path_data;
//     // Save Lambda opt
//     std::ofstream fileLambdaoptS(solutions_path + "/lambda_s_opt.csv");
//     if(fileLambdaoptS.is_open()){
//         fileLambdaoptS << std::setprecision(16) << best_lambda[0];
//         fileLambdaoptS.close();
//     }
//     std::ofstream fileLambdaoptT(solutions_path + "/lambda_t_opt.csv");
//     if (fileLambdaoptT.is_open()){
//         fileLambdaoptT << std::setprecision(16) << best_lambda[1];
//         fileLambdaoptT.close();
//     }

//     // Save GCV scores
//     std::ofstream fileGCV_scores(solutions_path + "/gcv_scores.csv");
//     for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
//         fileGCV_scores << std::setprecision(16) << GCV.gcvs()[i] << "\n" ; 



//     // ----RUN---

//     double lambda_D = read_csv<double>(solutions_path + "/lambda_s_opt.csv")(0,0); 
//     double lambda_T = read_csv<double>(solutions_path + "/lambda_t_opt.csv")(0,0);

//     STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::pointwise);
//     model.set_lambda_D(lambda_D);
//     model.set_lambda_T(lambda_T);
//     model.set_spatial_locations(space_locs);
//     model.set_temporal_locations(time_locs);
//     // set model's data
//     BlockFrame<double, int> df;
//     df.stack(OBSERVATIONS_BLK, y);
//     model.set_data(df);
//     // solve smoothing problem
//     model.init();
//     model.solve();

//     // Save C++ solution 
//     DMatrix<double> computedF = model.f();
//     const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//     std::ofstream filef(solutions_path + "/f_STRPDE.csv");
//     if (filef.is_open()){
//         filef << computedF.format(CSVFormatf);
//         filef.close();
//     }

// }






// -------------------------- PM 10 ------------------------

// // gcv time

// TEST(case_study_gcv, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

//     const std::string pollutant = "PM10";   // PM10 NO2
//     const std::string type_locs = "all";   // "reduced" "all"
//     const std::string periodicity = "by_week";   // "by_day" "by_week"  ---> MODIFICA ANCHE M!
//     const std::string data_threshold = "4"; 
//     const std::string outlier = "/no_out";   // "" "/no_out"
//     const std::string eps_string = "1e+1";   // "0"  "1e+0" "1e+1"
//     const bool esagoni = true; 

//     std::string gcv_type = "stochastic";   // "exact" "stochastic"  ---> MODIFICA ANCHE GIU'!
//     std::size_t seed = 438172;
//     unsigned int MC_run = 100; 

//     std::string est_type = "mean";    // mean quantile
//     std::vector<double> alphas = {0.9}; 

//     // Marco 
//     std::string path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia"; 
//     std::string path_data = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant; 
//     std::string solutions_path; 

//     // // Ilenia 
//     // std::string path_data = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant; 
//     // std::string solutions_path; 
//     // if(est_type == "mean")
//     //     solutions_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/STRPDE/" + pollutant + "/" + periodicity + "/gcv_Cpp_" + gcv_type; 
//     // if(est_type == "quantile")
//     //     solutions_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/QSTRPDE/" + pollutant + "/" + periodicity + "/gcv_Cpp_" + gcv_type + "/eps_" + eps_string; 

//     if(est_type == "mean")
//         solutions_path = path + "/STRPDE/" + pollutant + "/" + periodicity + "/gcv_Cpp_" + gcv_type; 
//     if(est_type == "quantile")
//         solutions_path = path + "/QSTRPDE/" + pollutant + "/" + periodicity + "/gcv_Cpp_" + gcv_type + "/eps_" + eps_string; 

//     // lambdas sequence 
//     std::vector<DVector<double>> lambdas_d_t; std::vector<double> lambdas_d; std::vector<double> lambdas_t;
//     std::vector<double> lambdas50_d; std::vector<double> lambdas50_t; std::vector<DVector<double>> lambdas50_d_t; 
//     std::vector<double> lambdas90_d; std::vector<double> lambdas90_t; std::vector<DVector<double>> lambdas90_d_t; 
//     std::vector<double> lambdas95_d; std::vector<double> lambdas95_t; std::vector<DVector<double>> lambdas95_d_t; 

//     if(est_type == "mean"){
//         for(double xs = -6.4; xs <= -1.0; xs += 0.2)
//             lambdas_d.push_back(std::pow(10,xs));

//         for(double xt = -5.0; xt <= +3.0; xt += 0.5)
//             lambdas_t.push_back(std::pow(10,xt));    

//         for(auto i = 0; i < lambdas_d.size(); ++i)
//             for(auto j = 0; j < lambdas_t.size(); ++j) 
//                 lambdas_d_t.push_back(SVector<2>(lambdas_d[i], lambdas_t[j]));
//     }
//     if(est_type == "quantile"){

//         // 50% 
//         {
//             for(double xs = -12.0; xs <= -1.0; xs += 1.0)
//                 lambdas50_d.push_back(std::pow(10,xs));

//             for(double xt = -5.0; xt <= -3.9; xt += 2.0)
//                 lambdas50_t.push_back(std::pow(10,xt));    

//             for(auto i = 0; i < lambdas50_d.size(); ++i)
//                 for(auto j = 0; j < lambdas50_t.size(); ++j) 
//                     lambdas50_d_t.push_back(SVector<2>(lambdas50_d[i], lambdas50_t[j]));
//         }

//         // 90% 
//         {
//             for(double xs = -11.0; xs <= -5.0; xs += 1.0)
//                 lambdas90_d.push_back(std::pow(10,xs));

//             for(double xt = -5.0; xt <= -4.0; xt += 2.0)
//                 lambdas90_t.push_back(std::pow(10,xt));    

//             for(auto i = 0; i < lambdas90_d.size(); ++i)
//                 for(auto j = 0; j < lambdas90_t.size(); ++j) 
//                     lambdas90_d_t.push_back(SVector<2>(lambdas90_d[i], lambdas90_t[j]));
//         }

//         // 95% 
//         {
//             for(double xs = -8.0; xs <= -1.0; xs += 1.0)
//                 lambdas95_d.push_back(std::pow(10,xs));

//             for(double xt = -5.0; xt <= 1.0; xt += 2.0)
//                 lambdas95_t.push_back(std::pow(10,xt));    

//             for(auto i = 0; i < lambdas95_d.size(); ++i)
//                 for(auto j = 0; j < lambdas95_t.size(); ++j) 
//                     lambdas95_d_t.push_back(SVector<2>(lambdas95_d[i], lambdas95_t[j]));
//         }

//     }

//     // define temporal domain
//     unsigned int M = 22; 
//     std::string M_string = std::to_string(M);
//     double tf; // final time 
//     if(periodicity == "by_day")
//         tf = 364.0;    
//     if(periodicity == "by_week")
//         tf = 357.0;     
//     Mesh<1, 1> time_mesh(0, tf, M-1);
//     // define spatial domain and regularizing PDE
//     MeshLoader<Mesh2D> domain("mesh_lombardia_esagoni");

//     // import data and locs from files
//     DMatrix<double> y; DMatrix<double> space_locs; DMatrix<double> time_locs; 
//     if(type_locs == "reduced"){
//         if(pollutant == "NO2"){
//             y = read_csv<double>(path_data + "/NO2_max_2022_Cpp_reduced.csv");
//         } else{
//             y = read_csv<double>(path_data + "/PM10_2022_Cpp_reduced.csv");
//         }
//         time_locs = read_csv<double>(path_data + "/time_locations_reduced.csv");
//     } else{
//         if(pollutant == "NO2"){
//             y = read_csv<double>(path_data + "/NO2_max_2022_Cpp.csv");
//             time_locs = read_csv<double>(path_data + "/time_locations.csv");
//         } else{
//             y = read_csv<double>(path_data + "/" + periodicity + outlier + "/PM10_2022_threshold" + data_threshold + ".csv");  
//             time_locs = read_csv<double>(path_data + "/" + periodicity + "/time_locs.csv"); 
//         }     
//     }
//     if(!esagoni)
//         space_locs = read_csv<double>(path_data + "/locs.csv"); 
//     else
//         space_locs = read_csv<double>(path_data + "/locs_esagoni.csv");
    
//     // check dimensions
//     std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
//     std::cout << "dim time loc " << time_locs.rows() << " " << time_locs.cols() << std::endl;
//     std::cout << "dim y " << y.rows() << " " << y.cols() << std::endl;

//     BlockFrame<double, int> df;
//     df.stack(OBSERVATIONS_BLK, y);
   
//     // define regularizing PDE in space 
//     auto Ld = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
//     PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

//     // define regularizing PDE in time
//     auto Lt = -bilaplacian<SPLINE>();
//     PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

//     std::cout << "-----------------------------GCV STARTS------------------------" << std::endl; 

//     if(est_type == "mean"){

//         STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::pointwise);
    
//         // set model's data
//         model.set_spatial_locations(space_locs);
//         model.set_temporal_locations(time_locs);
        
//         model.set_data(df);
//         model.init();

//         // define GCV function and grid of \lambda_D values

//         // stochastic
//         auto GCV = model.gcv<StochasticEDF>(MC_run, seed);
//         // // exact
//         // auto GCV = model.gcv<ExactEDF>();

//         // optimize GCV
//         Grid<fdapde::Dynamic> opt;
//         opt.optimize(GCV, lambdas_d_t);
//         SVector<2> best_lambda = opt.optimum();

//         // Save lambda sequence 
//         std::ofstream fileLambda_S_Seq(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/mesh_new/lambdas_S_seq.csv");
//         for(std::size_t i = 0; i < lambdas_d.size(); ++i) 
//             fileLambda_S_Seq << std::setprecision(16) << lambdas_d[i] << "\n"; 
//         fileLambda_S_Seq.close();

//         for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
//             std::cout << lambdas_t[i] << "\n"; 

//         std::ofstream fileLambda_T_Seq(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/mesh_new/lambdas_T_seq.csv");
//         for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
//             fileLambda_T_Seq << std::setprecision(16) << lambdas_t[i] << "\n"; 
//         fileLambda_T_Seq.close();

//         // Save Lambda opt
//         std::ofstream fileLambdaoptS(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/mesh_new/lambda_s_opt.csv");
//         if(fileLambdaoptS.is_open()){
//             fileLambdaoptS << std::setprecision(16) << best_lambda[0];
//             fileLambdaoptS.close();
//         }
//         std::ofstream fileLambdaoptT(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/mesh_new/lambda_t_opt.csv");
//         if (fileLambdaoptT.is_open()){
//             fileLambdaoptT << std::setprecision(16) << best_lambda[1];
//             fileLambdaoptT.close();
//         }
//         // Save GCV scores
//         std::ofstream fileGCV_scores(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/mesh_new/gcv_scores.csv");
//         for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
//             fileGCV_scores << std::setprecision(16) << GCV.gcvs()[i] << "\n"; 
//         fileGCV_scores.close();

//         // Save edfs
//         std::ofstream fileEDF(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/mesh_new/edfs.csv");
//         for(std::size_t i = 0; i < GCV.edfs().size(); ++i) 
//             fileEDF << std::setprecision(16) << GCV.edfs()[i] << "\n"; 
//         fileEDF.close();

//     }

//     if(est_type == "quantile"){
        
//         for(double alpha : alphas){

//             unsigned int alpha_int = alpha*100; 
//             std::string alpha_string = std::to_string(alpha_int); 

//             std::cout << "----------------ALPHA = " << alpha_string << "----------------" << std::endl; 

//             QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);
    
//             // set model's data
//             model.set_spatial_locations(space_locs);
//             model.set_temporal_locations(time_locs);

//             model.set_exact_gcv(eps_string == "0");   // true: no lisciamento curva gcv
//             if(eps_string == "1e+1"){
//                 model.set_eps_power(1.0); 
//             }
//             if(eps_string == "1e+0"){
//                 model.set_eps_power(0.0); 
//             }
//             if(eps_string == "1e-1"){
//                 model.set_eps_power(-1.0); 
//             }
            
//             model.set_data(df);
//             model.init();

//             // define GCV function and grid of \lambda_D values

//             // stochastic
//             auto GCV = model.gcv<StochasticEDF>(MC_run, seed);
//             // // exact
//             // auto GCV = model.gcv<ExactEDF>();

//             if(alpha_string == "50"){
//                 lambdas_d = lambdas50_d; 
//                 lambdas_t = lambdas50_t;
//                 lambdas_d_t = lambdas50_d_t;
//             }
                 
//             if(alpha_string == "90"){
//                 lambdas_d = lambdas90_d; 
//                 lambdas_t = lambdas90_t;
//                 lambdas_d_t = lambdas90_d_t;
//             }
                 
//             if(alpha_string == "95"){
//                 lambdas_d = lambdas95_d; 
//                 lambdas_t = lambdas95_t;
//                 lambdas_d_t = lambdas95_d_t;
//             }
//             // optimize GCV
//             Grid<fdapde::Dynamic> opt;
//             opt.optimize(GCV, lambdas_d_t);
//             SVector<2> best_lambda = opt.optimum();

//             // Save lambda sequences 
//             std::ofstream fileLambda_S_Seq(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/alpha_" + alpha_string + "/mesh_new/lambdas_S_seq.csv");
//             for(std::size_t i = 0; i < lambdas_d.size(); ++i) 
//                 fileLambda_S_Seq << std::setprecision(16) << lambdas_d[i] << "\n"; 
//             fileLambda_S_Seq.close();
//             for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
//                 std::cout << lambdas_t[i] << "\n"; 
//             std::ofstream fileLambda_T_Seq(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/alpha_" + alpha_string + "/mesh_new/lambdas_T_seq.csv");
//             for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
//                 fileLambda_T_Seq << std::setprecision(16) << lambdas_t[i] << "\n"; 
//             fileLambda_T_Seq.close();

//             // Save Lambda opt
//             std::ofstream fileLambdaoptS(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/alpha_" + alpha_string + "/mesh_new/lambda_s_opt.csv");
//             if(fileLambdaoptS.is_open()){
//                 fileLambdaoptS << std::setprecision(16) << best_lambda[0];
//                 fileLambdaoptS.close();
//             }
//             std::ofstream fileLambdaoptT(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/alpha_" + alpha_string + "/mesh_new/lambda_t_opt.csv");
//             if (fileLambdaoptT.is_open()){
//                 fileLambdaoptT << std::setprecision(16) << best_lambda[1];
//                 fileLambdaoptT.close();
//             }
//             // Save GCV scores
//             std::ofstream fileGCV_scores(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/alpha_" + alpha_string + "/mesh_new/gcv_scores.csv");
//             for(std::size_t i = 0; i < GCV.gcvs().size(); ++i){
//                 fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.gcvs()[i]) << "\n"; 
//                 std::cout << "----- " << std::sqrt(GCV.gcvs()[i]) << std::endl ; 
//             }
                
//             fileGCV_scores.close();
//             // Save edfs
//             std::ofstream fileEDF(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/alpha_" + alpha_string + "/mesh_new/edfs.csv");
//             for(std::size_t i = 0; i < GCV.edfs().size(); ++i) 
//                 fileEDF << std::setprecision(16) << GCV.edfs()[i] << "\n"; 
//             fileEDF.close();  
//             // Save gcv numerator
//             std::cout << "model.n_obs() = " << model.n_obs() << std::endl; 
//             std::ofstream fileGCVnum(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/alpha_" + alpha_string + "/mesh_new/gcv_num.csv");
//             for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
//                 fileGCVnum << std::setprecision(16) << std::sqrt(GCV.gcvs()[i]*model.n_obs()) << "\n"; 
//             fileGCVnum.close();

//         }

//     }

      
// }

// // run time

// TEST(case_study_run, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

//     std::string pollutant = "PM10";   // PM10 NO2
//     std::string type_locs = "all";   // "reduced" "all"
//     std::string gcv_type = "stochastic";   // "exact" "stochastic"  
//     const std::string periodicity = "by_week";   // "by_day" "by_week"  ---> MODIFICA ANCHE M!
//     const std::string data_threshold = "4"; 
//     const std::string outlier = "/no_out";   // "" "/no_out"
//     const std::string eps_string = "1e+1";   // "0" "1e+0" "1e+1"
//     const bool esagoni = true;

//     std::string est_type = "mean";    // mean quantile
//     std::vector<double> alphas = {0.9}; 

//     // Marco 
//     std::string path_data = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant; 
//     std::string path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia"; 
//     std::string solutions_path; 
//     if(est_type == "mean")
//         solutions_path = path + "/STRPDE/" + pollutant + "/" + periodicity + "/gcv_Cpp_" + gcv_type; 
//     if(est_type == "quantile")
//         solutions_path = path + "/QSTRPDE/" + pollutant +"/" + periodicity + "/gcv_Cpp_" + gcv_type + "/eps_" + eps_string; 
   
//     // Ilenia 
//     //std::string path_data = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant ; 
//     // std::string solutions_path; 
//     // if(est_type == "mean")
//     //     solutions_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/STRPDE/" + pollutant + "/" + periodicity + "/gcv_Cpp_" + gcv_type; 
//     // if(est_type == "quantile")
//     //     solutions_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/QSTRPDE/" + pollutant +"/" + periodicity + "/gcv_Cpp_" + gcv_type + "/eps_" + eps_string; 


//     // define temporal domain
//     unsigned int M = 22; 
//     std::string M_string = std::to_string(M);
//     double tf; 
//     if(periodicity == "by_day")
//         tf = 364.0;    // final time 
//     if(periodicity == "by_week")
//         tf = 357.0;    // final time 
//     Mesh<1, 1> time_mesh(0, tf, M-1);
//     // define spatial domain and regularizing PDE
//     MeshLoader<Mesh2D> domain("mesh_lombardia_esagoni");

//     // import data and locs from files
//     DMatrix<double> y; DMatrix<double> space_locs; DMatrix<double> time_locs; 
//     if(type_locs == "reduced"){
//         if(pollutant == "NO2"){
//             y = read_csv<double>(path_data + "/" + periodicity + "/NO2_max_2022_Cpp_reduced.csv");
//         } else{
//             y = read_csv<double>(path_data + "/" + periodicity + "/PM10_2022_Cpp_reduced.csv");
//         }
//         time_locs = read_csv<double>(path_data + "/" + periodicity + "/time_locations_reduced.csv");
//     } else{
//         if(pollutant == "NO2"){
//             y = read_csv<double>(path_data + "/" + periodicity + "/NO2_max_2022_Cpp.csv");
//             time_locs = read_csv<double>(path_data + "/time_locations.csv");
//         } else{
//             y = read_csv<double>(path_data + "/" + periodicity + outlier + "/PM10_2022_threshold" + data_threshold + ".csv");  
//             time_locs = read_csv<double>(path_data + "/" + periodicity + "/time_locs.csv"); 
//         }     
//     }
//     if(!esagoni)
//         space_locs = read_csv<double>(path_data + "/locs.csv"); 
//     else
//         space_locs = read_csv<double>(path_data + "/locs_esagoni.csv");

//     // check dimensions
//     std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
//     std::cout << "dim time loc " << time_locs.rows() << " " << time_locs.cols() << std::endl;
//     std::cout << "dim y " << y.rows() << " " << y.cols() << std::endl;

//     BlockFrame<double, int> df;
//     df.stack(OBSERVATIONS_BLK, y);
   
//     // define regularizing PDE in space 
//     auto Ld = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
//     PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

//     // define regularizing PDE in time
//     auto Lt = -bilaplacian<SPLINE>();
//     PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

//     std::cout << "--------------------------------RUN STARTS--------------------------------" << std::endl; 
//     if(est_type == "mean"){

//         STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::pointwise);
    
//         // set model's data
//         model.set_spatial_locations(space_locs);
//         model.set_temporal_locations(time_locs);

//         // Read optima lambdas 
//         double lambda_T; double lambda_S; 
//         std::ifstream fileLambdaT_opt(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/mesh_new/lambda_t_opt.csv");
//         if(fileLambdaT_opt.is_open()){
//             fileLambdaT_opt >> lambda_T; 
//             fileLambdaT_opt.close();
//         }
//         std::ifstream fileLambdaS_opt(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/mesh_new/lambda_s_opt.csv");
//         if(fileLambdaS_opt.is_open()){
//             fileLambdaS_opt >> lambda_S; 
//             fileLambdaS_opt.close();
//         }

//         std::cout << "lambda S " << lambda_S << std::endl;
//         std::cout << "lambda T " << lambda_T << std::endl;

//         model.set_lambda_D(lambda_S);
//         model.set_lambda_T(lambda_T);
        
//         model.set_data(df);

//         model.init();
//         model.solve();

//         // Save C++ solution 
//         DMatrix<double> computedF = model.f();
//         const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filef(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/mesh_new/f.csv");
//         if (filef.is_open()){
//             filef << computedF.format(CSVFormatf);
//             filef.close();
//         }

//         DMatrix<double> computedFn = model.Psi()*model.f();
//         const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filefn(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/mesh_new/fn.csv");
//         if (filefn.is_open()){
//             filefn << computedFn.format(CSVFormatfn);
//             filefn.close();
//         }
//     }

//     if(est_type == "quantile"){
        
//         for(double alpha : alphas){

//             unsigned int alpha_int = alpha*100; 
//             std::string alpha_string = std::to_string(alpha_int); 

//             std::cout << "----------------ALPHA = " << alpha_string << "----------------" << std::endl; 

//             QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);
    
//             // set model's data
//             model.set_spatial_locations(space_locs);
//             model.set_temporal_locations(time_locs);

//             // Read optima lambdas 
//             double lambda_T; double lambda_S; 
//             std::ifstream fileLambdaT_opt(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/alpha_" + alpha_string + "/mesh_new/lambda_t_opt.csv");
//             if(fileLambdaT_opt.is_open()){
//                 fileLambdaT_opt >> lambda_T; 
//                 fileLambdaT_opt.close();
//             }
//             std::ifstream fileLambdaS_opt(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/alpha_" + alpha_string + "/mesh_new/lambda_s_opt.csv");
//             if(fileLambdaS_opt.is_open()){
//                 fileLambdaS_opt >> lambda_S; 
//                 fileLambdaS_opt.close();
//             }

//             std::cout << "lambda S " << lambda_S << std::endl;
//             std::cout << "lambda T " << lambda_T << std::endl;

//             model.set_exact_gcv(true);   // true: no lisciamento curva gcv
//             model.set_lambda_D(lambda_S);
//             model.set_lambda_T(lambda_T);
            
//             model.set_data(df);

//             model.init();
//             model.solve();

//             // Save C++ solution 
//             DMatrix<double> computedF = model.f();
//             const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//             std::ofstream filef(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/alpha_" + alpha_string + "/mesh_new/f.csv");
//             if (filef.is_open()){
//                 filef << computedF.format(CSVFormatf);
//                 filef.close();
//             }

//             DMatrix<double> computedFn = model.Psi()*model.f();
//             const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//             std::ofstream filefn(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/alpha_" + alpha_string + "/mesh_new/fn.csv");
//             if (filefn.is_open()){
//                 filefn << computedFn.format(CSVFormatfn);
//                 filefn.close();
//             }

//         }

//     }        
// }





// // Run to save A()

// TEST(case_study_run, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

//     std::string pollutant = "PM10";   // PM10 NO2
//     std::string type_locs = "all";   // "reduced" "all"
//     std::string gcv_type = "stochastic";   // "exact" "stochastic"  
//     const std::string periodicity = "by_week";   // "by_day" "by_week"  ---> MODIFICA ANCHE M!
//     const std::string data_threshold = "4"; 
//     const std::string outlier = "/no_out";   // "" "/no_out"
//     const std::string eps_string = "1e+0";   // "0" "1e+0" "1e+1"

//     std::string est_type = "mean";    // mean quantile
//     std::vector<double> alphas = {0.9}; 

//     // Marco 
//     // std::string path_data = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant; 
   
//     // Ilenia 
//     std::string path_data = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant ; 
//     std::string solutions_path; 
//     if(est_type == "mean")
//         solutions_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/STRPDE/" + pollutant + "/" + periodicity + "/gcv_Cpp_" + gcv_type; 
//     if(est_type == "quantile")
//         solutions_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/QSTRPDE/" + pollutant +"/" + periodicity + "/gcv_Cpp_" + gcv_type + "/eps_" + eps_string; 

//     std::cout << "Solution path : " << solutions_path << std::endl ; 

//     // define temporal domain
//     unsigned int M = 22; 
//     std::string M_string = std::to_string(M);
//     double tf; 
//     if(periodicity == "by_day")
//         tf = 364.0;    // final time 
//     if(periodicity == "by_week")
//         tf = 357.0;    // final time 
//     Mesh<1, 1> time_mesh(0, tf, M-1);
//     // define spatial domain and regularizing PDE
//     MeshLoader<Mesh2D> domain("mesh_lombardia_coarse");

//     // import data and locs from files
//     DMatrix<double> y; DMatrix<double> space_locs; DMatrix<double> time_locs; 
//     if(type_locs == "reduced"){
//         if(pollutant == "NO2"){
//             y = read_csv<double>(path_data + "/" + periodicity + "/NO2_max_2022_Cpp_reduced.csv");
//         } else{
//             y = read_csv<double>(path_data + "/" + periodicity + "/PM10_2022_Cpp_reduced.csv");
//         }
//         time_locs = read_csv<double>(path_data + "/" + periodicity + "/time_locations_reduced.csv");
//     } else{
//         if(pollutant == "NO2"){
//             y = read_csv<double>(path_data + "/" + periodicity + "/NO2_max_2022_Cpp.csv");
//             time_locs = read_csv<double>(path_data + "/time_locations.csv");
//         } else{
//             y = read_csv<double>(path_data + "/" + periodicity + outlier + "/PM10_2022_threshold" + data_threshold + ".csv");  
//             time_locs = read_csv<double>(path_data + "/" + periodicity + "/time_locs.csv"); 
//         }     
//     }
//     space_locs = read_csv<double>(path_data + "/locs.csv");

//     // check dimensions
//     std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
//     std::cout << "dim time loc " << time_locs.rows() << " " << time_locs.cols() << std::endl;
//     std::cout << "dim y " << y.rows() << " " << y.cols() << std::endl;

//     BlockFrame<double, int> df;
//     df.stack(OBSERVATIONS_BLK, y);
   
//     // define regularizing PDE in space 
//     auto Ld = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
//     PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

//     // define regularizing PDE in time
//     auto Lt = -bilaplacian<SPLINE>();
//     PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

//     std::cout << "--------------------------------RUN STARTS--------------------------------" << std::endl; 
//     if(est_type == "mean"){

//         STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::pointwise);
    
//         // set model's data
//         model.set_spatial_locations(space_locs);
//         model.set_temporal_locations(time_locs);

//         // // Read optima lambdas 
//         // double lambda_T; double lambda_S; 
//         // std::ifstream fileLambdaT_opt(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/lambda_t_opt.csv");
//         // if(fileLambdaT_opt.is_open()){
//         //     fileLambdaT_opt >> lambda_T; 
//         //     fileLambdaT_opt.close();
//         // }
//         // std::ifstream fileLambdaS_opt(solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/lambda_s_opt.csv");
//         // if(fileLambdaS_opt.is_open()){
//         //     fileLambdaS_opt >> lambda_S; 
//         //     fileLambdaS_opt.close();
//         // }

//         // --------- cond A

//         // Impose lambdas 
//         double lambda_T = 1e-5; 
//         double lambda_S ;

//         std::vector<double> seq_lambda_S;
//         for(double xs = -10.0; xs <= +2.0; xs += 0.5)
//             seq_lambda_S.push_back(std::pow(10,xs));

//         std::size_t count_lambda = 1;

//         for(double lambda_S_choice : seq_lambda_S ){

//             lambda_S = lambda_S_choice;

//             std::cout << "lambda S " << lambda_S << std::endl;
//             std::cout << "lambda T " << lambda_T << std::endl;

//             model.set_lambda_D(lambda_S);
//             model.set_lambda_T(lambda_T);
            
//             model.set_data(df);

//             model.init();
//             model.solve();

//             // Save matrix A 
//             SpMatrix<double> computedA = model.A();

//             Eigen::saveMarket(computedA, solutions_path + "/" + type_locs + "_threshold" + data_threshold + "/M_" + M_string + outlier + "/A/A_" + std::to_string(count_lambda) + ".mtx");

//             count_lambda = count_lambda +1;
        
//         }

        
//     }

    
// }



















// ------------------------------------------------------------OLD---------------------------------------------------

// // STRPDE 
// TEST(strpde_case_study, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {
    
//     const std::string path_data = "C:/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/our_marco/data"; 

//     std::cout << "here 1 " << std::endl; 
//     // define temporal domain
//     DVector<double> time_mesh;
//     unsigned int M = 12; 
//     time_mesh.resize(M);
//     for (std::size_t i = 0; i < M; ++i) time_mesh[i] = i+1;
//     // define spatial domain and regularizing PDE
//     MeshLoader<Mesh2D> domain("/case_study/mesh_reshaped");
//     // import data from files
//     DMatrix<double> time_locs  = read_csv<double>(path_data + "/time_locs.csv");
//     DMatrix<double> space_locs = read_csv<double>(path_data + "/space_locs.csv");
//     DMatrix<double> y          = read_csv<double>(path_data + "/y.csv");
//     std::cout << "here 2" << std::endl;
//     std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
//     std::cout << "dim time loc " << time_locs.rows() << " " << time_locs.cols() << std::endl;
//     std::cout << "dim y " << y.rows() << " " << y.cols() << std::endl;
//     // define regularizing PDE
//     auto L = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.rows(), 1);
//     PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    
//     // ----GCV---

//     std::cout << "here 3" << std::endl;
//     std::vector<SVector<2>> lambdas_d_t;
//     for(double xs = -10.0; xs <= -8.0; xs +=0.1)
//         for(double xt = -8.0; xt <= -5.0; xt +=0.5) 
//             lambdas_d_t.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));
    
//     // define model
//     STRPDE<decltype(problem), SpaceTimeSeparable, GeoStatLocations, MonolithicSolver> model_gcv(problem, time_mesh);

//     model_gcv.set_spatial_locations(space_locs);
//     model_gcv.set_temporal_locations(time_locs);
//     std::cout << "here 4" << std::endl;

//     BlockFrame<double, int> df_gcv;
//     df_gcv.stack(OBSERVATIONS_BLK, y);
//     model_gcv.set_data(df_gcv);
//     std::cout << "here 5" << std::endl;
//     model_gcv.init();
//     std::cout << "here 6" << std::endl;

//     // define GCV function and grid of \lambda_D values
//     GCV<decltype(model_gcv), ExactEDF<decltype(model_gcv)>> GCV(model_gcv);
//     ScalarField<2, decltype(GCV)> obj(GCV);  
//     // optimize GCV
//     Grid<2> opt;
//     opt.optimize(obj, lambdas_d_t);
//     std::cout << "here 7" << std::endl;
//     SVector<2> best_lambda = opt.optimum();
//     std::cout << "here 8" << std::endl;

//     std::string solutions_path = "C:/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/our_marco/results"; 
//     // Save Lambda opt
//     std::ofstream fileLambdaoptS(solutions_path + "/lambda_s_opt.csv");
//     if(fileLambdaoptS.is_open()){
//         fileLambdaoptS << std::setprecision(16) << best_lambda[0];
//         fileLambdaoptS.close();
//     }
//     std::ofstream fileLambdaoptT(solutions_path + "/lambda_t_opt.csv");
//     if (fileLambdaoptT.is_open()){
//         fileLambdaoptT << std::setprecision(16) << best_lambda[1];
//         fileLambdaoptT.close();
//     }


//     // ----RUN---

//     double lambda_D = read_csv<double>(solutions_path + "/lambda_s_opt.csv")(0,0); 
//     double lambda_T = read_csv<double>(solutions_path + "/lambda_t_opt.csv")(0,0);

//     STRPDE<decltype(problem), SpaceTimeSeparable, GeoStatLocations, MonolithicSolver> model(problem, time_mesh);
//     model.set_lambda_D(lambda_D);
//     model.set_lambda_T(lambda_T);
//     model.set_spatial_locations(space_locs);
//     model.set_temporal_locations(time_locs);
//     // set model's data
//     BlockFrame<double, int> df;
//     df.stack(OBSERVATIONS_BLK, y);
//     model.set_data(df);
//     // solve smoothing problem
//     model.init();
//     model.solve();

//     // Save C++ solution 
//     DMatrix<double> computedF = model.f();
//     const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//     std::ofstream filef(solutions_path + "/f_STRPDE.csv");
//     if (filef.is_open()){
//         filef << computedF.format(CSVFormatf);
//         filef.close();
//     }

// }








// // gcv  

// TEST(case_study_gcv_space_only, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

//     std::string pollutant = "NO2";   // PM10 NO2

//     const std::string type_locs = "all";   // "reduced" "all"

//     const std::string type_space_mesh = "coarse";   // coarse ---> CAMBIA MESH SOTTO !!! 
    
//     const std::string type_data = "/data_rescaled";  // "data_original"  ---> MODIFICA ANCHE GIU'lettura dati !!!

//     const std::string time_step = "by_week";   // "" "by_day" "by_week"

//     std::string gcv_type = "stochastic";   // "exact" "stochastic"  ---> MODIFICA ANCHE GIU'!
//     std::size_t seed = 438172;
//     unsigned int MC_run = 100; 

//     std::string est_type = "quantile";    // mean quantile
//     std::vector<double> alphas = {0.5};     //  {0.5, 0.9, 0.95};     

//     std::string lambda_selection_type = "gcv"; // _smooth_eps1e-1.5

//     // // Marco 
//     // std::string path_data = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant; 
    
//     // std::string solutions_path; 
//     // if(est_type == "mean")
//     //     solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/STRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 
//     // if(est_type == "quantile")
//     //     solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/QSTRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 

//     // Ilenia 
//     std::string path_data = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant + "/week_52"  ; 
//     std::string solutions_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/QSRPDE/" + pollutant ; 


//     // lambdas sequence 

//     std::vector<DVector<double>> lambdas_50;
//     std::vector<DVector<double>> lambdas_90;
//     std::vector<DVector<double>> lambdas_95;


//     // 50% 
//     for(double x = -5.0; x <= -1.0; x += 0.1) lambdas_50.push_back(SVector<1>(std::pow(10, x)));


//     // define spatial domain and regularizing PDE

//     MeshLoader<Mesh2D> domain("mesh_lombardia_coarse");        // ATT: da cambiare in base a type_space_mesh ma con if non funziona


    
//     // import data and locs from files
//     DMatrix<double> y; DMatrix<double> space_locs; DMatrix<double> time_locs; 
    
//     if(pollutant == "NO2"){
//         y = read_csv<double>(path_data + "/NO2_max.csv");   
//     } else{
//         y = read_csv<double>(path_data + "/PM10_2022.csv");  
//     }     

//     space_locs = read_csv<double>(path_data + "/locs.csv");
    
//     // check dimensions
//     std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
//     std::cout << "dim y " << y.rows() << " " << y.cols() << std::endl;

//     BlockFrame<double, int> df;
//     df.stack(OBSERVATIONS_BLK, y);
   
//     // define regularizing PDE in space 
//     auto Ld = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 , 1);
//     PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

        
//     for(double alpha : alphas){

//         unsigned int alpha_int = alpha*100; 
//         std::string alpha_string = std::to_string(alpha_int); 

//         std::cout << "----------------ALPHA = " << alpha_string << "----------------" << std::endl; 

//         QSRPDE<SpaceOnly> model(space_penalty, Sampling::pointwise, alpha);

//         // set model's data
//         model.set_spatial_locations(space_locs);

//         model.set_exact_gcv(lambda_selection_type == "gcv");
//         if(lambda_selection_type == "gcv_smooth_eps1e-3")
//             model.set_eps_power(-3.0);
//         if(lambda_selection_type == "gcv_smooth_eps1e-2")
//             model.set_eps_power(-2.0);
//         if(lambda_selection_type == "gcv_smooth_eps1e-1.5")
//             model.set_eps_power(-1.5);
//         if(lambda_selection_type == "gcv_smooth_eps1e-1")
//             model.set_eps_power(-1.0);
        
//         model.set_data(df);
//         model.init();

//         // define GCV function and grid of \lambda_D values

//         // stochastic
//         auto GCV = model.gcv<StochasticEDF>(MC_run, seed);
//         // // exact
//         // auto GCV = model.gcv<ExactEDF>();

//         // optimize GCV
//         Grid<fdapde::Dynamic> opt;
//         opt.optimize(GCV, lambdas_50);
//         SVector<1> best_lambda = opt.optimum();

//         // Save Lambda opt
//         std::ofstream fileLambdaoptS(solutions_path +  "/lambda_s_opt.csv");
//         if(fileLambdaoptS.is_open()){
//             fileLambdaoptS << std::setprecision(16) << best_lambda[0];
//             fileLambdaoptS.close();
//         }

//         // Save GCV scores
//         std::ofstream fileGCV_scores(solutions_path + "/gcv_scores.csv");
//         for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
//             fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.gcvs()[i]) << "\n"; 

//         fileGCV_scores.close();

//         // Save edfs
//         std::ofstream fileEDF(solutions_path +  "/edfs.csv");
//         for(std::size_t i = 0; i < GCV.edfs().size(); ++i) 
//             fileEDF << std::setprecision(16) << GCV.edfs()[i] << "\n"; 
//         fileEDF.close();
//     }

      
// }



// // run

// TEST(case_study_run_space_only, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

//     std::string pollutant = "NO2";   // PM10 NO2

//     const std::string type_locs = "all";   // "reduced" "all"

//     const std::string type_space_mesh = "coarse";   // coarse ---> CAMBIA MESH SOTTO !!! 
    
//     const std::string type_data = "/data_rescaled";  // "data_original"  ---> MODIFICA ANCHE GIU'lettura dati !!!

//     const std::string time_step = "by_week";   // "" "by_day" "by_week"

//     std::string gcv_type = "stochastic";   // "exact" "stochastic"  ---> MODIFICA ANCHE GIU'!
//     std::size_t seed = 438172;
//     unsigned int MC_run = 100; 

//     std::string est_type = "quantile";    // mean quantile
//     std::vector<double> alphas = {0.5};     //  {0.5, 0.9, 0.95};     

//     std::string lambda_selection_type = "gcv"; // _smooth_eps1e-1.5

//     // // Marco 
//     // std::string path_data = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant; 
    
//     // std::string solutions_path; 
//     // if(est_type == "mean")
//     //     solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/STRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 
//     // if(est_type == "quantile")
//     //     solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/QSTRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 

//     // Ilenia 
//     std::string path_data = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant + "/week_52"  ; 
//     std::string solutions_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/QSRPDE/" + pollutant ; 


//     // lambdas sequence 

//     std::vector<DVector<double>> lambdas_50;
//     std::vector<DVector<double>> lambdas_90;
//     std::vector<DVector<double>> lambdas_95;


//     // 50% 
//     for(double x = -11.0; x <= +2.0; x += 1.0) lambdas_50.push_back(SVector<1>(std::pow(10, x)));


//     // define spatial domain and regularizing PDE

//     MeshLoader<Mesh2D> domain("mesh_lombardia_coarse");        // ATT: da cambiare in base a type_space_mesh ma con if non funziona


    
//     // import data and locs from files
//     DMatrix<double> y; DMatrix<double> space_locs; DMatrix<double> time_locs; 
    
//     if(pollutant == "NO2"){
//         y = read_csv<double>(path_data + "/NO2_max.csv");   
//     } else{
//         y = read_csv<double>(path_data + "/PM10_2022.csv");  
//     }     

//     space_locs = read_csv<double>(path_data + "/locs.csv");
    
//     // check dimensions
//     std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
//     std::cout << "dim y " << y.rows() << " " << y.cols() << std::endl;

//     BlockFrame<double, int> df;
//     df.stack(OBSERVATIONS_BLK, y);
   
//     // define regularizing PDE in space 
//     auto Ld = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 , 1);
//     PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

        
//     for(double alpha : alphas){

//         unsigned int alpha_int = alpha*100; 
//         std::string alpha_string = std::to_string(alpha_int); 

//         std::cout << "----------------ALPHA = " << alpha_string << "----------------" << std::endl; 

//         QSRPDE<SpaceOnly> model(space_penalty, Sampling::pointwise, alpha);

//         // set model's data
//         model.set_spatial_locations(space_locs);

//         // Read lambda
//         double lambda = 1e-6; 
//         // std::ifstream fileLambdaS_opt(solutions_path + "/lambda_s_opt.csv");
//         // if(fileLambdaS_opt.is_open()){
//         //     fileLambdaS_opt >> lambda; 
//         //     fileLambdaS_opt.close();
//         // }

//         std::cout << "LambdaS = " << lambda << std::endl; 
//         model.set_lambda_D(lambda);
        
//         // set model's data
//         BlockFrame<double, int> df;
//         df.insert(OBSERVATIONS_BLK, y);
//         model.set_data(df);
//         // solve smoothing problem
//         model.init();
//         model.solve();

//         // Save C++ solution 
//         DMatrix<double> computedF = model.f();
//         const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filef(solutions_path +  "/f_1e-6.csv");
//         if (filef.is_open()){
//             filef << computedF.format(CSVFormatf);
//             filef.close();
//         }

//         DMatrix<double> computedFn = model.Psi()*model.f();
//         const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filefn(solutions_path + "/fn_1e-6.csv");
//         if (filefn.is_open()){
//             filefn << computedFn.format(CSVFormatfn);
//             filefn.close();
//         }

        
//     }

      
// }










// // -----------------
// // gcv  

// TEST(case_study_gcv_space_only_srpde, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

//     std::string pollutant = "NO2";   // PM10 NO2

//     const std::string type_locs = "all";   // "reduced" "all"

//     const std::string type_space_mesh = "coarse";   // coarse ---> CAMBIA MESH SOTTO !!! 
    
//     const std::string type_data = "/data_rescaled";  // "data_original"  ---> MODIFICA ANCHE GIU'lettura dati !!!

//     const std::string time_step = "by_week";   // "" "by_day" "by_week"

//     std::string gcv_type = "stochastic";   // "exact" "stochastic"  ---> MODIFICA ANCHE GIU'!
//     std::size_t seed = 438172;
//     unsigned int MC_run = 100; 

//     std::string est_type = "quantile";    // mean quantile
//     std::vector<double> alphas = {0.5};     //  {0.5, 0.9, 0.95};     

//     std::string lambda_selection_type = "gcv"; // _smooth_eps1e-1.5

//     // // Marco 
//     // std::string path_data = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant; 
    
//     // std::string solutions_path; 
//     // if(est_type == "mean")
//     //     solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/STRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 
//     // if(est_type == "quantile")
//     //     solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/QSTRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 

//     // Ilenia 
//     std::string path_data = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant + "/week_52"  ; 
//     std::string solutions_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/SRPDE/" + pollutant ; 


//     // lambdas sequence 

//     std::vector<DVector<double>> lambdas_50;
//     std::vector<DVector<double>> lambdas_90;
//     std::vector<DVector<double>> lambdas_95;


//     // 50% 
//     for(double x = -5.0; x <= -1.0; x += 0.1) lambdas_50.push_back(SVector<1>(std::pow(10, x)));


//     // define spatial domain and regularizing PDE

//     MeshLoader<Mesh2D> domain("mesh_lombardia_coarse");        // ATT: da cambiare in base a type_space_mesh ma con if non funziona


    
//     // import data and locs from files
//     DMatrix<double> y; DMatrix<double> space_locs; DMatrix<double> time_locs; 
    
//     if(pollutant == "NO2"){
//         y = read_csv<double>(path_data + "/NO2_max.csv");   
//     } else{
//         y = read_csv<double>(path_data + "/PM10_2022.csv");  
//     }     

//     space_locs = read_csv<double>(path_data + "/locs.csv");
    
//     // check dimensions
//     std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
//     std::cout << "dim y " << y.rows() << " " << y.cols() << std::endl;

//     BlockFrame<double, int> df;
//     df.stack(OBSERVATIONS_BLK, y);
   
//     // define regularizing PDE in space 
//     auto Ld = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 , 1);
//     PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

        

//         SRPDE model(space_penalty, Sampling::pointwise);

//         // set model's data
//         model.set_spatial_locations(space_locs);
        

//         model.set_data(df);
//         model.init();

//         // define GCV function and grid of \lambda_D values

//         // stochastic
//         auto GCV = model.gcv<StochasticEDF>(MC_run, seed);
//         // // exact
//         // auto GCV = model.gcv<ExactEDF>();

//         // optimize GCV
//         Grid<fdapde::Dynamic> opt;
//         opt.optimize(GCV, lambdas_50);
//         SVector<1> best_lambda = opt.optimum();

//         // Save Lambda opt
//         std::ofstream fileLambdaoptS(solutions_path +  "/lambda_s_opt.csv");
//         if(fileLambdaoptS.is_open()){
//             fileLambdaoptS << std::setprecision(16) << best_lambda[0];
//             fileLambdaoptS.close();
//         }

//         // Save GCV scores
//         std::ofstream fileGCV_scores(solutions_path + "/gcv_scores.csv");
//         for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
//             fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.gcvs()[i]) << "\n"; 

//         fileGCV_scores.close();

//         // Save edfs
//         std::ofstream fileEDF(solutions_path +  "/edfs.csv");
//         for(std::size_t i = 0; i < GCV.edfs().size(); ++i) 
//             fileEDF << std::setprecision(16) << GCV.edfs()[i] << "\n"; 
//         fileEDF.close();


      
// }



// // run

// TEST(case_study_run_space_only, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

//     std::string pollutant = "NO2";   // PM10 NO2

//     const std::string type_locs = "all";   // "reduced" "all"

//     const std::string type_space_mesh = "coarse";   // coarse ---> CAMBIA MESH SOTTO !!! 
    
//     const std::string type_data = "/data_rescaled";  // "data_original"  ---> MODIFICA ANCHE GIU'lettura dati !!!

//     const std::string time_step = "by_week";   // "" "by_day" "by_week"

//     std::string gcv_type = "stochastic";   // "exact" "stochastic"  ---> MODIFICA ANCHE GIU'!
//     std::size_t seed = 438172;
//     unsigned int MC_run = 100; 

//     std::string est_type = "quantile";    // mean quantile
//     std::vector<double> alphas = {0.5};     //  {0.5, 0.9, 0.95};     

//     std::string lambda_selection_type = "gcv"; // _smooth_eps1e-1.5

//     // // Marco 
//     // std::string path_data = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant; 
    
//     // std::string solutions_path; 
//     // if(est_type == "mean")
//     //     solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/STRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 
//     // if(est_type == "quantile")
//     //     solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/QSTRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 

//     // Ilenia 
//     std::string path_data = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant + "/week_52"  ; 
//     std::string solutions_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/SRPDE/" + pollutant ; 


//     // lambdas sequence 

//     std::vector<DVector<double>> lambdas_50;
//     std::vector<DVector<double>> lambdas_90;
//     std::vector<DVector<double>> lambdas_95;


//     // 50% 
//     for(double x = -11.0; x <= +2.0; x += 1.0) lambdas_50.push_back(SVector<1>(std::pow(10, x)));


//     // define spatial domain and regularizing PDE

//     MeshLoader<Mesh2D> domain("mesh_lombardia_coarse");        // ATT: da cambiare in base a type_space_mesh ma con if non funziona


    
//     // import data and locs from files
//     DMatrix<double> y; DMatrix<double> space_locs; DMatrix<double> time_locs; 
    
//     if(pollutant == "NO2"){
//         y = read_csv<double>(path_data + "/NO2_max.csv");   
//     } else{
//         y = read_csv<double>(path_data + "/PM10_2022.csv");  
//     }     

//     space_locs = read_csv<double>(path_data + "/locs.csv");
    
//     // check dimensions
//     std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
//     std::cout << "dim y " << y.rows() << " " << y.cols() << std::endl;

//     BlockFrame<double, int> df;
//     df.stack(OBSERVATIONS_BLK, y);
   
//     // define regularizing PDE in space 
//     auto Ld = -laplacian<FEM>();
//     DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 , 1);
//     PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

        

//         SRPDE model(space_penalty, Sampling::pointwise);

//         // set model's data
//         model.set_spatial_locations(space_locs);

//         // Read lambda
//         double lambda; 
//         std::ifstream fileLambdaS_opt(solutions_path + "/lambda_s_opt.csv");
//         if(fileLambdaS_opt.is_open()){
//             fileLambdaS_opt >> lambda; 
//             fileLambdaS_opt.close();
//         }

//         std::cout << "LambdaS = " << lambda << std::endl; 
//         model.set_lambda_D(lambda);
        
//         // set model's data
//         df.insert(OBSERVATIONS_BLK, y);
//         model.set_data(df);
//         // solve smoothing problem
//         model.init();
//         model.solve();

//         // Save C++ solution 
//         DMatrix<double> computedF = model.f();
//         const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filef(solutions_path +  "/f.csv");
//         if (filef.is_open()){
//             filef << computedF.format(CSVFormatf);
//             filef.close();
//         }

//         DMatrix<double> computedFn = model.Psi()*model.f();
//         const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filefn(solutions_path + "/fn.csv");
//         if (filefn.is_open()){
//             filefn << computedFn.format(CSVFormatfn);
//             filefn.close();
//         }


      
// }
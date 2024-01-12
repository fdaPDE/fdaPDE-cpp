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
using fdapde::core::MatrixDataWrapper;
using fdapde::core::PDE;
using fdapde::core::VectorDataWrapper;

#include "../../fdaPDE/models/sampling_design.h"
#include "../../fdaPDE/models/regression/strpde.h"
#include "../../fdaPDE/models/regression/qsrpde.h"
using fdapde::models::STRPDE;
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


// gcv  

TEST(case_study_gcv, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

    std::string pollutant = "PM10";   // PM10 NO2

    std::string gcv_type = "stochastic";   // "exact" "stochastic"  ---> MODIFICA ANCHE GIU'!
    std::size_t seed = 438172;
    unsigned int MC_run = 100; 

    std::string est_type = "quantile";    // mean quantile
    std::vector<double> alphas = {0.1, 0.5, 0.9}; 

    // Marco 
    std::string path_data = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant; 
    std::string solutions_path; 
    if(est_type == "mean")
        solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/STRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 
    if(est_type == "quantile")
        solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/QSTRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 

    // Ilenia 
    // std::string path_data = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/NO2";

    const std::string type_locs = "all";   // "reduced" "all"

    // lambdas sequence 
    std::vector<DVector<double>> lambdas_d_t;
    std::vector<double> lambdas_d;
    std::vector<double> lambdas_t;

    if(est_type == "mean"){
        for(double xs = -3.0; xs <= -1.5; xs += 0.05)
            lambdas_d.push_back(std::pow(10,xs));

        for(double xt = 0.; xt <= 1.5; xt += 0.5)
            lambdas_t.push_back(std::pow(10,xt));    

        for(auto i = 0; i < lambdas_d.size(); ++i)
            for(auto j = 0; j < lambdas_t.size(); ++j) 
                lambdas_d_t.push_back(SVector<2>(lambdas_d[i], lambdas_t[j]));
    }
    if(est_type == "quantile"){
        for(double xs = -3.0; xs <= -2.0; xs += 1.0)
            lambdas_d.push_back(std::pow(10,xs));

        for(double xt = 0.; xt <= 0.2; xt += 0.5)
            lambdas_t.push_back(std::pow(10,xt));    

        for(auto i = 0; i < lambdas_d.size(); ++i)
            for(auto j = 0; j < lambdas_t.size(); ++j) 
                lambdas_d_t.push_back(SVector<2>(lambdas_d[i], lambdas_t[j]));
    }


    // define temporal domain
    unsigned int M = 50; 
    std::string M_string = std::to_string(M);
    double tf = 364.0;    // final time 
    Mesh<1, 1> time_mesh(0, tf, M-1);
    // define spatial domain and regularizing PDE
    MeshLoader<Mesh2D> domain("mesh_lombardia_coarse");

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
            y = read_csv<double>(path_data + "/NO2_max_2022_Cpp.csv");
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
        std::ofstream fileLambda_S_Seq(solutions_path + "/" + type_locs + "/M_" + M_string + "/lambdas_S_seq.csv");
        for(std::size_t i = 0; i < lambdas_d.size(); ++i) 
            fileLambda_S_Seq << std::setprecision(16) << lambdas_d[i] << "\n"; 
        fileLambda_S_Seq.close();

        for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
            std::cout << lambdas_t[i] << "\n"; 

        std::ofstream fileLambda_T_Seq(solutions_path + "/" + type_locs + "/M_" + M_string + "/lambdas_T_seq.csv");
        for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
            fileLambda_T_Seq << std::setprecision(16) << lambdas_t[i] << "\n"; 
        fileLambda_T_Seq.close();

        // Save Lambda opt
        std::ofstream fileLambdaoptS(solutions_path + "/" + type_locs + "/M_" + M_string + "/lambda_s_opt.csv");
        if(fileLambdaoptS.is_open()){
            fileLambdaoptS << std::setprecision(16) << best_lambda[0];
            fileLambdaoptS.close();
        }
        std::ofstream fileLambdaoptT(solutions_path + "/" + type_locs + "/M_" + M_string + "/lambda_t_opt.csv");
        if (fileLambdaoptT.is_open()){
            fileLambdaoptT << std::setprecision(16) << best_lambda[1];
            fileLambdaoptT.close();
        }
        // Save GCV scores
        std::ofstream fileGCV_scores(solutions_path + "/" + type_locs + "/M_" + M_string + "/gcv_scores.csv");
        for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
            fileGCV_scores << std::setprecision(16) << GCV.gcvs()[i] << "\n"; 

        fileGCV_scores.close();

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

            model.set_exact_gcv(true);   // true: no lisciamento curva gcv
            
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
            std::ofstream fileLambda_S_Seq(solutions_path + "/" + type_locs + "/M_" + M_string + "/alpha_" + alpha_string + "/lambdas_S_seq.csv");
            for(std::size_t i = 0; i < lambdas_d.size(); ++i) 
                fileLambda_S_Seq << std::setprecision(16) << lambdas_d[i] << "\n"; 
            fileLambda_S_Seq.close();

            for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
                std::cout << lambdas_t[i] << "\n"; 

            std::ofstream fileLambda_T_Seq(solutions_path + "/" + type_locs + "/M_" + M_string + "/alpha_" + alpha_string + "/lambdas_T_seq.csv");
            for(std::size_t i = 0; i < lambdas_t.size(); ++i) 
                fileLambda_T_Seq << std::setprecision(16) << lambdas_t[i] << "\n"; 
            fileLambda_T_Seq.close();

            // Save Lambda opt
            std::ofstream fileLambdaoptS(solutions_path + "/" + type_locs + "/M_" + M_string + "/alpha_" + alpha_string + "/lambda_s_opt.csv");
            if(fileLambdaoptS.is_open()){
                fileLambdaoptS << std::setprecision(16) << best_lambda[0];
                fileLambdaoptS.close();
            }
            std::ofstream fileLambdaoptT(solutions_path + "/" + type_locs + "/M_" + M_string + "/alpha_" + alpha_string + "/lambda_t_opt.csv");
            if (fileLambdaoptT.is_open()){
                fileLambdaoptT << std::setprecision(16) << best_lambda[1];
                fileLambdaoptT.close();
            }
            // Save GCV scores
            std::ofstream fileGCV_scores(solutions_path + "/" + type_locs + "/M_" + M_string + "/alpha_" + alpha_string + "/gcv_scores.csv");
            for(std::size_t i = 0; i < GCV.gcvs().size(); ++i) 
                fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.gcvs()[i]) << "\n"; 

            fileGCV_scores.close();

        }

        }

      
}



// // run  

// TEST(case_study_run, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {

//     std::string pollutant = "PM10";   // PM10 NO2

//     std::string gcv_type = "stochastic";   // "exact" "stochastic"  

//     std::string est_type = "mean";    // mean quantile
//     std::vector<double> alphas = {0.1, 0.5, 0.9}; 

//     // Marco 
//     std::string path_data = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/" + pollutant; 
//     std::string solutions_path; 
//     if(est_type == "mean")
//         solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/STRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 
//     if(est_type == "quantile")
//         solutions_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/ARPA/Lombardia/QSTRPDE/" + pollutant + "/gcv_Cpp_" + gcv_type; 

//     // Ilenia 
//     // std::string path_data = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared/case_study/ARPA/Lombardia/dati_Cpp/NO2";

//     const std::string type_locs = "all";   // "reduced" "all"

//     // define temporal domain
//     unsigned int M = 50; 
//     std::string M_string = std::to_string(M);
//     double tf = 364.0;    // final time 
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
//             y = read_csv<double>(path_data + "/NO2_max_2022_Cpp.csv");
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

//     // Read optima lambdas 
//     double lambda_T; double lambda_S; 
//     std::ifstream fileLambdaT_opt(solutions_path + "/" + type_locs + "/M_" + M_string + "/lambda_t_opt.csv");
//     if(fileLambdaT_opt.is_open()){
//         fileLambdaT_opt >> lambda_T; 
//         fileLambdaT_opt.close();
//     }
//     std::ifstream fileLambdaS_opt(solutions_path + "/" + type_locs + "/M_" + M_string + "/lambda_s_opt.csv");
//     if(fileLambdaS_opt.is_open()){
//         fileLambdaS_opt >> lambda_S; 
//         fileLambdaS_opt.close();
//     }

//     std::cout << "lambda S " << lambda_S << std::endl;
//     std::cout << "lambda T " << lambda_T << std::endl;

//     if(est_type == "mean"){

//         STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::pointwise);
    
//         // set model's data
//         model.set_spatial_locations(space_locs);
//         model.set_temporal_locations(time_locs);

//         model.set_lambda_D(lambda_S);
//         model.set_lambda_T(lambda_T);
        
//         model.set_data(df);

//         model.init();
//         model.solve();

//         // Save C++ solution 
//         DMatrix<double> computedF = model.f();
//         const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filef(solutions_path + "/" + type_locs + "/M_" + M_string + "/f.csv");
//         if (filef.is_open()){
//             filef << computedF.format(CSVFormatf);
//             filef.close();
//         }

//         DMatrix<double> computedFn = model.Psi()*model.f();
//         const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filefn(solutions_path + "/" + type_locs + "/M_" + M_string + "/fn.csv");
//         if (filefn.is_open()){
//             filefn << computedFn.format(CSVFormatfn);
//             filefn.close();
//         }

 

//     }

//     if(est_type == "quantile"){
        
//         for(double alpha : alphas){

//             unsigned int alpha_int = alpha*100; 
//             std::string alpha_string = std::to_string(alpha_int); 

//             QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);
    
//             // set model's data
//             model.set_spatial_locations(space_locs);
//             model.set_temporal_locations(time_locs);

//             model.set_exact_gcv(true);   // true: no lisciamento curva gcv
//             model.set_lambda_D(lambda_S);
//             model.set_lambda_T(lambda_T);
            
//             model.set_data(df);

//             model.init();
//             model.solve();

//             // Save C++ solution 
//             DMatrix<double> computedF = model.f();
//             const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//             std::ofstream filef(solutions_path + "/" + type_locs + "/M_" + M_string + "/f_" + alpha_string + ".csv");
//             if (filef.is_open()){
//                 filef << computedF.format(CSVFormatf);
//                 filef.close();
//             }

//             DMatrix<double> computedFn = model.Psi()*model.f();
//             const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//             std::ofstream filefn(solutions_path + "/" + type_locs + "/M_" + M_string + "/fn_" + alpha_string + ".csv");
//             if (filefn.is_open()){
//                 filefn << computedFn.format(CSVFormatfn);
//                 filefn.close();
//             }



//         }

//     }



        
// }






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




//------------------------------------------------------------OLD---------------------------------------------------

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
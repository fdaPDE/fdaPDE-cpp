#include <cstddef>
#include <gtest/gtest.h>   // testing framework

#include <fdaPDE/core.h>
using fdapde::core::advection;
using fdapde::core::diffusion;
using fdapde::core::dt;
using fdapde::core::FEM;
using fdapde::core::fem_order;
using fdapde::core::laplacian;
using fdapde::core::MatrixDataWrapper;
using fdapde::core::PDE;
using fdapde::core::VectorDataWrapper;

#include "../../fdaPDE/models/regression/strpde.h"
#include "../../fdaPDE/models/sampling_design.h"
using fdapde::models::STRPDE;
using fdapde::models::SpaceTimeSeparable;
using fdapde::models::SpaceTimeParabolic;
using fdapde::models::GeoStatMeshNodes;
using fdapde::models::GeoStatLocations;
using fdapde::models::Areal;
using fdapde::models::MonolithicSolver;
using fdapde::models::IterativeSolver;

#include "../../fdaPDE/calibration/gcv.h"
using fdapde::calibration::ExactEDF;
using fdapde::calibration::ExactGCV;
using fdapde::calibration::GCV;
using fdapde::calibration::StochasticEDF;
using fdapde::core::Grid; 

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_mtx;
using fdapde::testing::read_csv;


// STRPDE 
TEST(strpde_case_study, laplacian_nonparametric_samplingatlocations_timelocations_separable_monolithic_missingdata) {
    
    const std::string path_data = "C:/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/our_marco/data"; 

    std::cout << "here 1 " << std::endl; 
    // define temporal domain
    DVector<double> time_mesh;
    unsigned int M = 12; 
    time_mesh.resize(M);
    for (std::size_t i = 0; i < M; ++i) time_mesh[i] = i+1;
    // define spatial domain and regularizing PDE
    MeshLoader<Mesh2D> domain("/case_study/mesh_reshaped");
    // import data from files
    DMatrix<double> time_locs  = read_csv<double>(path_data + "/time_locs.csv");
    DMatrix<double> space_locs = read_csv<double>(path_data + "/space_locs.csv");
    DMatrix<double> y          = read_csv<double>(path_data + "/y.csv");
    std::cout << "here 2" << std::endl;
    std::cout << "dim space loc " << space_locs.rows() << " " << space_locs.cols() << std::endl;
    std::cout << "dim time loc " << time_locs.rows() << " " << time_locs.cols() << std::endl;
    std::cout << "dim y " << y.rows() << " " << y.cols() << std::endl;
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.rows(), 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    
    // ----GCV---

    std::cout << "here 3" << std::endl;
    std::vector<SVector<2>> lambdas_d_t;
    for(double xs = -10.0; xs <= -8.0; xs +=0.1)
        for(double xt = -8.0; xt <= -5.0; xt +=0.5) 
            lambdas_d_t.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));
    
    // define model
    STRPDE<decltype(problem), SpaceTimeSeparable, GeoStatLocations, MonolithicSolver> model_gcv(problem, time_mesh);

    model_gcv.set_spatial_locations(space_locs);
    model_gcv.set_temporal_locations(time_locs);
    std::cout << "here 4" << std::endl;

    BlockFrame<double, int> df_gcv;
    df_gcv.stack(OBSERVATIONS_BLK, y);
    model_gcv.set_data(df_gcv);
    std::cout << "here 5" << std::endl;
    model_gcv.init();
    std::cout << "here 6" << std::endl;

    // define GCV function and grid of \lambda_D values
    GCV<decltype(model_gcv), ExactEDF<decltype(model_gcv)>> GCV(model_gcv);
    ScalarField<2, decltype(GCV)> obj(GCV);  
    // optimize GCV
    Grid<2> opt;
    opt.optimize(obj, lambdas_d_t);
    std::cout << "here 7" << std::endl;
    SVector<2> best_lambda = opt.optimum();
    std::cout << "here 8" << std::endl;

    std::string solutions_path = "C:/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/Thesis_shared/case_study/our_marco/results"; 
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


    // ----RUN---

    double lambda_D = read_csv<double>(solutions_path + "/lambda_s_opt.csv")(0,0); 
    double lambda_T = read_csv<double>(solutions_path + "/lambda_t_opt.csv")(0,0);

    STRPDE<decltype(problem), SpaceTimeSeparable, GeoStatLocations, MonolithicSolver> model(problem, time_mesh);
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

    // Save C++ solution 
    DMatrix<double> computedF = model.f();
    const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
    std::ofstream filef(solutions_path + "/f_STRPDE.csv");
    if (filef.is_open()){
        filef << computedF.format(CSVFormatf);
        filef.close();
    }

}
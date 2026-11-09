// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.

#include <fdaPDE/fdapde.h>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include "my_fun.h"

#include <chrono>
#include <fstream>

using namespace fdapde;
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using sparse_matrix_t = Eigen::SparseMatrix<double>;

struct MethodSpec {
    std::string name;
    bool enabled;
};

struct CriterionSpec {
    std::string name;
    std::string seq_type;
    bool enabled;
    std::pair<double, double> lasco_bounds;
    double lasco_step;
    std::map<std::string, std::pair<double, double>> fine_bounds;
    double fine_step;
};

int test_2();

int main() {
    test_2();
    return 0;
}

int test_2() {

    // ==================== 1. PARAMETRI GENERALI SIMULAZIONE ====================

    const std::string trial_number = "3";
    const unsigned int L = 10;
    const double domain_area = static_cast<double>(L * L);

    const unsigned int sim_start = 1;
    const unsigned int n_sim = 50;

    const unsigned int n_folds = 30;

    const std::string R_path = "../../OneDrive - Politecnico di Milano/Corsi/PhD/Codice/models/COSP/space-only/Test_2/trial_" + trial_number + "/";


    // ==================== 2. PARAMETRI ALGORITMO FPIRLS ========================

    const bool init_with_residuals = false;

    ModelConfig cfg;
    cfg.domain_area = domain_area;
    cfg.random_seed = 37;
    cfg.random_seed_p = 17;
    cfg.random_seed_A = 23;
    cfg.tol = 1e-4;
    cfg.max_iter = 50;
    cfg.verbose = false;
    cfg.init_with_residuals = init_with_residuals;


    // ==================== 3. SELEZIONE MODELLI DA ESEGUIRE =====================

    const bool hetero_fit     = false;
    const bool hetero2_fit    = false;
    const bool hetero3_fit    = false;
    const bool hetero4_fit    = false;
    const bool hetero5_fit    = false;
    const bool hetero6_fit    = false;
    const bool hetero7_fit    = false;
    const bool hetero8_fit    = true;
    const bool naive_fit      = true;
    const bool only_point_fit = true;
    const bool only_area_fit  = true;


    // ==================== 4. SELEZIONE CRITERI E GRIGLIE LAMBDA ================

    const bool run_gcv   = true;  const std::string gcv_seq_type   = "fine"; // "lasco" o "fine"
    const bool run_cv    = false;  const std::string cv_seq_type    = "fine";
    const bool run_aicc  = false;  const std::string aicc_seq_type  = "fine";
    const bool run_gfaic = false;  const std::string gfaic_seq_type = "fine";
    const bool run_rcv   = false;  const std::string rcv_seq_type   = "fine";
    const bool run_fpe   = false;  const std::string fpe_seq_type   = "fine";

    const std::pair<double, double> lasco_bounds_default = {-8.0, 4.0};
    const double lasco_step_default = 0.50;
    const double fine_step_default  = 0.10;
    const double fine_half_window_xu = 1.0; // per hetero8

    std::map<std::string, std::pair<double, double>> fine_bounds_gcv = {
        {"hetero",     {-8.0,  2.0}},
        {"hetero2",    {-8.0,  2.0}},
        {"hetero3",    {-8.0,  2.0}},
        {"hetero4",    {-8.0,  2.0}},
        {"hetero5",    {-8.0,  2.0}},
        {"hetero6",    {-8.0,  2.0}},
        {"hetero7",    {-8.0,  2.0}},
        {"hetero8",    {-8.0,  2.0}},
        {"naive",      {-8.0,  2.0}},
        {"only_point", {-8.0,  2.0}},
        {"only_area",  {-8.0,  2.0}}
    };

    std::map<std::string, std::pair<double, double>> fine_bounds_cv = {
        {"hetero",     {-8.0,  2.0}},
        {"hetero2",    {-8.0,  2.0}},
        {"hetero3",    {-8.0,  2.0}},
        {"hetero4",    {-8.0,  2.0}},
        {"hetero5",    {-8.0,  2.0}},
        {"hetero6",    {-8.0,  2.0}},
        {"hetero7",    {-8.0,  2.0}},
        {"hetero8",    {-8.0,  2.0}},
        {"naive",      {-8.0,  2.0}},
        {"only_point", {-8.0,  2.0}},
        {"only_area",  {-8.0,  2.0}}
    };

    std::map<std::string, std::pair<double, double>> fine_bounds_aicc  = fine_bounds_gcv;
    std::map<std::string, std::pair<double, double>> fine_bounds_gfaic = fine_bounds_gcv;
    std::map<std::string, std::pair<double, double>> fine_bounds_rcv   = fine_bounds_gcv;
    std::map<std::string, std::pair<double, double>> fine_bounds_fpe   = fine_bounds_gcv;

    // ==================== POPOLAMENTO AUTOMATICO TABELLE =======================

    std::vector<MethodSpec> methods = {
        {"hetero",     hetero_fit},
        {"hetero2",    hetero2_fit},
        {"hetero3",    hetero3_fit},
        {"hetero4",    hetero4_fit},
        {"hetero5",    hetero5_fit},
        {"hetero6",    hetero6_fit},
        {"hetero7",    hetero7_fit},
        {"hetero8",    hetero8_fit},
        {"naive",      naive_fit},
        {"only_point", only_point_fit},
        {"only_area",  only_area_fit}
    };

    std::vector<CriterionSpec> criteria = {
        {"gcv",   gcv_seq_type,   run_gcv,   lasco_bounds_default, lasco_step_default, fine_bounds_gcv,   fine_step_default},
        {"cv",    cv_seq_type,    run_cv,    lasco_bounds_default, lasco_step_default, fine_bounds_cv,    fine_step_default},
        {"aicc",  aicc_seq_type,  run_aicc,  lasco_bounds_default, lasco_step_default, fine_bounds_aicc,  fine_step_default},
        {"gfaic", gfaic_seq_type, run_gfaic, lasco_bounds_default, lasco_step_default, fine_bounds_gfaic, fine_step_default},
        {"rcv",   rcv_seq_type,   run_rcv,   lasco_bounds_default, lasco_step_default, fine_bounds_rcv,   fine_step_default},
        {"fpe",   fpe_seq_type,   run_fpe,   lasco_bounds_default, lasco_step_default, fine_bounds_fpe,   fine_step_default}
    };


    // ==================== SETUP MESH E FISICA ==================================

    const std::string mesh_dir = "my_data/mesh/square_" + std::to_string(L) + "x" + std::to_string(L) + "/";
    Eigen::Matrix<double, Dynamic, Dynamic> points = read_csv<double>(mesh_dir + "points.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> elements = read_csv<int>(mesh_dir + "elements.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> boundary = read_csv<int>(mesh_dir + "boundary.csv").as_matrix();
    elements.array() -= 1;
    Triangulation<2, 2> D(points, elements, boundary);

    std::cout << "Define the physics..." << std::endl;
    FeSpace Vh(D, P1<1>);
    TrialFunction f_trial(Vh);
    TestFunction v_test(Vh);
    auto a = integral(D)(dot(grad(f_trial), grad(v_test)));
    ZeroField<2> u_fun;
    auto F = integral(D)(u_fun * v_test);

    sparse_matrix_t R1 = a.assemble();
    sparse_matrix_t R0 = integral(D)(f_trial * v_test).assemble();
    vector_t u = F.assemble();
    unsigned int n_dofs = Vh.n_dofs();


    // ==================== CICLO SIMULAZIONI ====================================
    
    for (unsigned int sim = sim_start; sim <= n_sim; ++sim) {
        
        std::cout << "\nSIMULATION " << sim << " / " << n_sim << std::endl;
        std::string sim_dir = R_path + "results/sim_" + std::to_string(sim) + "/";

        GeoFrame geo_data(D);
        auto& l_point = geo_data.insert_scalar_layer<POINT>("l1", R_path + "locs.csv");
        auto& l_areal = geo_data.insert_scalar_layer<POLYGON>("l2", R_path + "incidence_mat.csv");
        l_point.load_csv<double>(sim_dir + "response_point.csv");
        l_areal.load_csv<double>(sim_dir + "response_areal.csv");

        unsigned int n_points = l_point.rows();
        unsigned int n_areal = l_areal.rows();
        unsigned int n = n_points + n_areal;

        vector_t y_point = l_point.col<double>(0).as_matrix();
        vector_t y_areal = l_areal.col<double>(0).as_matrix();
        vector_t y(n);
        y << y_point, y_areal;

        std::vector<int> fold_ids(n);
        std::ifstream file_folds(sim_dir + "fold_ids.csv");
        for (unsigned int idx = 0; idx < n; ++idx) file_folds >> fold_ids[idx];
        file_folds.close();

        sparse_matrix_t Psi_point = internals::point_basis_eval(Vh, l_point.geometry<0>());
        auto [Psi_areal, D_areal] = internals::areal_basis_eval(Vh, l_areal.geometry<0>());

        std::vector<Eigen::Triplet<double>> triplets_Psi;
        for (int k = 0; k < Psi_point.outerSize(); ++k) {
            for (sparse_matrix_t::InnerIterator it(Psi_point, k); it; ++it) {
                triplets_Psi.emplace_back(it.row(), it.col(), it.value());
            }
        }
        for (int k = 0; k < Psi_areal.outerSize(); ++k) {
            for (sparse_matrix_t::InnerIterator it(Psi_areal, k); it; ++it) {
                triplets_Psi.emplace_back(it.row() + n_points, it.col(), it.value());
            }
        }
        sparse_matrix_t Psi(n, n_dofs);
        Psi.setFromTriplets(triplets_Psi.begin(), triplets_Psi.end());

        std::vector<Eigen::Triplet<double>> triplets_D;
        for (unsigned int i = 0; i < n_points; ++i) triplets_D.emplace_back(i, i, 1.0);
        for (unsigned int i = 0; i < n_areal; ++i) triplets_D.emplace_back(i + n_points, i + n_points, D_areal[i]);
        sparse_matrix_t D_matrix(n, n);
        D_matrix.setFromTriplets(triplets_D.begin(), triplets_D.end());


        double elapsed_s; 

        for (const auto& m : methods) {

            if (!m.enabled) continue;
            std::cout << "Running " << m.name << " selection..." << std::endl;
            std::string out_dir = sim_dir + m.name + "/";

            if (m.name == "hetero8") {

                auto t_start = std::chrono::high_resolution_clock::now();
                std::vector<double> final_lambda_seq;
                auto res_xu = run_hetero8_srpde(n_points, n_areal, n_dofs, D_matrix, Psi, R1, R0, y, u,
                                                lasco_bounds_default, lasco_step_default,
                                                fine_half_window_xu, fine_step_default,
                                                cfg, final_lambda_seq);

                auto t_end = std::chrono::high_resolution_clock::now();
                elapsed_s = std::chrono::duration<double>(t_end - t_start).count();                                

                std::string tag = "gcv_" + gcv_seq_type;
                write_csv(out_dir + "lambdas_seq_" + tag + ".csv", final_lambda_seq);
                write_csv(out_dir + "gcv_scores_" + gcv_seq_type + ".csv", res_xu.gcv_scores);
                save_scalar(out_dir + "lambda_opt_" + tag + ".csv", res_xu.lambda_opt);
                write_csv(out_dir + "lambda_history_gcv.csv", res_xu.lambda_history);
                save_scalar(out_dir + "time_" + tag + ".csv", elapsed_s);  // per hetero8 questo e' il tempo totale.

                save_fit_results(out_dir, "gcv", res_xu, false);

            } else {
                for(const auto& crit : criteria) {

                    if (!crit.enabled) continue;
                    std::cout << "  --- Criterion: " << crit.name << " (" << crit.seq_type << ")..." << std::endl;
                    
                    auto t_start = std::chrono::high_resolution_clock::now();

                    std::vector<double> lambda_seq;
                    if (crit.seq_type == "lasco" || crit.fine_bounds.find(m.name) == crit.fine_bounds.end()) {
                        lambda_seq = build_log10_seq(crit.lasco_bounds.first, crit.lasco_bounds.second, crit.lasco_step);
                    } else {
                        auto bounds = crit.fine_bounds.at(m.name);
                        lambda_seq = build_log10_seq(bounds.first, bounds.second, crit.fine_step);
                    }

                    vector_t scores(lambda_seq.size());
                    for (size_t i = 0; i < lambda_seq.size(); ++i) {
                        scores[i] = evaluate_criterion(crit.name, m.name, n_points, n_areal, n_dofs,
                                                    D_matrix, Psi, R1, R0, y, u, lambda_seq[i],
                                                    fold_ids, n_folds, cfg);
                    }

                    Eigen::Index min_idx;
                    scores.minCoeff(&min_idx);
                    double best_lambda = lambda_seq[min_idx];

                    auto t_end = std::chrono::high_resolution_clock::now();
                    elapsed_s = std::chrono::duration<double>(t_end - t_start).count(); 
  

                    std::string tag = crit.name + "_" + crit.seq_type;
                    write_csv(out_dir + "lambdas_seq_" + tag + ".csv", lambda_seq);
                    write_csv(out_dir + crit.name + "_scores_" + crit.seq_type + ".csv", scores);
                    save_scalar(out_dir + "lambda_opt_" + tag + ".csv", best_lambda);
                    save_scalar(out_dir + "time_selection_" + tag + ".csv", elapsed_s);
                }

            }

                      
        }

    }
    return 0;
}
// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.

#include <fdaPDE/fdapde.h>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <cmath>
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

int run_pm10_cosp();

int main() {
    run_pm10_cosp();
    return 0;
}

int run_pm10_cosp() {

    const std::string resolution_CAMS = "_res04";  // "_res01" "_res04" 


    // ==================== 1. PERCORSI E PARAMETRI GENERALI =====================

    const std::string R_path = "../../OneDrive - Politecnico di Milano/Corsi/PhD/Codice/case_studies/cosp/pm10/";
    const std::string mesh_dir = R_path + "mesh/convex_hull/canotto_fine/";
    const std::string data_dir = R_path + "data" + resolution_CAMS + "/";
    const std::string results_dir = R_path + "results" + resolution_CAMS + "/";
    const unsigned int n_folds = read_scalar(data_dir + "n_folds.csv"); 
    const double domain_area = read_scalar(data_dir + "domain_area.csv");

    // ==================== 2. SETUP MESH E FISICA ===============================

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

    // ==================== 3. PARAMETRI ALGORITMO FPIRLS ========================

    ModelConfig cfg;
    cfg.domain_area = domain_area;
    cfg.random_seed = 37;
    cfg.random_seed_p = 17;
    cfg.random_seed_A = 23;
    cfg.tol = 1e-4;
    cfg.max_iter = 50;
    cfg.verbose = false;
    cfg.init_with_residuals = false;

    // ==================== 4. SELEZIONE MODELLI E GRIGLIA GCV ===================

    const bool hetero8_fit    = true;
    const bool naive_fit      = true;
    const bool only_point_fit = true;
    const bool only_area_fit  = true;

    const bool compute_kfold_error = true; 

    const std::string gcv_seq_type = "fine";
    const std::pair<double, double> lasco_bounds_default = {-9.0, 2.0};
    const double lasco_step_default = 0.50;
    const double fine_step_default  = 0.05;
    const double fine_half_window_xu = 1.0;

    std::map<std::string, std::pair<double, double>> fine_bounds_gcv = {
        {"naive",      {-7.0, -5.0}},
        {"only_point", {-6.0, -4.0}},
        {"only_area",  {-13.0, -10.0}}
    };

    std::vector<MethodSpec> methods = {
        {"hetero8",    hetero8_fit},
        {"naive",      naive_fit},
        {"only_point", only_point_fit},
        {"only_area",  only_area_fit}
    };

    // ==================== 5. CARICAMENTO DATI PM10 =============================

    GeoFrame geo_data(D);
    auto& l_point = geo_data.insert_scalar_layer<POINT>("l1", data_dir + "locs.csv");
    auto& l_areal = geo_data.insert_scalar_layer<POLYGON>("l2", data_dir + "incidence_mat.csv");
    l_point.load_csv<double>(data_dir + "response_point.csv");
    l_areal.load_csv<double>(data_dir + "response_areal.csv");

    unsigned int n_points = l_point.rows();
    unsigned int n_areal = l_areal.rows();
    unsigned int n = n_points + n_areal;

    vector_t y_point = l_point.col<double>(0).as_matrix();
    vector_t y_areal = l_areal.col<double>(0).as_matrix();
    vector_t y(n);
    y << y_point, y_areal;

    std::vector<int> fold_ids(n);
    std::ifstream file_folds(data_dir + "fold_ids.csv");
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

    // ==================== 6. SELEZIONE GCV, FIT E K-FOLD CV ====================

    std::string tag = "gcv_" + gcv_seq_type;

    for (const auto& m : methods) {

         /**************************** SELECTION AND RUN *******************************/

        if (!m.enabled) continue;
        std::cout << "\nRunning model: " << m.name << " (GCV)..." << std::endl;
        std::string out_dir = results_dir + m.name + "/";
        bool is_single_variance = (m.name == "naive" || m.name == "only_point" || m.name == "only_area");

        double best_lambda = 0.0;
        std::vector<double> lambda_seq;

        if (m.name == "hetero8") {
            auto t_start = std::chrono::high_resolution_clock::now();
            std::vector<double> final_lambda_seq;
            auto res_xu = run_hetero8_srpde(n_points, n_areal, n_dofs, D_matrix, Psi, R1, R0, y, u,
                                            lasco_bounds_default, lasco_step_default,
                                            fine_half_window_xu, fine_step_default,
                                            cfg, final_lambda_seq);

            auto t_end = std::chrono::high_resolution_clock::now();
            double elapsed_s = std::chrono::duration<double>(t_end - t_start).count();

            best_lambda = res_xu.lambda_opt;
            write_csv(out_dir + "lambdas_seq_" + tag + ".csv", final_lambda_seq);
            write_csv(out_dir + "gcv_scores_" + gcv_seq_type + ".csv", res_xu.gcv_scores);
            save_scalar(out_dir + "lambda_opt_" + tag + ".csv", best_lambda);
            write_csv(out_dir + "lambda_history_gcv.csv", res_xu.lambda_history);
            save_scalar(out_dir + "time_" + tag + ".csv", elapsed_s);

            save_fit_results(out_dir, "gcv", res_xu, false);
        } else {
            auto t_sel_start = std::chrono::high_resolution_clock::now();

            if (gcv_seq_type == "lasco" || fine_bounds_gcv.find(m.name) == fine_bounds_gcv.end()) {
                lambda_seq = build_log10_seq(lasco_bounds_default.first, lasco_bounds_default.second, lasco_step_default);
            } else {
                auto bounds = fine_bounds_gcv.at(m.name);
                lambda_seq = build_log10_seq(bounds.first, bounds.second, fine_step_default);
            }

            vector_t scores(lambda_seq.size());
            for (size_t i = 0; i < lambda_seq.size(); ++i) {
                scores[i] = evaluate_criterion("gcv", m.name, n_points, n_areal, n_dofs,
                                               D_matrix, Psi, R1, R0, y, u, lambda_seq[i],
                                               fold_ids, n_folds, cfg);
            }

            Eigen::Index min_idx;
            scores.minCoeff(&min_idx);
            best_lambda = lambda_seq[min_idx];

            auto t_sel_end = std::chrono::high_resolution_clock::now();
            double elapsed_sel_s = std::chrono::duration<double>(t_sel_end - t_sel_start).count();

            write_csv(out_dir + "lambdas_seq_" + tag + ".csv", lambda_seq);
            write_csv(out_dir + "gcv_scores_" + gcv_seq_type + ".csv", scores);
            save_scalar(out_dir + "lambda_opt_" + tag + ".csv", best_lambda);
            save_scalar(out_dir + "time_selection_" + tag + ".csv", elapsed_sel_s);

            auto t_run_start = std::chrono::high_resolution_clock::now();
            auto res = dispatch_model(m.name, n_points, n_areal, n_dofs, D_matrix, Psi, R1, R0, y, u, best_lambda, cfg);
            auto t_run_end = std::chrono::high_resolution_clock::now();
            double elapsed_run_s = std::chrono::duration<double>(t_run_end - t_run_start).count();

            save_fit_results(out_dir, "gcv", res, is_single_variance);
            save_scalar(out_dir + "time_run_" + tag + ".csv", elapsed_run_s);
        }


        /**************************** CROSS VALIDATION *******************************/

        // Calcolo K-fold CV error per il confronto tra modelli
        if(compute_kfold_error){

            std::cout << "Running spatial K-fold CV for " << m.name << "..." << std::endl;
            vector_t cv_rmse_point(n_folds);
            vector_t cv_rmse_area(n_folds);
            vector_t cv_rmse_total(n_folds);

            for(unsigned int k = 1; k <= n_folds; ++k) {

                std::cout << "  Running fold number " << k << " of " << n_folds << " for " << m.name << "..." << std::endl;

                std::vector<unsigned int> idx_p_train, idx_p_test;
                for (unsigned int i = 0; i < n_points; ++i) {
                    if (fold_ids[i] == static_cast<int>(k)) idx_p_test.push_back(i);
                    else idx_p_train.push_back(i);
                }

                std::vector<unsigned int> idx_a_train, idx_a_test;
                for (unsigned int i = 0; i < n_areal; ++i) {
                    if (fold_ids[n_points + i] == static_cast<int>(k)) idx_a_test.push_back(i);
                    else idx_a_train.push_back(i);
                }

                unsigned int n_p_tr = idx_p_train.size();
                unsigned int n_a_tr = idx_a_train.size();
                unsigned int n_tr = n_p_tr + n_a_tr;

                vector_t y_tr(n_tr);
                std::vector<Eigen::Triplet<double>> triplets_D_tr;
                for (unsigned int r = 0; r < n_p_tr; ++r) {
                    y_tr[r] = y_point[idx_p_train[r]];
                    triplets_D_tr.emplace_back(r, r, 1.0);
                }
                for (unsigned int r = 0; r < n_a_tr; ++r) {
                    y_tr[n_p_tr + r] = y_areal[idx_a_train[r]];
                    triplets_D_tr.emplace_back(n_p_tr + r, n_p_tr + r, D_areal[idx_a_train[r]]);
                }
                sparse_matrix_t D_tr(n_tr, n_tr);
                D_tr.setFromTriplets(triplets_D_tr.begin(), triplets_D_tr.end());

                std::vector<int> row_map_p(n_points, -1);
                for (unsigned int r = 0; r < n_p_tr; ++r) row_map_p[idx_p_train[r]] = static_cast<int>(r);
                std::vector<int> row_map_a(n_areal, -1);
                for (unsigned int r = 0; r < n_a_tr; ++r) row_map_a[idx_a_train[r]] = static_cast<int>(n_p_tr + r);

                std::vector<Eigen::Triplet<double>> triplets_Psi_tr;
                for (int col = 0; col < Psi_point.outerSize(); ++col) {
                    for (sparse_matrix_t::InnerIterator it(Psi_point, col); it; ++it) {
                        int new_r = row_map_p[it.row()];
                        if (new_r >= 0) triplets_Psi_tr.emplace_back(new_r, it.col(), it.value());
                    }
                }
                for (int col = 0; col < Psi_areal.outerSize(); ++col) {
                    for (sparse_matrix_t::InnerIterator it(Psi_areal, col); it; ++it) {
                        int new_r = row_map_a[it.row()];
                        if (new_r >= 0) triplets_Psi_tr.emplace_back(new_r, it.col(), it.value());
                    }
                }
                sparse_matrix_t Psi_tr(n_tr, n_dofs);
                Psi_tr.setFromTriplets(triplets_Psi_tr.begin(), triplets_Psi_tr.end());

                vector_t f_fold;
                if (m.name == "hetero8") {
                    std::vector<double> dummy_seq;
                    auto res_fold = run_hetero8_srpde(n_p_tr, n_a_tr, n_dofs, D_tr, Psi_tr, R1, R0, y_tr, u,
                                                    lasco_bounds_default, lasco_step_default,
                                                    fine_half_window_xu, fine_step_default,
                                                    cfg, dummy_seq);
                    f_fold = res_fold.f;
                } else {
                    std::vector<int> dummy_folds(n_tr, 1);
                    vector_t fold_gcv_scores(lambda_seq.size());
                    for (size_t i = 0; i < lambda_seq.size(); ++i) {
                        fold_gcv_scores[i] = evaluate_criterion("gcv", m.name, n_p_tr, n_a_tr, n_dofs,
                                                                D_tr, Psi_tr, R1, R0, y_tr, u, lambda_seq[i],
                                                                dummy_folds, n_folds, cfg);
                    }
                    Eigen::Index min_fold_idx;
                    fold_gcv_scores.minCoeff(&min_fold_idx);
                    double lambda_fold_opt = lambda_seq[min_fold_idx];

                    auto res_fold = dispatch_model(m.name, n_p_tr, n_a_tr, n_dofs, D_tr, Psi_tr, R1, R0, y_tr, u, lambda_fold_opt, cfg);
                    f_fold = res_fold.f;
                }

                vector_t fn_pred_full = Psi * f_fold;

                double sse_p = 0.0;
                for (unsigned int idx : idx_p_test) {
                    double diff = y_point[idx] - fn_pred_full[idx];
                    sse_p += diff * diff;
                }

                double sse_a = 0.0;
                for (unsigned int idx : idx_a_test) {
                    double diff = y_areal[idx] - fn_pred_full[n_points + idx];
                    sse_a += diff * diff;
                }

                cv_rmse_point[k - 1] = idx_p_test.empty() ? 0.0 : std::sqrt(sse_p / idx_p_test.size());  // metto zero, tanto poi prendero' la media riscalando per il numero corretto
                cv_rmse_area[k - 1]  = idx_a_test.empty() ? 0.0 : std::sqrt(sse_a / idx_a_test.size());
                cv_rmse_total[k - 1] = (idx_p_test.empty() && idx_a_test.empty()) ? 0.0 :
                                    std::sqrt((sse_p + sse_a) / (idx_p_test.size() + idx_a_test.size()));
            }

            write_csv(out_dir + "cv_rmse_point.csv", cv_rmse_point);
            write_csv(out_dir + "cv_rmse_area.csv", cv_rmse_area);
            write_csv(out_dir + "cv_rmse_total.csv", cv_rmse_total);            


        }

    }

    return 0;
}
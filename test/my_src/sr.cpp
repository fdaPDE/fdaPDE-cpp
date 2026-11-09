// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.

#include <fdaPDE/fdapde.h>
#include <iostream>
#include <string>
#include <vector>
#include "my_fun.h"

#include <chrono>
#include <fstream>

using namespace fdapde;
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using sparse_matrix_t = Eigen::SparseMatrix<double>;

struct FitMethodSpec {
    std::string name;
    bool enabled;
    bool run_with_forced_lambda;
};

struct FitCriterionSpec {
    std::string name;
    std::string seq_type;
    bool enabled;
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

    // ==================== 3. SELEZIONE MODELLI DA FITTARE ======================

    const bool hetero_fit     = false;
    const bool hetero2_fit    = false;
    const bool hetero3_fit    = false;
    const bool hetero4_fit    = false;
    const bool hetero5_fit    = false;
    const bool hetero6_fit    = false;
    const bool hetero7_fit    = false;
    const bool naive_fit      = true;
    const bool only_point_fit = true;
    const bool only_area_fit  = true;

    // ==================== 4. OPZIONI LAMBDA FORZATO ============================

    const bool force_lambda = false;
    const double lambda_forced = 0.01;

    const bool force_hetero  = false;
    const bool force_hetero2 = false;
    const bool force_hetero3 = false;
    const bool force_hetero4 = false;
    const bool force_hetero5 = false;
    const bool force_hetero6 = false;
    const bool force_hetero7 = false;
    const bool force_naive   = false;

    // ==================== 5. SELEZIONE CRITERI PER IL FIT ======================

    const bool run_gcv   = true;  const std::string gcv_seq_type   = "fine";
    const bool run_cv    = false;  const std::string cv_seq_type    = "fine";
    const bool run_aicc  = false;  const std::string aicc_seq_type  = "fine";
    const bool run_gfaic = false;  const std::string gfaic_seq_type = "fine";
    const bool run_rcv   = false;  const std::string rcv_seq_type   = "fine";
    const bool run_fpe   = false;  const std::string fpe_seq_type   = "fine";

    // ==================== POPOLAMENTO AUTOMATICO TABELLE =======================

    std::vector<FitMethodSpec> methods = {
        {"hetero",     hetero_fit,     force_hetero},
        {"hetero2",    hetero2_fit,    force_hetero2},
        {"hetero3",    hetero3_fit,    force_hetero3},
        {"hetero4",    hetero4_fit,    force_hetero4},
        {"hetero5",    hetero5_fit,    force_hetero5},
        {"hetero6",    hetero6_fit,    force_hetero6},
        {"hetero7",    hetero7_fit,    force_hetero7},
        {"naive",      naive_fit,      force_naive},
        {"only_point", only_point_fit, false},
        {"only_area",  only_area_fit,  false}
    };

    std::vector<FitCriterionSpec> criteria = {
        {"gcv",   gcv_seq_type,   run_gcv},
        {"cv",    cv_seq_type,    run_cv},
        {"aicc",  aicc_seq_type,  run_aicc},
        {"gfaic", gfaic_seq_type, run_gfaic},
        {"rcv",   rcv_seq_type,   run_rcv},
        {"fpe",   fpe_seq_type,   run_fpe}
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
            std::cout << "Running " << m.name << " fit..." << std::endl;

            std::string out_dir = sim_dir + m.name + "/";
            bool is_single_variance = (m.name == "naive" || m.name == "only_point" || m.name == "only_area");

            for (const auto& crit : criteria) {

                if (!crit.enabled) continue;

                std::string lambda_file = out_dir + "lambda_opt_" + crit.name + "_" + crit.seq_type + ".csv";
                double lambda_opt = read_scalar(lambda_file);

                auto t_start = std::chrono::high_resolution_clock::now();

                auto res = dispatch_model(m.name, n_points, n_areal, n_dofs, D_matrix, Psi, R1, R0, y, u, lambda_opt, cfg);

                auto t_end = std::chrono::high_resolution_clock::now();
                elapsed_s = std::chrono::duration<double>(t_end - t_start).count(); 

                save_fit_results(out_dir, crit.name, res, is_single_variance);

                std::string tag = crit.name + "_" + crit.seq_type;
                save_scalar(out_dir + "time_run_" + tag + ".csv", elapsed_s);

            }

            if (force_lambda && m.run_with_forced_lambda) {
                std::string forced_dir = sim_dir + m.name + "_forced/";
                auto res_forced = dispatch_model(m.name, n_points, n_areal, n_dofs, D_matrix, Psi, R1, R0, y, u, lambda_forced, cfg);
                save_fit_results(forced_dir, "", res_forced, is_single_variance);
            }
        }
    }
    return 0;
}
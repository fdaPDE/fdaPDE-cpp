#ifndef MY_FUN_H
#define MY_FUN_H

#include <fdaPDE/fdapde.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace fdapde;
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using sparse_matrix_t = Eigen::SparseMatrix<double>;
using fdapde::GridSearch;

enum class EdfTarget { Global, Points, Areal };

/********************************************* HELPER FUNCTIONS *********************************************/

template <typename D_matrixType, typename W_type, typename InvAType>
double edf(const D_matrixType& D_matrix, const W_type& W_, const sparse_matrix_t& Psi, int n_dofs,
           const InvAType& invA_, EdfTarget target = EdfTarget::Global, int n_points = 0, int r = 100, int seed = 109) {

    int n_locs = Psi.rows();
    matrix_t Us = matrix_t::Zero(n_locs, r);
    int seed_ = (seed == random_seed) ? std::random_device()() : seed;
    std::mt19937 rng(seed_);
    rademacher_distribution rademacher;

    int start_idx = 0;
    int end_idx = n_locs;

    if (target == EdfTarget::Points) {
        end_idx = n_points;
    } else if (target == EdfTarget::Areal) {
        start_idx = n_points;
    }

    for (int i = start_idx; i < end_idx; ++i) {
        for (int j = 0; j < r; ++j) {
            Us(i, j) = rademacher(rng);
        }
    }

    matrix_t Ys = Us.transpose() * Psi;
    matrix_t Bs = matrix_t::Zero(2 * n_dofs, r);
    Bs.topRows(n_dofs) = -Psi.transpose() * D_matrix * W_ * Us;

    matrix_t x = invA_.solve(Bs);
    double trS = 0.0;
    for (int i = 0; i < r; ++i) {
        trS += Ys.row(i).dot(x.col(i).head(n_dofs));
    }
    return trS / r;
}

inline double ftPf(double lambda, const vector_t& g_, const sparse_matrix_t& R0_) {
    return lambda * g_.dot(R0_ * g_);
}

// Compute empirical variance for initialization 
inline double empirical_point_variance(const vector_t& values) {
    if (values.size() < 2) {
        throw std::invalid_argument(
            "At least two observations are required to compute a variance."
        );
    }

    const double mean = values.mean();
    const double ss = (values.array() - mean).square().sum();
    const double variance = ss / static_cast<double>(values.size() - 1);

    if (!std::isfinite(variance) || variance <= 0.0) {
        throw std::invalid_argument(
            "The empirical variance must be finite and strictly positive."
        );
    }

    return variance;
}
// Areal case. NOTE: since Var(Y_a) = sigma_a^2 / D_a, allora stimatore di sigma_a^2 = (1/(n-1)) * sum_i (Di_a*(Y_a_i - mean(Y_a))^2) 
inline double empirical_areal_variance(const vector_t& y_A, const sparse_matrix_t& D_matrix, unsigned int n_points) {
    const Eigen::Index n_A = y_A.size();

    if (n_A < 2) {
        throw std::invalid_argument(
            "At least two areal observations are required."
        );
    }

    if (D_matrix.rows() != D_matrix.cols() ||
        D_matrix.rows() < static_cast<Eigen::Index>(n_points) + n_A) {
        throw std::invalid_argument(
            "Incompatible dimensions for the areal variance initialization."
        );
    }

    double sum_weights = 0.0;
    double weighted_sum = 0.0;

    for (Eigen::Index i = 0; i < n_A; ++i) {
        const Eigen::Index idx =
            static_cast<Eigen::Index>(n_points) + i;

        const double d_i = D_matrix.coeff(idx, idx);

        if (!std::isfinite(d_i) || d_i <= 0.0) {
            throw std::invalid_argument(
                "Areal weights must be finite and strictly positive."
            );
        }

        sum_weights += d_i;
        weighted_sum += d_i * y_A(i);
    }

    const double weighted_mean = weighted_sum / sum_weights;

    double weighted_ss = 0.0;

    for (Eigen::Index i = 0; i < n_A; ++i) {
        const Eigen::Index idx =
            static_cast<Eigen::Index>(n_points) + i;

        const double d_i = D_matrix.coeff(idx, idx);
        const double residual = y_A(i) - weighted_mean;

        weighted_ss += d_i * residual * residual;
    }

    const double variance =
        weighted_ss / static_cast<double>(n_A - 1);

    if (!std::isfinite(variance) || variance <= 0.0) {
        throw std::invalid_argument(
            "The empirical areal variance scale must be finite and strictly positive."
        );
    }

    return variance;
}

inline sparse_matrix_t build_block_W(unsigned int n_points, unsigned int n_areal, double domain_area, double w_p, double w_a) {
    double scale = 1.0 / (static_cast<double>(n_points + n_areal) * domain_area); // normalizzazione per dimensione campionaria E area totale
    std::vector<Eigen::Triplet<double>> triplets_W;
    triplets_W.reserve(n_points + n_areal);
    for (unsigned int i = 0; i < n_points; ++i) {
        triplets_W.emplace_back(i, i, scale * w_p);
    }
    for (unsigned int i = 0; i < n_areal; ++i) {
        triplets_W.emplace_back(i + n_points, i + n_points, scale * w_a);
    }
    sparse_matrix_t W(n_points + n_areal, n_points + n_areal);
    W.setFromTriplets(triplets_W.begin(), triplets_W.end());
    return W;
}

inline std::vector<double> build_log10_seq(double exp_min, double exp_max, double step) {
    std::vector<double> seq;
    for (double xs = exp_min; xs <= exp_max + 1e-9; xs += step) {
        seq.push_back(std::pow(10.0, xs));
    }
    return seq;
}

inline void save_scalar(const std::string& filepath, double val) {
    std::ofstream f(filepath);
    if (f.is_open()) {
        f << std::setprecision(16) << val << "\n";
        f.close();
    }
}

inline double read_scalar(const std::string& filepath) {
    double val = 0.0;
    std::ifstream f(filepath);
    if (f.is_open()) {
        f >> val;
        f.close();
    }
    return val;
}

// Calcola la somma delle 4 tracce presenti in f_i(lambda) nella formula (39) di Xu:
// tr(N_lambda^{-1} N_i) + lambda * tr(C N_lambda^{-1} N_i)
// + lambda^2 * tr(C^2 N_lambda^{-1} N_i) + lambda^3 * tr(C^3 N_lambda^{-1} N_i)
// dove C = N_lambda^{-1} P, P = R1^T R0^{-1} R1 e N_i = Psi_i^T D_i W_i Psi_i.
template <typename D_matrixType, typename W_type, typename InvAType>
double edf_xu_block(const D_matrixType& D_matrix, const W_type& W, const sparse_matrix_t& Psi,
                    const sparse_matrix_t& R1, int n_dofs, const InvAType& invA,
                    EdfTarget target, int n_points, int r = 100, int seed = 109) {

    int n_locs = Psi.rows();
    matrix_t Us = matrix_t::Zero(n_locs, r);
    int seed_ = (seed == random_seed) ? std::random_device()() : seed;
    std::mt19937 rng(seed_);
    rademacher_distribution rademacher;

    int start_idx = (target == EdfTarget::Areal) ? n_points : 0;
    int end_idx   = (target == EdfTarget::Points) ? n_points : n_locs;

    for (int i = start_idx; i < end_idx; ++i) {
        for (int j = 0; j < r; ++j) {
            Us(i, j) = rademacher(rng);
        }
    }

    matrix_t Ys = Us.transpose() * Psi; // r x n_dofs

    // Termine di ordine 0: F0 = N_lambda^{-1} Psi_i^T D_i W_i Us_i
    matrix_t Bs = matrix_t::Zero(2 * n_dofs, r);
    Bs.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * Us;
    matrix_t F0 = invA.solve(Bs).topRows(n_dofs);

    // Per ogni matrice V, risolvere A_sp * [F_next; G_next] = [0; R1 * V]
    // restituisce esattamente F_next = (lambda * C) * V = lambda * N_lambda^{-1} R1^T R0^{-1} R1 * V.
    Bs.topRows(n_dofs).setZero();

    // Termine di ordine 1: F1 = (lambda * C) * F0
    Bs.bottomRows(n_dofs) = R1 * F0;
    matrix_t F1 = invA.solve(Bs).topRows(n_dofs);

    // Termine di ordine 2: F2 = (lambda * C)^2 * F0
    Bs.bottomRows(n_dofs) = R1 * F1;
    matrix_t F2 = invA.solve(Bs).topRows(n_dofs);

    // Termine di ordine 3: F3 = (lambda * C)^3 * F0
    Bs.bottomRows(n_dofs) = R1 * F2;
    matrix_t F3 = invA.solve(Bs).topRows(n_dofs);

    matrix_t F_tot = F0 + F1 + F2 + F3;

    double tr_sum = 0.0;
    for (int j = 0; j < r; ++j) {
        tr_sum += Ys.row(j).dot(F_tot.col(j));
    }
    return tr_sum / static_cast<double>(r);
}

// Valuta la funzione obiettivo del GCV pesato V(lambda, W_rho) secondo l'eq. (28a) di Xu
inline double compute_xu_weighted_gcv(unsigned int n_points, unsigned int n_areal, unsigned int n_dofs,
                                      double domain_area, const sparse_matrix_t& D_matrix,
                                      const sparse_matrix_t& Psi, const sparse_matrix_t& R1,
                                      const sparse_matrix_t& R0, const vector_t& y, const vector_t& u,
                                      double lambda, double rho, unsigned int random_seed) {

    const unsigned int n = n_points + n_areal;
    sparse_matrix_t W = build_block_W(n_points, n_areal, domain_area, 1.0, rho);

    SparseBlockMatrix<double, 2, 2> A(-Psi.transpose() * D_matrix * W * Psi, lambda * R1.transpose(), R1, R0);
    Eigen::SparseLU<sparse_matrix_t> invA(A);

    vector_t rhs = vector_t::Zero(2 * n_dofs);
    rhs.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * y;
    rhs.bottomRows(n_dofs) = u;

    vector_t f = invA.solve(rhs).topRows(n_dofs);

    // Residui regolarizzati r(lambda, W_rho) = (I - H(lambda, W_rho)) * y (eq. 28c)
    vector_t r_vec = y - Psi * f;

    // Numeratore pesato: (1/n) * r^T * W_rho * r, con W_rho = diag(I_np, rho * D_a)
    double r_W_r = 0.0;
    if (n_points > 0) {
        r_W_r += r_vec.head(n_points).squaredNorm();
    }
    if (n_areal > 0) {
        r_W_r += rho * (r_vec.tail(n_areal).array().square() * D_matrix.diagonal().tail(n_areal).array()).sum();
    }

    // Traccia di H(lambda, W_rho) (eq. 28b)
    double tr_H = edf(D_matrix, W, Psi, n_dofs, invA, EdfTarget::Global, n_points, 100, random_seed);
    double tr_I_minus_H = static_cast<double>(n) - tr_H;

    if (tr_I_minus_H <= 1e-12) {
        return std::numeric_limits<double>::max();
    }

    // V(lambda, W_rho) = [(1/n) * r^T W_rho r] / [(1/n) * tr(I - H)]^2 (eq. 28a)
    double num = r_W_r / static_cast<double>(n);
    double den = (tr_I_minus_H / static_cast<double>(n)) * (tr_I_minus_H / static_cast<double>(n));
    return num / den;
}



/************************************************** STRUCTS **************************************************/

struct HeteroSRPDEResults {
    vector_t f;
    vector_t g;
    vector_t fitted;
    double sigma_sq_p = 0.0;
    double sigma_sq_A = 0.0;
    double edfs = 0.0;
    unsigned int n_iter = 0;
    std::vector<double> obj_history = {};
    std::vector<double> loss_p_history = {};
    std::vector<double> loss_A_history = {};
    std::vector<double> loss_global_history = {};
    std::vector<double> sigma_sq_p_history = {};
    std::vector<double> sigma_sq_A_history = {};
    double lambda_opt = 0.0;
    vector_t gcv_scores = {};
    std::vector<double> lambda_history = {};
};

struct NaiveSRPDEResults {
    vector_t f;
    vector_t g;
    vector_t fitted;
    double sigma_sq = 0.0;
    double edfs = 0.0;
};

struct ModelConfig {
    double domain_area = 1.0;
    unsigned int random_seed_p = 17;
    unsigned int random_seed_A = 23;
    unsigned int random_seed = 37;
    double tol = 1e-6;
    unsigned int max_iter = 50;
    bool verbose = false;
    bool init_with_residuals = true;
};


// 2) Dichiarazioni anticipate:
inline HeteroSRPDEResults dispatch_model(const std::string& method,
                                         unsigned int n_points, unsigned int n_areal, unsigned int n_dofs,
                                         const sparse_matrix_t& D_matrix, const sparse_matrix_t& Psi,
                                         const sparse_matrix_t& R1, const sparse_matrix_t& R0,
                                         const vector_t& y, const vector_t& u, double lambda,
                                         const ModelConfig& cfg);

inline double ic_score_fun(const std::string& criterion, const std::string& method,
                           unsigned int n_points, unsigned int n_areal, unsigned int n_dofs,
                           const sparse_matrix_t& D_matrix, const sparse_matrix_t& Psi,
                           const sparse_matrix_t& R1, const sparse_matrix_t& R0,
                           const vector_t& y, const vector_t& u, double lambda,
                           const ModelConfig& cfg);


/****************************************** MODELS IMPLEMENTATION ******************************************/

// Hetero (v1): W = 1/(n * |Omega|) * diag(1, w_A), edf_p = tr(S_pp), edf_A = tr(S_aa)
HeteroSRPDEResults run_hetero_srpde(double init_wA, unsigned int n_points, unsigned int n_areal, unsigned int n_dofs,
                                    double domain_area, const sparse_matrix_t& D_matrix, const sparse_matrix_t& Psi,
                                    const sparse_matrix_t& R1, const sparse_matrix_t& R0,
                                    const vector_t& y, const vector_t& u, double lambda,
                                    unsigned int random_seed_p, unsigned int random_seed_A, unsigned int random_seed,
                                    const double tol_ = 1e-6, const unsigned int max_iter_ = 50, bool verbose = false) {

    double w_A = init_wA;
    double w_A_new = w_A;
    double sigma_sq_p = 0.0, sigma_sq_A = 0.0;
    vector_t f, g, fitted;
    double edf_p = 0.0, edf_A = 0.0, edfs = 0.0;

    double Jold = std::numeric_limits<double>::max(), Jnew = 0.0;
    unsigned int n_iter_ = 0;
    std::vector<double> obj_history, loss_p_history, loss_A_history, loss_global_history, sigma_sq_p_history, sigma_sq_A_history;

    while (n_iter_ < max_iter_ && std::abs(Jnew - Jold) > tol_) {
        w_A = w_A_new;
        sparse_matrix_t W = build_block_W(n_points, n_areal, domain_area, 1.0, w_A);

        SparseBlockMatrix<double, 2, 2> A(-Psi.transpose() * D_matrix * W * Psi, lambda * R1.transpose(), R1, R0);
        Eigen::SparseLU<sparse_matrix_t> invA(A);

        vector_t rhs = vector_t::Zero(2 * n_dofs, 1);
        rhs.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * y;
        rhs.bottomRows(n_dofs) = u;

        vector_t x = invA.solve(rhs);
        f = x.topRows(n_dofs);
        g = x.bottomRows(n_dofs);
        fitted = Psi * f;

        edf_p = edf(D_matrix, W, Psi, n_dofs, invA, EdfTarget::Points, n_points, 100, random_seed_p);
        edf_A = edf(D_matrix, W, Psi, n_dofs, invA, EdfTarget::Areal, n_points, 100, random_seed_A);

        vector_t res = y - fitted;
        double Q_p = res.head(n_points).squaredNorm();
        double Q_a = (res.tail(n_areal).array().square() * D_matrix.diagonal().tail(n_areal).array()).sum();

        sigma_sq_p = Q_p / (n_points - edf_p);
        sigma_sq_A = Q_a / (n_areal - edf_A);
        w_A_new = sigma_sq_p / sigma_sq_A;

        double data_loss_p = Q_p;
        double data_loss_A = w_A * Q_a;
        double data_loss = (data_loss_p + data_loss_A) / ((n_points + n_areal) * domain_area);

        Jold = Jnew;
        Jnew = data_loss + ftPf(lambda, g, R0);
        n_iter_++;

        if (verbose) {
            std::cout << "Iter " << n_iter_ << " | diff: " << std::abs(Jnew - Jold) << std::endl;
        }

        obj_history.push_back(Jnew);
        loss_p_history.push_back(data_loss_p / (n_points * domain_area));
        loss_A_history.push_back(data_loss_A / (n_areal * domain_area));
        loss_global_history.push_back(data_loss);
        sigma_sq_p_history.push_back(sigma_sq_p);
        sigma_sq_A_history.push_back(sigma_sq_A);
    }

    edfs = edf_p + edf_A;
    return {f, g, fitted, sigma_sq_p, sigma_sq_A, edfs, n_iter_, obj_history,
            loss_p_history, loss_A_history, loss_global_history, sigma_sq_p_history, sigma_sq_A_history};
}

// Hetero (v2): W = 1/(n * |Omega|) * diag(1, w_A), denominatore comune n - tr(S)
HeteroSRPDEResults run_hetero2_srpde(double init_wA, unsigned int n_points, unsigned int n_areal, unsigned int n_dofs,
                                     double domain_area, const sparse_matrix_t& D_matrix, const sparse_matrix_t& Psi,
                                     const sparse_matrix_t& R1, const sparse_matrix_t& R0,
                                     const vector_t& y, const vector_t& u, double lambda,
                                     unsigned int random_seed_p, unsigned int random_seed_A, unsigned int random_seed,
                                     const double tol_ = 1e-6, const unsigned int max_iter_ = 50, bool verbose = false) {

    double w_A = init_wA;
    sparse_matrix_t W = build_block_W(n_points, n_areal, domain_area, 1.0, w_A);
    double sigma_sq_p = 0.0, sigma_sq_A = 0.0, edf_glob = 0.0;
    vector_t f, g, fitted;

    double Jold = std::numeric_limits<double>::max(), Jnew = 0.0;
    unsigned int n_iter_ = 0;
    std::vector<double> obj_history, loss_p_history, loss_A_history, loss_global_history, sigma_sq_p_history, sigma_sq_A_history;

    while (n_iter_ < max_iter_ && std::abs(Jnew - Jold) > tol_) {
        SparseBlockMatrix<double, 2, 2> A(-Psi.transpose() * D_matrix * W * Psi, lambda * R1.transpose(), R1, R0);
        Eigen::SparseLU<sparse_matrix_t> invA(A);

        vector_t rhs = vector_t::Zero(2 * n_dofs, 1);
        rhs.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * y;
        rhs.bottomRows(n_dofs) = u;

        vector_t x = invA.solve(rhs);
        f = x.topRows(n_dofs);
        g = x.bottomRows(n_dofs);
        fitted = Psi * f;

        edf_glob = edf(D_matrix, W, Psi, n_dofs, invA, EdfTarget::Global, n_points + n_areal, 100, random_seed);

        vector_t res = y - fitted;
        double Q_p = res.head(n_points).squaredNorm();
        double Q_a = (res.tail(n_areal).array().square() * D_matrix.diagonal().tail(n_areal).array()).sum();

        double data_loss_p = Q_p;
        double data_loss_A = w_A * Q_a;
        double data_loss = (data_loss_p + data_loss_A) / ((n_points + n_areal) * domain_area);

        sigma_sq_p = Q_p / ((n_points + n_areal) - edf_glob);
        sigma_sq_A = Q_a / ((n_points + n_areal) - edf_glob);
        w_A = sigma_sq_p / sigma_sq_A;
        W = build_block_W(n_points, n_areal, domain_area, 1.0, w_A);

        Jold = Jnew;
        Jnew = data_loss + ftPf(lambda, g, R0);
        n_iter_++;

        obj_history.push_back(Jnew);
        loss_p_history.push_back(data_loss_p / (n_points * domain_area));
        loss_A_history.push_back(data_loss_A / (n_areal * domain_area));
        loss_global_history.push_back(data_loss);
        sigma_sq_p_history.push_back(sigma_sq_p);
        sigma_sq_A_history.push_back(sigma_sq_A);
    }

    return {f, g, fitted, sigma_sq_p, sigma_sq_A, edf_glob, n_iter_, obj_history,
            loss_p_history, loss_A_history, loss_global_history, sigma_sq_p_history, sigma_sq_A_history};
}

// Hetero (v3): W = 1/(n * |Omega|) * diag(1, w_A), sistema 2x2 su expected RSS
HeteroSRPDEResults run_hetero3_srpde(double init_wA, unsigned int n_points, unsigned int n_areal, unsigned int n_dofs,
                                     double domain_area, const sparse_matrix_t& D_matrix, const sparse_matrix_t& Psi,
                                     const sparse_matrix_t& R1, const sparse_matrix_t& R0,
                                     const vector_t& y, const vector_t& u, double lambda,
                                     unsigned int random_seed_p, unsigned int random_seed_A, unsigned int random_seed,
                                     const double tol_ = 1e-6, const unsigned int max_iter_ = 50, bool verbose = false) {

    double w_A = init_wA;
    sparse_matrix_t W = build_block_W(n_points, n_areal, domain_area, 1.0, w_A);
    double sigma_sq_p = 0.0, sigma_sq_A = 0.0, edfs_total = 0.0;
    vector_t f, g, fitted;

    double Jold = std::numeric_limits<double>::max(), Jnew = 0.0;
    unsigned int n_iter_ = 0;
    std::vector<double> obj_history, loss_p_history, loss_A_history, loss_global_history, sigma_sq_p_history, sigma_sq_A_history;

    while (n_iter_ < max_iter_ && std::abs(Jnew - Jold) > tol_) {
        SparseBlockMatrix<double, 2, 2> A(-Psi.transpose() * D_matrix * W * Psi, lambda * R1.transpose(), R1, R0);
        Eigen::SparseLU<sparse_matrix_t> invA(A);

        vector_t rhs = vector_t::Zero(2 * n_dofs, 1);
        rhs.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * y;
        rhs.bottomRows(n_dofs) = u;

        vector_t x = invA.solve(rhs);
        f = x.topRows(n_dofs);
        g = x.bottomRows(n_dofs);
        fitted = Psi * f;

        int r = 100;
        matrix_t X_p = matrix_t::Zero(n_points + n_areal, r);
        matrix_t X_a = matrix_t::Zero(n_points + n_areal, r);

        std::mt19937 rng_p(random_seed_p);
        std::mt19937 rng_A(random_seed_A);
        rademacher_distribution rademacher;

        for (unsigned int i = 0; i < n_points; ++i) {
            for (int j = 0; j < r; ++j) X_p(i, j) = rademacher(rng_p);
        }
        for (unsigned int i = 0; i < n_areal; ++i) {
            double d_a = D_matrix.coeff(n_points + i, n_points + i);
            double inv_sqrt_d = 1.0 / std::sqrt(d_a);
            for (int j = 0; j < r; ++j) X_a(n_points + i, j) = rademacher(rng_A) * inv_sqrt_d;
        }

        matrix_t Bs_p = matrix_t::Zero(2 * n_dofs, r);
        Bs_p.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * X_p;
        matrix_t V_p = Psi * invA.solve(Bs_p).topRows(n_dofs);

        matrix_t Bs_a = matrix_t::Zero(2 * n_dofs, r);
        Bs_a.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * X_a;
        matrix_t V_a = Psi * invA.solve(Bs_a).topRows(n_dofs);

        double tr_S_pp = 0.0, tr_S_pp_S_pp_T = 0.0, tr_S_ap_T_Da_S_ap = 0.0;
        double tr_S_aa = 0.0, tr_S_pa_Da_inv_S_pa_T = 0.0, tr_S_aa_T_Da_S_aa_Da_inv = 0.0;

        for (int j = 0; j < r; ++j) {
            vector_t S_pp_Up = V_p.col(j).head(n_points);
            vector_t S_ap_Up = V_p.col(j).tail(n_areal);
            vector_t Up = X_p.col(j).head(n_points);

            tr_S_pp += Up.dot(S_pp_Up);
            tr_S_pp_S_pp_T += S_pp_Up.squaredNorm();
            for (unsigned int i = 0; i < n_areal; ++i) {
                tr_S_ap_T_Da_S_ap += S_ap_Up(i) * S_ap_Up(i) * D_matrix.coeff(n_points + i, n_points + i);
            }

            vector_t S_pa_Za = V_a.col(j).head(n_points);
            vector_t S_aa_Za = V_a.col(j).tail(n_areal);
            vector_t Za = X_a.col(j).tail(n_areal);

            tr_S_pa_Da_inv_S_pa_T += S_pa_Za.squaredNorm();
            for (unsigned int i = 0; i < n_areal; ++i) {
                double d_a = D_matrix.coeff(n_points + i, n_points + i);
                tr_S_aa_T_Da_S_aa_Da_inv += S_aa_Za(i) * S_aa_Za(i) * d_a;
                tr_S_aa += Za(i) * S_aa_Za(i) * d_a;
            }
        }

        tr_S_pp /= r; tr_S_pp_S_pp_T /= r; tr_S_ap_T_Da_S_ap /= r;
        tr_S_aa /= r; tr_S_pa_Da_inv_S_pa_T /= r; tr_S_aa_T_Da_S_aa_Da_inv /= r;
        edfs_total = tr_S_pp + tr_S_aa;

        double c_pp = n_points - 2.0 * tr_S_pp + tr_S_pp_S_pp_T;
        double c_pa = tr_S_pa_Da_inv_S_pa_T;
        double c_ap = tr_S_ap_T_Da_S_ap;
        double c_aa = n_areal - 2.0 * tr_S_aa + tr_S_aa_T_Da_S_aa_Da_inv;

        vector_t res = y - fitted;
        double Q_p = res.head(n_points).squaredNorm();
        double Q_a = (res.tail(n_areal).array().square() * D_matrix.diagonal().tail(n_areal).array()).sum();

        double data_loss_p = Q_p;
        double data_loss_A = w_A * Q_a;
        double data_loss = (data_loss_p + data_loss_A) / ((n_points + n_areal) * domain_area);

        double den = (c_pp * c_aa - c_pa * c_ap);
        if (std::abs(den) < 1e-12) den = 1e-8;

        sigma_sq_p = (c_aa * Q_p - c_pa * Q_a) / den;
        sigma_sq_A = (c_pp * Q_a - c_ap * Q_p) / den;
        w_A = sigma_sq_p / sigma_sq_A;
        W = build_block_W(n_points, n_areal, domain_area, 1.0, w_A);

        Jold = Jnew;
        Jnew = data_loss + ftPf(lambda, g, R0);
        n_iter_++;

        obj_history.push_back(Jnew);
        loss_p_history.push_back(data_loss_p / (n_points * domain_area));
        loss_A_history.push_back(data_loss_A / (n_areal * domain_area));
        loss_global_history.push_back(data_loss);
        sigma_sq_p_history.push_back(sigma_sq_p);
        sigma_sq_A_history.push_back(sigma_sq_A);
    }

    return {f, g, fitted, sigma_sq_p, sigma_sq_A, edfs_total, n_iter_, obj_history,
            loss_p_history, loss_A_history, loss_global_history, sigma_sq_p_history, sigma_sq_A_history};
}

// Hetero (v4): W = 1/(n * |Omega|) * diag(1/sigma2_p, 1/sigma2_A), edf_p = tr(S_pp), edf_A = tr(S_aa)
HeteroSRPDEResults run_hetero4_srpde(double init_sigma_sq_p, double init_sigma_sq_A,
                                     unsigned int n_points, unsigned int n_areal, unsigned int n_dofs,
                                     double domain_area, const sparse_matrix_t& D_matrix, const sparse_matrix_t& Psi,
                                     const sparse_matrix_t& R1, const sparse_matrix_t& R0,
                                     const vector_t& y, const vector_t& u, double lambda,
                                     unsigned int random_seed_p, unsigned int random_seed_A, unsigned int random_seed,
                                     const double tol_ = 1e-6, const unsigned int max_iter_ = 50, bool verbose = false) {

    double sigma_sq_p = init_sigma_sq_p;
    double sigma_sq_A = init_sigma_sq_A;
    double sigma_sq_p_new = sigma_sq_p, sigma_sq_A_new = sigma_sq_A;

    vector_t f, g, fitted;
    double edf_p = 0.0, edf_A = 0.0, edfs = 0.0;
    double Jold = std::numeric_limits<double>::max(), Jnew = 0.0;
    unsigned int n_iter_ = 0;

    std::vector<double> obj_history, loss_p_history, loss_A_history, loss_global_history, sigma_sq_p_history, sigma_sq_A_history;

    while (n_iter_ < max_iter_ && std::abs(Jnew - Jold) > tol_) {
        sparse_matrix_t W = build_block_W(n_points, n_areal, domain_area, 1.0 / sigma_sq_p, 1.0 / sigma_sq_A);

        SparseBlockMatrix<double, 2, 2> A(-Psi.transpose() * D_matrix * W * Psi, lambda * R1.transpose(), R1, R0);
        Eigen::SparseLU<sparse_matrix_t> invA(A);

        vector_t rhs = vector_t::Zero(2 * n_dofs, 1);
        rhs.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * y;
        rhs.bottomRows(n_dofs) = u;

        vector_t x = invA.solve(rhs);
        f = x.topRows(n_dofs);
        g = x.bottomRows(n_dofs);
        fitted = Psi * f;

        edf_p = edf(D_matrix, W, Psi, n_dofs, invA, EdfTarget::Points, n_points, 100, random_seed_p);
        edf_A = edf(D_matrix, W, Psi, n_dofs, invA, EdfTarget::Areal, n_points, 100, random_seed_A);

        vector_t res = y - fitted;
        double Q_p = res.head(n_points).squaredNorm();
        double Q_a = (res.tail(n_areal).array().square() * D_matrix.diagonal().tail(n_areal).array()).sum();

        sigma_sq_p_new = Q_p / (n_points - edf_p);
        sigma_sq_A_new = Q_a / (n_areal - edf_A);

        double data_loss_p = Q_p / sigma_sq_p;
        double data_loss_A = Q_a / sigma_sq_A;
        double data_loss = (data_loss_p + data_loss_A) / ((n_points + n_areal) * domain_area);

        Jold = Jnew;
        Jnew = data_loss + ftPf(lambda, g, R0);
        ++n_iter_;

        obj_history.push_back(Jnew);
        loss_p_history.push_back(data_loss_p / (n_points * domain_area));
        loss_A_history.push_back(data_loss_A / (n_areal * domain_area));
        loss_global_history.push_back(data_loss);
        sigma_sq_p_history.push_back(sigma_sq_p_new);
        sigma_sq_A_history.push_back(sigma_sq_A_new);

        sigma_sq_p = sigma_sq_p_new;
        sigma_sq_A = sigma_sq_A_new;
    }

    edfs = edf_p + edf_A;
    return {f, g, fitted, sigma_sq_p, sigma_sq_A, edfs, n_iter_, obj_history,
            loss_p_history, loss_A_history, loss_global_history, sigma_sq_p_history, sigma_sq_A_history};
}

// Hetero (v5): W = 1/(n * |Omega|) * diag(1/sigma2_p, 1/sigma2_A), sistema 2x2 su expected RSS
HeteroSRPDEResults run_hetero5_srpde(double init_sigma_sq_p, double init_sigma_sq_A,
                                     unsigned int n_points, unsigned int n_areal, unsigned int n_dofs,
                                     double domain_area, const sparse_matrix_t& D_matrix, const sparse_matrix_t& Psi,
                                     const sparse_matrix_t& R1, const sparse_matrix_t& R0,
                                     const vector_t& y, const vector_t& u, double lambda,
                                     unsigned int random_seed_p, unsigned int random_seed_A, unsigned int random_seed,
                                     const double tol_ = 1e-6, const unsigned int max_iter_ = 50, bool verbose = false) {

    double sigma_sq_p = init_sigma_sq_p;
    double sigma_sq_A = init_sigma_sq_A;
    double sigma_sq_p_new = sigma_sq_p, sigma_sq_A_new = sigma_sq_A;

    vector_t f, g, fitted;
    double edfs = 0.0;
    double Jold = std::numeric_limits<double>::max(), Jnew = 0.0;
    unsigned int n_iter_ = 0;

    std::vector<double> obj_history, loss_p_history, loss_A_history, loss_global_history, sigma_sq_p_history, sigma_sq_A_history;

    while (n_iter_ < max_iter_ && std::abs(Jnew - Jold) > tol_) {
        sparse_matrix_t W = build_block_W(n_points, n_areal, domain_area, 1.0 / sigma_sq_p, 1.0 / sigma_sq_A);

        SparseBlockMatrix<double, 2, 2> A(-Psi.transpose() * D_matrix * W * Psi, lambda * R1.transpose(), R1, R0);
        Eigen::SparseLU<sparse_matrix_t> invA(A);

        vector_t rhs = vector_t::Zero(2 * n_dofs, 1);
        rhs.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * y;
        rhs.bottomRows(n_dofs) = u;

        vector_t x = invA.solve(rhs);
        f = x.topRows(n_dofs);
        g = x.bottomRows(n_dofs);
        fitted = Psi * f;

        vector_t res = y - fitted;
        double Q_p = res.head(n_points).squaredNorm();
        double Q_a = (res.tail(n_areal).array().square() * D_matrix.diagonal().tail(n_areal).array()).sum();

        int r = 100;
        matrix_t X_p = matrix_t::Zero(n_points + n_areal, r);
        matrix_t X_a = matrix_t::Zero(n_points + n_areal, r);

        std::mt19937 rng_p(random_seed_p);
        std::mt19937 rng_A(random_seed_A);
        rademacher_distribution rademacher;

        for (unsigned int i = 0; i < n_points; ++i) {
            for (int j = 0; j < r; ++j) X_p(i, j) = rademacher(rng_p);
        }
        for (unsigned int i = 0; i < n_areal; ++i) {
            double d_a = D_matrix.coeff(n_points + i, n_points + i);
            double inv_sqrt_d = 1.0 / std::sqrt(d_a);
            for (int j = 0; j < r; ++j) X_a(n_points + i, j) = rademacher(rng_A) * inv_sqrt_d;
        }

        matrix_t Bs_p = matrix_t::Zero(2 * n_dofs, r);
        Bs_p.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * X_p;
        matrix_t V_p = Psi * invA.solve(Bs_p).topRows(n_dofs);

        matrix_t Bs_a = matrix_t::Zero(2 * n_dofs, r);
        Bs_a.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * X_a;
        matrix_t V_a = Psi * invA.solve(Bs_a).topRows(n_dofs);

        double tr_S_pp = 0.0, tr_S_pp_S_pp_T = 0.0, tr_S_ap_T_Da_S_ap = 0.0;
        double tr_S_aa = 0.0, tr_S_pa_Da_inv_S_pa_T = 0.0, tr_S_aa_T_Da_S_aa_Da_inv = 0.0;

        for (int j = 0; j < r; ++j) {
            vector_t S_pp_Up = V_p.col(j).head(n_points);
            vector_t S_ap_Up = V_p.col(j).tail(n_areal);
            vector_t Up = X_p.col(j).head(n_points);

            tr_S_pp += Up.dot(S_pp_Up);
            tr_S_pp_S_pp_T += S_pp_Up.squaredNorm();
            for (unsigned int i = 0; i < n_areal; ++i) {
                tr_S_ap_T_Da_S_ap += S_ap_Up(i) * S_ap_Up(i) * D_matrix.coeff(n_points + i, n_points + i);
            }

            vector_t S_pa_Za = V_a.col(j).head(n_points);
            vector_t S_aa_Za = V_a.col(j).tail(n_areal);
            vector_t Za = X_a.col(j).tail(n_areal);

            tr_S_pa_Da_inv_S_pa_T += S_pa_Za.squaredNorm();
            for (unsigned int i = 0; i < n_areal; ++i) {
                double d_a = D_matrix.coeff(n_points + i, n_points + i);
                tr_S_aa_T_Da_S_aa_Da_inv += S_aa_Za(i) * S_aa_Za(i) * d_a;
                tr_S_aa += Za(i) * S_aa_Za(i) * d_a;
            }
        }

        tr_S_pp /= r; tr_S_pp_S_pp_T /= r; tr_S_ap_T_Da_S_ap /= r;
        tr_S_aa /= r; tr_S_pa_Da_inv_S_pa_T /= r; tr_S_aa_T_Da_S_aa_Da_inv /= r;
        edfs = tr_S_pp + tr_S_aa;

        double c_pp = n_points - 2.0 * tr_S_pp + tr_S_pp_S_pp_T;
        double c_pa = tr_S_pa_Da_inv_S_pa_T;
        double c_ap = tr_S_ap_T_Da_S_ap;
        double c_aa = n_areal - 2.0 * tr_S_aa + tr_S_aa_T_Da_S_aa_Da_inv;

        double den = c_pp * c_aa - c_pa * c_ap;
        if (std::abs(den) < 1e-12) den = 1e-8;

        sigma_sq_p_new = (c_aa * Q_p - c_pa * Q_a) / den;
        sigma_sq_A_new = (c_pp * Q_a - c_ap * Q_p) / den;

        double data_loss_p = Q_p / sigma_sq_p;
        double data_loss_A = Q_a / sigma_sq_A;
        double data_loss = (data_loss_p + data_loss_A) / ((n_points + n_areal) * domain_area);

        Jold = Jnew;
        Jnew = data_loss + ftPf(lambda, g, R0);
        ++n_iter_;

        obj_history.push_back(Jnew);
        loss_p_history.push_back(data_loss_p / (n_points * domain_area));
        loss_A_history.push_back(data_loss_A / (n_areal * domain_area));
        loss_global_history.push_back(data_loss);
        sigma_sq_p_history.push_back(sigma_sq_p_new);
        sigma_sq_A_history.push_back(sigma_sq_A_new);

        sigma_sq_p = sigma_sq_p_new;
        sigma_sq_A = sigma_sq_A_new;
    }

    return {f, g, fitted, sigma_sq_p, sigma_sq_A, edfs, n_iter_, obj_history,
            loss_p_history, loss_A_history, loss_global_history, sigma_sq_p_history, sigma_sq_A_history};
}


// Hetero (v6): come hetero3, ma edfs = tr(2 S_tilde - S_tilde^T S_tilde)
HeteroSRPDEResults run_hetero6_srpde(double init_wA, unsigned int n_points, unsigned int n_areal, unsigned int n_dofs,
                                     double domain_area, const sparse_matrix_t& D_matrix, const sparse_matrix_t& Psi,
                                     const sparse_matrix_t& R1, const sparse_matrix_t& R0,
                                     const vector_t& y, const vector_t& u, double lambda,
                                     unsigned int random_seed_p, unsigned int random_seed_A, unsigned int random_seed,
                                     const double tol_ = 1e-6, const unsigned int max_iter_ = 50, bool verbose = false) {

    double w_A = init_wA;
    sparse_matrix_t W = build_block_W(n_points, n_areal, domain_area, 1.0, w_A);
    double sigma_sq_p = 0.0, sigma_sq_A = 0.0, edfs_total = 0.0;
    vector_t f, g, fitted;

    double Jold = std::numeric_limits<double>::max(), Jnew = 0.0;
    unsigned int n_iter_ = 0;
    std::vector<double> obj_history, loss_p_history, loss_A_history, loss_global_history, sigma_sq_p_history, sigma_sq_A_history;

    while (n_iter_ < max_iter_ && std::abs(Jnew - Jold) > tol_) {
        SparseBlockMatrix<double, 2, 2> A(-Psi.transpose() * D_matrix * W * Psi, lambda * R1.transpose(), R1, R0);
        Eigen::SparseLU<sparse_matrix_t> invA(A);

        vector_t rhs = vector_t::Zero(2 * n_dofs, 1);
        rhs.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * y;
        rhs.bottomRows(n_dofs) = u;

        vector_t x = invA.solve(rhs);
        f = x.topRows(n_dofs);
        g = x.bottomRows(n_dofs);
        fitted = Psi * f;

        int r = 100;
        matrix_t X_p = matrix_t::Zero(n_points + n_areal, r);
        matrix_t X_a = matrix_t::Zero(n_points + n_areal, r);

        std::mt19937 rng_p(random_seed_p);
        std::mt19937 rng_A(random_seed_A);
        rademacher_distribution rademacher;

        for (unsigned int i = 0; i < n_points; ++i) {
            for (int j = 0; j < r; ++j) X_p(i, j) = rademacher(rng_p);
        }
        for (unsigned int i = 0; i < n_areal; ++i) {
            double d_a = D_matrix.coeff(n_points + i, n_points + i);
            double inv_sqrt_d = 1.0 / std::sqrt(d_a);
            for (int j = 0; j < r; ++j) X_a(n_points + i, j) = rademacher(rng_A) * inv_sqrt_d;
        }

        matrix_t Bs_p = matrix_t::Zero(2 * n_dofs, r);
        Bs_p.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * X_p;
        matrix_t V_p = Psi * invA.solve(Bs_p).topRows(n_dofs);

        matrix_t Bs_a = matrix_t::Zero(2 * n_dofs, r);
        Bs_a.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * X_a;
        matrix_t V_a = Psi * invA.solve(Bs_a).topRows(n_dofs);

        double tr_S_pp = 0.0, tr_S_pp_S_pp_T = 0.0, tr_S_ap_T_Da_S_ap = 0.0;
        double tr_S_aa = 0.0, tr_S_pa_Da_inv_S_pa_T = 0.0, tr_S_aa_T_Da_S_aa_Da_inv = 0.0;

        for (int j = 0; j < r; ++j) {
            vector_t S_pp_Up = V_p.col(j).head(n_points);
            vector_t S_ap_Up = V_p.col(j).tail(n_areal);
            vector_t Up = X_p.col(j).head(n_points);

            tr_S_pp += Up.dot(S_pp_Up);
            tr_S_pp_S_pp_T += S_pp_Up.squaredNorm();
            for (unsigned int i = 0; i < n_areal; ++i) {
                tr_S_ap_T_Da_S_ap += S_ap_Up(i) * S_ap_Up(i) * D_matrix.coeff(n_points + i, n_points + i);
            }

            vector_t S_pa_Za = V_a.col(j).head(n_points);
            vector_t S_aa_Za = V_a.col(j).tail(n_areal);
            vector_t Za = X_a.col(j).tail(n_areal);

            tr_S_pa_Da_inv_S_pa_T += S_pa_Za.squaredNorm();
            for (unsigned int i = 0; i < n_areal; ++i) {
                double d_a = D_matrix.coeff(n_points + i, n_points + i);
                tr_S_aa_T_Da_S_aa_Da_inv += S_aa_Za(i) * S_aa_Za(i) * d_a;
                tr_S_aa += Za(i) * S_aa_Za(i) * d_a;
            }
        }

        tr_S_pp /= r; tr_S_pp_S_pp_T /= r; tr_S_ap_T_Da_S_ap /= r;
        tr_S_aa /= r; tr_S_pa_Da_inv_S_pa_T /= r; tr_S_aa_T_Da_S_aa_Da_inv /= r;

        double c_pp = n_points - 2.0 * tr_S_pp + tr_S_pp_S_pp_T;
        double c_pa = tr_S_pa_Da_inv_S_pa_T;
        double c_ap = tr_S_ap_T_Da_S_ap;
        double c_aa = n_areal - 2.0 * tr_S_aa + tr_S_aa_T_Da_S_aa_Da_inv;

        edfs_total = static_cast<double>(n_points + n_areal) - (c_pp + c_ap + c_pa + c_aa);

        vector_t res = y - fitted;
        double Q_p = res.head(n_points).squaredNorm();
        double Q_a = (res.tail(n_areal).array().square() * D_matrix.diagonal().tail(n_areal).array()).sum();

        double data_loss_p = Q_p;
        double data_loss_A = w_A * Q_a;
        double data_loss = (data_loss_p + data_loss_A) / ((n_points + n_areal) * domain_area);

        double den = (c_pp * c_aa - c_pa * c_ap);
        if (std::abs(den) < 1e-12) den = 1e-8;

        sigma_sq_p = (c_aa * Q_p - c_pa * Q_a) / den;
        sigma_sq_A = (c_pp * Q_a - c_ap * Q_p) / den;
        w_A = sigma_sq_p / sigma_sq_A;
        W = build_block_W(n_points, n_areal, domain_area, 1.0, w_A);

        Jold = Jnew;
        Jnew = data_loss + ftPf(lambda, g, R0);
        n_iter_++;

        obj_history.push_back(Jnew);
        loss_p_history.push_back(data_loss_p / (n_points * domain_area));
        loss_A_history.push_back(data_loss_A / (n_areal * domain_area));
        loss_global_history.push_back(data_loss);
        sigma_sq_p_history.push_back(sigma_sq_p);
        sigma_sq_A_history.push_back(sigma_sq_A);
    }

    return {f, g, fitted, sigma_sq_p, sigma_sq_A, edfs_total, n_iter_, obj_history,
            loss_p_history, loss_A_history, loss_global_history, sigma_sq_p_history, sigma_sq_A_history};
}

// Hetero (v7): come hetero5, ma edfs = tr(2 S_tilde - S_tilde^T S_tilde)
HeteroSRPDEResults run_hetero7_srpde(double init_sigma_sq_p, double init_sigma_sq_A,
                                     unsigned int n_points, unsigned int n_areal, unsigned int n_dofs,
                                     double domain_area, const sparse_matrix_t& D_matrix, const sparse_matrix_t& Psi,
                                     const sparse_matrix_t& R1, const sparse_matrix_t& R0,
                                     const vector_t& y, const vector_t& u, double lambda,
                                     unsigned int random_seed_p, unsigned int random_seed_A, unsigned int random_seed,
                                     const double tol_ = 1e-6, const unsigned int max_iter_ = 50, bool verbose = false) {

    double sigma_sq_p = init_sigma_sq_p;
    double sigma_sq_A = init_sigma_sq_A;
    double sigma_sq_p_new = sigma_sq_p, sigma_sq_A_new = sigma_sq_A;

    vector_t f, g, fitted;
    double edfs = 0.0;
    double Jold = std::numeric_limits<double>::max(), Jnew = 0.0;
    unsigned int n_iter_ = 0;

    std::vector<double> obj_history, loss_p_history, loss_A_history, loss_global_history, sigma_sq_p_history, sigma_sq_A_history;

    while (n_iter_ < max_iter_ && std::abs(Jnew - Jold) > tol_) {
        sparse_matrix_t W = build_block_W(n_points, n_areal, domain_area, 1.0 / sigma_sq_p, 1.0 / sigma_sq_A);

        SparseBlockMatrix<double, 2, 2> A(-Psi.transpose() * D_matrix * W * Psi, lambda * R1.transpose(), R1, R0);
        Eigen::SparseLU<sparse_matrix_t> invA(A);

        vector_t rhs = vector_t::Zero(2 * n_dofs, 1);
        rhs.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * y;
        rhs.bottomRows(n_dofs) = u;

        vector_t x = invA.solve(rhs);
        f = x.topRows(n_dofs);
        g = x.bottomRows(n_dofs);
        fitted = Psi * f;

        vector_t res = y - fitted;
        double Q_p = res.head(n_points).squaredNorm();
        double Q_a = (res.tail(n_areal).array().square() * D_matrix.diagonal().tail(n_areal).array()).sum();

        int r = 100;
        matrix_t X_p = matrix_t::Zero(n_points + n_areal, r);
        matrix_t X_a = matrix_t::Zero(n_points + n_areal, r);

        std::mt19937 rng_p(random_seed_p);
        std::mt19937 rng_A(random_seed_A);
        rademacher_distribution rademacher;

        for (unsigned int i = 0; i < n_points; ++i) {
            for (int j = 0; j < r; ++j) X_p(i, j) = rademacher(rng_p);
        }
        for (unsigned int i = 0; i < n_areal; ++i) {
            double d_a = D_matrix.coeff(n_points + i, n_points + i);
            double inv_sqrt_d = 1.0 / std::sqrt(d_a);
            for (int j = 0; j < r; ++j) X_a(n_points + i, j) = rademacher(rng_A) * inv_sqrt_d;
        }

        matrix_t Bs_p = matrix_t::Zero(2 * n_dofs, r);
        Bs_p.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * X_p;
        matrix_t V_p = Psi * invA.solve(Bs_p).topRows(n_dofs);

        matrix_t Bs_a = matrix_t::Zero(2 * n_dofs, r);
        Bs_a.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * X_a;
        matrix_t V_a = Psi * invA.solve(Bs_a).topRows(n_dofs);

        double tr_S_pp = 0.0, tr_S_pp_S_pp_T = 0.0, tr_S_ap_T_Da_S_ap = 0.0;
        double tr_S_aa = 0.0, tr_S_pa_Da_inv_S_pa_T = 0.0, tr_S_aa_T_Da_S_aa_Da_inv = 0.0;

        for (int j = 0; j < r; ++j) {
            vector_t S_pp_Up = V_p.col(j).head(n_points);
            vector_t S_ap_Up = V_p.col(j).tail(n_areal);
            vector_t Up = X_p.col(j).head(n_points);

            tr_S_pp += Up.dot(S_pp_Up);
            tr_S_pp_S_pp_T += S_pp_Up.squaredNorm();
            for (unsigned int i = 0; i < n_areal; ++i) {
                tr_S_ap_T_Da_S_ap += S_ap_Up(i) * S_ap_Up(i) * D_matrix.coeff(n_points + i, n_points + i);
            }

            vector_t S_pa_Za = V_a.col(j).head(n_points);
            vector_t S_aa_Za = V_a.col(j).tail(n_areal);
            vector_t Za = X_a.col(j).tail(n_areal);

            tr_S_pa_Da_inv_S_pa_T += S_pa_Za.squaredNorm();
            for (unsigned int i = 0; i < n_areal; ++i) {
                double d_a = D_matrix.coeff(n_points + i, n_points + i);
                tr_S_aa_T_Da_S_aa_Da_inv += S_aa_Za(i) * S_aa_Za(i) * d_a;
                tr_S_aa += Za(i) * S_aa_Za(i) * d_a;
            }
        }

        tr_S_pp /= r; tr_S_pp_S_pp_T /= r; tr_S_ap_T_Da_S_ap /= r;
        tr_S_aa /= r; tr_S_pa_Da_inv_S_pa_T /= r; tr_S_aa_T_Da_S_aa_Da_inv /= r;

        double c_pp = n_points - 2.0 * tr_S_pp + tr_S_pp_S_pp_T;
        double c_pa = tr_S_pa_Da_inv_S_pa_T;
        double c_ap = tr_S_ap_T_Da_S_ap;
        double c_aa = n_areal - 2.0 * tr_S_aa + tr_S_aa_T_Da_S_aa_Da_inv;

        edfs = static_cast<double>(n_points + n_areal) - (c_pp + c_ap + c_pa + c_aa);

        double den = c_pp * c_aa - c_pa * c_ap;
        if (std::abs(den) < 1e-12) den = 1e-8;

        sigma_sq_p_new = (c_aa * Q_p - c_pa * Q_a) / den;
        sigma_sq_A_new = (c_pp * Q_a - c_ap * Q_p) / den;

        double data_loss_p = Q_p / sigma_sq_p;
        double data_loss_A = Q_a / sigma_sq_A;
        double data_loss = (data_loss_p + data_loss_A) / ((n_points + n_areal) * domain_area);

        Jold = Jnew;
        Jnew = data_loss + ftPf(lambda, g, R0);
        ++n_iter_;

        obj_history.push_back(Jnew);
        loss_p_history.push_back(data_loss_p / (n_points * domain_area));
        loss_A_history.push_back(data_loss_A / (n_areal * domain_area));
        loss_global_history.push_back(data_loss);
        sigma_sq_p_history.push_back(sigma_sq_p_new);
        sigma_sq_A_history.push_back(sigma_sq_A_new);

        sigma_sq_p = sigma_sq_p_new;
        sigma_sq_A = sigma_sq_A_new;
    }

    return {f, g, fitted, sigma_sq_p, sigma_sq_A, edfs, n_iter_, obj_history,
            loss_p_history, loss_A_history, loss_global_history, sigma_sq_p_history, sigma_sq_A_history};
}


// Hetero (v8): Iterative GCV estimation (Xu, Geophys. J. Int., 2009, Section 3)
HeteroSRPDEResults run_hetero8_srpde(unsigned int n_points, unsigned int n_areal, unsigned int n_dofs,
                                     const sparse_matrix_t& D_matrix, const sparse_matrix_t& Psi,
                                     const sparse_matrix_t& R1, const sparse_matrix_t& R0,
                                     const vector_t& y, const vector_t& u,
                                     const std::pair<double, double>& lasco_bounds, double lasco_step,
                                     double fine_half_window, double fine_step,
                                     const ModelConfig& cfg,
                                     std::vector<double>& final_lambda_seq_out) {

    // lasco_bounds: bounds for the coarse grid search of lambda
    // lasco_step: step size for the coarse grid search of lambda
    // fine_half_window: half window size for the fine grid search of lambda
    // fine_step: step size for the fine grid search of lambda
    // final_lambda_seq_out: output vector containing the final sequence of lambda values used in the optimization

    const unsigned int n = n_points + n_areal;
    std::vector<double> lasco_seq = build_log10_seq(lasco_bounds.first, lasco_bounds.second, lasco_step);

    // STEP 0.1: Inizializzazione delle varianze tramite stima empirica 
    double sigma_sq_p = 0.0;
    double sigma_sq_A = 0.0;

    if (cfg.init_with_residuals) {
        // Varianza emprica dei residui. Per ottenere i residui, fittiamo i modelli 
        // "only_point" e "only_area" separatamente, con lambda scelto tramite GCV
        double best_gcv_p = std::numeric_limits<double>::max();
        double lambda_p_init = lasco_seq.front();
        for (double lam : lasco_seq) {
            double score_p = ic_score_fun("gcv", "only_point", n_points, n_areal, n_dofs,
                                          D_matrix, Psi, R1, R0, y, u, lam, cfg);
            if (score_p < best_gcv_p) {
                best_gcv_p = score_p;
                lambda_p_init = lam;
            }
        }
        auto res_p = dispatch_model("only_point", n_points, n_areal, n_dofs,
                                    D_matrix, Psi, R1, R0, y, u, lambda_p_init, cfg);
        sigma_sq_p = empirical_point_variance(y.head(n_points) - res_p.fitted);

        double best_gcv_a = std::numeric_limits<double>::max();
        double lambda_a_init = lasco_seq.front();
        for (double lam : lasco_seq) {
            double score_a = ic_score_fun("gcv", "only_area", n_points, n_areal, n_dofs,
                                          D_matrix, Psi, R1, R0, y, u, lam, cfg);
            if (score_a < best_gcv_a) {
                best_gcv_a = score_a;
                lambda_a_init = lam;
            }
        }
        auto res_a = dispatch_model("only_area", n_points, n_areal, n_dofs,
                                    D_matrix, Psi, R1, R0, y, u, lambda_a_init, cfg);
        sigma_sq_A = empirical_areal_variance(y.tail(n_areal) - res_a.fitted, D_matrix, n_points);
    } else {
        // Varianza empirica della variabile risposta y 
        sigma_sq_p = empirical_point_variance(y.head(n_points));
        sigma_sq_A = empirical_areal_variance(y.tail(n_areal), D_matrix, n_points);
    }

    // inital rho guess
    double rho = sigma_sq_p / sigma_sq_A;  

    // STEP 0.2: inizializzazione del lambda tramite (weighted) GCV 
    vector_t gcv_scores(lasco_seq.size());
    for (size_t idx = 0; idx < lasco_seq.size(); ++idx) {
        gcv_scores[idx] = compute_xu_weighted_gcv(n_points, n_areal, n_dofs, cfg.domain_area,
                                                  D_matrix, Psi, R1, R0, y, u,
                                                  lasco_seq[idx], rho, cfg.random_seed);
    }
    Eigen::Index min_idx = 0;
    gcv_scores.minCoeff(&min_idx);
    double lambda = lasco_seq[min_idx];
    final_lambda_seq_out = lasco_seq;


    vector_t f, g, fitted;
    double edfs = 0.0;
    unsigned int n_iter = 0;
    bool converged = false;

    std::vector<double> obj_history, loss_p_history, loss_A_history, loss_global_history;
    std::vector<double> sigma_sq_p_history, sigma_sq_A_history, lambda_history;

    // CICLO ITERATIVO: STEP 1, STEP 2 (su griglia fine centrata nella stima corrente del lambda)
    while (n_iter < cfg.max_iter && !converged) {

        std::cout << "  ---Iteration #" << n_iter + 1 << "..." << std::endl;

        double sigma_sq_p_old = sigma_sq_p;
        double sigma_sq_A_old = sigma_sq_A;
        double rho_old = rho;
        double lambda_old = lambda;


        // ------STEP 1: Fitta il modello pesato, e stima varianze tramite formula (39)------

        // Calcolo dei pesi 
        sparse_matrix_t W = build_block_W(n_points, n_areal, cfg.domain_area, 1.0, rho);

        // Fit del modello SRPDE pesato 
        SparseBlockMatrix<double, 2, 2> A(-Psi.transpose() * D_matrix * W * Psi, lambda * R1.transpose(), R1, R0);
        Eigen::SparseLU<sparse_matrix_t> invA(A);
        vector_t rhs = vector_t::Zero(2 * n_dofs);
        rhs.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * y;
        rhs.bottomRows(n_dofs) = u;
        vector_t sol = invA.solve(rhs);
        f = sol.topRows(n_dofs);
        g = sol.bottomRows(n_dofs);
        fitted = Psi * f;
        
        // Calcolo bias-corrected residulas 
        vector_t r_lambda = y - fitted;   // residuals 
        vector_t rhs_bias = vector_t::Zero(2 * n_dofs);
        rhs_bias.bottomRows(n_dofs) = R1 * f;
        vector_t lambda_C_f = invA.solve(rhs_bias).topRows(n_dofs);  // calcolo di lambda_C_f con trick computazionale (vedi GoodNotes)
        vector_t r_bar = r_lambda - Psi * lambda_C_f;  // basis-corrected residuals 

        // Aggiornamento delle varianze (e rho)
        double num_p = r_bar.head(n_points).squaredNorm();
        double num_a = (r_bar.tail(n_areal).array().square() * D_matrix.diagonal().tail(n_areal).array()).sum();

        double trace_corr_p = edf_xu_block(D_matrix, W, Psi, R1, n_dofs, invA,
                                           EdfTarget::Points, n_points, 100, cfg.random_seed_p);
        double trace_corr_a = edf_xu_block(D_matrix, W, Psi, R1, n_dofs, invA,
                                           EdfTarget::Areal, n_points, 100, cfg.random_seed_A);

        double den_p = static_cast<double>(n_points) - trace_corr_p;
        double den_a = static_cast<double>(n_areal)  - trace_corr_a;

        sigma_sq_p = num_p / den_p;
        sigma_sq_A = num_a / den_a;
        rho = sigma_sq_p / sigma_sq_A;


        // ------STEP 2: aggiornamento del lambda------
        
        // Costruzione griglia fine centrata nel lambda corrente
        double log10_center = std::log10(lambda);
        double exp_min = std::max(lasco_bounds.first,  log10_center - fine_half_window);
        double exp_max = std::min(lasco_bounds.second, log10_center + fine_half_window);
        std::vector<double> fine_seq = build_log10_seq(exp_min, exp_max, fine_step);

        // Calcolo weighted GCV (eq. 28a)
        gcv_scores.resize(fine_seq.size());
        for (size_t idx = 0; idx < fine_seq.size(); ++idx) {
            gcv_scores[idx] = compute_xu_weighted_gcv(n_points, n_areal, n_dofs, cfg.domain_area,
                                                      D_matrix, Psi, R1, R0, y, u,
                                                      fine_seq[idx], rho, cfg.random_seed);
        }
        gcv_scores.minCoeff(&min_idx);
        // update lambda
        lambda = fine_seq[min_idx];
        final_lambda_seq_out = fine_seq;

        // Calcolo della funzione obiettivo
        double RSS_p = r_lambda.head(n_points).squaredNorm();
        double RSS_a = (r_lambda.tail(n_areal).array().square() * D_matrix.diagonal().tail(n_areal).array()).sum();
        double data_loss = (RSS_p + rho_old * RSS_a) / (static_cast<double>(n) * cfg.domain_area);
        double J_curr = data_loss + ftPf(lambda_old, g, R0);

        obj_history.push_back(J_curr);
        loss_p_history.push_back(RSS_p / (static_cast<double>(n_points) * cfg.domain_area));
        loss_A_history.push_back((rho_old * RSS_a) / (static_cast<double>(n_areal) * cfg.domain_area));
        loss_global_history.push_back(data_loss);
        sigma_sq_p_history.push_back(sigma_sq_p);
        sigma_sq_A_history.push_back(sigma_sq_A);
        lambda_history.push_back(lambda);

        // Controllo convergenza
        double diff_sigma_p = std::abs(sigma_sq_p - sigma_sq_p_old);
        double diff_sigma_A = std::abs(sigma_sq_A - sigma_sq_A_old);
        double diff_lambda  = std::abs(lambda - lambda_old);

        n_iter++;

        // criterio sia su sigma^2 che su lambda
        // if (std::max(diff_sigma_p, diff_sigma_A) <= cfg.tol && diff_lambda <= cfg.tol) {
        //     converged = true;
        // }
        // criterio solo su sigma^2
        if (std::max(diff_sigma_p, diff_sigma_A) <= cfg.tol) {
            converged = true;
        }
    
    }

    // At convergence: compute final results 
    sparse_matrix_t W_final = build_block_W(n_points, n_areal, cfg.domain_area, 1.0, rho);
    SparseBlockMatrix<double, 2, 2> A_final(-Psi.transpose() * D_matrix * W_final * Psi, lambda * R1.transpose(), R1, R0);
    Eigen::SparseLU<sparse_matrix_t> invA_final(A_final);

    vector_t rhs_final = vector_t::Zero(2 * n_dofs);
    rhs_final.topRows(n_dofs) = -Psi.transpose() * D_matrix * W_final * y;
    rhs_final.bottomRows(n_dofs) = u;

    vector_t sol_final = invA_final.solve(rhs_final);
    f = sol_final.topRows(n_dofs);
    g = sol_final.bottomRows(n_dofs);
    fitted = Psi * f;
    edfs = edf(D_matrix, W_final, Psi, n_dofs, invA_final, EdfTarget::Global, n_points, 100, cfg.random_seed);

    return {f, g, fitted, sigma_sq_p, sigma_sq_A, edfs, n_iter, obj_history,
            loss_p_history, loss_A_history, loss_global_history,
            sigma_sq_p_history, sigma_sq_A_history, lambda, gcv_scores, lambda_history};
}

// Naive: W = 1/(n * |Omega|) * I
NaiveSRPDEResults run_naive_srpde(unsigned int n_points, unsigned int n_areal, unsigned int n_dofs,
                                  double domain_area, const sparse_matrix_t& D_matrix, const sparse_matrix_t& Psi,
                                  const sparse_matrix_t& R1, const sparse_matrix_t& R0,
                                  const vector_t& y, const vector_t& u, double lambda,
                                  unsigned int random_seed, bool verbose = false) {

    sparse_matrix_t W = build_block_W(n_points, n_areal, domain_area, 1.0, 1.0);

    SparseBlockMatrix<double, 2, 2> A(-Psi.transpose() * D_matrix * W * Psi, lambda * R1.transpose(), R1, R0);
    Eigen::SparseLU<sparse_matrix_t> invA(A);

    vector_t rhs = vector_t::Zero(2 * n_dofs, 1);
    rhs.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * y;
    rhs.bottomRows(n_dofs) = u;

    vector_t x = invA.solve(rhs);
    vector_t f = x.topRows(n_dofs);
    vector_t g = x.bottomRows(n_dofs);
    vector_t fitted = Psi * f;

    double edfs = edf(D_matrix, W, Psi, n_dofs, invA, EdfTarget::Global, n_points, 100, random_seed);
    double sigma_sq = (fitted - y).squaredNorm() / (n_points + n_areal - edfs);

    return {f, g, fitted, sigma_sq, edfs};
}


/****************************************** UNIFIED DISPATCHER & I/O ******************************************/

inline HeteroSRPDEResults dispatch_model(const std::string& method,
                                         unsigned int n_points, unsigned int n_areal, unsigned int n_dofs,
                                         const sparse_matrix_t& D_matrix, const sparse_matrix_t& Psi,
                                         const sparse_matrix_t& R1, const sparse_matrix_t& R0,
                                         const vector_t& y, const vector_t& u, double lambda,
                                         const ModelConfig& cfg) {

    const bool needs_initialization =
        method == "hetero"  ||
        method == "hetero2" ||
        method == "hetero3" ||
        method == "hetero4" ||
        method == "hetero5" ||
        method == "hetero6" ||
        method == "hetero7";

    double init_sigma2p = 0.0;
    double init_sigma2A = 0.0;
    double init_wA = 0.0;

    if (needs_initialization) {
        if (cfg.init_with_residuals) {
            sparse_matrix_t D_p = D_matrix.topLeftCorner(n_points, n_points);
            sparse_matrix_t Psi_p = Psi.topRows(n_points);
            vector_t y_p = y.head(n_points);
            auto res_p = run_naive_srpde(n_points, 0, n_dofs, cfg.domain_area, D_p, Psi_p, R1, R0, y_p, u, lambda, cfg.random_seed, cfg.verbose);
            vector_t residuals_p = y_p - res_p.fitted;

            sparse_matrix_t D_a = D_matrix.bottomRightCorner(n_areal, n_areal);
            sparse_matrix_t Psi_a = Psi.bottomRows(n_areal);
            vector_t y_a = y.tail(n_areal);
            auto res_a = run_naive_srpde(0, n_areal, n_dofs, cfg.domain_area, D_a, Psi_a, R1, R0, y_a, u, lambda, cfg.random_seed, cfg.verbose);
            vector_t residuals_a = y_a - res_a.fitted;

            init_sigma2p = empirical_point_variance(residuals_p);
            init_sigma2A = empirical_areal_variance(residuals_a, D_matrix, n_points);
        } else {
            init_sigma2p = empirical_point_variance(y.head(n_points));
            init_sigma2A = empirical_areal_variance(y.tail(n_areal), D_matrix, n_points);
        }

        init_wA = init_sigma2p / init_sigma2A;
    }

    if (method == "hetero") {
        return run_hetero_srpde(init_wA, n_points, n_areal, n_dofs, cfg.domain_area, D_matrix, Psi, R1, R0, y, u, lambda,
                                cfg.random_seed_p, cfg.random_seed_A, cfg.random_seed, cfg.tol, cfg.max_iter, cfg.verbose);
    } else if (method == "hetero2") {
        return run_hetero2_srpde(init_wA, n_points, n_areal, n_dofs, cfg.domain_area, D_matrix, Psi, R1, R0, y, u, lambda,
                                 cfg.random_seed_p, cfg.random_seed_A, cfg.random_seed, cfg.tol, cfg.max_iter, cfg.verbose);
    } else if (method == "hetero3") {
        return run_hetero3_srpde(init_wA, n_points, n_areal, n_dofs, cfg.domain_area, D_matrix, Psi, R1, R0, y, u, lambda,
                                 cfg.random_seed_p, cfg.random_seed_A, cfg.random_seed, cfg.tol, cfg.max_iter, cfg.verbose);
    } else if (method == "hetero4") {
        return run_hetero4_srpde(init_sigma2p, init_sigma2A, n_points, n_areal, n_dofs, cfg.domain_area, D_matrix, Psi, R1, R0, y, u, lambda,
                                 cfg.random_seed_p, cfg.random_seed_A, cfg.random_seed, cfg.tol, cfg.max_iter, cfg.verbose);
    } else if (method == "hetero5") {
        return run_hetero5_srpde(init_sigma2p, init_sigma2A, n_points, n_areal, n_dofs, cfg.domain_area, D_matrix, Psi, R1, R0, y, u, lambda,
                                 cfg.random_seed_p, cfg.random_seed_A, cfg.random_seed, cfg.tol, cfg.max_iter, cfg.verbose);
    } else if (method == "hetero6") {
        return run_hetero6_srpde(init_wA, n_points, n_areal, n_dofs, cfg.domain_area, D_matrix, Psi, R1, R0, y, u, lambda,
                                 cfg.random_seed_p, cfg.random_seed_A, cfg.random_seed, cfg.tol, cfg.max_iter, cfg.verbose);
    } else if (method == "hetero7") {
        return run_hetero7_srpde(init_sigma2p, init_sigma2A, n_points, n_areal, n_dofs, cfg.domain_area, D_matrix, Psi, R1, R0, y, u, lambda,
                                 cfg.random_seed_p, cfg.random_seed_A, cfg.random_seed, cfg.tol, cfg.max_iter, cfg.verbose);
    } else if (method == "naive") {
        auto res = run_naive_srpde(n_points, n_areal, n_dofs, cfg.domain_area, D_matrix, Psi, R1, R0, y, u, lambda, cfg.random_seed, cfg.verbose);
        return {res.f, res.g, res.fitted, res.sigma_sq, res.sigma_sq, res.edfs, 1, {}, {}, {}, {}, {}, {}};
    } else if (method == "only_point") {
        sparse_matrix_t D_p = D_matrix.topLeftCorner(n_points, n_points);
        sparse_matrix_t Psi_p = Psi.topRows(n_points);
        vector_t y_p = y.head(n_points);
        auto res = run_naive_srpde(n_points, 0, n_dofs, cfg.domain_area, D_p, Psi_p, R1, R0, y_p, u, lambda, cfg.random_seed, cfg.verbose);
        return {res.f, res.g, res.fitted, res.sigma_sq, res.sigma_sq, res.edfs, 1, {}, {}, {}, {}, {}, {}};
    } else if (method == "only_area") {
        sparse_matrix_t D_a = D_matrix.bottomRightCorner(n_areal, n_areal);
        sparse_matrix_t Psi_a = Psi.bottomRows(n_areal);
        vector_t y_a = y.tail(n_areal);
        auto res = run_naive_srpde(0, n_areal, n_dofs, cfg.domain_area, D_a, Psi_a, R1, R0, y_a, u, lambda, cfg.random_seed, cfg.verbose);
        return {res.f, res.g, res.fitted, res.sigma_sq, res.sigma_sq, res.edfs, 1, {}, {}, {}, {}, {}, {}};
    }
    throw std::invalid_argument("Unknown method: " + method);
}




inline void save_fit_results(const std::string& dir_path, const std::string& suffix,
                             const HeteroSRPDEResults& res, bool is_single_variance) {
    std::string sfx = suffix.empty() ? "" : "_" + suffix;
    write_csv(dir_path + "f" + sfx + ".csv", res.f);
    write_csv(dir_path + "fn" + sfx + ".csv", res.fitted);

    if (is_single_variance) {
        save_scalar(dir_path + "sigma_sq" + sfx + ".csv", res.sigma_sq_p);
    } else {
        save_scalar(dir_path + "sigma_sq_p" + sfx + ".csv", res.sigma_sq_p);
        save_scalar(dir_path + "sigma_sq_a" + sfx + ".csv", res.sigma_sq_A);
        if (!res.obj_history.empty()) {
            write_csv(dir_path + "obj_history" + sfx + ".csv", res.obj_history);
            write_csv(dir_path + "loss_p_history" + sfx + ".csv", res.loss_p_history);
            write_csv(dir_path + "loss_a_history" + sfx + ".csv", res.loss_A_history);
            write_csv(dir_path + "loss_global_history" + sfx + ".csv", res.loss_global_history);
            write_csv(dir_path + "sigma_sq_p_history" + sfx + ".csv", res.sigma_sq_p_history);
            write_csv(dir_path + "sigma_sq_a_history" + sfx + ".csv", res.sigma_sq_A_history);
        }
    }
}


/****************************************** LAMBDA SELECTION CRITERIA ******************************************/

inline double compute_ic_metric(const std::string& criterion, unsigned int n_obs, double rss, double edfs) {
    double n = static_cast<double>(n_obs);
    double dor = n - edfs;

    if (criterion == "gcv") {
        return (n / (dor * dor)) * rss;
    } else if (criterion == "aicc") {
        double den = n - edfs - 2.0;
        return std::log(rss / n) + 1.0 + (2.0 * (edfs + 1.0)) / den;
    } else if (criterion == "gfaic") {
        return (rss / n) * std::exp(2.0 * edfs / n);
    } else if (criterion == "rcv") {
        double tr_n = edfs / n;
        double num_factor = (1.0 + tr_n) * (1.0 + tr_n);
        double den_factor = 1.0 - tr_n;
        return (rss / n) * (num_factor / (den_factor * den_factor));
    } else if (criterion == "fpe") {
        return (rss / n) * ((n + edfs) / dor);
    }
    throw std::invalid_argument("Unknown IC criterion: " + criterion);
}


inline double ic_score_fun(const std::string& criterion, const std::string& method,
                           unsigned int n_points, unsigned int n_areal, unsigned int n_dofs,
                           const sparse_matrix_t& D_matrix, const sparse_matrix_t& Psi,
                           const sparse_matrix_t& R1, const sparse_matrix_t& R0,
                           const vector_t& y, const vector_t& u, double lambda,
                           const ModelConfig& cfg) {

    auto res = dispatch_model(method, n_points, n_areal, n_dofs, D_matrix, Psi, R1, R0, y, u, lambda, cfg);

    unsigned int n_active = n_points + n_areal;
    double rss = 0.0;

    if (method == "only_point") {
        n_active = n_points;
        rss = (res.fitted - y.head(n_points)).squaredNorm();
    } else if (method == "only_area") {
        n_active = n_areal;
        vector_t res_a = res.fitted - y.tail(n_areal);
        rss = (res_a.array().square() * D_matrix.diagonal().tail(n_areal).array()).sum();
    } else if (method == "naive") {
        vector_t res_all = res.fitted - y;
        rss = (res_all.array().square() * D_matrix.diagonal().array()).sum();
    } else {
        vector_t res_all = res.fitted - y;
        double rss_p = res_all.head(n_points).squaredNorm();
        double rss_a = (res_all.tail(n_areal).array().square() * D_matrix.diagonal().tail(n_areal).array()).sum();
        double rho = res.sigma_sq_p / res.sigma_sq_A;
        rss = rss_p + rho * rss_a;
    }

    return compute_ic_metric(criterion, n_active, rss, res.edfs);
}


inline double k_fold_cv_score_fun(const std::string& method,
                                  unsigned int n_points, unsigned int n_areal, unsigned int n_dofs,
                                  const sparse_matrix_t& D_matrix, const sparse_matrix_t& Psi,
                                  const sparse_matrix_t& R1, const sparse_matrix_t& R0,
                                  const vector_t& y, const vector_t& u, double lambda,
                                  const std::vector<int>& fold_ids, int K,
                                  const ModelConfig& cfg) {

    unsigned int n_p_eff = (method == "only_area") ? 0 : n_points;
    unsigned int n_a_eff = (method == "only_point") ? 0 : n_areal;
    std::string fit_method = (method == "only_point" || method == "only_area") ? "naive" : method;

    double total_cv_error = 0.0;

    for (int k = 1; k <= K; ++k) {
        std::vector<int> train_idxs_p, train_idxs_a, test_idxs_p, test_idxs_a;
        for (unsigned int i = 0; i < n_p_eff; ++i) {
            if (fold_ids[i] != k) train_idxs_p.push_back(i);
            else test_idxs_p.push_back(i);
        }
        for (unsigned int i = 0; i < n_a_eff; ++i) {
            unsigned int global_idx = n_points + i;
            if (fold_ids[global_idx] != k) train_idxs_a.push_back(global_idx);
            else test_idxs_a.push_back(global_idx);
        }

        unsigned int n_p_train = train_idxs_p.size();
        unsigned int n_a_train = train_idxs_a.size();
        unsigned int n_train = n_p_train + n_a_train;

        vector_t y_train(n_train);
        sparse_matrix_t D_train(n_train, n_train);
        sparse_matrix_t Psi_train(n_train, Psi.cols());
        std::vector<Eigen::Triplet<double>> D_triplets, Psi_triplets;

        int row_idx = 0;
        for (int idx : train_idxs_p) {
            y_train(row_idx) = y(idx);
            D_triplets.emplace_back(row_idx, row_idx, D_matrix.coeff(idx, idx));
            for (sparse_matrix_t::InnerIterator it(Psi, idx); it; ++it) {
                Psi_triplets.emplace_back(row_idx, it.col(), it.value());
            }
            row_idx++;
        }
        for (int idx : train_idxs_a) {
            y_train(row_idx) = y(idx);
            D_triplets.emplace_back(row_idx, row_idx, D_matrix.coeff(idx, idx));
            for (sparse_matrix_t::InnerIterator it(Psi, idx); it; ++it) {
                Psi_triplets.emplace_back(row_idx, it.col(), it.value());
            }
            row_idx++;
        }

        D_train.setFromTriplets(D_triplets.begin(), D_triplets.end());
        Psi_train.setFromTriplets(Psi_triplets.begin(), Psi_triplets.end());

        vector_t f_train = dispatch_model(fit_method, n_p_train, n_a_train, n_dofs,
                                          D_train, Psi_train, R1, R0, y_train, u, lambda, cfg).f;

        vector_t y_pred = Psi * f_train;
        for (int idx : test_idxs_p) {
            total_cv_error += std::pow(y(idx) - y_pred(idx), 2);
        }
        for (int idx : test_idxs_a) {
            total_cv_error += std::pow(y(idx) - y_pred(idx), 2) * D_matrix.coeff(idx, idx);
        }
    }

    return total_cv_error / (n_p_eff + n_a_eff);
}

inline double evaluate_criterion(const std::string& criterion, const std::string& method,
                                 unsigned int n_points, unsigned int n_areal, unsigned int n_dofs,
                                 const sparse_matrix_t& D_matrix, const sparse_matrix_t& Psi,
                                 const sparse_matrix_t& R1, const sparse_matrix_t& R0,
                                 const vector_t& y, const vector_t& u, double lambda,
                                 const std::vector<int>& fold_ids, int K,
                                 const ModelConfig& cfg) {
    if (criterion == "cv") {
        return k_fold_cv_score_fun(method, n_points, n_areal, n_dofs, D_matrix, Psi, R1, R0, y, u, lambda, fold_ids, K, cfg);
    }
    return ic_score_fun(criterion, method, n_points, n_areal, n_dofs, D_matrix, Psi, R1, R0, y, u, lambda, cfg);
}






#endif
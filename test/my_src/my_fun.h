#ifndef MY_FUN_H
#define MY_FUN_H

#include <fdaPDE/fdapde.h>

#include <iostream>
#include <string>
#include <vector>

using namespace fdapde;
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using sparse_matrix_t = Eigen::SparseMatrix<double>;
using fdapde::GridSearch; 


enum class EdfTarget { Global, Points, Areal };

// NOTA: per ora e' solo NONparametrico !! 
template <typename D_matrixType, typename W_type, typename InvAType>
double edf(const D_matrixType& D_matrix, const W_type& W_, const sparse_matrix_t& Psi, int n_dofs, 
           const InvAType& invA_, EdfTarget target = EdfTarget::Global, int n_points = 0, int r = 100, int seed = 109) {

    // n_points: number of pointise locations         

    int n_locs = Psi.rows();
    matrix_t Us = matrix_t::Zero(n_locs, r);
    int seed_ = (seed == random_seed) ? std::random_device()() : seed;
    std::mt19937 rng(seed_);
    rademacher_distribution rademacher;

    int start_idx = 0;
    int end_idx = n_locs;

    // ATT assuming points stacked ABOVE areas !! 
    if (target == EdfTarget::Points) {
        end_idx = n_points;
    } else if (target == EdfTarget::Areal) {
        start_idx = n_points;
    }

    // fill the randemacher vector only where needed 
    for (int i = start_idx; i < end_idx; ++i) {
        for (int j = 0; j < r; ++j) { Us(i, j) = rademacher(rng); }
    }

    matrix_t Ys = Us.transpose() * Psi;
    matrix_t Bs = matrix_t::Zero(2 * n_dofs, r);
    Bs.topRows(n_dofs) = -Psi.transpose() * D_matrix * W_ * Us;
        
    matrix_t x = invA_.solve(Bs);
    double trS = 0;
    for (int i = 0; i < r; ++i) { trS += Ys.row(i).dot(x.col(i).head(n_dofs)); }
    return trS / r;
}

double ftPf(double lambda, const vector_t& g_, const sparse_matrix_t& R0_) {   
    return lambda * g_.dot(R0_ * g_);
}


struct HeteroSRPDEResults {
    vector_t f;
    vector_t g;
    vector_t fitted;
    double sigma_sq_p;
    double sigma_sq_A;
    double edfs; 
    unsigned int n_iter;
    std::vector<double> obj_history;
    std::vector<double> loss_p_history;
    std::vector<double> loss_A_history;
    std::vector<double> loss_global_history;
    std::vector<double> sigma_sq_p_history;
    std::vector<double> sigma_sq_A_history;
};

struct NaiveSRPDEResults {
    vector_t f;
    vector_t g;
    vector_t fitted;
    double sigma_sq;
    double edfs; 
};

HeteroSRPDEResults run_hetero_srpde(unsigned int n_points, unsigned int n_areal, unsigned int n_dofs, 
                         const sparse_matrix_t& D_matrix, const sparse_matrix_t& Psi, 
                         const sparse_matrix_t& R1, const sparse_matrix_t& R0, 
                         const vector_t& y, const vector_t& u, double lambda,
                         unsigned int random_seed_p, unsigned int random_seed_A, unsigned int random_seed,
                         const double tol_ = 1e-6, const unsigned int max_iter_ = 50, bool verbose = false) {

    // Weight matrix 
    // function to construcut the weight matrix: takes w_A as input and returns W = diag([1, ..., 1, w_A, ..., w_A])
    // w_A: ratio sigma2_p / sigma2_A
    auto construct_W = [&](double w_A) {
        std::vector<Eigen::Triplet<double>> triplets_W;
        for (int i = 0; i < n_points; ++i) {
            triplets_W.push_back(Eigen::Triplet<double>(i, i, 1.0));
        }
        for (int i = 0; i < n_areal; ++i) {
            triplets_W.push_back(Eigen::Triplet<double>(i + n_points, i + n_points, w_A));
        }
        sparse_matrix_t W(n_points + n_areal, n_points + n_areal);
        W.setFromTriplets(triplets_W.begin(), triplets_W.end());
        return W;
    };
    double w_A = 1.0;     // initial guess for the variance ratio
    sparse_matrix_t W = construct_W(w_A);  

    // store room for results 
    double sigma_sq_p, sigma_sq_A;  
    vector_t f, g, fitted;

    // iterative loop
    double Jold = std::numeric_limits<double>::max(), Jnew = 0;
    unsigned int n_iter_ = 0;
    std::vector<double> obj_history;
    std::vector<double> loss_p_history;
    std::vector<double> loss_A_history;
    std::vector<double> loss_global_history;
    std::vector<double> sigma_sq_p_history;
    std::vector<double> sigma_sq_A_history;
    while (n_iter_ < max_iter_ && std::abs(Jnew - Jold) > tol_) {

        if (verbose) {
            std::cout << "Iteration " << n_iter_ << std::endl;
        }

        // Step 1: solve weighted SRPDE
        SparseBlockMatrix<double, 2, 2> A(-Psi.transpose() * D_matrix * W * Psi, lambda * R1.transpose(), lambda * R1, lambda * R0);
        Eigen::SparseLU<sparse_matrix_t> invA(A);

        vector_t rhs = vector_t::Zero(2 * n_dofs, 1);
        rhs.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * y;
        rhs.bottomRows(n_dofs) = lambda * u;

        vector_t x = invA.solve(rhs); // expansion coefficient vector
        f = x.topRows(n_dofs);
        g = x.bottomRows(n_dofs);
        fitted = Psi * f;

        
        // Step 2: update weights: 
        
        // - 2.1 compute edf_p and edf_A. Note: the function edf() computes the global effective degrees of freedom. 
        // We need to mask areal data for edf_p and point data for edf_A.
        double edf_p = edf(D_matrix, W, Psi, n_dofs, invA, EdfTarget::Points, n_points, 100, random_seed_p*n_iter_);
        double edf_A = edf(D_matrix, W, Psi, n_dofs, invA, EdfTarget::Areal, n_points, 100, random_seed_A*n_iter_);
        
        // - 2.2 compute sigma_sq_p an sigma_sq_A  --> ATT: assumes that points are stacked ABOVE areas in the data vector y !!
        sigma_sq_p = (fitted.head(n_points) - y.head(n_points)).squaredNorm() / (n_points-edf_p);
        // sigma_sq_A = (fitted.tail(n_areal) - y.tail(n_areal)).squaredNorm() / (n_areal-edf_A);
        sigma_sq_A = ( (fitted.tail(n_areal) - y.tail(n_areal)).array().square() * D_matrix.diagonal().tail(n_areal).array() ).sum() / (n_areal - edf_A);
        // --> ATT: metto *D in modo che sia una somma di residui pesati per le aree, cosicché questa sia effettivamente
        // quella che nella teoria sto chiamando sigma_sq_a...


        // - 2.3 compute w_A = sigma_sq_p / sigma_sq_A
        w_A = sigma_sq_p / sigma_sq_A;
        
        // - 2.4 construct new weight matrix W
        W = construct_W(w_A); 

        // compute data loss
        vector_t res = y - fitted;
        double data_loss_p = res.head(n_points).squaredNorm();
        double data_loss_A = w_A * (res.tail(n_areal).array().square() * D_matrix.diagonal().tail(n_areal).array()).sum();
        double data_loss = (data_loss_p + data_loss_A) / (n_points + n_areal);  

        // prepare for next iteration
        Jold = Jnew;
        Jnew = data_loss + ftPf(lambda, g, R0); 
        n_iter_++;

        if (verbose) {
            std::cout << "Absolute difference in objective: " << std::abs(Jnew - Jold) << std::endl;
            std::cout << "--------------" << std::endl;
        }

        obj_history.push_back(Jnew);
        loss_p_history.push_back(data_loss_p);
        loss_A_history.push_back(data_loss_A);
        loss_global_history.push_back(data_loss);
        sigma_sq_p_history.push_back(sigma_sq_p);
        sigma_sq_A_history.push_back(sigma_sq_A);
    }

    // compute global edfs --> outside the loop for efficiency (global edf are not needed inside the loop)
    SparseBlockMatrix<double, 2, 2> A_conv(-Psi.transpose() * D_matrix * W * Psi, lambda * R1.transpose(), lambda * R1, lambda * R0);
    Eigen::SparseLU<sparse_matrix_t> invA_conv(A_conv);
    double edfs = edf(D_matrix, W, Psi, n_dofs, invA_conv, EdfTarget::Global, n_points, 100, random_seed);

    // returns results 
    return{f, g, fitted, sigma_sq_p, sigma_sq_A, edfs, n_iter_, obj_history, 
           loss_p_history, loss_A_history, loss_global_history, 
           sigma_sq_p_history, sigma_sq_A_history};

}



NaiveSRPDEResults run_naive_srpde(unsigned int n_points, unsigned int n_areal, unsigned int n_dofs, 
                         const sparse_matrix_t& D_matrix, const sparse_matrix_t& Psi, 
                         const sparse_matrix_t& R1, const sparse_matrix_t& R0, 
                         const vector_t& y, const vector_t& u, double lambda, unsigned int random_seed, bool verbose = true) {

    // set W = I 
    std::vector<Eigen::Triplet<double>> triplets_W;
    for (int i = 0; i < n_points; ++i) {
        triplets_W.push_back(Eigen::Triplet<double>(i, i, 1.0));
    }
    for (int i = 0; i < n_areal; ++i) {
        triplets_W.push_back(Eigen::Triplet<double>(i + n_points, i + n_points, 1.0));
    }
    sparse_matrix_t W(n_points + n_areal, n_points + n_areal);
    W.setFromTriplets(triplets_W.begin(), triplets_W.end());


    // store room for results 
    double sigma_sq_p, sigma_sq_A;  
    vector_t f, g, fitted;

    // solve weighted SRPDE
    SparseBlockMatrix<double, 2, 2> A(-Psi.transpose() * D_matrix * W * Psi, lambda * R1.transpose(), lambda * R1, lambda * R0);
    Eigen::SparseLU<sparse_matrix_t> invA(A);

    vector_t rhs = vector_t::Zero(2 * n_dofs, 1);
    rhs.topRows(n_dofs) = -Psi.transpose() * D_matrix * W * y;
    rhs.bottomRows(n_dofs) = lambda * u;

    vector_t x = invA.solve(rhs); // expansion coefficient vector
    f = x.topRows(n_dofs);
    g = x.bottomRows(n_dofs);
    fitted = Psi * f;

    double edfs = edf(D_matrix, W, Psi, n_dofs, invA, EdfTarget::Global, n_points, 100, random_seed); 

    // estimated variance
    double sigma_sq = (fitted - y).squaredNorm() / (n_points + n_areal - edfs);

    // return results
    return {f, g, fitted, sigma_sq, edfs};

}




// define the GCV score objective to minimize to find the optimal lambda 
double gcv_score_fun(bool run_hetero, unsigned int n_points, unsigned int n_areal, unsigned int n_dofs, 
                                const sparse_matrix_t& D_matrix, const sparse_matrix_t& Psi, 
                                const sparse_matrix_t& R1, const sparse_matrix_t& R0, 
                                const vector_t& y, const vector_t& u, double lambda,
                                unsigned int random_seed_p, unsigned int random_seed_A, unsigned int random_seed, bool verbose = true){ 

    // store room for results 
    double sigma_sq_p, sigma_sq_A;  
    vector_t f, g, fitted;
    double edfs; 

    sparse_matrix_t W(n_points + n_areal, n_points + n_areal);

    if(run_hetero){

        auto hetero_results = run_hetero_srpde(n_points, n_areal, n_dofs, 
                                        D_matrix, Psi, R1, R0, y, u, lambda, 
                                        random_seed_p, random_seed_A, random_seed, verbose); 


        edfs = hetero_results.edfs;   
        fitted = hetero_results.fitted;
   

    } else{

        auto naive_results = run_naive_srpde(n_points, n_areal, n_dofs, 
                                            D_matrix, Psi, R1, R0, y, u, lambda, 
                                            random_seed, verbose); 

        edfs = naive_results.edfs; 
        fitted = naive_results.fitted;

    }

    // Compute the GCV score: n/(n-trS) * (fitted - y)^2
    double gcv_score = (n_points + n_areal) / ( n_points + n_areal - edfs ) * (fitted - y).squaredNorm();
    
    return(gcv_score);

};



#endif
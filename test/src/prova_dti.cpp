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

#include "logger.h"
std::ofstream file;

using namespace fdapde;
using fdapde::test::almost_equal;
using vector_t = Eigen::Matrix<double, Dynamic, 1>;
using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;

double loss(double S0, const vector_t& S, const matrix_t& gradients, matrix_t L) {
    double l = 0;
    for (int i = 0; i < gradients.cols(); ++i) {
        vector_t gi = gradients.col(i);
        double Si = S[i];

        matrix_t exp_L = expm(L);
        double diff = std::log(S0 / Si) - gi.transpose() * exp_L * gi;

        l += diff * diff;

        vector_t dL = fdapde::internals::dG_exp(1, gi, L);
    }
    return l;
}

vector_t gradient(double S0, const vector_t& S, const matrix_t& gradients, matrix_t L) {
    vector_t grad_vec(3);
    for (int i = 0; i < gradients.cols(); ++i) {
        vector_t gi = gradients.col(i);
        double Si = S[i];

        matrix_t exp_L = expm(L);
        double diff = std::log(S0 / Si) - gi.transpose() * exp_L * gi;

        vector_t dL = fdapde::internals::dG_exp(1, gi, L);
        grad_vec -= 2.0 * diff * dL;
    }
    return grad_vec;
}

struct opt_functor_t {
    double operator()(const vector_t& L) { return loss(S0_, S_, gradients_, matrix_view(L)); }
    std::function<vector_t(const vector_t&)> derive() {
        return [this](const vector_t& L) { return gradient(S0_, S_, gradients_, matrix_view(L)); };
    }
    opt_functor_t(double S0, const vector_t& S, const matrix_t& gradients) : S0_(S0), S_(S), gradients_(gradients) { }
    double S0_;
    vector_t S_;
    matrix_t gradients_;
};

TEST(fdti_it, test_01) {
    // Diffusion tensor
    matrix_t D(2, 2);
    D << 0.970, 0.0, 0, 1.751;
    matrix_t L = logm(D);

    std::cout << "Diffusion tensor" << std::endl;
    std::cout << D << std::endl;
    std::cout << std::endl;
    std::cout << "Log Diffusion tensor" << std::endl;
    std::cout << L << std::endl;

    // gradient directions
    vector_t g0(2);
    g0 << 1, 0;
    vector_t g1(2);
    g1 << 0, 1;
    vector_t g2(2);
    g2 << 1, 1;
    vector_t g3(2);
    g3 << -1, 1;
    matrix_t gradients(2, 4);
    gradients.col(0) = g0.normalized();
    gradients.col(1) = g1.normalized();
    gradients.col(2) = g2.normalized();
    gradients.col(3) = g3.normalized();
    std::cout << std::endl;
    std::cout << "Gradient directions" << std::endl;
    std::cout << gradients << std::endl;

    // Simulated DWI images
    double S0 = 10.0;
    vector_t S(gradients.cols());
    for (int i = 0; i < gradients.cols(); ++i) {
        vector_t gi = gradients.col(i);
        S[i] = S0 * std::exp(-gi.dot(D * gi));
    }
    std::cout << std::endl;
    std::cout << "Simulated DWI images" << std::endl;
    std::cout << S.transpose() << std::endl;

    dwi_data data {vector_t::Ones(4), gradients, S0 * vector_t::Ones(1), matrix_t {S.transpose()}};

    // Grid resolution
    int res = 50;
    int total = res * res * res;
    int index = 0;

    // Output containers
    vector_t obj = vector_t::Zero(total);
    matrix_t grad = matrix_t::Zero(total, 3);
    matrix_t L_vec(total, 3);

    // Sampling
    for (double lxx = -1.0; lxx <= 1.0; lxx += 2.0 / res) {
        for (double lyy = -1.0; lyy <= 1.0; lyy += 2.0 / res) {
            for (double lxy = -1.0 / sqrt(2); lxy <= 1.0 / sqrt(2); lxy += sqrt(2) / res) {
                L_vec.row(index) << lxx, lyy, sqrt(2) * lxy;

                // Build symmetric matrix L from vector
                matrix_t L(2, 2);
                L << lxx, lxy, lxy, lyy;

                // Compute loss and gradient
                LossFunctor loss_fn = riccian_loss;
                obj[index] = loss_fn.loss(data, matrix_t {L_vec.row(index)});
                grad.row(index) = loss_fn.grad_loss(data, matrix_t {L_vec.row(index)}).row(0);
                index++;
            }
        }
    }

    std::string filename = "../data/vsr/tensors/RESULTS/descent.csv";
    {
        std::ofstream clearFile(filename, std::ios::trunc);
        if (!clearFile.is_open()) { std::cerr << "Error clearing file.\n"; }
        // File is cleared here
    }

    file.open(filename, std::ios::app);

    vector_t initialization(3);   // = vector_t::Zero(3);
    initialization << -0.5, -0.8, 0.8;
    // initialization.normalized();

    std::cout << std::endl;
    BFGS<Dynamic, BacktrackingLineSearch> opt {50000, 1e-5, 1e-3};   // BacktrackingLineSearch
    vector_t L_opt = opt.optimize(
      opt_functor_t(S0, S, gradients),   //
      initialization                     // ,                                           //
                                         // [](auto value) { std::cout << value << ", " << std ::endl; }   //
    );
    std::cout << std::endl;
    std::cout << "Optimizer result" << std::endl;
    // std::cout << matrix_view(L_opt) << std::endl;
    // std::cout << "loss value achieved: " << loss(S0, S, gradients, L_opt) << std::endl;

    write_csv("../data/vsr/tensors/RESULTS/L_vec.csv", L_vec);
    write_csv("../data/vsr/tensors/RESULTS/obj.csv", obj);
    write_csv("../data/vsr/tensors/RESULTS/grad.csv", grad);
    file.close();
}
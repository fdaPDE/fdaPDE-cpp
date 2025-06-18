
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

#ifndef __DTI_UTILITY_H__
#define __DTI_UTILITY_H__

namespace fdapde {
namespace internals {

// data structure to store Diffusion-Weighted Imaging (DWI) measurements
struct dwi_data {
   protected:
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
   private:
    vector_t b_;    // b-values (diffusion weightings), size: n_gradients
    matrix_t g_;    // gradient directions, size: n_gradients x n_dim
    vector_t S0_;   // baseline signal (b = 0), size: n_voxels
    matrix_t S_;    // diffusion signals, size: n_voxels x n_gradients
   public:
    // constructor
    dwi_data() = default;
    dwi_data(const vector_t& b, const matrix_t& g, const vector_t& S0, const matrix_t& S) :
        b_(b), g_(g), S0_(S0), S_(S) { }
    // getters
    const vector_t& b() const { return b_; }
    const matrix_t& g() const { return g_; }
    const vector_t& S0() const { return S0_; }
    const matrix_t& S() const { return S_; }
};

// Compute the directional derivative ∂_G exp(L) for symmetric matrices G = g g^top and L = log(D)
// using the spectral decomposition of L.
// This corresponds to evaluating the Fréchet derivative of the matrix exponential at L,
// applied to G, leveraging the divided difference formulation in the eigenbasis of L.

// compute M = R^T G R with entries scaled by divided differences of exp at eigenvalues of L
inline Eigen::Matrix<double, Dynamic, Dynamic> computeM(
  const Eigen::Matrix<double, Dynamic, Dynamic>& G, const Eigen::Matrix<double, Dynamic, Dynamic>& R,
  const Eigen::Matrix<double, Dynamic, 1>& S_diag) {
    int d = S_diag.size();
    Eigen::Matrix<double, Dynamic, Dynamic> M = R.transpose() * G * R;

    for (int l = 0; l < d; ++l) {
        for (int m = 0; m < d; ++m) {
            double s_l = S_diag(l);
            double s_m = S_diag(m);
            // apply divided difference of exp(s) to M(l,m)
            if (std::abs(s_l - s_m) > 1e-9) {
                M(l, m) *= (std::exp(s_m) - std::exp(s_l)) / (s_m - s_l);
            } else {
                M(l, m) *= std::exp(s_l);
            }
        }
    }
    return M;
}

// compute ∂_G exp(L) where G = g g^T and L = log(D), using spectral decomposition of L
inline Eigen::Matrix<double, Dynamic, 1>
dG_exp(double b, const Eigen::Matrix<double, Dynamic, 1>& g, const Eigen::Matrix<double, Dynamic, Dynamic>& L) {
    Eigen::Matrix<double, Dynamic, Dynamic> G = g * g.transpose();   // rank-one symmetric matrix G = g g^T
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, Dynamic, Dynamic>> es_L(L);
    Eigen::Matrix<double, Dynamic, Dynamic> R = es_L.eigenvectors();   // R^T L R = S
    Eigen::Matrix<double, Dynamic, 1> S_diag = es_L.eigenvalues();     // S = diag(s₁, ..., s_d)

    Eigen::Matrix<double, Dynamic, Dynamic> M = computeM(G, R, S_diag);         // M = ∂_{Rᵀ G R} exp(S)
    Eigen::Matrix<double, Dynamic, Dynamic> dG_exp_L = R * M * R.transpose();   // ∂_G exp(L) = R M R^T

    return vector_view(dG_exp_L);   // flatten result to vector form
}

}   // namespace internals

}   // namespace fdapde

#endif   // __DTI_UTILITY__

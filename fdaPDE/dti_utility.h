
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

// utilities for riccian noise distribution
double bessel_I0(double x) {
    double I0 = 0.0;
    double term = 1.0;
    double x2_over_4 = x * x / 4.0;
    for (int k = 1; k < 50; ++k) {
        term *= x2_over_4 / (k * k);   // (x^2/4)^k / (k!)^2
        I0 += term;
        if (term < 1e-12) break;
    }
    return I0;
}
double log_bessel_I0(double x) {
    // For large x, use asymptotic expansion
    if (x > 20.0) { return x - 0.5 * std::log(2 * M_PI * x); }

    // For small x, use log of series
    double sum = 1.0;
    double term = 1.0;
    double x2_over_4 = x * x / 4.0;

    for (int k = 1; k < 50; ++k) {
        term *= x2_over_4 / (k * k);
        // log(sum + term) - log(sum) ≈ log(1 + term/sum) = log1p(term/sum)
        // std::log1p(y) computes log(1 + y) in a numerically stable way even if y is very small
        if (std::log1p(term / sum) < 1e-12) break;
        sum += term;
    }

    return std::log(sum);
}
double bessel_I0_ratio(double x) {
    if (x > 10) {
        // Taylor expansion
        std::array<double, 12> coeff {
          1.0,           -0.5,     -0.125,          -0.125,        -0.1953125,        -0.40625,
          -1.0478515625, -3.21875, -11.46646118164, -46.478515625, -211.276149749755, -1064.67822265625};
        double x_pow = 1.0;
        double result = 0.0;
        for (int i = 0; i < 12; ++i) {
            result += coeff[i] * x_pow;
            x_pow /= x;
        }
        return result;
    } else {
        // Series definition
        double x2_4 = (x * x) / 4.0;
        double num = 0.0;
        double den = 0.0;
        double num_term = 1.0;
        double den_term = 1.0;
        for (int k = 0; k < 20; ++k) {
            if (k > 0) {
                num_term *= x2_4 / (k * (k + 1));
                den_term *= x2_4 / (k * k);
            }
            num += num_term;
            den += den_term;
            if (num_term < 1e-12 && den_term < 1e-12) break;
        }
        return (x / 2.0) * (num / den);
    }
}

double conditional_riccian_density(double x, double y, double sigma_sq) {
    // p(x | y; sigma^2)
    return x / sigma_sq * std::exp(-(x * x + y * y) / (2 * sigma_sq)) * bessel_I0(x * y / sigma_sq);
}
double log_conditional_riccian_density(double Shat, double S, double sigma_sq) {
    if (Shat <= 0.0 || S <= 0.0) return -1e10;   // log(0) guard

    double argument = Shat * S / sigma_sq;
    return std::log(Shat) - std::log(sigma_sq) - (Shat * Shat + S * S) / (2.0 * sigma_sq) + log_bessel_I0(argument);
}

}   // namespace internals

// data structure to store Diffusion-Weighted Imaging (DWI) measurements
struct dwi_data {
   protected:
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
   private:
    vector_t b_;            // b-values (diffusion weightings), size: n_gradients
    matrix_t g_;            // gradient directions, size: n_gradients x n_dim
    vector_t S0_;           // baseline signal (b = 0), size: n_voxels
    matrix_t S_;            // diffusion signals, size: n_voxels x n_gradients
    double sigma_ = 1e-2;   // variance of the Riccian noise
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
    double sigma() const { return sigma_; }
    int n_locs() const { return S_.rows(); }
    int n_obs() const { return S_.cols(); }
};

// non so dove metterli ma mi servono
using vector_t = Eigen::Matrix<double, Dynamic, 1>;
using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;

struct LossFunctor {
    using loss_fun_t = std::function<double(const dwi_data&, const matrix_t&)>;
    using grad_loss_fun_t = std::function<matrix_t(const dwi_data&, const matrix_t&)>;
    LossFunctor(loss_fun_t loss, grad_loss_fun_t grad_loss) : loss(std::move(loss)), grad_loss(std::move(grad_loss)) { }
    loss_fun_t loss;
    grad_loss_fun_t grad_loss;
};

// linearized gaussian loss
inline const LossFunctor linearized_gaussian_loss {
  [](const dwi_data& data, const matrix_t& L) -> double {
      double loss = 0.;
      int n_locs = data.n_locs();
      int n_obs = data.n_obs();
      vector_t S0 = data.S0();
      for (int j = 0; j < n_locs; ++j) {
          matrix_t expLj = expm(matrix_view(L.row(j)));
          for (int i = 0; i < n_obs; ++i) {
              double bi = data.b()[i];
              vector_t gi = data.g().col(i);
              double Sij = data.S()(j, i);
              double diff = std::log(S0[j] / Sij) - bi * gi.dot(expLj * gi);
              loss += diff * diff;
          }
      }
      loss /= n_locs * n_obs;
      return loss;
  },
  [](const dwi_data& data, const matrix_t& L) -> matrix_t {
      int n_locs = data.n_locs();
      int n_obs = data.n_obs();
      int n_cols = L.cols();
      matrix_t gradient = matrix_t::Zero(n_locs, n_cols);
      vector_t S0 = data.S0();
      for (int j = 0; j < n_locs; ++j) {
          matrix_t Lj = matrix_view(L.row(j));
          matrix_t expLj = expm(Lj);
          for (int i = 0; i < n_obs; ++i) {
              double bi = data.b()[i];
              vector_t gi = data.g().col(i);
              vector_t dG_exp_L = internals::dG_exp(bi, gi, Lj);
              double Sij = data.S()(j, i);
              double diff = std::log(S0[j] / Sij) - bi * gi.dot(expLj * gi);
              gradient.row(j) -= bi * diff * dG_exp_L;
          }
      }
      gradient *= 2.0 / (n_locs * n_obs);
      return gradient;
  }};

// linearized gaussian loss
inline const LossFunctor gaussian_loss {
  [](const dwi_data& data, const matrix_t& L) -> double {
      double loss = 0.;
      int n_locs = data.n_locs();
      int n_obs = data.n_obs();
      vector_t S0 = data.S0();
      for (int j = 0; j < n_locs; ++j) {
          matrix_t expLj = expm(matrix_view(L.row(j)));
          for (int i = 0; i < n_obs; ++i) {
              double bi = data.b()[i];
              vector_t gi = data.g().col(i);
              double Sij = data.S()(j, i);
              double Sij_hat = S0[j] * std::exp(-bi * gi.dot(expLj * gi));
              double diff = Sij - Sij_hat;
              loss += diff * diff;
          }
      }
      loss /= n_locs * n_obs;
      return loss;
  },
  [](const dwi_data& data, const matrix_t& L) -> matrix_t {
      int n_locs = data.n_locs();
      int n_obs = data.n_obs();
      int n_cols = L.cols();
      matrix_t gradient = matrix_t::Zero(n_locs, n_cols);
      vector_t S0 = data.S0();
      for (int j = 0; j < n_locs; ++j) {
          matrix_t Lj = matrix_view(L.row(j));
          matrix_t expLj = expm(Lj);
          for (int i = 0; i < n_obs; ++i) {
              double bi = data.b()[i];
              vector_t gi = data.g().col(i);
              vector_t dG_exp_L = internals::dG_exp(bi, gi, Lj);
              double Sij = data.S()(j, i);
              double Sij_hat = S0[j] * std::exp(-bi * gi.dot(expLj * gi));
              double diff = Sij - Sij_hat;
              gradient.row(j) += bi * diff * Sij_hat * dG_exp_L;
          }
      }
      gradient *= 2.0 / (n_locs * n_obs);
      return gradient;
  }};

// linearized riccian loss
inline const LossFunctor riccian_loss {
  [](const dwi_data& data, const matrix_t& L) -> double {
      double loss = 0.;
      int n_locs = data.n_locs();
      int n_obs = data.n_obs();
      vector_t S0 = data.S0();
      double sigma_sq = data.sigma() * data.sigma();
      for (int j = 0; j < n_locs; ++j) {
          matrix_t expLj = expm(matrix_view(L.row(j)));
          for (int i = 0; i < n_obs; ++i) {
              double bi = data.b()[i];
              vector_t gi = data.g().col(i);
              double Sij = data.S()(j, i);
              double Sij_hat = S0[j] * std::exp(-bi * gi.dot(expLj * gi));
              loss -= internals::log_conditional_riccian_density(Sij, Sij_hat, sigma_sq);
          }
      }
      loss /= n_locs * n_obs;
      return loss;
  },
  [](const dwi_data& data, const matrix_t& L) -> matrix_t {
      int n_locs = data.n_locs();
      int n_obs = data.n_obs();
      int n_cols = L.cols();
      matrix_t gradient = matrix_t::Zero(n_locs, n_cols);
      vector_t S0 = data.S0();
      double sigma_sq = data.sigma() * data.sigma();
      for (int j = 0; j < n_locs; ++j) {
          matrix_t Lj = matrix_view(L.row(j));
          matrix_t expLj = expm(Lj);
          for (int i = 0; i < n_obs; ++i) {
              double bi = data.b()[i];
              vector_t gi = data.g().col(i);
              vector_t dG_exp_L = internals::dG_exp(bi, gi, Lj);
              double Sij = data.S()(j, i);
              double Sij_hat = S0[j] * std::exp(-bi * gi.dot(expLj * gi));
              double alpha = internals::bessel_I0_ratio(Sij * Sij_hat / sigma_sq);
              double diff = Sij_hat - alpha * Sij;
              gradient.row(j) -= bi / sigma_sq * diff * Sij_hat * dG_exp_L;
          }
      }
      gradient *= 1.0 / (n_locs * n_obs);
      return gradient;
  }};

}   // namespace fdapde

#endif   // __DTI_UTILITY__

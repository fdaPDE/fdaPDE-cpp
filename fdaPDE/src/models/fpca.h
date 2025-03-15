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

#ifndef __FPCA_H__
#define __FPCA_H__

#include "header_check.h"

namespace fdapde {

[[maybe_unused]] constexpr int ComputeRandSVD = 0x1;
[[maybe_unused]] constexpr int ComputeXactSVD = 0x0;
  
template <typename VariationalSolver> class fPCA {
   private:
    using solver_t = std::decay_t<VariationalSolver>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using data_t   = Eigen::Map<const Eigen::Matrix<double, Dynamic, Dynamic, Eigen::ColMajor>>;
    static constexpr int n_lambda = solver_t::n_lambda;
   public:
    fPCA() noexcept = default;
    template <typename GeoFrame, typename Penalty>
    fPCA(const std::string& colname, const GeoFrame& gf, Penalty&& penalty) noexcept :
        solver_(), data_(gf[0].data().template col<double>(colname).as_matrix()) {
        fdapde_assert(gf.n_layers() == 1);
        n_locs_ = data_.rows();
	n_units_ = data_.cols();
        if constexpr (requires(Penalty p) { p.get(); }) {
            solver_ = solver_t(gf, penalty.get());
        } else {
            solver_ = solver_t(gf, penalty(gf.template triangulation<0>()).get());
        }
    }

    // power method
    template <typename... LambdaT>
        requires((std::is_convertible_v<LambdaT, double> && ...) && (sizeof...(LambdaT) == n_lambda))
    void fit(int rank, int flag, LambdaT... lambda) {
        matrix_t X = data_.transpose();   // create temporary of mapped data
        // first guess of PCs set to a multivariate PCA (SVD)
        matrix_t V;
        if (flag & ComputeRandSVD) {
            RSI<matrix_t> svd(X, rank);
            V = std::move(svd.matrixV());
        } else {
            Eigen::JacobiSVD<matrix_t> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
	    V = std::move(svd.matrixV());
        }
        f_.resize(solver_.n_dofs(), rank);
        s_.resize(n_units_, rank);
        f_norm_.resize(rank);
        lambda_.resize(sizeof...(lambda), rank);
	
        for (int i = 0; i < rank; ++i) {
            // initialization
            vector_t fn = V.col(i);
	    vector_t sn(n_units_);
            double Jold = std::numeric_limits<double>::max(), Jnew = 1.0;
            n_iter_ = 0;
            while (!almost_equal(Jnew, Jold, tol_) && n_iter_ < max_iter_) {
                // s = X * fn / \norm(X * fn)
                sn = X * fn;
                sn = sn / sn.norm();		
                // compute loadings as \argmin_f \sum_i (y_i - f(p_i))^2 + \int_D (\Delta f)^2, with y = X^\top * s
                solver_.update_response(X.transpose() * sn);
                solver_.fit(lambda...);
                // prepare for next iteration
                n_iter_++;
                fn = solver_.Psi() * solver_.f();
                Jold = Jnew;
                Jnew = (X - sn * fn.transpose()).squaredNorm() + solver_.ftPf(lambda...);
            }
            // store results
            f_norm_[i] = std::sqrt(solver_.f().dot(solver_.mass() * solver_.f()));
            f_.col(i) = solver_.f() / f_norm_[i];
            s_.col(i) = sn;
            internals::apply_index_pack<sizeof...(lambda)>([&]<int... Ns_>() { ((lambda_(Ns_, i) = lambda), ...); });
            // deflate
            X = X - sn * (solver_.Psi() * f_.col(i)).transpose() * f_norm_[i];
        }
        return;
    }
    // observers
    const matrix_t& scores() const { return s_; }
    const matrix_t& loading() const { return f_; }
    const std::vector<double>& loadings_norm() const { return f_norm_; }
    const matrix_t& lambda() const { return lambda_; }
    int n_iter() const { return n_iter_; }
   private:
    data_t data_;       // mapped geoframe data
    solver_t solver_;   // variational solver used in the smoothing step

    int n_locs_ = 0, n_units_ = 0;
    matrix_t f_;                   // PCs expansion coefficient vector
    matrix_t s_;                   // PCs scores
    std::vector<double> f_norm_;   // L^2 norm of estimated components
    matrix_t lambda_;              // selected level of smoothing for each component

    // power method algorithm parameters
    int n_iter_ = 0;
    double tol_ = 1e-6;
    int max_iter_ = 20;
};

// deduction guide
template <typename GeoFrame, typename Penalty>
fPCA(const std::string& colname, const GeoFrame& gf, Penalty&& solver) -> fPCA<typename Penalty::solver_t>;

}   // namespace fdapde

#endif   // __FPCA_H__

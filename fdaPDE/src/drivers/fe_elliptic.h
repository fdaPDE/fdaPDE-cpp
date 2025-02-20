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

#ifndef __FE_ELLIPTIC_DRIVER_H__
#define __FE_ELLIPTIC_DRIVER_H__

#include "header_check.h"

namespace fdapde {
namespace internals {

// solves \min_{f, \beta} \| W^{1/2} * (y_i - x_i^\top * \beta - f(p_i)) \|_2^2 + \int_D (Lf - u)^2, L elliptic operator
struct fe_elliptic_driver {
   private:
    using vector_t        = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t        = Eigen::Matrix<double, Dynamic, Dynamic>;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using diag_matrix_t   = Eigen::DiagonalMatrix<double, Dynamic, Dynamic>;
    using sparse_solver_t = Eigen::SparseLU<sparse_matrix_t>;
    using dense_solver_t  = Eigen::PartialPivLU<matrix_t>;

    template <typename Penalty, typename GeoFrame, typename WeightMatrix>
    void init_(const std::string& formula, const GeoFrame& gf, Penalty&& penalty, const WeightMatrix& W) {
        using BilinearForm = std::tuple_element_t<0, std::decay_t<Penalty>>;
        using LinearForm = std::tuple_element_t<1, std::decay_t<Penalty>>;
        using FeSpace = typename BilinearForm::TrialSpace;
	
        const BilinearForm& bilinear_form = std::get<0>(penalty);
        const LinearForm& linear_form     = std::get<1>(penalty);
        n_dofs_ = bilinear_form.n_dofs();   // number of basis functions over physical domain
        internals::fe_mass_assembly_loop<FeSpace> mass_assembler(bilinear_form.trial_space());
        R0_ = mass_assembler.assemble();
	R1_ = bilinear_form.assemble();
	u_  = linear_form.assemble();
	
        // evaluate basis system on physical domain
        switch (gf.category(0)[0]) {
        case ltype::point: {
            const auto& layer = geo_cast<POINT>(gf[0]).template geometry<0>();
            if (layer.points_at_dofs()) {
                Psi_.resize(n_dofs_, n_dofs_);
                Psi_.setIdentity();   // \psi_i(p_j) = 1 \iff i == j, otherwise \psi_i(p_j) = 0
            } else {
                Psi_ = internals::point_basis_eval(bilinear_form.trial_space(), layer.coordinates());
            }
            D_ = vector_t::Ones(n_obs_).asDiagonal();
            break;
        }
        case ltype::areal: {
            const auto& layer = geo_cast<POLYGON>(gf[0]).template geometry<0>();
            const auto& [psi, measure_vect] =
              internals::areal_basis_eval(bilinear_form.trial_space(), layer.incidence_matrix());
            Psi_ = psi;
            D_ = measure_vect.asDiagonal();   // regions' measure
            break;
        }
        }
      
        // parse formula, extract data from geoframe
        Formula formula_(formula);
        std::vector<std::string> covs;
        for (const std::string& token : formula_.rhs()) {
            if (gf.contains(token)) { covs.push_back(token); }
        }
        n_covs_ = covs.size();
        y_.resize(n_obs_);
        {
            const auto& y_data = gf[0].data().template col<double>(formula_.lhs());
            y_data.assign_to(y_);
            if (y_data.has_nan()) {   // correct \Psi for missing observations
                Psi_ = y_data.nan().as_matrix().repeat(1, n_dofs_).select(Psi_);
            }
        }
        if (n_covs_ != 0) {
            // assemble design matrix
            X_.resize(n_obs_, n_covs_);
            for (int i = 0; i < n_covs_; ++i) { gf[0].data().template col<double>(covs[i]).assign_to(X_.col(i)); }
            XtWX_ = X_.transpose() * W * X_;
            invXtWX_ = XtWX_.partialPivLu();
            invXtWXXtW_ = invXtWX_.solve(X_.transpose() * W);   // (X^\top * X)^{-1} * (X^\top * W)
            // woodbury decomposition matrices
            U_ = matrix_t::Zero(2 * n_dofs_, n_covs_);
            U_.block(0, 0, n_dofs_, n_covs_) = Psi_.transpose() * D_ * W * X_;
            V_ = matrix_t::Zero(n_covs_, 2 * n_dofs_);
            V_.block(0, 0, n_covs_, n_dofs_) = X_.transpose() * W * Psi_;
        }
	return;
    }
   public:
    fe_elliptic_driver() noexcept = default;
    template <typename Penalty, typename GeoFrame>
        requires(internals::is_pair_v<Penalty>)
    fe_elliptic_driver(const std::string& formula, const GeoFrame& gf, Penalty&& penalty) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        fdapde_assert(gf.n_layers() == 1);
        n_obs_ = gf[0].rows();   // number of data locations on physical domain
        init_(formula, gf, penalty, Eigen::Matrix<double, Dynamic, 1>::Ones(n_obs_).asDiagonal());
    }
    template <typename Penalty, typename GeoFrame, typename WeightMatrix>
        requires(internals::is_pair_v<Penalty>)
    fe_elliptic_driver(const std::string& formula, const GeoFrame& gf, Penalty&& penalty, const WeightMatrix& W) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        fdapde_assert(gf.n_layers() == 1);
	n_obs_ = gf[0].rows();   // number of data locations on physical domain
        init_(formula, gf, penalty, W);
    }

    void operator()(double lambda) {
        // assemble system matrix for nonparameteric part
        SparseBlockMatrix<double, 2, 2> A_(
          -Psi_.transpose() * D_ * Psi_, lambda * R1_.transpose(), lambda * R1_, lambda * R0_);
        invA_.compute(A_);
        // linear system rhs
        vector_t b_(2 * n_dofs_);
        b_.block(n_dofs_, 0, n_dofs_, 1) = lambda * u_;

        vector_t x;
        if (n_covs_ == 0) {   // nonparametric case
            b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * y_;
            x = invA_.solve(b_);
            f_ = x.head(n_dofs_);
        } else {   // parametric case
            b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * internals::lmbQ(X_, invXtWX_, y_);
            x = woodbury_system_solve(invA_, U_, XtWX_, V_, b_);
            f_ = x.head(n_dofs_);
            beta_ = invXtWXXtW_ * (y_ - Psi_ * f_);
        } 
        g_ = x.tail(n_dofs_);   // PDE misfit
        return;
    }
    template <typename WeightMatrix> void operator()(double lambda, WeightMatrix&& W) {
        // assemble system matrix for nonparameteric part
        SparseBlockMatrix<double, 2, 2> A_(
          -Psi_.transpose() * D_ * W * Psi_, lambda * R1_.transpose(), lambda * R1_, lambda * R0_);
        invA_.compute(A_);
        // linear system rhs
        vector_t b_(2 * n_dofs_);
        b_.block(n_dofs_, 0, n_dofs_, 1) = lambda * u_;

        vector_t x;
        if (n_covs_ == 0) {   // nonparametric case
            b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * W * y_;
            x = invA_.solve(b_);
            f_ = x.head(n_dofs_);
        } else {   // parametric case
            XtWX_ = X_.transpose() * W * X_;
            invXtWX_ = XtWX_.partialPivLu();
            b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * internals::lmbQ(W, X_, invXtWX_, y_);
            // woodbury matrices
            U_.block(0, 0, n_dofs_, n_covs_) = Psi_.transpose() * D_ * W * X_;
            V_.block(0, 0, n_covs_, n_dofs_) = X_.transpose() * W * Psi_;
            // solve A * x = (A_ + U_ * (X^\top*W*X) * V_) * x = b
            x = woodbury_system_solve(invA_, U_, XtWX_, V_, b_);
            f_ = x.head(n_dofs_);
            beta_ = invXtWX_.solve(X_.transpose() * W) * (y_ - Psi_ * f_);
        }
        g_ = x.tail(n_dofs_);   // PDE misfit
        return;
    }

    // observers
    int n_dofs() const { return n_dofs_; }
    const sparse_matrix_t& R0() const { return R0_; }
    const sparse_matrix_t& R1() const { return R1_; }
    const sparse_matrix_t& Psi() const { return Psi_; }
    const vector_t& u() const { return u_; }
    const vector_t& f() const { return f_; }
    const vector_t& beta() const { return beta_; }
    const vector_t& g() const { return g_; }
    // penalty matrix: \lambda * R1^\top * (R0)^{-1} * R1
    matrix_t P(double lambda) const {
        if (!invR0_.has_value()) { invR0_->compute(R0_); }
        return lambda * R1_.transpose() * invR0_->solve(R1_);
    }
    matrix_t P() const { return P(1.0); }
   protected:
    int n_dofs_ = 0, n_obs_ = 0, n_covs_ = 0;
    sparse_matrix_t R0_;    // n_dofs x n_dofs matrix [R0]_{ij} = \int_D \psi_i * \psi_j
    sparse_matrix_t R1_;    // n_dofs x n_dofs matrix [R1]_{ij} = \int_D a(\psi_i, \psi_j)
    sparse_matrix_t Psi_;   // n_obs x n_dofs matrix [Psi]_{ij} = \psi_j(p_i)
    vector_t u_;            // n_dofs x 1 vector u_i = \int_D u * \psi_i
    diag_matrix_t D_;       // vector of regions' measures (areal sampling)
    mutable std::optional<sparse_solver_t> invR0_;
    vector_t f_, beta_, g_;
    sparse_solver_t invA_;   // factorization of (2 * n_dofs) x (2 * n_dofs) nonparametric matrix

    matrix_t X_;                // n_obs x n_covs design matrix
    vector_t y_;                // n_obs x 1 observation vector
    matrix_t U_, V_;            // (2 * n_dofs) x n_covs matrices [\Psi^\top * D * W * y, 0] and [X^\top * W * \Psi, 0]
    matrix_t XtWX_;             // n_covs x n_covs matrix X^\top * W * X
    dense_solver_t invXtWX_;    // factorization of n_covs x n_covs matrix X^\top * W * X
    matrix_t invXtWXXtW_;       // n_covs x n_obs matrix (X^\top * X)^{-1} * (X^\top W)
};

}   // namespace internals
}   // namespace fdapde

#endif // __FE_ELLIPTIC_DRIVER_H__

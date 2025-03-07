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

#ifndef __FE_ELLIPTIC_SOLVER_H__
#define __FE_ELLIPTIC_SOLVER_H__

#include "header_check.h"

namespace fdapde {
namespace internals {

// solves \min_{f, \beta} \| W^{1/2} * (y_i - x_i^\top * \beta - f(p_i)) \|_2^2 + \int_D (Lf - u)^2, L elliptic operator
struct fe_elliptic_solver {
   private:
    using vector_t        = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t        = Eigen::Matrix<double, Dynamic, Dynamic>;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using diag_matrix_t   = Eigen::DiagonalMatrix<double, Dynamic, Dynamic>;
    using sparse_solver_t = eigen_sparse_solver_movable_wrap<Eigen::SparseLU<sparse_matrix_t>>;
    using dense_solver_t  = Eigen::PartialPivLU<matrix_t>;

    template <typename GeoFrame, typename Penalty>
    void init_(const std::string& formula, const GeoFrame& gf, Penalty&& penalty) {
        fdapde_static_assert(internals::is_valid_penalty_pair_v<Penalty>, INVALID_PENALTY_DESCRIPTION);
        using BilinearForm = std::tuple_element_t<0, std::decay_t<Penalty>>;
        using LinearForm = std::tuple_element_t<1, std::decay_t<Penalty>>;
        using FeSpace = typename BilinearForm::TrialSpace;
	// discretization
        const BilinearForm& bilinear_form = std::get<0>(penalty);
        const LinearForm& linear_form = std::get<1>(penalty);
        n_dofs_ = bilinear_form.n_dofs();   // number of basis functions over physical domain
        internals::fe_mass_assembly_loop<FeSpace> mass_assembler(bilinear_form.trial_space());
        R0_ = mass_assembler.assemble();
        R1_ = bilinear_form.assemble();
        u_  = linear_form.assemble();

        // basis system evaluation
        switch (gf.category(0)[0]) {
        case ltype::point: {
            const auto& spatial_index = geo_index_cast<0, POINT>(gf[0]);
            // evaluate basis at locations
            Psi_ = internals::point_basis_eval(bilinear_form.trial_space(), spatial_index);
            D_ = vector_t::Ones(n_obs_).asDiagonal();
            break;
        }
        case ltype::areal: {
            const auto& spatial_index = geo_index_cast<0, POLYGON>(gf[0]);
            const auto& [psi, measure_vect] =
              internals::areal_basis_eval(bilinear_form.trial_space(), spatial_index);
            Psi_ = psi;
            D_ = measure_vect.asDiagonal();   // regions' measure
            break;
        }
        }
	analyze_data_(formula, gf);
	return;
    }
    template <typename GeoFrame> void analyze_data_(const std::string& formula, const GeoFrame& gf) {
        // data extraction
        Formula formula_(formula);
        std::vector<std::string> covs;
        for (const std::string& token : formula_.rhs()) {
            if (gf.contains(token)) { covs.push_back(token); }
        }
        n_covs_ = covs.size();
        const auto& y_data = gf[0].data().template col<double>(formula_.lhs());
        fdapde_assert((y_data.blk_sz() == 1 || !y_data.has_nan()) && y_data.rows() > 0);
        y_.resize(n_obs_, y_data.blk_sz());
        y_data.assign_to(y_);
        if (y_data.has_nan()) {   // correct \Psi for missing observations
            Psi_ = y_data.nan().as_matrix().col(0).repeat(1, n_dofs_).select(Psi_);
        }
        b_.resize(2 * n_dofs_, y_.cols());
        if (n_covs_ == 0) {   // prepare linear system rhs (as it depends on y_)
            b_.block(0, 0, n_dofs_, y_.cols()) = -Psi_.transpose() * D_ * y_;
        } else {
            // assemble design matrix
            X_.resize(n_obs_, n_covs_);
            for (int i = 0; i < n_covs_; ++i) { gf[0].data().template col<double>(covs[i]).assign_to(X_.col(i)); }
            XtWX_ = X_.transpose() * W_ * X_;
            invXtWX_ = XtWX_.partialPivLu();
            invXtWXXtW_ = invXtWX_.solve(X_.transpose() * W_);   // (X^\top * X)^{-1} * (X^\top * W)
            // woodbury decomposition matrices
            U_ = matrix_t::Zero(2 * n_dofs_, n_covs_);
            U_.block(0, 0, n_dofs_, n_covs_) = Psi_.transpose() * D_ * W_ * X_;
            V_ = matrix_t::Zero(n_covs_, 2 * n_dofs_);
            V_.block(0, 0, n_covs_, n_dofs_) = X_.transpose() * W_ * Psi_;
	    // prepare linear system rhs (as it depends on y and X)
            b_.block(0, 0, n_dofs_, y_.cols()) = -Psi_.transpose() * D_ * internals::lmbQ(W_, X_, invXtWX_, y_);
        }
	return;
    }
   public:
    static constexpr int n_lambda = 1;

    fe_elliptic_solver() noexcept = default;
    template <typename GeoFrame, typename Penalty>
        requires(internals::is_pair_v<Penalty>)
    fe_elliptic_solver(const std::string& formula, const GeoFrame& gf, Penalty&& penalty) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        fdapde_assert(gf.n_layers() == 1);
        n_obs_ = gf[0].rows();   // number of data locations on physical domain
        W_.resize(n_obs_, n_obs_);
        W_.setIdentity();
        W_ /= n_obs_;   // data loss normalization
        init_(formula, gf, penalty);
    }
    template <typename GeoFrame, typename Penalty, typename WeightMatrix>
        requires(internals::is_pair_v<Penalty> && std::is_convertible_v<WeightMatrix, diag_matrix_t>)
    fe_elliptic_solver(const std::string& formula, const GeoFrame& gf, Penalty&& penalty, const WeightMatrix& W) :
        W_(W) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        fdapde_assert(gf.n_layers() == 1);
        n_obs_ = gf[0].rows();   // number of data locations on physical domain
        W_ /= n_obs_;            // data loss normalization
        init_(formula, gf, penalty);
    }
    template <typename GeoFrame>
    fe_elliptic_solver(const std::string& formula, const GeoFrame& gf, const fe_elliptic_solver& other) :
        n_dofs_(other.n_dofs_),
        R0_(other.R0_),
        R1_(other.R1_),
        Psi_(other.Psi_),
        u_(other.u_),
        D_(other.D_),
        invR0_(other.invR0_) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        fdapde_assert(gf.n_layers() == 1);
        n_obs_ = gf[0].rows();   // number of data locations on physical domain
        W_.resize(n_obs_, n_obs_);
        W_.setIdentity();
        W_ /= n_obs_;   // data loss normalization
        analyze_data_(formula, gf);
    }

    std::pair<matrix_t, matrix_t> operator()(double lambda) {
        if (!lambda_saved_.has_value() || lambda_saved_.value() != lambda) {
            // assemble and factorize system matrix for nonparameteric part
            SparseBlockMatrix<double, 2, 2> A_(
              -Psi_.transpose() * D_ * W_ * Psi_, lambda * R1_.transpose(), lambda * R1_, lambda * R0_);
            invA_.compute(A_);
            // linear system rhs
            b_.block(n_dofs_, 0, n_dofs_, y_.cols()) = lambda * u_.replicate(1, y_.cols());
            lambda_saved_ = lambda;
        }
        matrix_t x;
        if (n_covs_ == 0) {
            x = invA_.solve(b_);
            f_ = x.topRows(n_dofs_);
        } else {
            x = woodbury_system_solve(invA_, U_, XtWX_, V_, b_);
            f_ = x.topRows(n_dofs_);
            beta_ = invXtWXXtW_ * (y_ - Psi_ * f_);
        }
        g_ = x.bottomRows(n_dofs_);   // PDE misfit
        return std::make_pair(f_, beta_);
    }
    template <typename ResponseMatrix, typename WeightMatrix>
    std::pair<matrix_t, matrix_t> operator()(double lambda, ResponseMatrix&& y, WeightMatrix&& W) {
        // assemble and factorize system matrix for nonparameteric part
        SparseBlockMatrix<double, 2, 2> A_(
          -Psi_.transpose() * D_ * W * Psi_ / n_obs_, lambda * R1_.transpose(), lambda * R1_, lambda * R0_);
        invA_.compute(A_);
        // linear system rhs
        if (!lambda_saved_.has_value() || lambda_saved_.value() != lambda) {
            b_.block(n_dofs_, 0, n_dofs_, y.cols()) = lambda * u_.replicate(1, y.cols());
            lambda_saved_ = lambda;
        }
        vector_t x;
        if (n_covs_ == 0) {
            b_.block(0, 0, n_dofs_, y.cols()) = -Psi_.transpose() * D_ * W * y / n_obs_;
            x = invA_.solve(b_);
            f_ = x.topRows(n_dofs_);
        } else {
            XtWX_ = X_.transpose() * W * X_ / n_obs_;
            invXtWX_ = XtWX_.partialPivLu();
            b_.block(0, 0, n_dofs_, y.cols()) = -Psi_.transpose() * D_ * internals::lmbQ(W, X_, invXtWX_, y) / n_obs_;
            // woodbury matrices
            U_.block(0, 0, n_dofs_, n_covs_) = Psi_.transpose() * D_ * W * X_ / n_obs_;
            V_.block(0, 0, n_covs_, n_dofs_) = X_.transpose() * W * Psi_ / n_obs_;
            // solve A * x = (A_ + U_ * (X^\top * W * X) * V_) * x = b
            x = woodbury_system_solve(invA_, U_, XtWX_, V_, b_);
            f_ = x.topRows(n_dofs_);
            beta_ = invXtWX_.solve(X_.transpose() * W / n_obs_) * (y - Psi_ * f_);
        }
        g_ = x.bottomRows(n_dofs_);   // PDE misfit
        return std::make_pair(f_, beta_);
    }
    template <typename WeightMatrix> std::pair<matrix_t, matrix_t> operator()(double lambda, WeightMatrix&& W) {
        return operator()(lambda, y_, W);
    }
    // hutchinson approximation for Tr[S]
    double edf(int r = 100, int seed = random_seed) {
        fdapde_assert(lambda_saved_.has_value());
        if (!Ys_.has_value() || !Bs_.has_value()) {
            int seed_ = (seed == random_seed) ? std::random_device()() : seed;
            std::mt19937 rng(seed_);
            rademacher_distribution rademacher;
            matrix_t Us(n_obs_, r);
            for (int i = 0; i < n_obs_; ++i) {
                for (int j = 0; j < r; ++j) { Us(i, j) = rademacher(rng); }
            }
            Ys_ = Us.transpose() * Psi_;
            Bs_ = matrix_t::Zero(2 * n_dofs_, r);   // implicitly enforce homogeneous forcing
            if (n_covs_ == 0) {
                Bs_->topRows(n_dofs_) = -Psi_.transpose() * D_ * W_ * Us;
            } else {
                Bs_->topRows(n_dofs_) = -Psi_.transpose() * D_ * internals::lmbQ(W_, X_, invXtWX_, Us);
            }
        }
        matrix_t x = n_covs_ == 0 ? invA_.solve(*Bs_) : woodbury_system_solve(invA_, U_, XtWX_, V_, *Bs_);
        double trS = 0;   // monte carlo Tr[S] approximation
        for (int i = 0; i < r; ++i) { trS += Ys_->row(i).dot(x.col(i).head(n_dofs_)); }
        return trS / r;
    }
    vector_t ftPf() {
        vector_t ftPf_(y_.cols());
        for (int i = 0; i < y_.cols(); ++i) { ftPf_[i] = (*lambda_saved_) * g_.col(i).dot(R0_ * g_.col(i)); }
        return ftPf_;
    }

    // observers
    int n_dofs() const { return n_dofs_; }
    const sparse_matrix_t& mass() const { return R0_; }
    const sparse_matrix_t& stiff() const { return R1_; }
    const sparse_matrix_t& Psi() const { return Psi_; }
    const vector_t& force() const { return u_; }
    const matrix_t& f() const { return f_; }
    const matrix_t& beta() const { return beta_; }
    const matrix_t& misfit() const { return g_; }
    const matrix_t& design_matrix() const { return X_; }
    const matrix_t& response() const { return y_; }
    // penalty matrix: \lambda * R1^\top * (R0)^{-1} * R1
    matrix_t P(double lambda) const {
        if (!invR0_.has_value()) { invR0_->compute(R0_); }
        return lambda * R1_.transpose() * invR0_->solve(R1_);
    }
    matrix_t P() const { return P(1.0); }
   protected:
    std::optional<double> lambda_saved_;
    sparse_solver_t invA_;
    matrix_t b_;
    // Tr[S] hutchinson stochastic approximation matrices
    std::optional<matrix_t> Ys_;
    std::optional<matrix_t> Bs_;

    int n_dofs_ = 0, n_obs_ = 0, n_covs_ = 0;
    sparse_matrix_t R0_;    // n_dofs x n_dofs matrix [R0]_{ij} = \int_D \psi_i * \psi_j
    sparse_matrix_t R1_;    // n_dofs x n_dofs matrix [R1]_{ij} = \int_D a(\psi_i, \psi_j)
    sparse_matrix_t Psi_;   // n_obs x n_dofs matrix [Psi]_{ij} = \psi_j(p_i)
    vector_t u_;            // n_dofs x 1 vector u_i = \int_D u * \psi_i
    diag_matrix_t D_;       // vector of regions' measures (areal sampling)
    mutable std::optional<sparse_solver_t> invR0_;
    matrix_t f_, beta_, g_;

    matrix_t X_;               // n_obs x n_covs design matrix
    matrix_t y_;               // n_obs x 1 observation vector
    sparse_matrix_t W_;        // n_obs x n_obs matrix of observation weights
    matrix_t U_, V_;           // (2 * n_dofs) x n_covs matrices [\Psi^\top * D * W * y, 0] and [X^\top * W * \Psi, 0]
    matrix_t XtWX_;            // n_covs x n_covs matrix X^\top * W * X
    dense_solver_t invXtWX_;   // factorization of n_covs x n_covs matrix X^\top * W * X
    matrix_t invXtWXXtW_;      // n_covs x n_obs matrix (X^\top * X)^{-1} * (X^\top W)
};

}   // namespace internals

// general non-parametrized elliptic solver factory method
template <typename BilinearForm, typename LinearForm> struct fe_elliptic_penalty {
    fdapde_static_assert(
      std::is_same_v<typename BilinearForm::discretization_category FDAPDE_COMMA finite_element_tag>&&
        std::is_same_v<typename LinearForm::discretization_category FDAPDE_COMMA finite_element_tag>,
      FE_ELLIPTIC_PENALTY_IS_FOR_FINITE_ELEMENT_DISCRETIZATIONS_ONLY);
    using solver_t = internals::fe_elliptic_solver;

    fe_elliptic_penalty(const BilinearForm& bilinear_form, const LinearForm& linear_form) :
        penalty_(std::make_pair(bilinear_form, linear_form)) { }
    const std::tuple<BilinearForm, LinearForm>& get() const { return penalty_; }
   private:
    std::tuple<BilinearForm, LinearForm> penalty_;
};
template <typename BilinearForm, typename LinearForm>
    requires(internals::is_bilinear_form_v<BilinearForm> && internals::is_linear_form_v<LinearForm>)
fe_elliptic_penalty<BilinearForm, LinearForm>
fe_elliptic(const BilinearForm& bilinear_form, const LinearForm& linear_form) {
    return fe_elliptic_penalty(bilinear_form, linear_form);
}
template <typename BilinearForm>
    requires(internals::is_bilinear_form_v<BilinearForm>)
auto fe_elliptic(const BilinearForm& bilinear_form) {   // implicit homogeneous forcing
    using FeSpace = typename BilinearForm::TrialSpace;
    const FeSpace& Vh = bilinear_form.trial_space();
    static constexpr int embed_dim = FeSpace::embed_dim;

    ScalarField<embed_dim, decltype([](const Eigen::Matrix<double, embed_dim, 1>&) { return 0; })> u;
    TestFunction v(Vh);
    auto linear_form = integral(Vh.triangulation())(u * v);
    return fe_elliptic(bilinear_form, linear_form);
}

// catalogue of standard elliptic penalizations
template <typename Functor> struct fe_elliptic_factory {
    using solver_t = internals::fe_elliptic_solver;
    fe_elliptic_factory(const Functor& f) : f_(f) { }
    template <typename Triangulation> auto operator()(const Triangulation& D) const { return f_(D); }
  private:
    Functor f_;
};
// the laplace equation: -\Delta f = 0
auto fe_laplace() {
    return fe_elliptic_factory([]<typename Triangulation>(const Triangulation& D) {
        auto Vh = std::make_shared<FeSpace<Triangulation, FeP<1, 1>>>(D, P1<1>);
        TrialFunction f(Vh);
        TestFunction  v(Vh);
        auto a = integral(D)(dot(grad(f), grad(v)));
        return fe_elliptic(a);
    });
}
// the poisson equation: -\Delta f = u
template <typename Force> auto fe_poisson(Force&& u) {
    return fe_elliptic_factory([u]<typename Triangulation>(const Triangulation& D) {
        auto Vh = std::make_shared<FeSpace<Triangulation, FeP<1, 1>>>(D, P1<1>);
        TrialFunction f(Vh);
        TestFunction  v(Vh);
        auto a = integral(D)(dot(grad(f), grad(v)));
	auto F = integral(D)(u * v);
        return fe_elliptic(a, F);
    });
}
// the general homogeneous diffusion-transport-reaction equation: -div[K + grad(f)] + b \cdot grad(f) + c * f = 0
template <typename Diffusion, typename Transport, typename Reaction>
auto fe_diffusion_transport_reaction(Diffusion&& K, Transport&& b, Reaction&& c) {
    return fe_elliptic_factory([K, b, c]<typename Triangulation>(const Triangulation& D) {
        auto Vh = std::make_shared<FeSpace<Triangulation, FeP<1, 1>>>(D, P1<1>);
        TrialFunction f(Vh);
        TestFunction  v(Vh);
        auto a = integral(D)(dot(K * grad(f), grad(v)) + dot(b, grad(f)) * v + c * f * v);
        return fe_elliptic(a);
    });
}
// the general non-homogeneous diffusion-transport-reaction equation: -div[K + grad(f)] + b \cdot grad(f) + c * f = u
template <typename Diffusion, typename Transport, typename Reaction, typename Force>
auto fe_diffusion_transport_reaction(Diffusion&& K, Transport&& b, Reaction&& c, Force&& u) {
  return fe_elliptic_factory([K, b, c, u]<typename Triangulation>(const Triangulation& D) {
        auto Vh = std::make_shared<FeSpace<Triangulation, FeP<1, 1>>>(D, P1<1>);
        TrialFunction f(Vh);
        TestFunction  v(Vh);
        auto a = integral(D)(dot(K * grad(f), grad(v)) + dot(b, grad(f)) * v + c * f * v);
        auto F = integral(D)(u * v);
        return fe_elliptic(a, F);
    });
}

}   // namespace fdapde

#endif // __FE_ELLIPTIC_SOLVER_H__

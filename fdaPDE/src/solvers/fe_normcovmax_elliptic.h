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

#ifndef __FE_NORMCOVMAX_ELLIPTIC_SOLVER_H__
#define __FE_NORMCOVMAX_ELLIPTIC_SOLVER_H__

#include "fdaPDE/src/models/sr.h"
#include "header_check.h"

namespace fdapde {
namespace internals {

// solves \max_{f} f^\top \Psi^\top z   s.t.   f^\top \Psi^\top W \Psi f + \int_D (Lf - u)^2 = 1, L elliptic operator
struct fe_normcovmax_elliptic {
   private:
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using binarz_t = BinaryMatrix<Dynamic, Dynamic>;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using diag_matrix_t   = Eigen::DiagonalMatrix<double, Dynamic, Dynamic>;
    using sparse_solver_t = eigen_sparse_solver_movable_wrap<Eigen::SparseLU<sparse_matrix_t>>;
    using dense_solver_t  = Eigen::PartialPivLU<matrix_t>;
    using size_t = std::size_t;
    template <typename DataLocs>
    static constexpr bool is_valid_data_locs_descriptor_v =
      std::is_same_v<DataLocs, matrix_t> || std::is_same_v<DataLocs, binarz_t>;
    template <typename Penalty> struct is_valid_penalty {
        static constexpr bool value = requires(Penalty penalty) {
            penalty.bilinear_form();
            penalty.linear_form();
        };
    };
    template <typename Penalty> static constexpr bool is_valid_penalty_v = is_valid_penalty<Penalty>::value;

    // evaluation of basis system at spatial locations
    template <typename DataLocs>
    requires(is_valid_data_locs_descriptor_v<DataLocs>)
    void eval_basis_at_(const DataLocs& locs) {
        fdapde_assert(n_locs_ == locs.rows());
        if constexpr (std::is_same_v<DataLocs, matrix_t>) {   // pointwise sampling
            Psi_ = point_eval_(locs);
            D_ = vector_t::Ones(n_locs_).asDiagonal();
        } else {   // areal sampling
            const auto& [psi, measure_vect] = areal_eval_(locs);
            Psi_ = psi;
            D_ = measure_vect.asDiagonal();
        }
        return;
    }
    // optimized basis evaluation at geoframe
    template <typename GeoFrame> void eval_basis_at_(const GeoFrame& gf) {
        switch (gf.category(0)[0]) {
        case ltype::point: {
            const auto& spatial_index = geo_index_cast<0, POINT>(gf[0]);
            if (spatial_index.points_at_dofs()) {
                Psi_.resize(n_locs_, n_dofs_);
                Psi_.setIdentity();
            } else {
                Psi_ = point_eval_(spatial_index.coordinates());
            }
            D_ = vector_t::Ones(n_locs_).asDiagonal();
            break;
        }
        case ltype::areal: {
            const auto& spatial_index = geo_index_cast<0, POLYGON>(gf[0]);
            const auto& [psi, measure_vect] = areal_eval_(spatial_index.incidence_matrix());
            Psi_ = psi;
            D_ = measure_vect.asDiagonal();
            break;
        }
        }
	return;
    }
    void enforce_lhs_dirichlet_bc_(SparseBlockMatrix<double, 2, 2>& A) {
        if (dirichlet_dofs_.size() == 0) { return; }
        for (size_t i = 0; i < dirichlet_dofs_.size(); ++i) {
	      // zero out row and column in correspondance of Dirichlet-type dofs
	      A.row(dirichlet_dofs_[i]) *= 0;
	      A.col(dirichlet_dofs_[i]) *= 0;
	      A.row(n_dofs_ + dirichlet_dofs_[i]) *= 0;
	      A.col(n_dofs_ + dirichlet_dofs_[i]) *= 0;
	      // set diagonal elements to 1
	      A.coeffRef(dirichlet_dofs_[i], dirichlet_dofs_[i]) = 1;
	      A.coeffRef(n_dofs_ + dirichlet_dofs_[i], n_dofs_ + dirichlet_dofs_[i]) = 1;
        }
	    return;
    }
   public:
    static constexpr int n_lambda = 1;
    using solver_category = normcovmax_solver;

    fe_normcovmax_elliptic() noexcept = default;
    template <typename GeoFrame, typename Penalty, typename WeightMatrix>
    requires(is_valid_penalty_v<Penalty>)
    fe_normcovmax_elliptic(const GeoFrame& gf, Penalty&& penalty, const WeightMatrix& W) : W_(W) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        discretize(penalty);
        eval_basis_at_(gf);
    }
    template <typename GeoFrame, typename Penalty>
    requires(is_valid_penalty_v<Penalty>)
    fe_normcovmax_elliptic(const GeoFrame& gf, Penalty&& penalty) :
        fe_normcovmax_elliptic(gf, penalty, vector_t::Ones(gf[0].rows()).asDiagonal()) { }

    // perform finite element based numerical discretization
    template <typename Penalty> void discretize(Penalty&& penalty) {
        using BilinearForm = typename std::decay_t<Penalty>::BilinearForm;
        using LinearForm = typename std::decay_t<Penalty>::LinearForm;
        fdapde_static_assert(internals::is_valid_penalty_pair_v<BilinearForm FDAPDE_COMMA LinearForm>, INVALID_PENALTY_DESCRIPTION);
        using FeSpace = typename BilinearForm::TrialSpace;
	    // discretization
        const BilinearForm& bilinear_form = penalty.bilinear_form();
        const LinearForm& linear_form = penalty.linear_form();
        n_dofs_ = bilinear_form.n_dofs();   // number of basis functions over physical domain
        internals::fe_mass_assembly_loop<FeSpace> mass_assembler(bilinear_form.trial_space());
        R0_ = mass_assembler.assemble();
        R1_ = bilinear_form.assemble();
        u_  = linear_form.assemble();
	    // store handles for basis system evaluation at locations
        point_eval_ = [fe_space = bilinear_form.trial_space()](const matrix_t& locs) -> decltype(auto) {
            return internals::point_basis_eval(fe_space, locs);
        };
        areal_eval_ = [fe_space = bilinear_form.trial_space()](const binarz_t& locs) -> decltype(auto) {
            return internals::areal_basis_eval(fe_space, locs);
        };
	    b_.resize(2 * n_dofs_, 1);
	    // store Dirichlet boundary condition
	    auto& dof_handler = bilinear_form.trial_space().dof_handler();
	    dirichlet_dofs_ = dof_handler.dirichlet_dofs();
	    dirichlet_vals_ = dof_handler.dirichlet_values();
        // store boundary dofs
        boundary_dofs_ = bilinear_form.trial_space().triangulation().boundary_nodes().which(true);
        return;
    }
    // non-parametric fit
    // \sum_i w_i * (z_i - f(p_i))^2 + \int_D (Lf - u)^2
    template <typename DataLocs, typename WeightMatrix>
    requires(std::is_same_v<DataLocs, matrix_t> || std::is_same_v<DataLocs, binarz_t>)
    void analyze_data(const DataLocs& locs, const matrix_t& y, const WeightMatrix& W) {
        fdapde_assert(
          locs.rows() > 0 && y.rows() == locs.rows() && y.cols() == 1 && W.rows() == locs.rows() &&
          W.rows() == W.cols());
        n_obs_  = locs.rows();
	    n_locs_ = n_obs_;
        eval_basis_at_(locs);   // update \Psi matrix
        update_z_and_weights(y, W);
        return;
    }
    // evaluates basis system at physical locations
    template <typename GeoFrame, typename WeightMatrix> void analyze_data(const GeoFrame& gf, const WeightMatrix& W) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        fdapde_assert(gf.n_layers() == 1);
        n_obs_ = gf[0].rows();
        n_locs_ = n_obs_;
        eval_basis_at_(gf);   // update \Psi matrix
        W_ = W;
        return;
    }

    // modifiers
    void update_z(const vector_t& z) {
        fdapde_assert(Psi_.rows() > 0 && z.rows() == n_locs_ && z.cols() == 1);
        z_ = z;
        // correct \Psi for missing observations
        auto nan_pattern = na_matrix(z_);
	    int old_n_obs = n_obs_;
        if (nan_pattern.any()) {
            n_obs_ = n_locs_ - nan_pattern.count();
            B_ = (~nan_pattern).repeat(1, n_dofs_).select(Psi_, 0);
            z_ = (~nan_pattern).select(z_, 0);
        }
        // if (old_n_obs != n_obs_) { W_ *= (double)old_n_obs / n_obs_; }
        b_.block(0, 0, n_dofs_, 1) = -PsiNA().transpose() * z_;
	    // enforce dirichlet bc, if any
        for (size_t i = 0; i < dirichlet_dofs_.size(); ++i) {
            b_.row(dirichlet_dofs_[i]).setConstant(dirichlet_vals_[i]);
        }
        return;
    }
    template <typename WeightMatrix> void update_weights(const WeightMatrix& W) {
        fdapde_assert(Psi_.rows() > 0 && W.rows() == n_locs_ && W.rows() == W.cols());
        W_ = W;
	    // W_ /= n_obs_;
        b_.block(0, 0, n_dofs_, 1) = -PsiNA().transpose() * z_;
        // enforce dirichlet bc, if any
        for (size_t i = 0; i < dirichlet_dofs_.size(); ++i) {
            b_.row(dirichlet_dofs_[i]).setConstant(dirichlet_vals_[i]);
        }
        W_changed_ = true;
        return;
    }
    template <typename WeightMatrix> void update_z_and_weights(const vector_t& z, const WeightMatrix& W) {
        fdapde_assert(Psi_.rows() > 0 && z.rows() == n_locs_ && z.cols() == 1 && W.rows() == W.cols() && W.rows() == n_locs_);
        z_ = z;
        // correct \Psi for missing observations
        auto nan_pattern = na_matrix(z);
        if (nan_pattern.any()) {	  
            n_obs_ = n_locs_ - nan_pattern.count();
            B_ = (~nan_pattern).repeat(1, n_dofs_).select(Psi_, 0);
            z_ = (~nan_pattern).select(z_, 0);
        }
        update_weights(W);
        return;
    }

    // main fit entry point
    vector_t fit(double lambda) {
        fdapde_assert(lambda > 0 && n_dofs_ > 0 && n_obs_ > 0);
        if (lambda_saved_.value() != lambda || W_changed_) {
            // assemble and factorize system matrix
            SparseBlockMatrix<double, 2, 2> A(-PsiNA().transpose() * D_ * W_ * PsiNA(), lambda * R1_.transpose(), lambda * R1_, lambda * R0_);
	        enforce_lhs_dirichlet_bc_(A);
            invA_.compute(A);
	        W_changed_ = false;
        }
        if (lambda_saved_.value() != lambda) {
            // update linear system rhs
            b_.block(n_dofs_, 0, n_dofs_, 1) = lambda * u_;
            for (size_t i = 0; i < dirichlet_dofs_.size(); ++i) { b_.row(n_dofs_ + dirichlet_dofs_[i]).setZero(); }
        }
        lambda_saved_ = lambda;
        vector_t x;
        x = invA_.solve(b_);
        f_ = x.topRows(n_dofs_);
        g_ = x.bottomRows(n_dofs_);

        // Normalization
        double rho = fn().dot(W_ * fn()) + ftPf(lambda);
        if (rho <= 0.0) rho = 1.;
        rho = std::sqrt(rho);
        f_ /= rho;
        g_ /= rho;

        return f_;
    }
    template <typename LambdaT>
    requires(internals::is_vector_like_v<LambdaT>)
    std::pair<vector_t, vector_t> fit(LambdaT&& lambda) {
        fdapde_assert(lambda.size() == n_lambda);
        return fit(lambda[0]);
    }

    // penalty matrix: \lambda * R1^\top * (R0)^{-1} * R1
    matrix_t P(double lambda) const {
        if (!invR0_.has_value()) { invR0_.compute(R0_); }
        return lambda * R1_.transpose() * invR0_.solve(R1_);
    }
    template <typename LambdaT>
    requires(internals::is_vector_like_v<LambdaT>)
    matrix_t P(const LambdaT& lambda) const {
        fdapde_assert(lambda.size() == n_lambda);
        return P(lambda[0]);
    }
    matrix_t P() const { return P(1.0); }
    template <typename MassFactorization> matrix_t P(double lambda, const MassFactorization& invR0) const {
        return lambda * R1_.transpose() * invR0.solve(R1_);
    }
    sparse_matrix_t P_lumped(const double lambda = 1.) const {
        sparse_matrix_t R0_lumped = lump(R0_);
        vector_t d = R0_lumped.diagonal().eval();
        vector_t inv_d = d.cwiseInverse();
        sparse_matrix_t invR0_R1 = R1_;
        for (int k = 0; k < invR0_R1.outerSize(); ++k) {
            for (typename sparse_matrix_t::InnerIterator it(invR0_R1, k); it; ++it) {
                it.valueRef() *= inv_d[it.row()];
            }
        }
        return lambda * R1_.transpose() * invR0_R1;
    }
    // efficient evaluation of f^\top * P * f = g^\top * R0 * g
    double ftPf(double lambda) {
        if (lambda_saved_.value() != lambda || W_changed_) { fit(lambda); }
        return lambda * g_.dot(R0_ * g_);
    }
    template <typename LambdaT>
    requires(internals::is_vector_like_v<LambdaT>)
    double ftPf(const LambdaT& lambda) {
        fdapde_assert(lambda.size() == n_lambda);
        return ftPf(lambda[0]);
    }
    // left multiplication by \Psi
    vector_t lmbPsi(const vector_t& rhs) const { return Psi_ * rhs; }
    vector_t fn() const { return Psi_ * f_; }

    // observers
    int n_obs() const { return n_obs_; }
    int n_covs() const { return n_covs_; }
    int n_dofs() const { return n_dofs_; }
    const sparse_matrix_t& mass() const { return R0_; }
    const sparse_matrix_t& stiff() const { return R1_; }
    const sparse_matrix_t& Psi() const { return Psi_; }
    const sparse_matrix_t& PsiNA() const { return B_.has_value() ? *B_ : Psi_; }
    const vector_t& force() const { return u_; }
    const vector_t& f() const { return f_; }
    const vector_t& misfit() const { return g_; }
    const vector_t& response() const { return z_; }
    const sparse_matrix_t& weights() const { return W_; }
    double lambda() const { return *lambda_saved_; }
    const std::vector<int>& dirichlet_dofs() const  { return dirichlet_dofs_; }
    const std::vector<int>& boundary_dofs() const  { return boundary_dofs_; }

   protected:
    std::optional<double> lambda_saved_ = -1;
    sparse_solver_t invA_;
    matrix_t b_;
  
    int n_dofs_ = 0, n_locs_ = 0, n_obs_ = 0, n_covs_ = 0;
    sparse_matrix_t R0_;    // n_dofs x n_dofs matrix [R0]_{ij} = \int_D \psi_i * \psi_j
    sparse_matrix_t R1_;    // n_dofs x n_dofs matrix [R1]_{ij} = \int_D a(\psi_i, \psi_j)
    sparse_matrix_t Psi_;   // n_obs x n_dofs matrix [Psi]_{ij} = \psi_j(p_i)
    vector_t u_;            // n_dofs x 1 vector u_i = \int_D u * \psi_i
    diag_matrix_t D_;       // vector of regions' measures (areal sampling)
    mutable sparse_solver_t invR0_;
    std::optional<sparse_matrix_t> B_;   // \Psi matrix corrected for missing observations
    vector_t f_, g_;

    // basis system evaluation handles
    std::function<sparse_matrix_t(const matrix_t& locs)> point_eval_;
    std::function<std::pair<sparse_matrix_t, vector_t>(const binarz_t& locs)> areal_eval_;
    std::vector<int> dirichlet_dofs_;      // dofs where Dirichlet boundary conditions are imposed
    std::vector<double> dirichlet_vals_;   // values imposed at Dirichlet dofs
    std::vector<int> boundary_dofs_;

    vector_t z_;               // n_obs x 1 observation vector
    sparse_matrix_t W_;        // n_obs x n_obs matrix of observation weights
    bool W_changed_ {true};
};

}   // namespace internals

// elliptic solver API
template <typename BilinearForm_, typename LinearForm_> struct fe_normcovmax_elliptic {
    using solver_t = internals::fe_normcovmax_elliptic;
   private:
    struct penalty_packet {
        using BilinearForm = std::decay_t<BilinearForm_>;
        using LinearForm = std::decay_t<LinearForm_>;
       private:
        BilinearForm bilinear_form_;
        LinearForm linear_form_;
       public:
        penalty_packet(const BilinearForm_& bilinear_form, const LinearForm_& linear_form) :
            bilinear_form_(bilinear_form), linear_form_(linear_form) { }
        // observers
        const BilinearForm& bilinear_form() const { return bilinear_form_; }
        const LinearForm& linear_form() const { return linear_form_; }
    };
   public:
    fe_normcovmax_elliptic(const BilinearForm_& bilinear_form, const LinearForm_& linear_form) :
        penalty_(bilinear_form, linear_form) { }
    const penalty_packet& get() const { return penalty_; }
   private:
    penalty_packet penalty_;
};

}   // namespace fdapde

#endif // __FE_NORMCOVMAX_ELLIPTIC_SOLVER_H__

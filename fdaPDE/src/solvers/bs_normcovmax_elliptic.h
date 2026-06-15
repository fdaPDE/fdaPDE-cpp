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

#ifndef BS_NORMCOVMAX_ELLIPTIC_H
#define BS_NORMCOVMAX_ELLIPTIC_H

#include "fdaPDE/src/models/sr.h"
#include "header_check.h"

namespace fdapde {
namespace internals {

// solves \max_{f} f^\top \Psi^\top z   s.t.   f^\top \Psi^\top W \Psi f + \int_D (Lf - u)^2 = 1, L elliptic operator
struct bs_normcovmax_elliptic {
   private:
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using binary_t = BinaryMatrix<Dynamic, Dynamic>;
    using diag_matrix_t   = Eigen::DiagonalMatrix<double, Dynamic, Dynamic>;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using sparse_solver_t = eigen_sparse_solver_movable_wrap<Eigen::SparseLU<sparse_matrix_t>>;
    using dense_solver_t  = Eigen::PartialPivLU<matrix_t>;

    template <typename DataLocs>
    static constexpr bool is_valid_data_locs_descriptor_v =
      std::is_same_v<DataLocs, matrix_t> || std::is_same_v<DataLocs, binary_t>;

    template <typename Penalty> struct is_valid_penalty {
        static constexpr bool value = requires(Penalty penalty) {
            penalty.bilinear_form();
            penalty.linear_form();
        };
    };
    template <typename Penalty> static constexpr bool is_valid_penalty_v = is_valid_penalty<Penalty>::value;

    // evaluate basis at locations (point or areal)
    template <typename DataLocs>
        requires(is_valid_data_locs_descriptor_v<DataLocs>)
    void eval_basis_at_(const DataLocs& locs) {
        fdapde_assert(n_locs_ == locs.rows());
        if constexpr (std::is_same_v<DataLocs, matrix_t>) {
            Psi_ = point_eval_(locs);
            D_ = vector_t::Ones(n_locs_).asDiagonal();
        } else {
            auto pr = areal_eval_(locs);
            Psi_ = pr.first;
            D_ = pr.second.asDiagonal();
        }
    }
    // optimized basis evaluation at geoframe
    template <typename GeoFrame> void eval_basis_at_(const GeoFrame& gf) {
        switch (gf.category(0)[0]) {
            case ltype::point: {
                const auto& spatial_index = geo_index_cast<0, POINT>(gf[0]);
                Psi_ = point_eval_(spatial_index.coordinates());
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

    bool update_missing_pattern_(const vector_t& z) {
        auto new_pattern = na_matrix(z);

        bool pattern_changed = true;

        if (na_pattern_ready_ && new_pattern.size() == nan_pattern_.size()) {
            pattern_changed = new_pattern != nan_pattern_;
        }

        if (!pattern_changed) {
            return false;
        }

        nan_pattern_ = new_pattern;
        na_pattern_ready_ = true;

        if (nan_pattern_.any()) {
            n_obs_ = n_locs_ - nan_pattern_.count();
            B_ = (~nan_pattern_).repeat(1, n_dofs_).select(Psi_, 0);
        } else {
            n_obs_ = n_locs_;
            B_.reset();
        }

        return true;
    }

    void enforce_dof_(sparse_matrix_t& A, vector_t& rhs, const int dof, const double value) const {
        rhs -= A.col(dof) * value;

        A.row(dof) *= 0.0;
        A.col(dof) *= 0.0;

        A.coeffRef(dof, dof) = 1.0;
        rhs[dof] = value;
    }

    void enforce_lhs_dirichlet_bc_(sparse_matrix_t& A, vector_t& rhs) const {
        for (size_t i = 0; i < dirichlet_dofs_.size(); ++i) {
            enforce_dof_(A, rhs, dirichlet_dofs_[i], dirichlet_vals_[i]);
        }
    }

    void normalize_solution_(vector_t& f, const double lambda) {
        double rho = f.dot(Omega(lambda) * f);
        if (rho <= 0.0 || !std::isfinite(rho)) rho = 1.0;

        rho = std::sqrt(rho);
        f /= rho;
    }

   public:
    static constexpr int n_lambda = 1;
    using solver_category = normcovmax_solver;

    bs_normcovmax_elliptic() noexcept = default;
    template <typename GeoFrame, typename Penalty, typename WeightMatrix>
        requires(is_valid_penalty_v<Penalty>)
    bs_normcovmax_elliptic(const GeoFrame& gf, Penalty&& penalty, const WeightMatrix& W) : W_(W) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        // fdapde_assert(gf.n_layers() == 1);
        // n_obs_  = gf[0].rows();
        // n_locs_ = n_obs_;
        discretize(penalty);
        eval_basis_at_(gf);
    }
    template <typename GeoFrame, typename Penalty>
        requires(is_valid_penalty_v<Penalty>)
    bs_normcovmax_elliptic(const GeoFrame& gf, Penalty&& penalty) :
        bs_normcovmax_elliptic(gf, penalty, vector_t::Ones(gf[0].rows()).asDiagonal()) { }

    bs_normcovmax_elliptic(const bs_normcovmax_elliptic& other)
    : lambda_saved_(-1),
      c_(other.c_),
      n_dofs_(other.n_dofs_),
      n_locs_(other.n_locs_),
      n_obs_(other.n_obs_),
      R0_(other.R0_),
      R1_(other.R1_),
      Psi_(other.Psi_),
      u_(other.u_),
      D_(other.D_),
      B_(other.B_),
      nan_pattern_(other.nan_pattern_),
      na_pattern_ready_(other.na_pattern_ready_),
      f_(other.f_),
      point_eval_(other.point_eval_),
      areal_eval_(other.areal_eval_),
      dirichlet_vals_(other.dirichlet_vals_),
      dirichlet_dofs_(other.dirichlet_dofs_),
      boundary_dofs_(other.boundary_dofs_),
      z_(other.z_),
      W_(other.W_),
      Omega_(other.Omega_),
      Omega_changed_(false) {
        R0_.makeCompressed();
        R1_.makeCompressed();
        Psi_.makeCompressed();
        W_.makeCompressed();
        if (B_.has_value()) B_->makeCompressed();
    }

    bs_normcovmax_elliptic& operator=(const bs_normcovmax_elliptic&) = delete;


    // discretization: assemble R1 and basis evaluation handles
    template <typename Penalty> void discretize(Penalty&& penalty) {
        using BilinearForm = typename std::decay_t<Penalty>::BilinearForm;
        using LinearForm = typename std::decay_t<Penalty>::LinearForm;
        fdapde_static_assert(internals::is_valid_penalty_pair_v<BilinearForm FDAPDE_COMMA LinearForm>, INVALID_PENALTY_DESCRIPTION);

        // discretization
        using BsSpace  = typename BilinearForm::TrialSpace;
        const BilinearForm& bilinear_form = penalty.bilinear_form();
        const LinearForm& linear_form = penalty.linear_form();
        n_dofs_ = bilinear_form.n_dofs();

        // assemble penalty (stiffness) R1 and forcing u (if provided)
        auto& space = bilinear_form.trial_space();
        TrialFunction u(space);
        TestFunction v(space);
        R0_ = integral(space.triangulation())(u * v).assemble();
        R1_ = bilinear_form.assemble();
        u_  = linear_form.assemble();

        // store handles to evaluate the basis at locations
        point_eval_ = [bs_space = bilinear_form.trial_space()](const matrix_t& locs) -> decltype(auto) {
            return internals::point_basis_eval(bs_space, locs);
        };
        areal_eval_ = [bs_space = bilinear_form.trial_space()](const binary_t& locs) -> decltype(auto) {
            return internals::areal_basis_eval(bs_space, locs);
        };

        // preallocate
        c_.resize(n_dofs_, 1);
        f_.resize(n_dofs_);

        // auto& dof_handler = bilinear_form.trial_space().dof_handler();
        // dirichlet_dofs_ = dof_handler.dirichlet_dofs();
        // dirichlet_vals_ = dof_handler.dirichlet_values();
        boundary_dofs_ = bilinear_form.trial_space().triangulation().boundary_nodes().which(true);
    }

    // analyze data (raw data)
    template <typename DataLocs, typename WeightMatrix>
        requires(is_valid_data_locs_descriptor_v<DataLocs>)
    void analyze_data(const DataLocs& locs, const matrix_t& z, const WeightMatrix& W) {
        fdapde_assert(
            locs.rows() > 0 &&
            z.rows() == locs.rows() &&
            z.cols() == 1 &&
            W.rows() == locs.rows() &&
            W.rows() == W.cols()
        );

        n_obs_ = locs.rows();
        n_locs_ = n_obs_;

        eval_basis_at_(locs);
        update_z_and_weights(z, W);
    }

    // analyze data (geoframe only)
    template <typename GeoFrame, typename WeightMatrix> void analyze_data(const GeoFrame& gf, const WeightMatrix& W) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        fdapde_assert(gf.n_layers() == 1);

        n_obs_ = gf[0].rows();
        n_locs_ = n_obs_;

        eval_basis_at_(gf);   // update \Psi matrix
        W_ = W;
        Omega_changed_ = true;
    }

    // modifiers
    void update_z(const vector_t& z) {
        fdapde_assert(Psi_.rows() > 0 && z.rows() == n_locs_ && z.cols() == 1);

        const bool pattern_changed = update_missing_pattern_(z);

        if (pattern_changed) {
            Omega_changed_ = true;
        }

        if (na_pattern_ready_ && nan_pattern_.any()) {
            z_ = (~nan_pattern_).select(z, 0);
        } else {
            z_ = z;
        }

        c_ = PsiNA().transpose() * z_;
    }

    template <typename WeightMatrix>
    void update_weights(const WeightMatrix& W) {
        fdapde_assert(Psi_.rows() > 0 && W.rows() == n_locs_ && W.rows() == W.cols());

        W_ = W;
        Omega_changed_ = true;
    }

    template <typename WeightMatrix>
    void update_z_and_weights(const vector_t& z, const WeightMatrix& W) {
        fdapde_assert(
            Psi_.rows() > 0 &&
            z.rows() == n_locs_ &&
            z.cols() == 1 &&
            W.rows() == W.cols() &&
            W.rows() == n_locs_
        );

        const bool pattern_changed = update_missing_pattern_(z);

        if (na_pattern_ready_ && nan_pattern_.any()) {
            z_ = (~nan_pattern_).select(z, 0);
        } else {
            z_ = z;
        }

        update_weights(W);

        if (pattern_changed) {
            Omega_changed_ = true;
        }

        c_ = PsiNA().transpose() * z_;
    }

    const sparse_matrix_t& Omega(const double lambda) {
        if (lambda_saved_.value() != lambda || Omega_changed_) {
            Omega_ = PsiNA().transpose() * D_ * W_ * PsiNA() + lambda * P();
            Omega_.makeCompressed();
            Omega_changed_ = false;
            lambda_saved_ = lambda;
        }
        return Omega_;
    }
    const vector_t& c() const { return c_; }

    vector_t fit(double lambda) {
        fdapde_assert(lambda > 0 && n_dofs_ > 0 && n_obs_ > 0);

        auto A = Omega(lambda);
        vector_t rhs = c();

        enforce_lhs_dirichlet_bc_(A, rhs);

        invA_.compute(A);

        f_ = invA_.solve(rhs);
        normalize_solution_(f_, lambda);

        return f_;
    }
    template <typename LambdaT>
    requires(internals::is_vector_like_v<LambdaT>)
    vector_t fit(LambdaT&& lambda) {
        fdapde_assert(lambda.size() == n_lambda);
        return fit(lambda[0]);
    }

    // penalty matrix: \lambda * R1
    sparse_matrix_t P() const { return R1_; }
    sparse_matrix_t P(const double lambda) const { return lambda * P(); }
    double ftPf(const double lambda) const { return f_.dot(P(lambda) * f_); }

    // left multiplication by \Psi
    vector_t fn() const { return Psi_ * f_; }

    // observers
    int n_dofs() const { return n_dofs_; }
    int n_obs() const { return n_obs_;}

    const binary_t& nan_pattern() const { return nan_pattern_; }

    const sparse_matrix_t mass() const { return R0_; }
    const sparse_matrix_t& stiff() const { return R1_; }
    const sparse_matrix_t& Psi() const { return Psi_; }
    const sparse_matrix_t& PsiNA() const { return B_.has_value() ? *B_ : Psi_; }

    const vector_t& force() const { return u_; }
    const vector_t& f() const { return f_; }
    const vector_t& response() const { return z_; }

    const sparse_matrix_t& weights() const { return W_; }

    double lambda() const { return *lambda_saved_; }

    const std::vector<int>& dirichlet_dofs() const  { return dirichlet_dofs_; }
    const std::vector<int>& boundary_dofs() const  { return boundary_dofs_; }

protected:
    std::optional<double> lambda_saved_ = -1;

    // state
    int n_dofs_ = 0, n_locs_ = 0, n_obs_ = 0;
    sparse_matrix_t R0_;    // n_dofs x n_dofs matrix [R0]_{ij} = \int_D \psi_i * \psi_j
    sparse_matrix_t R1_;    // n_dofs x n_dofs matrix [R1]_{ij} = \int_D a(\psi_i, \psi_j)
    sparse_matrix_t Psi_;   // n_obs x n_dofs matrix [Psi]_{ij} = \psi_j(p_i)
    vector_t u_;            // n_dofs x 1 vector u_i = \int_D u * \psi_i
    vector_t z_;            // n_obs x 1 observation vector
    sparse_matrix_t W_;     // n_obs x n_obs matrix of observation weights
    diag_matrix_t D_;       // vector of regions' measures (areal sampling)

    // omega
    sparse_matrix_t Omega_;
    bool Omega_changed_ {true};

    // penalty
    sparse_matrix_t P_;

    // solvers
    sparse_solver_t invA_;
    mutable sparse_solver_t invR0_;

    // results
    vector_t c_;
    vector_t f_;

    // missingness
    binary_t nan_pattern_;
    bool na_pattern_ready_ {false};
    std::optional<sparse_matrix_t> B_;   // \Psi matrix corrected for missing observations

    // basis system evaluation handles
    std::function<sparse_matrix_t(const matrix_t& locs)> point_eval_;
    std::function<std::pair<sparse_matrix_t, vector_t>(const binary_t& locs)> areal_eval_;

    // boundary
    std::vector<int> boundary_dofs_;
    std::vector<int> dirichlet_dofs_ {};      // dofs where Dirichlet boundary conditions are imposed
    std::vector<double> dirichlet_vals_ {};   // values imposed at Dirichlet dofs
};

} // namespace internals

// solver factory
template <typename BilinearForm_, typename LinearForm_> struct bs_normcovmax_elliptic {
using solver_t = internals::bs_normcovmax_elliptic;
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
    bs_normcovmax_elliptic(const BilinearForm_& bilinear_form, const LinearForm_& linear_form) :
        penalty_(bilinear_form, linear_form) { }
    const penalty_packet& get() const { return penalty_; }
private:
    penalty_packet penalty_;
};
} // namespace fdapde

#endif // BS_NORMCOVMAX_ELLIPTIC_H

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

#ifndef __FDAPDE_RGCCA_SAMPLING_H__
#define __FDAPDE_RGCCA_SAMPLING_H__

#include "types.h"

namespace fdapde {
namespace internals {

struct empty_t {
    template<class... Args>
    explicit empty_t(Args&&...) noexcept {}
};

struct identity_ls {
    using vector_t = Eigen::VectorXd;
    using matrix_t = Eigen::MatrixXd;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;

    explicit identity_ls(const int n_dofs) : n_dofs_(n_dofs) {}

    [[nodiscard]] int n_dofs()  const { return n_dofs_; }
    [[nodiscard]] int n_obs()   const { return static_cast<int>(y_.size()); }
    [[nodiscard]] int n_covs()  const { return 0; }
    [[nodiscard]] double edf(int = 0, int = 0) const { return 0; }

    void analyze_data() {}
    void update_response_and_weights(const vector_t& y, const sparse_matrix_t&) { y_ = y; }
    void fit(double) { f_ = y_; }

    [[nodiscard]] const vector_t& response() const { return y_; }
    [[nodiscard]] vector_t fn() const { return f_; }
    [[nodiscard]] const vector_t& f() const { return f_; }
    [[nodiscard]] const sparse_matrix_t& Psi() const {
        if (Psi_.rows() == 0 && n_dofs_ > 0) {
            Psi_.resize(n_dofs_, n_dofs_);
            Psi_.setIdentity();
        }
        return Psi_;
    }
    [[nodiscard]] double ftPf(double) { return 0.; }

private:
    int n_dofs_{0};
    vector_t y_;
    vector_t f_;
    mutable sparse_matrix_t Psi_;
};

} // namespace internals

namespace rgcca {

struct IndependentSampling {
    using solver_t = ::fdapde::internals::identity_ls;
};

struct TimeDependentSampling {
    using solver_t = ::fdapde::internals::bs_ls_elliptic;
    using Matrix = Eigen::MatrixXd;
    using SparseMatrix = Eigen::SparseMatrix<double>;

    static void discretize(const Triangulation<1, 1>& T, solver_t& solver_) {
        BsSpace Bh(T, 3);
        TrialFunction f_T(Bh);
        TestFunction  v_T(Bh);
        auto a_T = integral(T)(dxx(f_T) * dxx(v_T));
        ZeroField<1> u_T;
        auto F_T = integral(T)(u_T * v_T);
        auto penalty = ::fdapde::bs_ls_elliptic(a_T, F_T);
        solver_.discretize(penalty.get());
    }

    static void compute_Psi(const Triangulation<1, 1>& T, const Matrix& times, SparseMatrix& Psi) {
        BsSpace Bh(T, 3);
        TrialFunction f_T(Bh);
        TestFunction  v_T(Bh);
        auto a_T = integral(T)(dx(f_T) * dx(v_T));
        ZeroField<1> u_T;
        auto F_T = integral(T)(u_T * v_T);
        auto penalty = ::fdapde::bs_ls_elliptic(a_T, F_T);
        using BilinearForm = typename std::decay_t<decltype(penalty.get())>::BilinearForm;
        const BilinearForm& bilinear_form = penalty.get().bilinear_form();
        Psi = ::fdapde::internals::point_basis_eval(bilinear_form.trial_space(), times);
    }
};

} // namespace rgcca
} // namespace fdapde

#endif // __FDAPDE_RGCCA_SAMPLING_H__

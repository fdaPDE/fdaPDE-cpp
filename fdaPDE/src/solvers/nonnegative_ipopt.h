#include <cassert>
#include <cmath>
#include <iostream>

#include <Eigen/Dense>

#include "IpIpoptApplication.hpp"
#include "IpTNLP.hpp"

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

class NonNegativeWeightProblem : public Ipopt::TNLP {
public:
    NonNegativeWeightProblem(const Matrix& Psi,
                            const Matrix& Omega,
                            const Vector& z,
                            const Vector& x0,
                            const std::vector<int>& boundary_dofs = {})
        : Psi_(Psi),
          Omega_(0.5 * (Omega + Omega.transpose())),
          z_(z),
          x0_(x0),
          boundary_dofs_(boundary_dofs)
    {
        c_ = Psi_.transpose() * z_;
        n_ = static_cast<Ipopt::Index>(c_.size());

        // build mask for boundary HBC
        is_boundary_.assign(n_, false);
        for (int idx : boundary_dofs_) {
            if (idx >= 0 && idx < n_) {
                is_boundary_[idx] = true;
            }
        }

        // apply dirichlet HBC to the starting point
        assert(x0_.size() == n_);
        for (Ipopt::Index i = 0; i < n_; ++i) {
            x0_[i] = std::max(0.0, x0_[i]);
            if (is_boundary_[i]) {
                x0_[i] = 0.0;
            }
        }
        const double norm2 = x0_.dot(Omega_ * x0_);
        if (norm2 > 0) x0_ /= std::sqrt(norm2);
    }

    bool get_nlp_info(
        Ipopt::Index& n,
        Ipopt::Index& m,
        Ipopt::Index& nnz_jac_g,
        Ipopt::Index& nnz_h_lag,
        IndexStyleEnum& index_style
    ) override {
        n = n_;
        m = 1;
        nnz_jac_g = n_;
        nnz_h_lag = n_ * (n_ + 1) / 2;
        index_style = Ipopt::TNLP::C_STYLE;
        return true;
    }

    bool get_bounds_info(
        Ipopt::Index n,
        Ipopt::Number* x_l,
        Ipopt::Number* x_u,
        Ipopt::Index m,
        Ipopt::Number* g_l,
        Ipopt::Number* g_u
    ) override {
        for (Ipopt::Index i = 0; i < n; ++i) {
            if (is_boundary_[i]) {
                x_l[i] = 0.0;
                x_u[i] = 0.0;
            } else {
                x_l[i] = 0.0;
                x_u[i] = 2e19;
            }
        }

        g_l[0] = 1.0;
        g_u[0] = 1.0;

        return true;
    }

    bool get_starting_point(
        Ipopt::Index n,
        bool init_x,
        Ipopt::Number* x,
        bool init_z,
        Ipopt::Number* z_L,
        Ipopt::Number* z_U,
        Ipopt::Index m,
        bool init_lambda,
        Ipopt::Number* lambda
    ) override {
        assert(init_x);

        for (Ipopt::Index i = 0; i < n; ++i) {
            x[i] = x0_[i];
        }

        return true;
    }

    bool eval_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool new_x,
        Ipopt::Number& obj_value
    ) override {
        Eigen::Map<const Vector> a(x, n);

        double s = c_.dot(a);
        obj_value = -s * s;

        return true;
    }

    bool eval_grad_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool new_x,
        Ipopt::Number* grad_f
    ) override {
        Eigen::Map<const Vector> a(x, n);
        Eigen::Map<Vector> grad(grad_f, n);

        double s = c_.dot(a);
        grad = -2.0 * s * c_;

        return true;
    }

    bool eval_g(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool new_x,
        Ipopt::Index m,
        Ipopt::Number* g
    ) override {
        Eigen::Map<const Vector> a(x, n);
        g[0] = a.dot(Omega_ * a);
        return true;
    }

    bool eval_jac_g(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool new_x,
        Ipopt::Index m,
        Ipopt::Index nele_jac,
        Ipopt::Index* iRow,
        Ipopt::Index* jCol,
        Ipopt::Number* values
    ) override {
        if (values == nullptr) {
            for (Ipopt::Index j = 0; j < n; ++j) {
                iRow[j] = 0;
                jCol[j] = j;
            }
        } else {
            Eigen::Map<const Vector> a(x, n);
            Eigen::Map<Vector> jac(values, n);

            jac = 2.0 * Omega_ * a;
        }

        return true;
    }

    bool eval_h(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool new_x,
        Ipopt::Number obj_factor,
        Ipopt::Index m,
        const Ipopt::Number* lambda,
        bool new_lambda,
        Ipopt::Index nele_hess,
        Ipopt::Index* iRow,
        Ipopt::Index* jCol,
        Ipopt::Number* values
    ) override {
        if (values == nullptr) {
            Ipopt::Index idx = 0;
            for (Ipopt::Index i = 0; i < n; ++i) {
                for (Ipopt::Index j = 0; j <= i; ++j) {
                    iRow[idx] = i;
                    jCol[idx] = j;
                    ++idx;
                }
            }
        } else {
            Matrix hess_obj = -2.0 * (c_ * c_.transpose());
            Matrix hess_con =  2.0 * Omega_;

            Matrix hess_lag =
                obj_factor * hess_obj +
                lambda[0]  * hess_con;

            Ipopt::Index idx = 0;
            for (Ipopt::Index i = 0; i < n; ++i) {
                for (Ipopt::Index j = 0; j <= i; ++j) {
                    values[idx] = hess_lag(i, j);
                    ++idx;
                }
            }
        }

        return true;
    }

    void finalize_solution(
        Ipopt::SolverReturn status,
        Ipopt::Index n,
        const Ipopt::Number* x,
        const Ipopt::Number* z_L,
        const Ipopt::Number* z_U,
        Ipopt::Index m,
        const Ipopt::Number* g,
        const Ipopt::Number* lambda,
        Ipopt::Number obj_value,
        const Ipopt::IpoptData* ip_data,
        Ipopt::IpoptCalculatedQuantities* ip_cq
    ) override {
        solution_ = Eigen::Map<const Vector>(x, n);
        status_ = status;
        obj_value_ = obj_value;
    }

    const Vector& solution() const { return solution_; }
    Ipopt::SolverReturn status() const { return status_; }
    double obj_value() const { return obj_value_; }

private:
    Matrix Psi_;
    Matrix Omega_;
    Vector z_;
    Vector c_;
    Vector x0_;
    Vector solution_;
    std::vector<bool> is_boundary_;
    const std::vector<int>& boundary_dofs_;

    Ipopt::Index n_;
    Ipopt::SolverReturn status_;
    double obj_value_ = 0.0;
};
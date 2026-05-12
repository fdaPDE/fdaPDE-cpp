#include <cassert>
#include <cmath>

#include <Eigen/Dense>

#include "IpIpoptApplication.hpp"
#include "IpTNLP.hpp"

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

class NonNegativeWeightProblem : public Ipopt::TNLP {
public:
    NonNegativeWeightProblem(
        const Matrix& Psi,
        const Matrix& Omega,
        const Vector& z,
        const std::vector<int>& boundary_dofs = {}
    ) : Psi_(Psi), Omega_(Omega), z_(z), boundary_dofs_(boundary_dofs) {

        n_ = static_cast<Ipopt::Index>(Psi_.cols());

        // cache constant quantities
        c_ = Psi_.transpose() * z_;

        // enforce symmetry once
         Omega_ = 0.5 * (Omega_ + Omega_.transpose());

        // precompute constant Hessian pieces
        hess_obj_ = -2.0 * (c_ * c_.transpose());
        hess_con_ =  2.0 * Omega_;

        // build mask for boundary dofs
        is_boundary_.assign(n_, false);
        for (int idx : boundary_dofs_) {
            if (idx < 0 || idx >= n_) {
                throw std::out_of_range("NonNegativeWeightProblem: boundary dof out of range");
            }
            is_boundary_[idx] = true;
        }

        // starting point
        x0_ = Vector::Ones(n_);
        for (Ipopt::Index i = 0; i < n_; ++i) {
            if (is_boundary_[i]) x0_[i] = 0.0;
        }
        const double norm2 = x0_.dot(Omega_ * x0_);
        if (norm2 > 0.0) x0_ /= std::sqrt(norm2);
        else throw std::runtime_error("NonNegativeWeightProblem: invalid starting point");

    }

    // returns the size of the problem
    bool get_nlp_info(
        Ipopt::Index& n,                // (out) number of variables
        Ipopt::Index& m,                // (out) number of constraints
        Ipopt::Index& nnz_jac_g,        // (out) number of nonzero entries in the Jacobian
        Ipopt::Index& nnz_h_lag,        // (out) number of nonzero entries in the Hessian
        IndexStyleEnum& index_style     // (out) numbering style used for row/col entries in the sparse matrix format
    ) override {

        n = n_;
        m = 1;

        nnz_jac_g = n_;
        nnz_h_lag = n_ * (n_ + 1) / 2;

        index_style = Ipopt::TNLP::C_STYLE; // 0-based

        return true;
    }

    // returns the variable bounds
    bool get_bounds_info(
        Ipopt::Index n,                 // (in) the number of variables x in the problem
        Ipopt::Number* x_l,             // (out) the lower bounds x_L for the variables x
        Ipopt::Number* x_u,             // (out) the upper bounds x_U for the variables x
        Ipopt::Index m,                 // (in) the number of constraints g(x) in the problem
        Ipopt::Number* g_l,             // (out) the lower bounds g_L for the constraints g(x)
        Ipopt::Number* g_u              // (out) the upper bounds g_H for the constraints g(x)
    ) override {

        for (Ipopt::Index i = 0; i < n; ++i) {
            if (is_boundary_[i]) {
                // Dirichlet Homogeneous BC
                x_l[i] = 0.0;
                x_u[i] = 0.0;
            } else {
                // Non-Negativity constraint
                x_l[i] = 0.0;
                x_u[i] = 2e19;
            }
        }

        g_l[0] = 1.0;
        g_u[0] = 1.0;

        return true;
    }

    // returns the initial point for the problem
    bool get_starting_point(
        Ipopt::Index n,                 // (in) the number of variables x in the problem
        bool init_x,                    // (in) if true, this method must provide an initial value for x (true)
        Ipopt::Number* x,               // (out) the initial values for the primal variables x
        bool init_z,                    // (in) if true, this method must provide an initial value for the bound multipliers z_L and z_u (false)
        Ipopt::Number* z_L,             // (out) the initial values for the bound multipliers z_L
        Ipopt::Number* z_U,             // (out) the initial values for the bound multipliers z_U
        Ipopt::Index m,                 // (in) the number of constraints g(x) in the problem
        bool init_lambda,               // (in) if true, this method must provide an initial value for the constraint multipliers \lambda (false)
        Ipopt::Number* lambda           // (out) the initial values for the constraint multipliers, \lambda
    ) override {

        assert(init_x == true);
        assert(init_z == false);
        assert(init_lambda == false);

        for (Ipopt::Index i = 0; i < n; ++i) {
            x[i] = x0_[i];
        }

        return true;
    }

    // returns the value of the objective function
    bool eval_f(
        Ipopt::Index n,                 // (in) the number of variables x in the problem
        const Ipopt::Number* x,         // (in) the values for the primal variables x at which the objective function f(x) is to be evaluated
        bool new_x,                     // (in) false if any evaluation method (eval_*) was previously called with the same values in x, true otherwise
        Ipopt::Number& obj_value        // (out) storage for the value of the objective function f(x)
    ) override {

        Eigen::Map<const Vector> a(x, n);
        const double s = c_.dot(a);

        obj_value = -s * s;

        return true;
    }

    // return the gradient of the objective function grad_{x} f(x)
    bool eval_grad_f(
        Ipopt::Index n,                 // (in) the number of variables x in the problem
        const Ipopt::Number* x,         // (in) the values for the primal variables x at which the objective function f(x) is to be evaluated
        bool new_x,                     // (in) false if any evaluation method (eval_*) was previously called with the same values in x, true otherwise
        Ipopt::Number* grad_f           // (out) array to store values of the gradient of the objective function grad_{x} f(x)
    ) override {

        Eigen::Map<const Vector> a(x, n);
        Eigen::Map<Vector> grad(grad_f, n);
        const double s = c_.dot(a);

        grad.noalias() = -2.0 * s * c_;

        return true;
    }

    // return the value of the constraints: g(x)
    bool eval_g(
        Ipopt::Index n,                 // (in) the number of variables x in the problem
        const Ipopt::Number* x,         // (in) the values for the primal variables x at which the constraint functions g(x) are to be evaluated
        bool new_x,                     // (in) false if any evaluation method (eval_*) was previously called with the same values in x, true otherwise
        Ipopt::Index m,                 // (in) the number of constraints g(x) in the problem
        Ipopt::Number* g                // (out) array to store constraint function values g(x), do not add or subtract the bound values g_L or g_u.
    ) override {

        Eigen::Map<const Vector> a(x, n);

        g[0] = a.dot(Omega_ * a);

        return true;
    }

    // return the structure or values of the Jacobian
    bool eval_jac_g(
        Ipopt::Index n,                 // (in) the number of variables x in the problem
        const Ipopt::Number* x,         // (in) first call: NULL; later calls: the values for the primal variables x at which the constraint Jacobian grad_{x} g(x)^\top is to be evaluated
        bool new_x,                     // (in) false if any evaluation method (eval_*) was previously called with the same values in x, true otherwise
        Ipopt::Index m,                 // (in) the number of constraints g(x) in the problem
        Ipopt::Index nele_jac,          // (in) the number of nonzero elements in the Jacobian;
        Ipopt::Index* iRow,             // (out) first call: array of length nele_jac to store the row indices of entries in the Jacobian of the constraints; later calls: NULL
        Ipopt::Index* jCol,             // (out) first call: array of length nele_jac to store the column indices of entries in the Jacobian of the constraints; later calls: NULL
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
            jac.noalias() = 2.0 * (Omega_ * a);
        }

        return true;
    }

    // return the structure or values of the Hessian
    bool eval_h(
        Ipopt::Index n,                 // (in) the number of variables x in the problem
        const Ipopt::Number* x,         // (in) first call: NULL; later calls: the values for the primal variables x at which the Hessian is to be evaluated
        bool new_x,                     // (in) false if any evaluation method (eval_*) was previously called with the same values in x, true otherwise
        Ipopt::Number obj_factor,       // (in) factor \sigma_f in front of the objective term in the Hessian
        Ipopt::Index m,                 // (in) the number of constraints g(x) in the problem
        const Ipopt::Number* lambda,    // (in) the values for the constraint multipliers \lambda at which the Hessian is to be evaluated
        bool new_lambda,                // (in) false if any evaluation method was previously called with the same values in \lambda, true otherwise
        Ipopt::Index nele_hess,         // (in) the number of nonzero elements in the Hessian
        Ipopt::Index* iRow,             // (out) first call: array of length nele_hess to store the row indices of entries in the Hessian; later calls: NULL
        Ipopt::Index* jCol,             // (out) first call: array of length nele_hess to store the column indices of entries in the Hessian; later calls: NULL
        Ipopt::Number* values           // (out) first call: NULL; later calls: array of length nele_hess to store the values of the entries in the Hessian
    ) override {

        if (values == nullptr) {
            // return the structure. this is a symmetric matrix, fill the lower left triangle only
            Ipopt::Index idx = 0;
            for (Ipopt::Index i = 0; i < n; ++i) {
                for (Ipopt::Index j = 0; j <= i; ++j) {
                    iRow[idx] = i;
                    jCol[idx] = j;
                    ++idx;
                }
            }
        } else {
            // return the values. this is a symmetric matrix, fill the lower left triangle only
            Ipopt::Index idx = 0;
            for (Ipopt::Index i = 0; i < n; ++i) {
                for (Ipopt::Index j = 0; j <= i; ++j) {
                    values[idx] = obj_factor * hess_obj_(i, j) + lambda[0] * hess_con_(i, j);
                    ++idx;
                }
            }
        }

        return true;
    }

    void finalize_solution(
        Ipopt::SolverReturn status,     // (in) gives the status of the algorithm
        Ipopt::Index n,                 // (in) the number of variables x in the problem
        const Ipopt::Number* x,         // (in) the final values for the primal variables
        const Ipopt::Number* z_L,       // (in) the final values for the lower bound multipliers
        const Ipopt::Number* z_U,       // (in) the final values for the lower bound multipliers
        Ipopt::Index m,                 // (in) the number of constraints g(x) in the problem
        const Ipopt::Number* g,         // (in) the final values of the constraint functions
        const Ipopt::Number* lambda,    // (in) the final values of the constraint multipliers
        Ipopt::Number obj_value,        // (in) the final value of the objective function
        const Ipopt::IpoptData* ip_data,
        Ipopt::IpoptCalculatedQuantities* ip_cq
    ) override {

        solution_ = Eigen::Map<const Vector>(x, n);
        status_ = status;
        obj_value_ = obj_value;

        // clean numerical negativity
        if ( status == Ipopt::SUCCESS || status == Ipopt::STOP_AT_ACCEPTABLE_POINT ) {

            // clip negatives
            solution_ = solution_.cwiseMax(0.0);

            // renormalize
            const double norm2 = solution_.dot(Omega_ * solution_);

            if (norm2 > 0.0) {
                solution_ /= std::sqrt(norm2);
            }

        }
    }


    // getters
    const Vector& solution() const { return solution_; }
    Ipopt::SolverReturn status() const { return status_; }
    double obj_value() const { return obj_value_; }

private:

    // inputs
    Matrix Psi_;
    Matrix Omega_;
    Vector z_;
    Vector x0_;
    std::vector<int> boundary_dofs_;

    // dimensions
    Ipopt::Index n_;

    // cache for constant quantities
    Vector c_;                          // \Psi^\top * z;
    Matrix hess_obj_;                   // - 2 * (c * c^\top)
    Matrix hess_con_;                   //   2 * \Omega;
    std::vector<bool> is_boundary_;     // dofs indexes of the boundary elements

    // results
    Vector solution_;
    Ipopt::SolverReturn status_ = Ipopt::SolverReturn::SUCCESS;
    double obj_value_ = 0.0;

};
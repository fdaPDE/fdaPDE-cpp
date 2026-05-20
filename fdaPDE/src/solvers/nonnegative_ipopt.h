#include <cassert>
#include <cmath>
#include <optional>

#include <Eigen/Dense>

#include "IpIpoptApplication.hpp"
#include "IpTNLP.hpp"

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;
using SparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;

class NonNegativeWeightProblem : public Ipopt::TNLP {
public:
    NonNegativeWeightProblem(
        const SparseMatrix& Psi,
        const SparseMatrix& Omega,
        const Vector& c,
        const Vector& x0,
        const std::vector<bool>& is_boundary = {}
    ) : Omega_(Omega), c_(c), x0_(x0), is_boundary_(is_boundary) {

        // dimensions
        n_ = static_cast<Ipopt::Index>(x0_.size());

        // hessian sparsity structure
        for (int k = 0; k < Omega_.outerSize(); ++k) {
            for (SparseMatrix::InnerIterator it(Omega_, k); it; ++it) {
                const Ipopt::Index i = it.row();
                const Ipopt::Index j = it.col();

                if (i >= j) {
                    hess_pos_.emplace_back(i, j);
                }
            }
        }

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
        nnz_h_lag = static_cast<Ipopt::Index>(hess_pos_.size());

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

        const Eigen::Map<const Vector> a(x, n);
        obj_value = - c_.dot(a);

        return true;
    }

    // return the gradient of the objective function grad_{x} f(x)
    bool eval_grad_f(
        Ipopt::Index n,                 // (in) the number of variables x in the problem
        const Ipopt::Number* x,         // (in) the values for the primal variables x at which the objective function f(x) is to be evaluated
        bool new_x,                     // (in) false if any evaluation method (eval_*) was previously called with the same values in x, true otherwise
        Ipopt::Number* grad_f           // (out) array to store values of the gradient of the objective function grad_{x} f(x)
    ) override {

        Eigen::Map<Vector> grad(grad_f, n);
        grad.noalias() = - c_;

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
            for (Ipopt::Index k = 0; k < static_cast<Ipopt::Index>(hess_pos_.size()); ++k) {
                iRow[k] = hess_pos_[k].first;
                jCol[k] = hess_pos_[k].second;
            }
        } else {
            // return the values. this is a symmetric matrix, fill the lower left triangle only
            for (Ipopt::Index k = 0; k < static_cast<Ipopt::Index>(hess_pos_.size()); ++k) {
                const auto [i, j] = hess_pos_[k];
                values[k] = 2.0 * lambda[0] * Omega_.coeff(i, j);
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
        if (status == Ipopt::SUCCESS || status == Ipopt::STOP_AT_ACCEPTABLE_POINT) {
            solution_ = solution_.cwiseMax(0.0);
            const double norm2 = solution_.dot(Omega_ * solution_);
            if (norm2 > 0.0) solution_ /= std::sqrt(norm2);
        }
    }


    // getters
    const Vector& solution() const { return solution_; }
    Ipopt::SolverReturn status() const { return status_; }
    double obj_value() const { return obj_value_; }

private:

    // inputs
    SparseMatrix Omega_;
    Vector c_;
    Vector x0_;

    // dimensions
    Ipopt::Index n_;

    // cache for constant quantities
    std::vector<std::pair<Ipopt::Index, Ipopt::Index>> hess_pos_;
    std::vector<bool> is_boundary_;

    // results
    Vector solution_;
    Ipopt::SolverReturn status_ = Ipopt::SolverReturn::SUCCESS;
    double obj_value_ = 0.0;
};


class NonNegativeWeightSolver {
public:
    NonNegativeWeightSolver(
        const SparseMatrix& Psi,
        const SparseMatrix& Omega,
        const std::vector<int>& boundary_dofs = {}
    ) : Psi_(Psi), Omega_(Omega) {

        app_ = IpoptApplicationFactory();

        // dimensions
        const int n = static_cast<Ipopt::Index>(Psi_.cols());

        // boundary conditions
        is_boundary_.assign(n, false);
        for (const int idx : boundary_dofs) {
            if (idx < 0 || idx >= n) {
                throw std::out_of_range("NonNegativeWeightProblem: boundary dof out of range");
            }
            is_boundary_[idx] = true;
        }

        // starting point
        x0_ = Vector::Ones(n);
        for (Ipopt::Index i = 0; i < n; ++i) {
            if (is_boundary_[i]) x0_[i] = 0.0;
        }

        const double norm2 = x0_.dot(Omega_ * x0_);
        if (norm2 > 0.0) x0_ /= std::sqrt(norm2);
        else throw std::runtime_error("NonNegativeWeightProblem: invalid starting point");

        xopt_ = x0_;

        const auto status = app_->Initialize();
        if (status != Ipopt::Solve_Succeeded) {
            throw std::runtime_error("Ipopt initialization failed.");
        }
    }

    Vector solve(const Vector& z) {

        // scaling
        double s = abs(z.dot(Psi_ * x0_));
        if (s * s <= 0.0) s = 1.0;
        const Vector c = Psi_.transpose() * z / s;

        auto* raw_pos = new NonNegativeWeightProblem(Psi_, Omega_, c, x0_, is_boundary_);
        auto* raw_neg = new NonNegativeWeightProblem(Psi_, Omega_, -c, x0_, is_boundary_);

        Ipopt::SmartPtr<Ipopt::TNLP> problem_pos = raw_pos;
        Ipopt::SmartPtr<Ipopt::TNLP> problem_neg = raw_neg;

        app_->OptimizeTNLP(problem_pos);
        app_->OptimizeTNLP(problem_neg);

        const bool pos_ok = raw_pos->status() == Ipopt::SUCCESS || raw_pos->status() == Ipopt::STOP_AT_ACCEPTABLE_POINT;
        const bool neg_ok = raw_neg->status() == Ipopt::SUCCESS || raw_neg->status() == Ipopt::STOP_AT_ACCEPTABLE_POINT;
        const bool pos_is_better = raw_pos->obj_value() <= raw_neg->obj_value(); // minimization problem

        if (pos_ok && neg_ok) xopt_ = (pos_is_better) ? raw_pos->solution() : raw_neg->solution();
        else if (pos_ok || neg_ok) {
            if (pos_ok) xopt_ =  raw_pos->solution();
            if (neg_ok) xopt_ =  raw_neg->solution();
        } else {
            std::cerr << "NonNegativeWeightSolver: optimization failed, returning the last admissible solution" << std::endl;
        }

        return xopt_;
    }

private:
    SparseMatrix Psi_;
    SparseMatrix Omega_;

    std::vector<bool> is_boundary_;
    Vector x0_;
    Vector xopt_;

    Ipopt::SmartPtr<Ipopt::IpoptApplication> app_;
};
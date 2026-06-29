#include <cassert>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "IpIpoptApplication.hpp"
#include "IpTNLP.hpp"

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;
using SparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;

class NonNegativeWeightProblem : public Ipopt::TNLP {
public:
    NonNegativeWeightProblem(
        const SparseMatrix& Omega,
        const Vector& c,
        const Vector& x0
    ) : Omega_(Omega), c_(c), x0_(x0) {

        if (Omega_.rows() != Omega_.cols())
            throw std::invalid_argument("Omega must be square");

        if (Omega_.rows() != c_.size() || Omega_.rows() != x0_.size())
            throw std::invalid_argument("incompatible dimensions");

        // dimensions
        n_ = static_cast<Ipopt::Index>(x0_.size());

        // hessian sparsity structure
        for (int k = 0; k < Omega_.outerSize(); ++k) {
            for (SparseMatrix::InnerIterator it(Omega_, k); it; ++it) {
                const Ipopt::Index i = it.row();
                const Ipopt::Index j = it.col();

                if (i >= j) {
                    hess_pos_.emplace_back(i, j);
                    hess_values_.push_back(it.value());
                }
            }
        }

    }
    void reset(const Vector& c, const Vector& x0) {
        if (c.size() != n_ || x0.size() != n_)
            throw std::invalid_argument("incompatible dimensions");

        c_ = c;
        x0_ = x0;
        solution_.resize(n_);
        status_ = Ipopt::SolverReturn::SUCCESS;
        obj_value_ = 0.0;
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

        // non-negativity constraint
        for (Ipopt::Index i = 0; i < n; ++i) {
            x_l[i] = 0.0;
            x_u[i] = 2e19;
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
                values[k] = 2.0 * lambda[0] * hess_values_[k];
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
    std::vector<double> hess_values_;

    // results
    Vector solution_;
    Ipopt::SolverReturn status_ = Ipopt::SolverReturn::SUCCESS;
    double obj_value_ = 0.0;
};


class NonNegativeWeightSolver {
public:
    NonNegativeWeightSolver(
        const SparseMatrix& Psi,
        const SparseMatrix& Omega
    ) : Psi_(Psi), Omega_(Omega) {

        Psi_.makeCompressed();
        Omega_.makeCompressed();

        app_ = IpoptApplicationFactory();

        // dimensions
        const int n = static_cast<Ipopt::Index>(Psi_.cols());
        const int m = static_cast<Ipopt::Index>(Psi_.rows());

        // starting point
        x_init_ = Vector::Ones(n);
        x_init_ = normalize_convex_comb_(x_init_, x_init_, 0);

        last_solution_ = x_init_;
        last_solution_pos_ = x_init_;
        last_solution_neg_ = x_init_;

        last_z_ = Vector::Zero(m);
        problem_pos_raw_ = new NonNegativeWeightProblem(Omega_, Vector::Zero(n), x_init_);
        problem_neg_raw_ = new NonNegativeWeightProblem(Omega_, Vector::Zero(n), x_init_);
        problem_pos_ = problem_pos_raw_;
        problem_neg_ = problem_neg_raw_;

        const auto status = app_->Initialize();
        if (status != Ipopt::Solve_Succeeded) throw std::runtime_error("Ipopt initialization failed.");
    }

    NonNegativeWeightSolver(const NonNegativeWeightSolver& other)
    : Psi_(other.Psi_),
      Omega_(other.Omega_),
      x_init_(other.x_init_),
      last_solution_(other.last_solution_),
      last_solution_pos_(other.last_solution_pos_),
      last_solution_neg_(other.last_solution_neg_),
      last_z_(other.last_z_),
      has_last_z_(other.has_last_z_)
    {
        Psi_.makeCompressed();
        Omega_.makeCompressed();

        app_ = IpoptApplicationFactory();
        problem_pos_raw_ = new NonNegativeWeightProblem(Omega_, Vector::Zero(Omega_.rows()), x_init_);
        problem_neg_raw_ = new NonNegativeWeightProblem(Omega_, Vector::Zero(Omega_.rows()), x_init_);
        problem_pos_ = problem_pos_raw_;
        problem_neg_ = problem_neg_raw_;

        const auto status = app_->Initialize();
        if (status != Ipopt::Solve_Succeeded) throw std::runtime_error("Ipopt initialization failed in copy constructor.");
    }

    NonNegativeWeightSolver& operator=(const NonNegativeWeightSolver&) = delete;

    Vector solve(const Vector& z) {
        const double z_norm2 = z.squaredNorm();
        const double scale = static_cast<double>(z.size() > 0 ? z.size() : 1);
        const double zero_tol = std::numeric_limits<double>::epsilon() * scale;
        if (z_norm2 <= zero_tol * zero_tol || !std::isfinite(z_norm2)) {
            return Vector::Zero(Psi_.cols());
        }

        double alpha = 0.0;
        if (has_last_z_) {
            const double nz  = std::sqrt(z.dot(z));
            const double nlz = std::sqrt(last_z_.dot(last_z_));
            if (nz > 0.0 && nlz > 0.0) {
                alpha = z.dot(last_z_) / (nz * nlz);
                alpha = std::clamp(alpha, 0.0, 0.95);
            }
        }

        // positive
        const Vector p = Psi_.transpose() * z;
        const bool positive_side_only = p.minCoeff() >= 0.0 && p.maxCoeff() > 0.0;
        const bool negative_side_only = p.maxCoeff() <= 0.0 && p.minCoeff() < 0.0;

        bool pos_ok = false;
        bool neg_ok = false;

        auto solve_pos = [&]() {
            Vector x0_pos = normalize_convex_comb_(x_init_, last_solution_pos_, alpha);
            double s_pos = std::abs(p.dot(x0_pos));
            if (s_pos <= 0.0 || !std::isfinite(s_pos)) s_pos = 1.0;
            problem_pos_raw_->reset(p / s_pos, x0_pos);
            app_->OptimizeTNLP(problem_pos_);
            pos_ok = problem_pos_raw_->status() == Ipopt::SUCCESS ||
                problem_pos_raw_->status() == Ipopt::STOP_AT_ACCEPTABLE_POINT;
            if (pos_ok) last_solution_pos_ = problem_pos_raw_->solution();
        };

        auto solve_neg = [&]() {
            Vector x0_neg = normalize_convex_comb_(x_init_, last_solution_neg_, alpha);
            double s_neg = std::abs(p.dot(x0_neg));
            if (s_neg <= 0.0 || !std::isfinite(s_neg)) s_neg = 1.0;
            problem_neg_raw_->reset(-p / s_neg, x0_neg);
            app_->OptimizeTNLP(problem_neg_);
            neg_ok = problem_neg_raw_->status() == Ipopt::SUCCESS ||
                problem_neg_raw_->status() == Ipopt::STOP_AT_ACCEPTABLE_POINT;
            if (neg_ok) last_solution_neg_ = problem_neg_raw_->solution();
        };

        if (positive_side_only) {
            solve_pos();
            if (!pos_ok) solve_neg();
        } else if (negative_side_only) {
            solve_neg();
            if (!neg_ok) solve_pos();
        } else {
            solve_pos();
            solve_neg();
        }

        if (pos_ok && neg_ok) {
            const double score_pos =  p.dot(last_solution_pos_);
            const double score_neg = -p.dot(last_solution_neg_);
            const bool choose_pos = score_pos >= score_neg;
            last_solution_ = choose_pos ? last_solution_pos_ : last_solution_neg_;
        } else if (pos_ok) {
            last_solution_ = last_solution_pos_;
        } else if (neg_ok) {
            last_solution_ = last_solution_neg_;
        } else {
            std::cerr << "NonNegativeWeightSolver: optimization failed, returning the last admissible solution\n";
        }

        if (pos_ok || neg_ok) {
            last_z_ = z;
            has_last_z_ = true;
        }

        return last_solution_;
    }

private:
    Vector normalize_convex_comb_(const Vector& x1, const Vector& x2, const double alpha) const {
        Vector x = (1.0-alpha) * x1 + alpha * x2;
        const double norm2 = x.dot(Omega_ * x);
        if (norm2 > 0.0) x /= std::sqrt(norm2);
        else throw std::runtime_error("NonNegativeWeightProblem: invalid starting point");
        return x;
    }


    SparseMatrix Psi_;
    SparseMatrix Omega_;

    Vector last_z_;
    bool has_last_z_ = false;

    Vector x_init_;
    Vector last_solution_pos_;
    Vector last_solution_neg_;
    Vector last_solution_;

    NonNegativeWeightProblem* problem_pos_raw_ = nullptr;
    NonNegativeWeightProblem* problem_neg_raw_ = nullptr;
    Ipopt::SmartPtr<Ipopt::TNLP> problem_pos_;
    Ipopt::SmartPtr<Ipopt::TNLP> problem_neg_;
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app_;
};

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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.

#ifndef __FDAPDE_NONNEGATIVE_WEIGHT_H__
#define __FDAPDE_NONNEGATIVE_WEIGHT_H__

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "header_check.h"

namespace fdapde {
namespace internals {

using Vector = Eigen::VectorXd;
using SparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;

// Carries the exact convex NNQP state needed to reproduce a KKT failure.
class NonNegativeWeightKKTFailure : public std::runtime_error {
public:
    NonNegativeWeightKKTFailure(
        const int sweeps,
        const int max_sweeps,
        const int accelerated_iterations,
        const int max_accelerated_iterations,
        const int accelerated_restarts,
        const double violation,
        const double tolerance,
        const double relaxation,
        const double accelerated_lipschitz,
        const Vector& c,
        const Vector& warm_start,
        const Vector& final_iterate,
        const SparseMatrix& omega
    ) : std::runtime_error(message_(
            sweeps, accelerated_iterations, violation, tolerance
        )),
        sweeps_(sweeps),
        max_sweeps_(max_sweeps),
        accelerated_iterations_(accelerated_iterations),
        max_accelerated_iterations_(max_accelerated_iterations),
        accelerated_restarts_(accelerated_restarts),
        violation_(violation),
        tolerance_(tolerance),
        relaxation_(relaxation),
        accelerated_lipschitz_(accelerated_lipschitz),
        c_(c),
        warm_start_(warm_start),
        final_iterate_(final_iterate),
        omega_(omega) {}

    void set_block_context(
        const int component,
        std::string block_name,
        const double lambda
    ) {
        component_ = component;
        block_name_ = std::move(block_name);
        lambda_ = lambda;
    }

    [[nodiscard]] int sweeps() const { return sweeps_; }
    [[nodiscard]] int max_sweeps() const { return max_sweeps_; }
    [[nodiscard]] int accelerated_iterations() const { return accelerated_iterations_; }
    [[nodiscard]] int max_accelerated_iterations() const { return max_accelerated_iterations_; }
    [[nodiscard]] int accelerated_restarts() const { return accelerated_restarts_; }
    [[nodiscard]] double violation() const { return violation_; }
    [[nodiscard]] double tolerance() const { return tolerance_; }
    [[nodiscard]] double relaxation() const { return relaxation_; }
    [[nodiscard]] double accelerated_lipschitz() const { return accelerated_lipschitz_; }
    [[nodiscard]] int component() const { return component_; }
    [[nodiscard]] const std::string& block_name() const { return block_name_; }
    [[nodiscard]] double lambda() const { return lambda_; }
    [[nodiscard]] const char* method() const { return "projected_fista"; }
    [[nodiscard]] const Vector& c() const { return c_; }
    [[nodiscard]] const Vector& warm_start() const { return warm_start_; }
    [[nodiscard]] const Vector& final_iterate() const { return final_iterate_; }
    [[nodiscard]] const SparseMatrix& omega() const { return omega_; }

private:
    static std::string message_(
        const int sweeps,
        const int accelerated_iterations,
        const double violation,
        const double tolerance
    ) {
        std::ostringstream message;
        message << "NonNegativeWeightSolver: projected FISTA did not satisfy KKT after "
                << accelerated_iterations << " iterations";
        if (sweeps > 0) message << " following " << sweeps << " PSOR sweeps";
        message << " (violation="
                << std::scientific << std::setprecision(3)
                << violation << ", tolerance=" << tolerance << ')';
        return message.str();
    }

    int sweeps_;
    int max_sweeps_;
    int accelerated_iterations_;
    int max_accelerated_iterations_;
    int accelerated_restarts_;
    double violation_;
    double tolerance_;
    double relaxation_;
    double accelerated_lipschitz_;
    int component_ = -1;
    std::string block_name_;
    double lambda_ = std::numeric_limits<double>::quiet_NaN();
    Vector c_;
    Vector warm_start_;
    Vector final_iterate_;
    SparseMatrix omega_;
};

struct NonNegativeWeightSolveStats {
    std::uint64_t calls = 0;
    std::uint64_t zero_signals = 0;
    std::uint64_t closed_form_sides = 0;
    std::uint64_t direct_sides = 0;
    std::uint64_t coordinate_attempts = 0;
    std::uint64_t coordinate_converged = 0;
    std::uint64_t coordinate_sweeps = 0;
    std::uint64_t accelerated_attempts = 0;
    std::uint64_t accelerated_converged = 0;
    std::uint64_t accelerated_iterations = 0;
    std::uint64_t accelerated_restarts = 0;
    std::uint64_t horst_boundaries = 0;
    std::uint64_t would_fallback = 0;

    NonNegativeWeightSolveStats& operator+=(const NonNegativeWeightSolveStats& other) {
        calls += other.calls;
        zero_signals += other.zero_signals;
        closed_form_sides += other.closed_form_sides;
        direct_sides += other.direct_sides;
        coordinate_attempts += other.coordinate_attempts;
        coordinate_converged += other.coordinate_converged;
        coordinate_sweeps += other.coordinate_sweeps;
        accelerated_attempts += other.accelerated_attempts;
        accelerated_converged += other.accelerated_converged;
        accelerated_iterations += other.accelerated_iterations;
        accelerated_restarts += other.accelerated_restarts;
        horst_boundaries += other.horst_boundaries;
        would_fallback += other.would_fallback;
        return *this;
    }
};

inline NonNegativeWeightSolveStats operator-(
    const NonNegativeWeightSolveStats& lhs,
    const NonNegativeWeightSolveStats& rhs
) {
    return {
        lhs.calls - rhs.calls,
        lhs.zero_signals - rhs.zero_signals,
        lhs.closed_form_sides - rhs.closed_form_sides,
        lhs.direct_sides - rhs.direct_sides,
        lhs.coordinate_attempts - rhs.coordinate_attempts,
        lhs.coordinate_converged - rhs.coordinate_converged,
        lhs.coordinate_sweeps - rhs.coordinate_sweeps,
        lhs.accelerated_attempts - rhs.accelerated_attempts,
        lhs.accelerated_converged - rhs.accelerated_converged,
        lhs.accelerated_iterations - rhs.accelerated_iterations,
        lhs.accelerated_restarts - rhs.accelerated_restarts,
        lhs.horst_boundaries - rhs.horst_boundaries,
        lhs.would_fallback - rhs.would_fallback
    };
}

// Solves max c^T w subject to w >= 0 and w^T Omega w = 1. If Omega is SPD
// and c has a positive entry, the direction is the unique nonzero solution of
// min_{x >= 0} 0.5 x^T Omega x - c^T x. Cyclic coordinate descent solves that
// convex problem; normalization recovers the constrained maximum.
class NonNegativeWeightSolver {
public:
    NonNegativeWeightSolver(
        const SparseMatrix& Psi,
        const SparseMatrix& Omega,
        const bool objective_sign_invariant = true,
        const bool use_closed_form_solution = false
    ) : Psi_(Psi),
        Omega_(Omega),
        objective_sign_invariant_(objective_sign_invariant),
        use_closed_form_solution_(use_closed_form_solution) {

        Psi_.makeCompressed();
        Omega_.makeCompressed();
        validate_omega_();

        if (use_closed_form_solution_) {
            SparseMatrix identity(Omega_.rows(), Omega_.cols());
            identity.setIdentity();
            const double identity_tol = 100.0 * std::numeric_limits<double>::epsilon() *
                std::max(1.0, Omega_.norm());
            if ((Omega_ - identity).norm() > identity_tol)
                throw std::invalid_argument(
                    "NonNegativeWeightSolver: closed-form solve requires identity Omega"
                );
        } else {
            omega_solver_.compute(Omega_);
            if (omega_solver_.info() != Eigen::Success)
                throw std::invalid_argument("NonNegativeWeightSolver: Omega factorization failed");

            const Vector pivots = omega_solver_.vectorD();
            const double max_pivot = pivots.cwiseAbs().maxCoeff();
            const double pivot_tol = std::numeric_limits<double>::epsilon() *
                static_cast<double>(std::max<Eigen::Index>(1, pivots.size())) * max_pivot;
            if (!pivots.allFinite() || !(max_pivot > 0.0) || pivots.minCoeff() <= pivot_tol)
                throw std::invalid_argument("NonNegativeWeightSolver: Omega must be numerically positive definite");
        }
        initialize_workspace_();

        x_init_ = Vector::Ones(Psi_.cols());
        x_init_ = normalize_(x_init_);
        last_solution_ = x_init_;
        last_solution_pos_ = x_init_;
        last_solution_neg_ = x_init_;
    }

    NonNegativeWeightSolver(const NonNegativeWeightSolver& other)
    : Psi_(other.Psi_),
      Omega_(other.Omega_),
      objective_sign_invariant_(other.objective_sign_invariant_),
      use_closed_form_solution_(other.use_closed_form_solution_),
      x_init_(other.x_init_),
      last_solution_(other.last_solution_),
      last_solution_pos_(other.last_solution_pos_),
      last_solution_neg_(other.last_solution_neg_) {
        Psi_.makeCompressed();
        Omega_.makeCompressed();
        if (!use_closed_form_solution_) {
            omega_solver_.compute(Omega_);
            if (omega_solver_.info() != Eigen::Success)
                throw std::invalid_argument("NonNegativeWeightSolver: Omega factorization failed in copy constructor");
        }
        initialize_workspace_();
    }

    NonNegativeWeightSolver& operator=(const NonNegativeWeightSolver&) = delete;

    static NonNegativeWeightSolveStats thread_stats() { return thread_stats_; }

    void reset_warm_start(const Vector& solution) {
        if (solution.size() != x_init_.size())
            throw std::invalid_argument("NonNegativeWeightSolver: incompatible warm start");
        last_solution_ = solution;
        last_solution_pos_ = solution;
        last_solution_neg_ = solution;
    }

    Vector solve(const Vector& z) {
        ++thread_stats_.calls;
        if (z.size() != Psi_.rows())
            throw std::invalid_argument("NonNegativeWeightSolver: incompatible signal dimension");
        if (!z.allFinite())
            throw std::invalid_argument("NonNegativeWeightSolver: signal must be finite");

        const double scale = static_cast<double>(std::max<Eigen::Index>(1, z.size()));
        const double zero_tol = std::numeric_limits<double>::epsilon() * scale;
        if (z.squaredNorm() <= zero_tol * zero_tol) {
            ++thread_stats_.zero_signals;
            return Vector::Zero(Psi_.cols());
        }

        transformed_signal_.noalias() = Psi_.transpose() * z;
        const Vector& p = transformed_signal_;
        if (!p.allFinite())
            throw std::runtime_error("NonNegativeWeightSolver: transformed signal is not finite");
        if (p.lpNorm<Eigen::Infinity>() <= zero_tol) {
            ++thread_stats_.zero_signals;
            return Vector::Zero(Psi_.cols());
        }

        // Horst is orientation-sensitive. If every coefficient is nonpositive,
        // the convex NNQP has the zero solution; the sphere optimum is instead
        // the least-negative normalized coordinate direction.
        if (!objective_sign_invariant_ && p.maxCoeff() <= 0.0) {
            ++thread_stats_.horst_boundaries;
            last_solution_ = nonpositive_horst_solution_(p);
            return last_solution_;
        }

        bool pos_ok = false;
        bool neg_ok = false;
        auto solve_side = [&](const Vector& c, Vector& solution) {
            if (!(c.maxCoeff() > 0.0))
                throw std::logic_error("NonNegativeWeightSolver: requested an empty NNQP orientation");

            if (use_closed_form_solution_) {
                solution = closed_form_solution_(c);
                ++thread_stats_.closed_form_sides;
                return true;
            }
            bool projected_direct_valid = false;
            if (solution.minCoeff() > 1e-14 &&
                try_direct_solution_(c, solution, &projected_direct_, &projected_direct_valid)) {
                ++thread_stats_.direct_sides;
                return true;
            }

            ++thread_stats_.coordinate_attempts;
            int sweeps = 0;
            int accelerated_iterations = 0;
            int accelerated_restarts = 0;
            double violation = std::numeric_limits<double>::infinity();
            double tolerance = 0.0;
            const Vector& warm_start = projected_direct_valid &&
                c.dot(projected_direct_) > c.dot(solution) ? projected_direct_ : solution;
            const bool converged = try_coordinate_solution_(
                c, warm_start, &solution, sweeps, accelerated_iterations,
                accelerated_restarts, violation, tolerance
            );
            thread_stats_.coordinate_sweeps += static_cast<std::uint64_t>(sweeps);
            if (accelerated_iterations > 0) {
                ++thread_stats_.accelerated_attempts;
                thread_stats_.accelerated_iterations +=
                    static_cast<std::uint64_t>(accelerated_iterations);
                thread_stats_.accelerated_restarts +=
                    static_cast<std::uint64_t>(accelerated_restarts);
            }
            if (converged) {
                ++thread_stats_.coordinate_converged;
                if (accelerated_iterations > 0)
                    ++thread_stats_.accelerated_converged;
                return true;
            }

            ++thread_stats_.would_fallback;
            throw NonNegativeWeightKKTFailure(
                sweeps,
                coordinate_max_sweeps_,
                accelerated_iterations,
                accelerated_max_iterations_,
                accelerated_restarts,
                violation,
                tolerance,
                coordinate_relaxation_,
                accelerated_lipschitz_,
                c,
                warm_start,
                coordinate_x_,
                Omega_
            );
        };

        const bool positive_side_only = p.minCoeff() >= 0.0;
        const bool negative_side_only = p.maxCoeff() <= 0.0;
        if (!objective_sign_invariant_ || positive_side_only) {
            pos_ok = solve_side(p, last_solution_pos_);
        } else if (negative_side_only) {
            neg_ok = solve_side(-p, last_solution_neg_);
        } else {
            auto pos_dominates = [&](const double neg_bound) {
                return pos_ok && p.dot(last_solution_pos_) +
                    1e-10 * (1.0 + std::abs(neg_bound)) >= neg_bound;
            };
            auto neg_dominates = [&](const double pos_bound) {
                return neg_ok && -p.dot(last_solution_neg_) +
                    1e-10 * (1.0 + std::abs(pos_bound)) >= pos_bound;
            };

            const double pos_warm_objective = p.dot(last_solution_pos_);
            const double neg_warm_objective = -p.dot(last_solution_neg_);
            if (pos_warm_objective >= neg_warm_objective) {
                pos_ok = solve_side(p, last_solution_pos_);
                if (!pos_dominates(upper_bound_(-p)))
                    neg_ok = solve_side(-p, last_solution_neg_);
            } else {
                neg_ok = solve_side(-p, last_solution_neg_);
                if (!neg_dominates(upper_bound_(p)))
                    pos_ok = solve_side(p, last_solution_pos_);
            }
        }

        if (pos_ok && neg_ok) {
            last_solution_ = p.dot(last_solution_pos_) >= -p.dot(last_solution_neg_) ?
                last_solution_pos_ : last_solution_neg_;
        } else if (pos_ok) {
            last_solution_ = last_solution_pos_;
        } else if (neg_ok) {
            last_solution_ = last_solution_neg_;
        } else {
            throw std::logic_error("NonNegativeWeightSolver: no orientation was solved");
        }
        return last_solution_;
    }

private:
    void initialize_workspace_() {
        diagonal_ = Omega_.diagonal();
        sqrt_diagonal_ = diagonal_.array().sqrt();
        inv_sqrt_diagonal_ = sqrt_diagonal_.array().inverse();

        const Eigen::Index n = Psi_.cols();
        transformed_signal_.resize(n);
        direct_solution_.resize(n);
        projected_direct_.resize(n);
        coordinate_x_.resize(n);
        coordinate_u_.resize(n);
        coordinate_gradient_.resize(n);
        accelerated_y_.resize(n);
        accelerated_next_.resize(n);
        scaled_c_.resize(n);
        bound_c_.resize(n);
        bound_solution_.resize(n);

        Vector scaled_absolute_row_sum = Vector::Zero(n);
        for (int j = 0; j < Omega_.outerSize(); ++j) {
            for (SparseMatrix::InnerIterator it(Omega_, j); it; ++it) {
                scaled_absolute_row_sum[it.row()] += std::abs(it.value()) *
                    inv_sqrt_diagonal_[it.row()] * inv_sqrt_diagonal_[j];
            }
        }
        accelerated_lipschitz_ = scaled_absolute_row_sum.maxCoeff() *
            (1.0 + std::sqrt(std::numeric_limits<double>::epsilon()));
        if (!(accelerated_lipschitz_ > 0.0) || !std::isfinite(accelerated_lipschitz_))
            throw std::invalid_argument("NonNegativeWeightSolver: invalid accelerated Lipschitz bound");
    }

    void validate_omega_() const {
        if (Omega_.rows() == 0 || Omega_.rows() != Omega_.cols())
            throw std::invalid_argument("NonNegativeWeightSolver: Omega must be nonempty and square");
        if (Psi_.cols() != Omega_.rows())
            throw std::invalid_argument("NonNegativeWeightSolver: incompatible Psi and Omega dimensions");
        for (int k = 0; k < Omega_.outerSize(); ++k)
            for (SparseMatrix::InnerIterator it(Omega_, k); it; ++it)
                if (!std::isfinite(it.value()))
                    throw std::invalid_argument("NonNegativeWeightSolver: Omega must be finite");

        const SparseMatrix asymmetry = Omega_ - SparseMatrix(Omega_.transpose());
        const double symmetry_tol = 100.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, Omega_.norm());
        if (asymmetry.norm() > symmetry_tol)
            throw std::invalid_argument("NonNegativeWeightSolver: Omega must be symmetric");
    }

    Vector normalize_(Vector x) const {
        const double norm2 = x.dot(Omega_ * x);
        if (!(norm2 > 0.0) || !std::isfinite(norm2))
            throw std::invalid_argument("NonNegativeWeightSolver: invalid Omega norm");
        return x / std::sqrt(norm2);
    }

    Vector closed_form_solution_(const Vector& c) const {
        Vector x = c.cwiseMax(0.0);
        const double norm = x.norm();
        if (!(norm > 0.0) || !std::isfinite(norm))
            throw std::logic_error("NonNegativeWeightSolver: empty closed-form orientation");
        return x / norm;
    }

    Vector nonpositive_horst_solution_(const Vector& c) const {
        int best_i = 0;
        double best_score = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < c.size(); ++i) {
            const double score = c[i] * inv_sqrt_diagonal_[i];
            if (score > best_score) {
                best_score = score;
                best_i = i;
            }
        }
        Vector out = Vector::Zero(c.size());
        out[best_i] = inv_sqrt_diagonal_[best_i];
        return out;
    }

    bool try_direct_solution_(
        const Vector& c,
        Vector& out,
        Vector* projected_out = nullptr,
        bool* projected_valid = nullptr
    ) {
        if (projected_valid != nullptr) *projected_valid = false;
        direct_solution_ = omega_solver_.solve(c);
        if (omega_solver_.info() != Eigen::Success || !direct_solution_.allFinite())
            return false;

        const double scale = std::max(1.0, direct_solution_.cwiseAbs().maxCoeff());
        const double tol = 100.0 * std::numeric_limits<double>::epsilon() * scale;
        if (direct_solution_.minCoeff() < -tol) {
            if (projected_out != nullptr) {
                *projected_out = direct_solution_.cwiseMax(0.0);
                const double norm2 = projected_out->dot(Omega_ * *projected_out);
                if (norm2 > 0.0 && std::isfinite(norm2)) {
                    *projected_out /= std::sqrt(norm2);
                    if (projected_valid != nullptr) *projected_valid = true;
                }
            }
            return false;
        }

        direct_solution_ = direct_solution_.cwiseMax(0.0);
        const double norm2 = direct_solution_.dot(Omega_ * direct_solution_);
        if (!(norm2 > 0.0) || !std::isfinite(norm2))
            return false;
        out = direct_solution_ / std::sqrt(norm2);
        return true;
    }

    bool try_coordinate_solution_(
        const Vector& c,
        const Vector& warm_start,
        Vector* out,
        int& sweeps,
        int& accelerated_iterations,
        int& accelerated_restarts,
        double& violation,
        double& tolerance
    ) {
        const int n = static_cast<int>(c.size());
        coordinate_x_ = warm_start.cwiseMax(0.0);
        const double warm_norm2 = coordinate_x_.dot(Omega_ * coordinate_x_);
        const double warm_scale = warm_norm2 > 0.0 ?
            std::max(0.0, c.dot(coordinate_x_) / warm_norm2) : 0.0;
        coordinate_u_ = sqrt_diagonal_.cwiseProduct(coordinate_x_) * warm_scale;
        coordinate_x_ = inv_sqrt_diagonal_.cwiseProduct(coordinate_u_);

        if constexpr (coordinate_max_sweeps_ > 0) {
            coordinate_gradient_.noalias() = Omega_ * coordinate_x_;
            coordinate_gradient_ -= c;
        }
        scaled_c_ = inv_sqrt_diagonal_.cwiseProduct(c);
        tolerance = 1e-8 * std::max(1.0, scaled_c_.lpNorm<Eigen::Infinity>());
        constexpr int coordinate_violation_check_every = 5;
        constexpr int accelerated_violation_check_every = 20;
        const auto accept_solution = [&]() {
            coordinate_x_ = inv_sqrt_diagonal_.cwiseProduct(coordinate_u_);
            const double norm2 = coordinate_x_.dot(Omega_ * coordinate_x_);
            if (!(norm2 > 0.0) || !std::isfinite(norm2)) return false;
            *out = coordinate_x_ / std::sqrt(norm2);
            return true;
        };
        const auto accept_solution_from_exact_gradient = [&]() {
            const double norm2 = coordinate_x_.dot(coordinate_gradient_) +
                coordinate_x_.dot(c);
            if (!(norm2 > 0.0) || !std::isfinite(norm2)) return false;
            *out = coordinate_x_ / std::sqrt(norm2);
            return true;
        };

        for (sweeps = 1; sweeps <= coordinate_max_sweeps_; ++sweeps) {
            for (int i = 0; i < n; ++i) {
                const double old_value = coordinate_u_[i];
                const double gradient_i = coordinate_gradient_[i] * inv_sqrt_diagonal_[i];
                const double new_value = std::max(
                    0.0, old_value - coordinate_relaxation_ * gradient_i
                );
                const double delta = new_value - old_value;
                if (delta == 0.0) continue;

                coordinate_u_[i] = new_value;
                const double delta_x = delta * inv_sqrt_diagonal_[i];
                for (SparseMatrix::InnerIterator it(Omega_, i); it; ++it)
                    coordinate_gradient_[it.row()] += delta_x * it.value();
            }

            if (sweeps % 20 == 0) {
                coordinate_x_ = inv_sqrt_diagonal_.cwiseProduct(coordinate_u_);
                coordinate_gradient_.noalias() = Omega_ * coordinate_x_;
                coordinate_gradient_ -= c;
            }
            if (sweeps % coordinate_violation_check_every != 0) continue;

            violation = 0.0;
            for (int i = 0; i < n; ++i) {
                const double gradient_i = coordinate_gradient_[i] * inv_sqrt_diagonal_[i];
                const double current = coordinate_u_[i] > 1e-14 ?
                    std::abs(gradient_i) : std::max(0.0, -gradient_i);
                violation = std::max(violation, current);
            }
            if (violation > tolerance) continue;
            return accept_solution();
        }
        --sweeps;

        // Projected FISTA acts on the diagonally scaled NNQP, with a Gershgorin
        // upper bound on its Lipschitz constant and deterministic gradient restart.
        accelerated_y_ = coordinate_u_;
        double momentum = 1.0;
        const auto accelerated_kkt_satisfied = [&]() {
            coordinate_x_ = inv_sqrt_diagonal_.cwiseProduct(coordinate_u_);
            coordinate_gradient_.noalias() = Omega_ * coordinate_x_;
            coordinate_gradient_ -= c;
            violation = 0.0;
            for (int i = 0; i < n; ++i) {
                const double gradient_i = coordinate_gradient_[i] * inv_sqrt_diagonal_[i];
                const double current = coordinate_u_[i] > 1e-14 ?
                    std::abs(gradient_i) : std::max(0.0, -gradient_i);
                violation = std::max(violation, current);
            }
            return violation <= tolerance;
        };
        for (accelerated_iterations = 1;
             accelerated_iterations <= accelerated_max_iterations_;
             ++accelerated_iterations) {
            coordinate_x_ = inv_sqrt_diagonal_.cwiseProduct(accelerated_y_);
            coordinate_gradient_.noalias() = Omega_ * coordinate_x_;
            coordinate_gradient_ -= c;
            accelerated_next_ = accelerated_y_ -
                inv_sqrt_diagonal_.cwiseProduct(coordinate_gradient_) /
                    accelerated_lipschitz_;
            accelerated_next_ = accelerated_next_.cwiseMax(0.0);
            if (!accelerated_next_.allFinite()) {
                violation = std::numeric_limits<double>::infinity();
                coordinate_x_ = inv_sqrt_diagonal_.cwiseProduct(coordinate_u_);
                return false;
            }

            const double next_momentum =
                0.5 * (1.0 + std::sqrt(1.0 + 4.0 * momentum * momentum));
            const bool restart = (accelerated_y_ - accelerated_next_)
                .dot(accelerated_next_ - coordinate_u_) > 0.0;
            if (restart) {
                ++accelerated_restarts;
                accelerated_y_ = accelerated_next_;
                momentum = 1.0;
            } else {
                accelerated_y_ = accelerated_next_ +
                    ((momentum - 1.0) / next_momentum) *
                        (accelerated_next_ - coordinate_u_);
                momentum = next_momentum;
            }
            coordinate_u_.swap(accelerated_next_);
            if (!accelerated_y_.allFinite()) {
                violation = std::numeric_limits<double>::infinity();
                coordinate_x_ = inv_sqrt_diagonal_.cwiseProduct(coordinate_u_);
                return false;
            }

            if (accelerated_iterations % accelerated_violation_check_every != 0) continue;
            if (accelerated_kkt_satisfied()) return accept_solution_from_exact_gradient();
        }
        --accelerated_iterations;
        if (accelerated_iterations % accelerated_violation_check_every != 0 &&
            accelerated_kkt_satisfied())
            return accept_solution_from_exact_gradient();
        coordinate_x_ = inv_sqrt_diagonal_.cwiseProduct(coordinate_u_);
        return false;
    }

    double upper_bound_(const Vector& c) {
        if (use_closed_form_solution_)
            return c.cwiseMax(0.0).norm();

        bound_c_ = c.cwiseMax(0.0);
        if (bound_c_.squaredNorm() == 0.0) return 0.0;
        bound_solution_ = omega_solver_.solve(bound_c_);
        if (omega_solver_.info() != Eigen::Success || !bound_solution_.allFinite())
            return std::numeric_limits<double>::infinity();
        const double q = bound_c_.dot(bound_solution_);
        return q > 0.0 && std::isfinite(q) ?
            std::sqrt(q) : std::numeric_limits<double>::infinity();
    }

    SparseMatrix Psi_;
    SparseMatrix Omega_;
    bool objective_sign_invariant_ = true;
    bool use_closed_form_solution_ = false;
    Eigen::SimplicialLDLT<SparseMatrix> omega_solver_;

    Vector diagonal_;
    Vector sqrt_diagonal_;
    Vector inv_sqrt_diagonal_;
    Vector transformed_signal_;
    Vector direct_solution_;
    Vector projected_direct_;
    Vector coordinate_x_;
    Vector coordinate_u_;
    Vector coordinate_gradient_;
    Vector accelerated_y_;
    Vector accelerated_next_;
    Vector scaled_c_;
    Vector bound_c_;
    Vector bound_solution_;

    // Projected FISTA is the primary bounded solve; retaining PSOR sweeps here
    // is useful only as a measured compatibility/performance trade-off.
    static constexpr int coordinate_max_sweeps_ = 0;
    static constexpr double coordinate_relaxation_ = 1.8;
    static constexpr int accelerated_max_iterations_ = 25000;
    double accelerated_lipschitz_ = 0.0;

    Vector x_init_;
    Vector last_solution_pos_;
    Vector last_solution_neg_;
    Vector last_solution_;

    inline static thread_local NonNegativeWeightSolveStats thread_stats_ {};
};

} // namespace internals
} // namespace fdapde

#endif // __FDAPDE_NONNEGATIVE_WEIGHT_H__

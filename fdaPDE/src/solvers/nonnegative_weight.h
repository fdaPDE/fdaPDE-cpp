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

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "header_check.h"

namespace fdapde {
namespace internals {

using Vector = Eigen::VectorXd;
using SparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;

struct NonNegativeWeightSolveStats {
    std::uint64_t calls = 0;
    std::uint64_t zero_signals = 0;
    std::uint64_t closed_form_sides = 0;
    std::uint64_t direct_sides = 0;
    std::uint64_t coordinate_attempts = 0;
    std::uint64_t coordinate_converged = 0;
    std::uint64_t coordinate_sweeps = 0;
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

        const Vector p = Psi_.transpose() * z;
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
            Vector projected_direct;
            if (try_direct_solution_(c, solution, &projected_direct)) {
                ++thread_stats_.direct_sides;
                return true;
            }

            ++thread_stats_.coordinate_attempts;
            int sweeps = 0;
            double violation = std::numeric_limits<double>::infinity();
            double tolerance = 0.0;
            const Vector& warm_start = projected_direct.size() == solution.size() &&
                c.dot(projected_direct) > c.dot(solution) ? projected_direct : solution;
            const bool converged = try_coordinate_solution_(
                c, warm_start, &solution, sweeps, violation, tolerance
            );
            thread_stats_.coordinate_sweeps += static_cast<std::uint64_t>(sweeps);
            if (converged) {
                ++thread_stats_.coordinate_converged;
                return true;
            }

            ++thread_stats_.would_fallback;
            std::ostringstream message;
            message << "NonNegativeWeightSolver: fallback disabled; coordinate descent did not "
                    << "satisfy KKT after "
                    << sweeps << " sweeps (violation=" << std::scientific << std::setprecision(3)
                    << violation << ", tolerance=" << tolerance << ')';
            throw std::runtime_error(message.str());
        };

        const bool positive_side_only = p.minCoeff() >= 0.0;
        const bool negative_side_only = p.maxCoeff() <= 0.0;
        if (!objective_sign_invariant_ || positive_side_only) {
            pos_ok = solve_side(p, last_solution_pos_);
        } else if (negative_side_only) {
            neg_ok = solve_side(-p, last_solution_neg_);
        } else {
            const double pos_bound = upper_bound_(p);
            const double neg_bound = upper_bound_(-p);
            auto pos_dominates = [&]() {
                return pos_ok && p.dot(last_solution_pos_) +
                    1e-10 * (1.0 + std::abs(neg_bound)) >= neg_bound;
            };
            auto neg_dominates = [&]() {
                return neg_ok && -p.dot(last_solution_neg_) +
                    1e-10 * (1.0 + std::abs(pos_bound)) >= pos_bound;
            };

            if (pos_bound >= neg_bound) {
                pos_ok = solve_side(p, last_solution_pos_);
                if (!pos_dominates()) neg_ok = solve_side(-p, last_solution_neg_);
            } else {
                neg_ok = solve_side(-p, last_solution_neg_);
                if (!neg_dominates()) pos_ok = solve_side(p, last_solution_pos_);
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
        const Vector diagonal = Omega_.diagonal();
        int best_i = 0;
        double best_score = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < c.size(); ++i) {
            const double score = c[i] / std::sqrt(diagonal[i]);
            if (score > best_score) {
                best_score = score;
                best_i = i;
            }
        }
        Vector out = Vector::Zero(c.size());
        out[best_i] = 1.0 / std::sqrt(diagonal[best_i]);
        return out;
    }

    bool try_direct_solution_(const Vector& c, Vector& out, Vector* projected_out = nullptr) const {
        Vector y = omega_solver_.solve(c);
        if (omega_solver_.info() != Eigen::Success || !y.allFinite())
            return false;

        const double scale = std::max(1.0, y.cwiseAbs().maxCoeff());
        const double tol = 100.0 * std::numeric_limits<double>::epsilon() * scale;
        if (y.minCoeff() < -tol) {
            if (projected_out != nullptr) {
                Vector projected = y.cwiseMax(0.0);
                const double norm2 = projected.dot(Omega_ * projected);
                if (norm2 > 0.0 && std::isfinite(norm2))
                    *projected_out = projected / std::sqrt(norm2);
            }
            return false;
        }

        y = y.cwiseMax(0.0);
        const double norm2 = y.dot(Omega_ * y);
        if (!(norm2 > 0.0) || !std::isfinite(norm2))
            return false;
        out = y / std::sqrt(norm2);
        return true;
    }

    bool try_coordinate_solution_(
        const Vector& c,
        const Vector& warm_start,
        Vector* out,
        int& sweeps,
        double& violation,
        double& tolerance
    ) const {
        const int n = static_cast<int>(c.size());
        const Vector diagonal = Omega_.diagonal();
        const Vector sqrt_diagonal = diagonal.array().sqrt();
        const Vector inv_sqrt_diagonal = sqrt_diagonal.array().inverse();
        Vector x = warm_start.cwiseMax(0.0);
        const double warm_norm2 = x.dot(Omega_ * x);
        const double warm_scale = warm_norm2 > 0.0 ? std::max(0.0, c.dot(x) / warm_norm2) : 0.0;
        Vector u = sqrt_diagonal.cwiseProduct(x) * warm_scale;
        x = inv_sqrt_diagonal.cwiseProduct(u);

        Vector gradient = inv_sqrt_diagonal.cwiseProduct(Omega_ * x - c);
        const Vector scaled_c = inv_sqrt_diagonal.cwiseProduct(c);
        tolerance = 1e-8 * std::max(1.0, scaled_c.lpNorm<Eigen::Infinity>());
        constexpr int max_sweeps = 20000;
        constexpr double relaxation = 1.8;

        for (sweeps = 1; sweeps <= max_sweeps; ++sweeps) {
            for (int i = 0; i < n; ++i) {
                const double old_value = u[i];
                const double new_value = std::max(0.0, old_value - relaxation * gradient[i]);
                const double delta = new_value - old_value;
                if (delta == 0.0) continue;

                u[i] = new_value;
                for (SparseMatrix::InnerIterator it(Omega_, i); it; ++it)
                    gradient[it.row()] += delta * it.value() *
                        inv_sqrt_diagonal[it.row()] * inv_sqrt_diagonal[i];
            }

            if (sweeps % 20 == 0) {
                x = inv_sqrt_diagonal.cwiseProduct(u);
                gradient = inv_sqrt_diagonal.cwiseProduct(Omega_ * x - c);
            }

            violation = 0.0;
            for (int i = 0; i < n; ++i) {
                const double current = u[i] > 1e-14 ?
                    std::abs(gradient[i]) : std::max(0.0, -gradient[i]);
                violation = std::max(violation, current);
            }
            if (violation > tolerance) continue;

            x = inv_sqrt_diagonal.cwiseProduct(u);
            const double norm2 = x.dot(Omega_ * x);
            if (!(norm2 > 0.0) || !std::isfinite(norm2)) return false;
            *out = x / std::sqrt(norm2);
            return true;
        }
        --sweeps;
        return false;
    }

    double upper_bound_(const Vector& c) const {
        if (use_closed_form_solution_)
            return c.cwiseMax(0.0).norm();

        const Vector c_pos = c.cwiseMax(0.0);
        if (c_pos.squaredNorm() == 0.0) return 0.0;
        const Vector y = omega_solver_.solve(c_pos);
        if (omega_solver_.info() != Eigen::Success || !y.allFinite())
            return std::numeric_limits<double>::infinity();
        const double q = c_pos.dot(y);
        return q > 0.0 && std::isfinite(q) ?
            std::sqrt(q) : std::numeric_limits<double>::infinity();
    }

    SparseMatrix Psi_;
    SparseMatrix Omega_;
    bool objective_sign_invariant_ = true;
    bool use_closed_form_solution_ = false;
    Eigen::SimplicialLDLT<SparseMatrix> omega_solver_;

    Vector x_init_;
    Vector last_solution_pos_;
    Vector last_solution_neg_;
    Vector last_solution_;

    inline static thread_local NonNegativeWeightSolveStats thread_stats_ {};
};

} // namespace internals
} // namespace fdapde

#endif // __FDAPDE_NONNEGATIVE_WEIGHT_H__

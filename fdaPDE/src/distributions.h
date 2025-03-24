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

#ifndef __FDAPDE_DISTRIBUTIONS_H__
#define __FDAPDE_DISTRIBUTIONS_H__

namespace fdapde {
namespace internals {

template <typename Distribution_> class distribution_base {
   public:
    using Distribution = Distribution_;
    distribution_base() noexcept = default;
    template <typename... Args> distribution_base(Args&&... args) : distr_(std::forward<Args>(args)...) { }
   protected:
    Distribution distr_;

    template <typename T, typename F>
        requires(
          internals::is_vector_like_v<T> &&
          std::is_convertible_v<internals::subscript_result_of_t<T, int>, double>)
    std::vector<double> apply_(const T& x, F&& f) const {
        std::vector<double> res(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) { res[i] = f(x[i]); }
        return res;
    }
};

}   // namespace internals

struct bernoulli_distribution : public internals::distribution_base<std::bernoulli_distribution> {
    using result_type = double;
    using param_type  = double;
   private:
    using Base = internals::distribution_base<std::bernoulli_distribution>;
    using Base::distr_;
    param_type p_ = 0;
   public:
    bernoulli_distribution() noexcept = default;
    explicit bernoulli_distribution(param_type p) noexcept : Base(p), p_(p) { };
    // density function
    template <typename InputType>
        requires(std::is_convertible_v<InputType, bool>)
    constexpr result_type pdf(InputType x) const {
        return bool(x) == true ? 1 - p_ : p_;
    };
    constexpr result_type cdf(double x) const { return x < 0 ? 0 : ((x >= 0 && x < 1) ? 1 - p_ : 1); }
    constexpr result_type mean() const { return p_; }
    constexpr result_type variance() const { return p_ * (1 - p_); }
    constexpr result_type quantile(double alpha) const {
        fdapde_assert(alpha >= 0 && alpha <= 1);
        return alpha <= 1 - p_ ? 0 : 1;
    }
    // random sampling
    template <typename RandomNumberGenerator> result_type operator()(RandomNumberGenerator& rng) { return distr_(rng); }

    template <typename T> constexpr const T& mean(const T& data) const { return data; }
    template <typename T> requires(!internals::is_eigen_dense_xpr_v<T>) constexpr auto variance(const T& data) const {
        return Base::apply_(data, [](auto v) { return v * (1 - v); });
    }
    template <typename T> requires(!internals::is_eigen_dense_xpr_v<T>) constexpr auto link(const T& data) const {
        return Base::apply_(data, [](auto v) { return std::log(v / (1 - v)); });
    }
    template <typename T> requires(!internals::is_eigen_dense_xpr_v<T>) constexpr auto inv_link(const T& data) const {
        return Base::apply_(data, [](auto v) { return 1.0 / (1 + std::exp(-v)); });
    }
    template <typename T> requires(!internals::is_eigen_dense_xpr_v<T>) constexpr auto der_link(const T& data) const {
        return Base::apply_(data, [](auto v) { return 1.0 / (v * (1 - v)); });
    }
    template <typename T>
        requires(std::is_floating_point_v<T>)
    constexpr result_type deviance(T x, T y) const {
        return almost_zero(y) ? 2 * std::log(1.0 / (1.0 - x)) : 2.0 * std::log(1.0 / x);
    }
    // SIMD vectorized
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    matrix_t variance(const matrix_t& x) const { return x.array() * (1 - x.array()); }
    matrix_t link    (const matrix_t& x) const { return ((1 - x.array()).inverse() * x.array()).log(); }
    matrix_t inv_link(const matrix_t& x) const { return (1 + ((-x).array().exp())).inverse(); }
    matrix_t der_link(const matrix_t& x) const { return (x.array() * (1 - x.array())).inverse(); }

    template <typename T> constexpr auto transform(const T& data) const {
        if constexpr (internals::is_eigen_dense_xpr_v<T>) {
            return 0.5 * (data.array() + 0.5);
        } else {
            return Base::apply_(data, [](auto v) { return 0.5 * (v + 0.5); });
        }
    }
    void set_param(param_type p) { p_ = p; }
};

struct rademacher_distribution : public internals::distribution_base<std::bernoulli_distribution> {
    using result_type = double;
    using param_type  = double;
   private:
    using Base = internals::distribution_base<std::bernoulli_distribution>;
    using Base::distr_;
   public:
    rademacher_distribution() noexcept : Base(0.5) { }
    // density function
    template <typename InputType> constexpr result_type pdf(InputType x) const {
        return (x == 1 || x == -1) ? 0.5 : 0.0;
    }
    constexpr result_type cdf(double x) const { return x < -1 ? 0 : ((-1 <= x < 1) ? 0.5 : 1.0); }
    constexpr result_type mean() const { return 0.0; }
    constexpr result_type variance() const { return 1.0; }
    // random sampling
    template <typename RandomNumberGenerator> result_type operator()(RandomNumberGenerator& rng) {
        return distr_(rng) ? 1.0 : -1.0;
    }
};

struct poisson_distribution : public internals::distribution_base<std::poisson_distribution<int>> {
    using result_type = double;
    using param_type  = double;
   private:
    using Base = internals::distribution_base<std::poisson_distribution<int>>;
    using Base::distr_;
    param_type l_;
   public:
    poisson_distribution() noexcept = default;
    explicit poisson_distribution(param_type l) noexcept : Base(l), l_(l) { };
    // density function
    template <typename InputType>
        requires(std::is_convertible_v<InputType, std::size_t>)
    constexpr result_type pdf(InputType k) const {
        return std::pow(l_, k) * std::exp(-l_) / factorial(k);
    }
    constexpr result_type cdf(double x) const {
        double cdf_ = 0;
        for (int i = 0; i < std::floor(x); ++i) { cdf_ += std::pow(l_, i) / factorial(i); }
        return cdf_ * std::exp(-l_);
    }
    constexpr result_type mean() const { return l_; }
    constexpr result_type variance() const { return l_; }
    constexpr result_type quantile(double alpha) const {
        fdapde_assert(alpha >= 0 && alpha <= 1);
        int x = 0;
        double cdf_ = 0;
        while (cdf_ < alpha) {
            cdf_ += std::exp(-l_) * std::pow(l_, x) / factorial(x);
            x++;
        }
        return std::floor(x);
    }
    // random sampling
    template <typename RandomNumberGenerator> result_type operator()(RandomNumberGenerator& rng) { return distr_(rng); }

    template <typename T> constexpr const T& mean(const T& data) const { return data; }
    template <typename T> constexpr const T& variance(const T& data) const { return data; }
    template <typename T> requires(!internals::is_eigen_dense_xpr_v<T>) constexpr auto link(const T& data) const {
        return Base::apply_(data, [](auto v) { return std::log(v); });
    }
    template <typename T> requires(!internals::is_eigen_dense_xpr_v<T>) constexpr auto inv_link(const T& data) const {
        return Base::apply_(data, [](auto v) { return std::exp(v); });
    }
    template <typename T> requires(!internals::is_eigen_dense_xpr_v<T>) constexpr auto der_link(const T& data) const {
        return Base::apply_(data, [](auto v) { return 1.0 / v; });
    }
    template <typename T>
        requires(std::is_floating_point_v<T>)
    constexpr result_type deviance(T x, T y) const {
        return y > 0 ? y * std::log(y / x) - (y - x) : x;
    }
    // SIMD vectorized
    using matrix_t = Eigen::Matrix<double, Dynamic, 1>;
    matrix_t variance(const matrix_t& x) const { return x; }
    matrix_t link    (const matrix_t& x) const { return x.array().log(); }
    matrix_t inv_link(const matrix_t& x) const { return x.array().exp(); }
    matrix_t der_link(const matrix_t& x) const { return x.array().inverse(); }

    template <typename T> constexpr auto transform(const T& data) const {
        if constexpr (internals::is_eigen_dense_xpr_v<T>) {
            return matrix_t((data.array() <= 0).select(1.0, data));
        } else {
            return Base::apply_(data, [](auto v) { return v <= 0 ? 1.0 : v; });
        }
    }
    void set_param(param_type l) { l_ = l; }
};

struct exponential_distribution : public internals::distribution_base<std::exponential_distribution<double>> {
    using result_type = double;
    using param_type  = double;
   private:
    using Base = internals::distribution_base<std::exponential_distribution<double>>;
    using Base::distr_;
    param_type l_;
   public:
    exponential_distribution() noexcept = default;
    explicit exponential_distribution(param_type l) : Base(l), l_(l) { };
    // density function
    template <typename InputType>
        requires(std::is_convertible_v<InputType, double>)
    constexpr result_type pdf(InputType x) const {
        return l_ * std::exp(-l_ * x);
    };
    constexpr result_type cdf(double x) const { return 1.0 - std::exp(-l_ * x); }
    constexpr result_type mean() const { return 1.0 / l_; }
    constexpr result_type variance() const { return 1.0 / (l_ * l_); }
    constexpr result_type quantile(double alpha) const { return -std::log(1 - alpha) / l_; }
    // random sampling
    template <typename RandomNumberGenerator> result_type operator()(RandomNumberGenerator& rng) { return distr_(rng); }

    template <typename T> requires(!internals::is_eigen_dense_xpr_v<T>) constexpr const T& mean(const T& data) const {
        return Base::apply_(data, [](auto v) { return 1.0 / v; });
    }
    template <typename T>
        requires(!internals::is_eigen_dense_xpr_v<T>)
    constexpr const T& variance(const T& data) const {
        return Base::apply_(data, [](auto v) { return v * v; });
    }
    template <typename T> constexpr auto link(const T& data) const {
        return Base::apply_(data, [](auto v) { return -1.0 / v; });
    }
    template <typename T> constexpr auto inv_link(const T& data) const {
        return Base::apply_(data, [](auto v) { return -1.0 / v; });
    }
    template <typename T> constexpr auto der_link(const T& data) const {
        return Base::apply_(data, [](auto v) { return 1.0 / (v * v); });
    }
    template <typename T>
        requires(std::is_floating_point_v<T>)
    constexpr result_type deviance(T x, T y) const {
        return 2 * ((y - x) / x - std::log(y / x));
    }
#ifdef __FDAPDE_HAS_EIGEN__  // SIMD vectorized
    using matrix_t = Eigen::Matrix<double, Dynamic, 1>;
    matrix_t variance(const matrix_t& x) const { return x.array().pow(2); }
    matrix_t link    (const matrix_t& x) const { return (-x).array().inverse(); }
    matrix_t inv_link(const matrix_t& x) const { return (-x).array().inverse(); }
    matrix_t der_link(const matrix_t& x) const { return x.array().pow(2).inverse(); }
#endif
    void set_param(param_type l) { l_ = l; }
};

class gamma_distribution : public internals::distribution_base<std::gamma_distribution<double>> {
    using result_type = double;
    using param_type  = double;
   private:
    using Base = internals::distribution_base<std::gamma_distribution<double>>;
    using Base::distr_;
    param_type k_;       // shape parameter
    param_type theta_;   // scale parameter
   public:
    gamma_distribution() noexcept = default;
    gamma_distribution(double k, double theta) : Base(k, theta), k_(k), theta_(theta) { };
    // density function
    template <typename InputType>
        requires(std::is_convertible_v<InputType, double>)
    constexpr result_type pdf(InputType x) const {
        return 1 / (std::tgamma(k_) * std::pow(theta_, k_)) * std::pow(x, k_ - 1) * std::exp(-x / theta_);
    }
    constexpr result_type mean() const { return k_ * theta_; }
    constexpr result_type variance() const { return k_ * theta_ * theta_; }
    // random sampling
    template <typename RandomNumberGenerator> result_type operator()(RandomNumberGenerator& rng) { return distr_(rng); }
  
    template <typename T> constexpr const T& mean(const T& data) const {
        return Base::apply_(data, [](auto v) { return v; });
    }
    template <typename T> constexpr const T& variance(const T& data) const {
        return Base::apply_(data, [](auto v) { return v * v; });
    }
    template <typename T> constexpr auto link(const T& data) const {
        return Base::apply_(data, [](auto v) { return -1.0 / v; });
    }
    template <typename T> constexpr auto inv_link(const T& data) const {
        return Base::apply_(data, [](auto v) { return -1.0 / v; });
    }
    template <typename T> constexpr auto der_link(const T& data) const {
        return Base::apply_(data, [](auto v) { return 1.0 / (v * v); });
    }
    template <typename T>
        requires(std::is_floating_point_v<T>)
    constexpr result_type deviance(T x, T y) const {
        return 2 * ((y - x) / x - std::log(y / x));
    }
#ifdef __FDAPDE_HAS_EIGEN__  // SIMD vectorized
    using matrix_t = Eigen::Matrix<double, Dynamic, 1>;
    matrix_t variance(const matrix_t& x) const { return x.array().pow(2); }
    matrix_t link    (const matrix_t& x) const { return (-x).array().inverse(); }
    matrix_t inv_link(const matrix_t& x) const { return (-x).array().inverse(); }
    matrix_t der_link(const matrix_t& x) const { return x.array().pow(2).inverse(); }
#endif
    void set_param(param_type k, param_type theta) {
        k_ = k;
        theta_ = theta;
    }
};

[[maybe_unused]] bernoulli_distribution   Bernoulli {};
[[maybe_unused]] poisson_distribution     Poisson {};
[[maybe_unused]] exponential_distribution Exponential {};
[[maybe_unused]] gamma_distribution       Gamma {};

}   // namespace fdapde

#endif   // __FDAPDE_DISTRIBUTIONS_H__

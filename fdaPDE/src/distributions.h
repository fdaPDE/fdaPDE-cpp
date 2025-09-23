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

template <typename FunctorT>
    requires(requires(FunctorT f, double x) { f(x); })
constexpr double adaptive_simpson_integrate(FunctorT&& f, double a, double b, double eps = 1e-10) {
    struct interval_t {
        double a, b;         // intervals extremes
        double fa, fb, fm;   // value of the objective at a, b and midpoint m = (a + b) / 2
        double S;            // estimated integral at interval [a, b]
    };

    std::stack<interval_t> intervals;
    double value = 0.0;
    // initialization
    double m = 0.5 * (a + b);
    double fa = f(a);
    double fb = f(b);
    double fm = f(m);
    double S = (b - a) / 6.0 * (fa + 4 * fm + fb);   // simpson approximation for \int_a^b f(t) dt
    intervals.emplace(a, b, fa, fb, fm, S);

    while (!intervals.empty()) {
        interval_t i = intervals.top();
        intervals.pop();
        double a = i.a, b = i.b, fa = i.fa, fb = i.fb, fm = i.fm, S = i.S;

        double m = 0.5 * (a + b);
        // left subinterval
        double lm = 0.5 * (a + m);
        double flm = f(lm);
        double Sl = (m - a) / 6.0 * (fa + 4 * flm + fm);
        // right subinterval
        double rm = 0.5 * (m + b);
        double frm = f(rm);
        double Sr = (b - m) / 6.0 * (fm + 4 * frm + fb);

        double S_ = Sl + Sr;

        if (std::abs(S_ - S) < 15 * eps) {
            value += S_ + (S_ - S) / 15.0;   // Richardson extrapolation
        } else {
            // error too high, recurse
            intervals.emplace(m, b, fm, fb, frm, Sr);   // left  subinterval
            intervals.emplace(a, m, fa, fm, flm, Sl);   // right subinterval
        }
    }
    return value;
}

// computes the lower incomplete gamma function gamma(a, x) = \int_0^x (t^{a-1} * exp(-t))dt
inline double lower_incomplete_gamma(double a, double x) {
    if (almost_zero(x)) { return 0.0; }
    // lower incomplete gamma integrand
    auto f = [a](double t) {
        if (t == 0.0) {
            // evaluate at small epsilon to avoid singularity at t = 0 when a < 1
            double eps = std::numeric_limits<double>::epsilon();
            return std::pow(eps, a - 1) * std::exp(-eps);
        }
        return std::pow(t, a - 1) * std::exp(-t);
    };
    return adaptive_simpson_integrate(f, 0.0, x, 1e-10);   // integral approximation by adaptive Simpson rule
}
// normalized lower incomplete gamma function
inline double gamma_p(double a, double x) { return lower_incomplete_gamma(a, x) / std::tgamma(a); }
// inverse error function (based on Newton-Rapson root finder)
inline double inverse_erf(double x, double eps = 1e-10) {
    double y = 0.0;
    double delta;
    do {
        delta = (std::erf(y) - x) / (2.0 / std::sqrt(std::numbers::pi) * std::exp(-y * y));
        y -= delta;
    } while (std::fabs(delta) > eps);
    return y;
}
  
}   // namespace internals

class simd_distribution {
#ifdef __FDAPDE_HAS_EIGEN__   // SIMD vectorized
   public:
    using matrix_t = Eigen::Matrix<double, Dynamic, 1>;
    constexpr simd_distribution() = default;
  
    virtual matrix_t variance(const matrix_t& x) const = 0;
    virtual matrix_t link    (const matrix_t& x) const = 0;
    virtual matrix_t inv_link(const matrix_t& x) const = 0;
    virtual matrix_t der_link(const matrix_t& x) const = 0;
    virtual double   deviance(const matrix_t& x, const matrix_t& y) const = 0;
#endif
    virtual ~simd_distribution() = default;
};

namespace internals {
  
template <typename Distribution_> class distribution_base : public simd_distribution {
   public:
    using Distribution = Distribution_;
    constexpr distribution_base() noexcept { }
    template <typename... Args> constexpr distribution_base(Args&&... args) : distr_(std::forward<Args>(args)...) { }
   protected:
    Distribution distr_;

    template <typename T, typename F>
        requires(
          internals::is_vector_like_v<T> && std::is_convertible_v<internals::subscript_result_of_t<T, int>, double>)
    constexpr std::vector<double> apply_(const T& x, F&& f) const {
        std::vector<double> res(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) { res[i] = f(x[i]); }
        return res;
    }
};

}   // namespace internals

struct bernoulli_distribution : public internals::distribution_base<std::bernoulli_distribution> {
    using result_type = double;
    using param_type  = double;
   protected:
    using Base = internals::distribution_base<std::bernoulli_distribution>;
    using Base::distr_;
    param_type p_ = 0;
   public:
    constexpr bernoulli_distribution() noexcept = default;
    constexpr explicit bernoulli_distribution(param_type p) noexcept : Base(p), p_(p) { };
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
#ifdef __FDAPDE_HAS_EIGEN__   // SIMD vectorized
    using matrix_t = Eigen::Matrix<double, Dynamic, 1>;
    matrix_t variance(const matrix_t& x) const override { return x.array() * (1 - x.array()); }
    matrix_t link    (const matrix_t& x) const override { return ((1 - x.array()).inverse() * x.array()).log(); }
    matrix_t inv_link(const matrix_t& x) const override { return (1 + ((-x).array().exp())).inverse(); }
    matrix_t der_link(const matrix_t& x) const override { return (x.array() * (1 - x.array())).inverse(); }
    double deviance(const matrix_t& x, const matrix_t& y) const override {
        fdapde_assert(x.cols() == 1 && y.cols() == 1 && x.rows() == y.rows());
        double dev_ = 0;
        for (int i = 0, n = x.rows(); i < n; ++i) { dev_ += deviance(x(i, 0), y(i, 0)); }
        return dev_;
    }
#endif
    template <typename T> constexpr auto transform(const T& data) const {
        if constexpr (internals::is_eigen_dense_xpr_v<T>) {
            return 0.5 * (data.array() + 0.5);
        } else {
            return Base::apply_(data, [](auto v) { return 0.5 * (v + 0.5); });
        }
    }
    void set_param(param_type p) { p_ = p; }
};

struct rademacher_distribution : public bernoulli_distribution {
    using result_type = double;
    using param_type = double;
   private:
    using Base = bernoulli_distribution;
    using Base::distr_;
   public:
    constexpr rademacher_distribution() noexcept : Base(0.5) { }
    // density function
    template <typename InputType> constexpr result_type pdf(InputType x) const {
        return (x == 1 || x == -1) ? 0.5 : 0.0;
    }
    constexpr result_type cdf(double x) const { return x < -1 ? 0 : ((-1 <= x && x < 1) ? 0.5 : 1.0); }
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
    constexpr poisson_distribution() noexcept = default;
    constexpr explicit poisson_distribution(param_type l) noexcept : Base(l), l_(l) { };
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
#ifdef __FDAPDE_HAS_EIGEN__   // SIMD vectorized
    using matrix_t = Eigen::Matrix<double, Dynamic, 1>;
    matrix_t variance(const matrix_t& x) const override { return x; }
    matrix_t link    (const matrix_t& x) const override { return x.array().log(); }
    matrix_t inv_link(const matrix_t& x) const override { return x.array().exp(); }
    matrix_t der_link(const matrix_t& x) const override { return x.array().inverse(); }
    double deviance(const matrix_t& x, const matrix_t& y) const override {
        fdapde_assert(x.cols() == 1 && y.cols() == 1 && x.rows() == y.rows());
        return ((y.array() > 0).select(y.array() * ((y.array() / x.array()).log() - 1) + x.array(), x.array()))
          .sum();
    }
#endif
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
    constexpr exponential_distribution() noexcept = default;
    constexpr explicit exponential_distribution(param_type l) : Base(l), l_(l) { };
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
#ifdef __FDAPDE_HAS_EIGEN__   // SIMD vectorized
    using matrix_t = Eigen::Matrix<double, Dynamic, 1>;
    matrix_t variance(const matrix_t& x) const override { return x.array().pow(2); }
    matrix_t link    (const matrix_t& x) const override { return (-x).array().inverse(); }
    matrix_t inv_link(const matrix_t& x) const override { return (-x).array().inverse(); }
    matrix_t der_link(const matrix_t& x) const override { return x.array().pow(2).inverse(); }
    double deviance(const matrix_t& x, const matrix_t& y) const override {
        fdapde_assert(x.cols() == 1 && y.cols() == 1 && x.rows() == y.rows());
        return (2 * ((y.array() - x.array()) / x.array() - (y.array() / x.array()).log())).sum();
    }
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
    constexpr gamma_distribution() noexcept = default;
    constexpr gamma_distribution(double k, double theta) : Base(k, theta), k_(k), theta_(theta) { }
    // density function
    template <typename InputType>
        requires(std::is_convertible_v<InputType, double>)
    result_type pdf(InputType x) const {
        return 1 / (std::tgamma(k_) * std::pow(theta_, k_)) * std::pow(x, k_ - 1) * std::exp(-x / theta_);
    }
    result_type cdf(double x) const { return internals::gamma_p(k_, x / theta_); }
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
#ifdef __FDAPDE_HAS_EIGEN__   // SIMD vectorized
    using matrix_t = Eigen::Matrix<double, Dynamic, 1>;
    matrix_t variance(const matrix_t& x) const override { return x.array().pow(2); }
    matrix_t link    (const matrix_t& x) const override { return (-x).array().inverse(); }
    matrix_t inv_link(const matrix_t& x) const override { return (-x).array().inverse(); }
    matrix_t der_link(const matrix_t& x) const override { return x.array().pow(2).inverse(); }
    double deviance(const matrix_t& x, const matrix_t& y) const override {
        fdapde_assert(x.cols() == 1 && y.cols() == 1 && x.rows() == y.rows());
        return (2 * ((y.array() - x.array()) / x.array() - (y.array() / x.array()).log())).sum();
    }
#endif
    void set_param(param_type k, param_type theta) {
        k_ = k;
        theta_ = theta;
    }
};

class normal_distribution : public internals::distribution_base<std::normal_distribution<double>> {
    using result_type = double;
    using param_type  = double;
  private:
    using Base = internals::distribution_base<std::normal_distribution<double>>;
    using Base::distr_;
    param_type mu_;
    param_type sigma_;
   public:
    constexpr normal_distribution() noexcept = default;
    constexpr normal_distribution(param_type mu, param_type sigma) : Base(mu, sigma), mu_(mu), sigma_(sigma) { }
    // density function
    result_type pdf(double x) const {
        constexpr double pi = std::numbers::pi;
        return 1.0 / (std::sqrt(2 * pi) * sigma_) * std::exp(-std::pow(x - mu_, 2) / (2 * std::pow(sigma_, 2)));
    }
    result_type cdf(double x) const { return 0.5 * (1 + std::erf((x - mu_) / (sigma_ * std::sqrt(2)))); }
    constexpr param_type mean() const { return mu_; }
    constexpr param_type variance() const { return sigma_ * sigma_; }
    double quantile(double alpha) const { return std::sqrt(2.0) * internals::inverse_erf(2.0 * alpha - 1.0); }
    // random sampling
    template <typename RandomNumberGenerator> result_type operator()(RandomNumberGenerator& rng) { return distr_(rng); }

    template <typename T> constexpr const T& mean(const T& data) const {
        return Base::apply_(data, []([[maybe_unused]] auto v) { return 1.0; });
    }
    template <typename T> constexpr const T& variance(const T& data) const {
        return Base::apply_(data, []([[maybe_unused]] auto v) { return 1.0; });
    }
    template <typename T> constexpr auto link(const T& data) const {
        return Base::apply_(data, [](auto v) { return v; });
    }
    template <typename T> constexpr auto inv_link(const T& data) const {
        return Base::apply_(data, [](auto v) { return v; });
    }
    template <typename T> constexpr auto der_link(const T& data) const {
        return Base::apply_(data, []([[maybe_unused]] auto v) { return 1.0; });
    }
    template <typename T>
        requires(std::is_floating_point_v<T>)
    constexpr result_type deviance(T x, T y) const {
        return (x - y) * (x - y);
    }
#ifdef __FDAPDE_HAS_EIGEN__   // SIMD vectorized
    using matrix_t = Eigen::Matrix<double, Dynamic, 1>;
    matrix_t variance(const matrix_t& x) const override { return matrix_t::Ones(x.rows()); }
    matrix_t link    (const matrix_t& x) const override { return x; }
    matrix_t inv_link(const matrix_t& x) const override { return x; }
    matrix_t der_link(const matrix_t& x) const override { return matrix_t::Ones(x.rows()); }
    double deviance(const matrix_t& x, const matrix_t& y) const override {
        fdapde_assert(x.cols() == 1 && y.cols() == 1 && x.rows() == y.rows());
        return (x - y).squaredNorm();
    }
#endif
    void set_param(param_type mu, param_type sigma) {
        mu_ = mu;
        sigma_ = sigma;
    }
};

class chi_squared_distribution : public internals::distribution_base<std::chi_squared_distribution<double>> {
    using result_type = double;
    using param_type  = double;
   private:
    using Base = internals::distribution_base<std::chi_squared_distribution<double>>;
    using Base::distr_;
    param_type n_;   // degrees of freedom
    gamma_distribution gamma_;
   public:
    constexpr chi_squared_distribution() noexcept = default;
    constexpr chi_squared_distribution(param_type n) : n_(n), gamma_(n / 2, 2) { }
    constexpr chi_squared_distribution(param_type n, param_type s) :
        n_(n), gamma_(n / 2, 2 * std::pow(s, 2)) { }   // scaled constructor
    // density function
    result_type pdf(double x) const { return gamma_.pdf(x); }
    result_type cdf(double x) const { return gamma_.cdf(x); }
    constexpr result_type mean() const { return n_; }
    constexpr result_type variance() const { return 2 * n_; }
    // quantile function (implemented as a binary search loop)
    double quantile(double alpha, double tol = 1e-6) {
        // support range [ql, qh] where quantile is searched
        double ql = 0.0;
        double qh = 1000.0;
        // binary search loop
        while (qh - ql > tol) {
            double q = (ql + qh) / 2.0;
            double y = cdf(q);
            if (y < alpha) {
                ql = q;
            } else {
                qh = q;
            }
        }
        return (ql + qh) / 2.0;
    }
    constexpr param_type dofs() const { return n_; }
    // random sampling
    template <typename RandomNumberGenerator> result_type operator()(RandomNumberGenerator& rng) { return distr_(rng); }

    template <typename T> constexpr const T& mean(const T& data) const { return gamma_.mean(data); }
    template <typename T> constexpr const T& variance(const T& data) const { return gamma_.variance(data); }
    template <typename T> constexpr auto link(const T& data) const { return gamma_.link(data); }
    template <typename T> constexpr auto inv_link(const T& data) const { return gamma_.inv_link(data); }
    template <typename T> constexpr auto der_link(const T& data) const { return gamma_.der_link(data); }
    template <typename T>
        requires(std::is_floating_point_v<T>)
    constexpr result_type deviance(T x, T y) const {
        return gamma_.deviance(x, y);
    }
#ifdef __FDAPDE_HAS_EIGEN__  // SIMD vectorized
    using matrix_t = Eigen::Matrix<double, Dynamic, 1>;
    matrix_t variance(const matrix_t& x) const override { return gamma_.variance(x); }
    matrix_t link    (const matrix_t& x) const override { return gamma_.link(x); }
    matrix_t inv_link(const matrix_t& x) const override { return gamma_.inv_link(x); }
    matrix_t der_link(const matrix_t& x) const override { return gamma_.der_link(x); }
    double deviance(const matrix_t& x, const matrix_t& y) const override { return gamma_.deviance(x, y); }
#endif
    void set_param(param_type n) { n_ = n; }
};
  
}   // namespace fdapde

#endif   // __FDAPDE_DISTRIBUTIONS_H__

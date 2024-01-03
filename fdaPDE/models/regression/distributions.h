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

#ifndef __DISTRIBUTION_H__
#define __DISTRIBUTION_H__

#include <fdaPDE/utils.h>

#include <cmath>
#include <cstddef>

namespace fdapde {
namespace models {

  // distribution concept
struct IDistribution {
    template <typename D>
    using fn_ptrs =
      fdapde::mem_fn_ptrs<&D::preprocess, &D::variance, &D::link, &D::inv_link, &D::der_link, &D::deviance>;

    decltype(auto) preprocess(const DVector<double>& data) const {
        return fdapde::invoke<DVector<double>, 0>(*this, data);
    }
    decltype(auto) variance(const DMatrix<double>& x) const { return fdapde::invoke<DMatrix<double>, 1>(*this, x); }
    decltype(auto) link(const DMatrix<double>& x)     const { return fdapde::invoke<DMatrix<double>, 2>(*this, x); }
    decltype(auto) inv_link(const DMatrix<double>& x) const { return fdapde::invoke<DMatrix<double>, 3>(*this, x); }
    decltype(auto) der_link(const DMatrix<double>& x) const { return fdapde::invoke<DMatrix<double>, 4>(*this, x); }
    decltype(auto) deviance(double x, double y)       const { return fdapde::invoke<double         , 5>(*this, x, y); }
};
using Distribution = fdapde::erase<fdapde::heap_storage, IDistribution>;
  
class Bernulli {
   private:
    double p_;   // distribution parameter
   public:
    // constructor
    Bernulli() = default;
    Bernulli(double p) : p_(p) {};
    // density function
    double pdf(std::size_t x) const { return x == 0 ? 1 - p_ : p_; };
    double mean() const { return p_; }
    // preprocess data to work with this distribution
    DVector<double> preprocess(const DVector<double>& data) const { return 0.5 * (data.array() + 0.5); }
    // vectorized operations
    DMatrix<double> variance(const DMatrix<double>& x) const { return x.array() * (1 - x.array()); }
    DMatrix<double> link    (const DMatrix<double>& x) const { return ((1 - x.array()).inverse() * x.array()).log(); }
    DMatrix<double> inv_link(const DMatrix<double>& x) const { return (1 + ((-x).array().exp())).inverse(); }
    DMatrix<double> der_link(const DMatrix<double>& x) const { return (x.array() * (1 - x.array())).inverse(); }
    // deviance function
    double deviance(double x, double y) const { return y == 0 ? 2 * std::log(1 / (1 - x)) : 2 * std::log(1 / x); };
};

class Poisson {
   private:
    double l_;   // distribution parameter
    // a simple utility to compute the factorial of an integer number
    std::size_t factorial(std::size_t k) const {
        std::size_t result = 1;
        while (k > 0) {
            result *= k;
            --k;
        }
        return result;
    }
   public:
    // constructor
    Poisson() = default;
    Poisson(double l) : l_(l) {};
    // density function
    double pdf(std::size_t k) const { return std::pow(l_, k) * std::exp(-l_) / factorial(k); };
    double mean() const { return l_; }
    // preprocess data to work with this distribution
    DVector<double> preprocess(const DVector<double>& data) const { return (data.array() <= 0).select(1.0, data); }
    // vectorized operations
    DMatrix<double> variance(const DMatrix<double>& x) const { return x; }
    DMatrix<double> link    (const DMatrix<double>& x) const { return x.array().log(); }
    DMatrix<double> inv_link(const DMatrix<double>& x) const { return x.array().exp(); }
    DMatrix<double> der_link(const DMatrix<double>& x) const { return x.array().inverse(); }
    // deviance function
    double deviance(double x, double y) const { return y > 0 ? y * std::log(y / x) - (y - x) : x; };
};

class Exponential {
   private:
    double l_;   // distribution parameter
   public:
    // constructor
    Exponential() = default;
    Exponential(double l) : l_(l) {};
    // density function
    double pdf(double x) const { return l_ * std::exp(-l_ * x); };
    double mean() const { return 1 / l_; }
    DVector<double> preprocess(const DVector<double>& data) const { return data; }
    // vectorized operations
    DMatrix<double> variance(const DMatrix<double>& x) const { return x.array().pow(2); }
    DMatrix<double> link    (const DMatrix<double>& x) const { return (-x).array().inverse(); }
    DMatrix<double> inv_link(const DMatrix<double>& x) const { return (-x).array().inverse(); }
    DMatrix<double> der_link(const DMatrix<double>& x) const { return x.array().pow(2).inverse(); }
    // deviance function
    double deviance(double x, double y) const { return 2 * ((y - x) / x - std::log(y / x)); };
};

class Gamma {
   private:
    double k_;       // shape parameter
    double theta_;   // scale parameter
   public:
    // constructor
    Gamma() = default;
    Gamma(double k, double theta) : k_(k), theta_(theta) {};
    // density function
    double pdf(double x) const {
        return 1 / (std::tgamma(k_) * std::pow(theta_, k_)) * std::pow(x, k_ - 1) * std::exp(-x / theta_);
    };
    double mean() const { return k_ * theta_; }
    DVector<double> preprocess(const DVector<double>& data) const { return data; }
    // vectorized operations
    DMatrix<double> variance(const DMatrix<double>& x) const { return x.array().pow(2); }
    DMatrix<double> link    (const DMatrix<double>& x) const { return (-x).array().inverse(); }
    DMatrix<double> inv_link(const DMatrix<double>& x) const { return (-x).array().inverse(); }
    DMatrix<double> der_link(const DMatrix<double>& x) const { return x.array().pow(2).inverse(); }
    // deviance function
    double deviance(double x, double y) const { return 2 * ((y - x) / x - std::log(y / x)); };
};

class Gaussian {
   private:
    double mu_;
    double sigma_;
   public:
    // constructor
    Gaussian() = default;
    Gaussian(double mu, double sigma) : mu_(mu), sigma_(sigma) {};
    // density function
    double pdf(double x) const {
        return 1 / (std::sqrt(2 * M_PI) * sigma_) * std::exp(-(x - mu_) * (x - mu_) / (2 * sigma_ * sigma_));
    };
    double mean() const { return mu_; }
    DVector<double> preprocess(const DVector<double>& data) const { return data; }
    // vectorized operations
    DMatrix<double> variance(const DMatrix<double>& x) const { return DMatrix<double>::Ones(x.rows(), x.cols()); }
    DMatrix<double> link    (const DMatrix<double>& x) const { return x; }
    DMatrix<double> inv_link(const DMatrix<double>& x) const { return x; }
    DMatrix<double> der_link(const DMatrix<double>& x) const { return DMatrix<double>::Ones(x.rows(), x.cols()); }
    // deviance function
    double deviance(double x, double y) const { return (x - y) * (x - y); };
};

}   // namespace models
}   // namespace fdapde

#endif   // __DISTRIBUTION_H__

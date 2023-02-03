#ifndef __DISTRIBUTION_H__
#define __DISTRIBUTION_H__

#include "../../core/utils/Symbols.h"
#include <cmath>
#include <cstddef>

namespace fdaPDE {
namespace models {

  // definition of common distributions from the exponential family

  class Bernulli {
  private:
    double p_; // distribution parameter
  public:
    // constructor
    Bernulli() = default;
    Bernulli(double p) : p_(p) {};
    // density function
    double pdf(std::size_t x) const { return x == 0 ? 1-p_ : p_; };
    double mean() const { return p_; }
    void preprocess(DVector<double>& data) const {
      for(std::size_t i = 0; i < data.rows(); ++i){
	data[i] = 0.5 * (data[i] + 0.5);
      }
      return;
    }
    // vectorized operations
    DMatrix<double> variance(const DMatrix<double>& x) const { return x.array()*(1-x.array()); }
    DMatrix<double> link(const DMatrix<double>& x) const { return ((1 - x.array()).inverse()*x.array()).log();  }
    DMatrix<double> inv_link(const DMatrix<double>& x) const { return (1+((-x).array().exp())).inverse(); }
    DMatrix<double> der_link(const DMatrix<double>& x) const { return (x.array()*(1-x.array())).inverse(); }

    // deviance
    double deviance(std::size_t x) { return x == 0 ? 2*std::log(1/(1-p_)) : 2*std::log(1/p_); };
  };

  class Poisson {
  private:
    double l_; // distribution parameter
    // a simple utility to compute the factorial of an integer number
    std::size_t factorial(std::size_t k) const {
      std::size_t result = 1;
      while(k>0) { result*=k; --k; }
      return result;
    }
  public:
    // constructor
    Poisson() = default;
    Poisson(double l) : l_(l) {};
    // density function
    double pdf(std::size_t k) const { return std::pow(l_, k)*std::exp(-l_)/factorial(k); };
    double mean() const { return l_; }    
    // preprocess data to work with this distribution
    void preprocess(DVector<double>& data) const {
      for(std::size_t i = 0; i < data.rows(); ++i){
	if(data[i] <= 0) data[i] = 1;
      }
      return;
    }
    // vectorized operations
    DMatrix<double> variance(const DMatrix<double>& x) { return x; }    
    DMatrix<double> link(const DMatrix<double>& x) const { return x.array().log(); }
    DMatrix<double> inv_link(const DMatrix<double>& x) const { return x.array().exp(); }
    DMatrix<double> der_link(const DMatrix<double>& x) const { return x.array().inverse(); }
    
    // deviance
    double deviance(std::size_t x) { return x > 0 ? x*std::log(x/l_) - (x-l_) : l_; };
  };

  class Exponential {
  private:
    double l_; // distribution parameter
  public:
    // constructor
    Exponential() = default;
    Exponential(double l) : l_(l) {};
    // density function
    double pdf(double x) const { return l_*std::exp(-l_*x); };
    double mean() const { return 1/l_; }
    void preprocess(DVector<double>& data) const { return; }
    // vectorized operations
    DMatrix<double> variance(const DMatrix<double>& x) { return x.array().pow(2); }    
    DMatrix<double> link(const DMatrix<double>& x) const { return (-x).array().inverse(); }
    DMatrix<double> inv_link(const DMatrix<double>& x) const { return (-x).array().inverse(); }
    DMatrix<double> der_link(const DMatrix<double>& x) const { return x.array().pow(2).inverse(); }
   
    // deviance
    double deviance(std::size_t x) { return 2*((x-l_)/l_ - std::log(x/l_)); };    
  };

  class Gamma {
  private:
    double k_;     // shape parameter
    double theta_; // scale parameter
  public:
    // constructor
    Gamma() = default;
    Gamma(double k, double theta) : k_(k), theta_(theta) {};
    // density function
    double pdf(double x) const { return 1/(std::tgamma(k_)*std::pow(theta_, k_))*std::pow(x, k_-1)*std::exp(-x/theta_); };
    double mean() const { return k_*theta_; }
    void preprocess(DVector<double>& data) const { return; }
    // vectorized operations
    DMatrix<double> variance(const DMatrix<double>& x) { return x.array().pow(2); }    
    DMatrix<double> link(const DMatrix<double>& x) const { return (-x).array().inverse(); }
    DMatrix<double> inv_link(const DMatrix<double>& x) const { return (-x).array().inverse(); }
    DMatrix<double> der_link(const DMatrix<double>& x) const { return x.array().pow(2).inverse(); }

    // deviance
    double deviance(std::size_t x) { return 2*((x-theta_)/theta_ - std::log(x/theta_)); };    
  };
  
}}

#endif // __DISTRIBUTION_H__

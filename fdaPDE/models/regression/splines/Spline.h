#ifndef __SPLINE_H__
#define __SPLINE_H__

#include "../../../core/utils/fields/ScalarField.h"
#include <cstddef>
using fdaPDE::core::ScalarExpr;
using fdaPDE::core::ZeroField;
#include <iomanip>

namespace fdaPDE{
namespace models{

  // Let u_0, u_1, ..., u_n n distinct knots. Call U = [u_0, u_1, ..., u_n] knot vector.
  // By Cox-DeBoor formula the i-th spline basis of order j N_{ij} is recursively defined as
  //
  // N_i0(x) = 1 if x \in [u_i, u_i+1) 0 otherwise
  // N_ij(x) = [(x-u_i)/(u_i+j - u_i)]*N_i,j-1(x) + [(u_i+j+1 - x)/(u_i+j+1 - u_i+1)]*N_i+1,j-1(x)
  
  // a spline of order R centered in knot u_i
  template <unsigned int R>
  class Spline : public ScalarExpr<Spline<R>> {
  private:
    const DVector<double>& knots_;
    std::size_t i_; // knot index where this basis is centered

    // store constants a_ = (u_i+j - u_i), b_ = (u_i+j+1 - u_i+1)
    double a_ = 0, b_ = 0;
  public:
    Spline(const DVector<double>& knots, std::size_t i) : knots_(knots), i_(i) {
      // avoid possible divisions by zero
      a_ = knots_[i_+R] - knots_[i_] != 0 ? knots_[i_+R] - knots_[i_] : 0;
      b_ = knots_[i_+R+1] - knots_[i_+1] != 0 ? knots_[i_+R+1] - knots_[i_+1] : 0;
    };
    // evaluates the spline at a given point
    inline double operator()(SVector<1> x) const {
      // exploit local support of splines
      if(x[0] < knots_[i_] || knots_[i_+R+1] < x[0]) return 0;
      // develop Cox-DeBoor recursion
      return (a_ != 0 ? (x[0] - knots_[i_])/a_*Spline<R-1>(knots_, i_)(x) : 0) +
	     (b_ != 0 ? (knots_[i_+R+1] - x[0])/b_*Spline<R-1>(knots_, i_+1)(x) : 0);
    }
    // compute derivative of spline as a ScalarExpr
    // d/dx N_ij(x) = j/(u_i+j - u_i)*N_i,j-1(x)  - j/(u_i+j+1 - u_i+1)*N_i+1,j-1(x)
    auto gradient() const {
      return R/a_*Spline<R-1>(knots_, i_) - R/b_*Spline<R-1>(knots_, i_+1);
    }
  };

  // specialization for order 0 splines (end of recursion)
  template <>
  class Spline<0> : ScalarExpr<Spline<0>> {
  private:
    const DVector<double>& knots_;
    std::size_t i_; // knot index where this basis is centered
  public:
    // constructor
    Spline(const DVector<double>& knots, std::size_t i) : knots_(knots), i_(i) {};
    // implements the indicator function over the closed interval [u_i, u_{i+1}]
    inline double operator()(SVector<1> x) const {
      return (knots_[i_] <= x[0] && x[0] < knots_[i_+1]) ? 1 : 0;
    }
    // return derivative as a zero callable field
    auto gradient() const { return ZeroField<1>(); }
  };

}}
#endif // __SPLINE_H__

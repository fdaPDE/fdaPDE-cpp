#ifndef __SPLINE_BASIS_H__
#define __SPLINE_BASIS_H__

#include "../../../core/utils/Symbols.h"
#include "Spline.h"
using fdaPDE::models::Spline;

namespace fdaPDE {
namespace models {

  // a SplineBasis of order R. This class is compliant with the functional basis concept adopted in the library.
  template <unsigned int R>
  class SplineBasis{
  private:
    DVector<double> knots_{}; // vector of knots
    std::vector<Spline<R>> basis_{};
  public:
    using const_iterator = typename std::vector<Spline<R>>::const_iterator;
    // constructor
    SplineBasis(const DVector<double>& knots) {
      //reserve space
      std::size_t n = knots.size();
      knots_.resize(n+2*R);
      // pad the knot vector to obtain a full basis for the whole knot span [u_0, u_n]
      for(std::size_t i = 0; i < n+2*R; ++i){
	if(i < R) knots_[i] = knots[0];
	else{
	  if(i < n+R) knots_[i] = knots[i-R];
	  else knots_[i] = knots[n-1];
	}
      }
      // reserve space and compute spline basis
      basis_.reserve(knots_.rows()-R-1);
      for(std::size_t k = 0; k < knots_.size()-R-1; ++k){
	// create spline centered at k-th point of knots_
	basis_.emplace_back(knots_, k);
      }
    };
    // return i-th element of the basis
    const Spline<R>& operator[](std::size_t i) const { return basis_[i]; }
    // return the number of basis elements
    int size() const { return basis_.size(); }
    // allow range-for over basis elements
    const_iterator begin() const { return basis_.cbegin(); }
    const_iterator end() const { return basis_.cend(); }
  };

}}

#endif // __SPLINE_BASIS_H__

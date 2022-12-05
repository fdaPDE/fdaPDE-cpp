#ifndef __SPACE_TIME_H__
#define __SPACE_TIME_H__

#include "../../core/utils/Symbols.h"
#include "SplineBasis.h"
using fdaPDE::models::SplineBasis;
#include "../../core/FEM/integration/IntegratorTables.h"
using fdaPDE::core::FEM::IntegratorTable;
using fdaPDE::core::FEM::GaussLegendre;

#include <type_traits>

namespace fdaPDE {
namespace models {

  // base class providing time discretization matrices. B is the type of basis used for discretizing the time dimension
  template <typename B> class TimeAssembler;

  // specialization of TimeAssembler for a spline basis of order R
  template <unsigned int R>
  class TimeAssembler<SplineBasis<R>> {
  private:
    const DVector<double>& time_;

    typedef SplineBasis<R> B;
    B basis_;

    // only support for cubic splines
    IntegratorTable<1,3,GaussLegendre> integrationTable_{};
    
    // specialized routine for the computation of 1D integrals over interval [a,b]. We rely on this routine instead of one in the
    // core/FEM/integration module because it is enought simple to be specialized here and more importantly we don't have an
    // explicit 1D mesh object for the time dimension while the integration module relies explicitly on it
    template <typename E>
    double integrate(double a, double b, const E& f) const {
      double result = 0;
      for(std::size_t iq = 0; iq < integrationTable_.num_nodes; ++iq){
	// quadrature rule is relative to the [-1,1] interval, need to perform change of variable
	// \int_{[a,b]} f(x) -> (b-a)/2 * \sum_{iq} w_{iq} * f((b-a)/2*x + (b+a)/2)
	result += f(SVector<1>(((b-a)/2)*integrationTable_.nodes[iq][0] + (b+a)/2))*integrationTable_.weights[iq];
      }
      // correct for measure of interval
      return (b-a)/2*result;
    }
    
  public:
    // constructor
    TimeAssembler(const DVector<double>& time) : time_(time), basis_(time_) {};
    
    template <typename E>
    SpMatrix<double> assemble(const E& f) const;
  };

  template <unsigned int R>
  template <typename E>
  SpMatrix<double> TimeAssembler<SplineBasis<R>>::assemble(const E& f) const {
    // compute result dimensions
    std::size_t M = basis_.size();
    // resize result matrix
    SpMatrix<double> discretizationMatrix;
    discretizationMatrix.resize(M,M);
    // prepare triplet list for sparse matrix construction
    std::vector<fdaPDE::Triplet<double>> tripletList;
    tripletList.reserve(M*M);

    // start assembly loop (exploit local support of spline basis)
    std::size_t m = time_.rows();
    for(std::size_t i = 0; i < M; ++i){
      for(std::size_t j = i; j < std::min(M-1, i+R+1); ++j){
	double phi = 0;
	for(std::size_t k = j; k < std::min(m-1, j+R); ++k){
	  phi += integrate(time_[k], time_[k+1], f);
	}
	tripletList.emplace_back(i,j, phi);
      }
    }
    // finalize construction
    discretizationMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    discretizationMatrix.makeCompressed();
    return discretizationMatrix;
  };

  // general class for the assembly of matrix \Phi = [\Phi]_{ij} = \phi_i(t_j),
  // being \phi_1, \phi_2, ... \phi_M a basis system of type B
  template <typename B> class PhiAssembler;

  // specialization of \Phi for a spline basis of order R
  template <unsigned int R> struct PhiAssembler<SplineBasis<R>> {
    typedef SplineBasis<R> B;
    static SpMatrix<double> compute(const DVector<double>& time_);
  };
  
  template <unsigned int R>
  SpMatrix<double> PhiAssembler<SplineBasis<R>>::compute(const DVector<double>& time_) {
    // create time basis
    B basis(time_);
    // define Phi matrix dimensions
    int m = time_.rows();
    int M = basis.size();
    // resize result matrix
    SpMatrix<double> Phi;
    Phi.resize(m, M);
    
    // triplet list to fill sparse matrix \Phi
    std::vector<fdaPDE::Triplet<double>> tripletList;
    tripletList.reserve(m*M);

    for(int i = 0; i < M; ++i){
      // exploit local support of splines
      for(int j = ((int)(i-B::order) > 0 ? (int)(i-B::order) : 0); j < std::min(m-1, i+1); ++j){
	tripletList.emplace_back(j, i, basis[i](SVector<1>(time_[j])));
      }
    }
    // finalize construction
    Phi.setFromTriplets(tripletList.begin(), tripletList.end());
    Phi.makeCompressed();    
    return Phi;
  }

  // convenient function template for the construction of \Phi
  template <typename B>
  SpMatrix<double> Phi(const DVector<double>& time) {
    return PhiAssembler<B>::compute(time);
  }
  
}}

#endif // __SPACE_TIME_H__

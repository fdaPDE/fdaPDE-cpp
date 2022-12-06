#ifndef __SPACE_TIME_H__
#define __SPACE_TIME_H__

#include "../../core/utils/Symbols.h"
#include "SpaceTime/Spline.h"
#include "SplineBasis.h"
#include <functional>
#include <limits>
using fdaPDE::models::SplineBasis;
#include "../../core/FEM/integration/IntegratorTables.h"
using fdaPDE::core::FEM::IntegratorTable;
using fdaPDE::core::FEM::GaussLegendre;
#include "../../core/utils/fields/ScalarField.h"
using fdaPDE::core::ScalarBinOp;
#include "../../core/utils/fields/FieldPtrs.h"
using fdaPDE::core::ScalarPtr;

#include <type_traits>

namespace fdaPDE {
namespace models {

  // base class providing time discretization matrices. B is the type of basis used for discretizing the time dimension
  // this is a specialized assembly loop for 1D problems which do not require all the machinery put in place in the FEM assembler
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

    template <typename E>
    constexpr std::size_t compute_internal_loop_limit(std::size_t i, std::size_t M) const {
      if constexpr(E::is_symmetric) return i;
      else return M;
    }
    
  public:
    // constructor
    TimeAssembler(const DVector<double>& time) : time_(time), basis_(time_) {};
    
    template <typename E>
    SpMatrix<double> assemble(const E& f) const;
  };

  // functor for the computation of the (i,j)-th element of the time mass discretization matrix \phi_i*\phi_j
  template <typename B> class TimeMass;  
  // partial specialization for SplineBasis case.
  template <unsigned int R>
  struct TimeMass<SplineBasis<R>> {
    static constexpr bool is_symmetric = true;
    auto integrate(const ScalarPtr<Spline<R>>& phi_i, const ScalarPtr<Spline<R>>& phi_j) const {
      return phi_i*phi_j;
    };
  };

  // functor for the computation of the (i,j)-th element of the time penalty matrix (\phi_i)_tt * (\phi_j)_tt
  template <typename B> class TimePenalty;
  // partial specialization for SplineBasis case.
  template <unsigned int R>
  struct TimePenalty<SplineBasis<R>> {
    static constexpr bool is_symmetric = true;
    auto integrate(ScalarPtr<Spline<R>>& phi_i, ScalarPtr<Spline<R>>& phi_j) const {
      return (phi_i->template derive<2>())*(phi_j->template derive<2>());
    };
  };
  
  template <unsigned int R>
  template <typename E>
  SpMatrix<double> TimeAssembler<SplineBasis<R>>::assemble(const E& op) const {
    // compute result dimensions
    std::size_t M = basis_.size();
    // resize result matrix
    SpMatrix<double> discretizationMatrix;
    discretizationMatrix.resize(M,M);
    // prepare triplet list for sparse matrix construction
    std::vector<fdaPDE::Triplet<double>> tripletList;
    tripletList.reserve(M*M);

    // store space for operands
    using basis_type = typename SplineBasis<R>::element_type;
    basis_type buff_phi_i(basis_.knots()), buff_phi_j(basis_.knots());
    ScalarPtr<Spline<R>> phi_i(&buff_phi_i), phi_j(&buff_phi_j);
    // develop integand here once
    auto f = op.integrate(phi_i, phi_j);
    
    // start assembly loop (exploit local support of spline basis)
    for(std::size_t i = 0; i < M; ++i){
      buff_phi_i = basis_[i];
      for(std::size_t j = 0; j <= compute_internal_loop_limit<E>(i,M); ++j){
	buff_phi_j = basis_[j]; // update buffer content
	double phi = 0;
	for(std::size_t k = j; k <= i+R; ++k){
	  phi += integrate(basis_.knots()[k], basis_.knots()[k+1], f);
	}
	tripletList.emplace_back(i,j, phi);
      }
    }
    // finalize construction
    discretizationMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    discretizationMatrix.makeCompressed();
    if constexpr(E::is_symmetric)
      return discretizationMatrix.selfadjointView<Eigen::Lower>();
    else
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

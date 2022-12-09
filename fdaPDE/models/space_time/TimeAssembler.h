#ifndef __TIME_ASSEMBLER_H__
#define __TIME_ASSEMBLER_H__

#include "SplineBasis.h"
using fdaPDE::models::SplineBasis;
#include "../../core/FEM/integration/IntegratorTables.h"
using fdaPDE::core::FEM::IntegratorTable;
using fdaPDE::core::FEM::GaussLegendre;

namespace fdaPDE {
namespace models {

  // base class providing time discretization matrices. B is the type of basis used for discretizing the time dimension
  // this is a specialized assembly loop for 1D problems which do not require all the machinery put in place in the FEM assembler
  template <typename B> class TimeAssembler;
  template <typename B> struct traits;
  
  // specialization of TimeAssembler for a spline basis of order R
  template <unsigned int R>
  class TimeAssembler<SplineBasis<R>> {
  private:
    const DVector<double>& time_; // time mesh
    // basis used for discrization of time operators
    typename traits<TimeAssembler<SplineBasis<R>>>::Basis basis_{};
    
    // quadrature rule for the computation of 1D integrals over interval [a,b].
    typename traits<TimeAssembler<SplineBasis<R>>>::Integrator integrationTable_{};    
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

    // if the operator is symmetric returns i, otherwise the inner assembly loop must run over the whole set of basis funcitons
    template <typename E>
    constexpr std::size_t compute_inner_loop_limit(std::size_t i, std::size_t M) const {
      if constexpr(E::is_symmetric) return i;
      else return M;
    }
  public:
    // constructor
    TimeAssembler(const DVector<double>& time) : time_(time), basis_(time_) {};
    // computes the discretization of E
    template <typename E>
    SpMatrix<double> assemble(const E& f) const;
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
    
    // start assembly loop (exploit local support of spline basis)
    for(std::size_t i = 0; i < M; ++i){
      for(std::size_t j = 0; j <= compute_inner_loop_limit<E>(i,M); ++j){
	// develop integrand field
	auto f = op.integrate(basis_[i], basis_[j]);
	// perform integration of f over [knots[j], knots[i+R+1]]
        double value = 0;
	for(std::size_t k = j; k <= i+R; ++k){
	  value += integrate(basis_.knots()[k], basis_.knots()[k+1], f);
	}
	tripletList.emplace_back(i,j, value); // store computed integral
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

  // currently only cubic splines are supported
  template <>
  struct traits<TimeAssembler<SplineBasis<3>>> {
    typedef SplineBasis<3> Basis;
    typedef IntegratorTable<1,3,GaussLegendre> Integrator;
  };
  
  // functor for the computation of the (i,j)-th element of the time mass discretization matrix \phi_i*\phi_j
  template <typename B> class TimeMass;  
  // partial specialization for SplineBasis case.
  template <unsigned int R>
  struct TimeMass<SplineBasis<R>> {
    static constexpr bool is_symmetric = true;
    // provide \phi_i * \phi_j
    auto integrate(const Spline<R>& phi_i, const Spline<R>& phi_j) const {
      return phi_i * phi_j;
    };
  };

  // functor for the computation of the (i,j)-th element of the time penalty matrix (\phi_i)_tt * (\phi_j)_tt
  template <typename B> class TimePenalty;
  // partial specialization for SplineBasis case.
  template <unsigned int R>
  struct TimePenalty<SplineBasis<R>> {
    static constexpr bool is_symmetric = true;
    // provide (\phi_i)_tt * (\phi_j)_tt
    auto integrate(const Spline<R>& phi_i, const Spline<R>& phi_j) const {
      return phi_i.template derive<2>() * phi_j.template derive<2>();
    };
  };
  
}}

#endif // __TIME_ASSEMBLER_H__

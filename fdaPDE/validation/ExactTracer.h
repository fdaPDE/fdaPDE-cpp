#ifndef __EXACT_TRACER_H__
#define __EXACT_TRACER_H__

#include "../core/utils/Symbols.h"
#include "../core/utils/fields/ScalarField.h"
using fdaPDE::core::ScalarField;

// Evaluates exactly the trace of matrix S = \Psi*T^{-1}*\Psi^T*Q.
// This means to explicitly compute matrix S and thus extract its trace.
class ExactTracer {
private:
  // compute matrix S = \Psi*T^{-1}*\Psi^T*Q
  template <typename M>
  std::unique_ptr<DMatrix<double>> S(const M& model) const{
    // compute \Psi^T*Q
    // an optimization can be put in place if \Psi is a permutation matrix (locations are a subset of nodes)
    DMatrix<double> E{};
    if(model.Q() != nullptr)
      E = model.Psi()->transpose()*(*model.Q());
    else
      E = model.Psi()->transpose();
    // factorize matrix T
    Eigen::PartialPivLU<DMatrix<double>> invT = model.T()->partialPivLu();

    // compute V = invT*E = T^{-1}*\Psi^T*Q
    DMatrix<double> V = invT.solve(E);

    // S = \Psi*T^{-1}*\Psi^T*Q
    // if locations are equal to mesh nodes then premultiply by \Psi reduces to the permutation of the rows of T^{-1}*\Psi^T*Q = V
    std::unique_ptr<DMatrix<double>> S = std::make_unique<DMatrix<double>>(model.n(), model.n());
    for(std::size_t k = 0; k < model.Psi()->outerSize(); ++k){
      for (SpMatrix<double>::InnerIterator it(*model.Psi(),k); it; ++it){
	S->row(it.row()) = V.row(it.col());
      }
    }
    return S;
  };
  
public:
  // constructor
  ExactTracer() = default;
  // evaluate trace of S in an exact way
  template <typename M>
  double compute(const M& model) const {
    return S(model)->trace();
  }
};

#endif // __EXACT_TRACER_H__

#ifndef __GCV_H__
#define __GCV_H__

#include <cstddef>
#include <memory>
#include <vector>
#include <functional>
#include "../core/utils/Symbols.h"
#include <Eigen/LU>

// M: statistical model of which the GCV should be optimized
template <typename M>
class GCV {
private:
  // 3 std::function each modeling the expression of the GCV, its first derivative and its second derivative
  std::function<double(SVector<1>)> gcv{};

  M& model_;
  // utilities to compute S and its trace in an efficient way...
  double trS_;
  
  // compute matrix S = \Psi*T^{-1}*\Psi^T*Q
  std::unique_ptr<DMatrix<double>> S(const SVector<1>& lambda) const{
    // compute \Psi^T*Q
    // an optimization can be put in place if \Psi is a permutation matrix (locations are a subset of nodes)
    DMatrix<double> E{};
    if(model_.Q() != nullptr)
      E = model_.Psi()->transpose()*(*model_.Q());
    else
      E = model_.Psi()->transpose();
    // factorize matrix T
    Eigen::PartialPivLU<DMatrix<double>> invT = model_.T(lambda[0])->partialPivLu();
    // compute V = invT*E = T^{-1}*\Psi^T*Q
    DMatrix<double> V = invT.solve(E);
    // S = \Psi*T^{-1}*\Psi^T*Q
    // if locations are equal to mesh nodes then premultiply by \Psi reduces to the permutation of the rows of T^{-1}*\Psi^T*Q = V
    std::unique_ptr<DMatrix<double>> S = std::make_unique<DMatrix<double>>(model_.n(), model_.n());
    for(std::size_t k = 0; k < model_.Psi()->outerSize(); ++k){
      for (SpMatrix<double>::InnerIterator it(*model_.Psi(),k); it; ++it){
	S->row(it.row()) = V.row(it.col());
      }
    }
    return S;
  };
  
public:
  GCV() = default;
  GCV(M& model) : model_(model) {
    gcv = [*this](SVector<1> lambda) -> double {
      // fit the model with this value of lambda
      model_.setLambda(lambda[0]);
      model_.smooth();
      // compute trace of matrix S given current lambda
      double trS = S(lambda)->trace();

      // GCV(\lambda) = n/(n-(q+trS(\lambda))^2)*norm(z - \hat_z(\lambda))^2
      double q = model_.q();      // number of covariates
      std::size_t n = model_.n(); // number of locations
      return (n/std::pow(n - (q + trS), 2))*(model_.z() - model_.fitted()).squaredNorm();
    };
  };
  
  // optimizes the GCV in an exact way
  double exact();
  // optimizes the GCV in an approximate way using finite differences to approximate first and second derivative
  template <typename O, typename... Args>
  SVector<1> FDApprox(O& optimizer, Args&... args) {
    // wrap gcv in a ScalarField object.
    // This forces optimizers to employ a finite difference approximation of gcv derivatives when it is required
    ScalarField<1> obj(gcv);
    // optimize gcv field
    optimizer.findMinimum(obj, args...);
    // recover found solution
    SVector<1> x = optimizer.getSolution();
    return x;
  };
  // optimizes the GCV in an approximate way using a monte carlo approach
  double MCApprox();
};

#endif // __GCV_H__

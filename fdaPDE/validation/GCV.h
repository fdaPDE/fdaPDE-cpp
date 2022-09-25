#ifndef __GCV_H__
#define __GCV_H__

#include <cstddef>
#include <memory>
#include <random>
#include <type_traits>
#include <vector>
#include <functional>
#include "../core/utils/Symbols.h"
#include "../core/utils/fields/ScalarField.h"
using fdaPDE::core::ScalarField;
#include "../core/NLA/SMW.h"
using fdaPDE::core::NLA::SMW;

class ExactTraceEvaluator {
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
  ExactTraceEvaluator() = default;
  // evaluate trace of S in an exact way
  template <typename M>
  double trace(const M& model) const { return S(model)->trace(); }
};

class StochasticTraceEvaluator {
private:
  std::size_t r_; // number of realizations

  template <typename M>
  std::shared_ptr<DMatrix<double>> Us(const M& model) const {
    // define bernulli distribution with parameter p = 0.5 to simulate Rademacher distribution
    std::random_device rng;
    std::bernoulli_distribution Be(0.5);
    // preallocate memory for matrix Us
    std::shared_ptr<DMatrix<double>> Us = std::make_shared<DMatrix<double>>(model.n(), r_);
    // fill matrix
    for(std::size_t i = 0; i < model.n(); ++i){
      for(std::size_t j = 0; j < r_; ++j){
	if(Be(rng)) Us->coeffRef(i,j) =  1.0;
	else        Us->coeffRef(i,j) = -1.0;
      }
    }
    return Us;
  }
  
public:
  // constructor
  StochasticTraceEvaluator(std::size_t r) : r_(r) { }

  // evaluate trace of S exploiting a monte carlo approximation
  template<typename M>
  double trace(const M& model) const {
    std::size_t n = model.n(); // number of observations
    // apply SMW decomposition
    std::size_t q_ = model.W()->cols(); // number of covariates
    // definition of matrices U and V
    DMatrix<double> U = DMatrix<double>::Zero(model.A()->rows(), q_);
    U.block(0,0, model.A()->rows()/2, q_) = model.Psi()->transpose()*(*model.W());
  
    DMatrix<double> V = DMatrix<double>::Zero(q_, model.A()->rows());
    V.block(0,0, q_, model.A()->rows()/2) = model.W()->transpose()*(*model.Psi());

    // Define system solver. Use SMW solver from NLA module
    SMW<> solver{};
    solver.compute(*model.A());
    
    // define matrix Bs
    DMatrix<double> Bs = DMatrix<double>::Zero(2*model.obs(), r_);
    DMatrix<double> UU = *Us(model);
    Bs.topRows(model.obs()) = model.Psi()->transpose() * lmbQ(model, UU);
    
    // solve system Mx = b
    DMatrix<double> solution;
    solution = solver.solve(U, *(model.WTW()), V, Bs);
    // compute approximated Tr[S] using monte carlo mean
    DMatrix<double> Y = UU.transpose() * (*model.Psi());
    double MCmean = 0;
    for(std::size_t i = 0; i < r_; ++i){
      MCmean += Y.row(i).dot(solution.col(i).head(n));
    }

    return MCmean/r_;
  }
};

// M: statistical model of which the GCV should be optimized
// trace evaluator, defaulted to stochastic for efficiency reasons. Use tag dispathcing to select the trace computation policy
template <typename M, typename TraceEvaluator = StochasticTraceEvaluator>
class GCV {
private:
  // 3 std::function each modeling the expression of the GCV, its first derivative and its second derivative
  std::function<double(SVector<1>)> gcv;
  
  M& model_;
  // utilities to compute S and its trace in an efficient way...
  TraceEvaluator t;

  // initialize gcv functors
  void init() {
    gcv = [*this](SVector<1> lambda) -> double {
      // fit the model with this value of lambda
      model_.setLambda(lambda[0]);
      model_.smooth();
      // compute trace of matrix S given current lambda
      double trS = t.trace(model_);

      // GCV(\lambda) = n/(n-(q+trS(\lambda))^2)*norm(z - \hat_z(\lambda))^2
      double q = model_.q();      // number of covariates
      std::size_t n = model_.n(); // number of locations

      std::cout << (n/std::pow(n - (q + trS), 2))*(model_.z() - model_.fitted()).squaredNorm() << std::endl;
      
      return (n/std::pow(n - (q + trS), 2))*(model_.z() - model_.fitted()).squaredNorm();
    };
  }
  
public:
  // SFINAE selection of constructor depending on trace evaluation strategy
  template <typename U = TraceEvaluator, // fake type to enable substitution
      typename std::enable_if<
	!std::is_same<U, StochasticTraceEvaluator>::value,
	int>::type = 0> 
  GCV(M& model)
    : model_(model) { init(); };

  template <typename U = TraceEvaluator,
      typename std::enable_if<
	std::is_same<U, StochasticTraceEvaluator>::value,
	int>::type = 0>
  GCV(M& model, std::size_t r)
    : model_(model), t(r) { init(); };
  
  // optimizes the GCV in an exact way
  double exact();
  // optimizes the GCV in an approximate way using finite differences to approximate first and second derivative
  template <typename O, typename... Args>
  SVector<1> FDApprox(O& optimizer, double gcvFDStep, Args&... args) {
    // wrap gcv in a ScalarField object.
    // This forces optimizers to employ a finite difference approximation of gcv derivatives when it is required
    ScalarField<1> obj(gcv);
    obj.setStep(gcvFDStep);
    // optimize gcv field
    optimizer.findMinimum(obj, args...);
    // recover found solution
    SVector<1> x = optimizer.getSolution();
    return x;
  };
};

#endif // __GCV_H__

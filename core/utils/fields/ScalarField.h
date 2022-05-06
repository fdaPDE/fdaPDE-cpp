#ifndef __SCALAR_FIELD_H__
#define __SCALAR_FIELD_H__

#include <functional>
#include "../Symbols.h"
#include "ScalarFieldExpressions.h"

namespace fdaPDE{
namespace core{

  // a template class for handling general scalar fields. A field wrapped
  // by this template doesn't guarantee any regularity condition.
  // N is the dimension of the domain.
  template <unsigned int N>
  class ScalarField : public FieldExpr<ScalarField<N>>{
  protected:
    // the function this class wraps
    std::function<double(SVector<N>)> f;
  
  public:
    // constructor
    ScalarField(std::function<double(SVector<N>)> f_) : f(f_) {};

    // preserve std::function syntax for evaluating a function at point
    inline double operator()(const SVector<N>& x) const { return f(x); };
    // method for explicitly evaluate the function at point
    inline double evaluateAtPoint(const SVector<N>& x) const { return f(x); };
  
    // gradient approximation of function around a given point using central differences.
    SVector<N> getGradientApprox (const SVector<N>& x, double step) const;

    // hessian approximation of function around a given point using central differences.
    SMatrix<N> getHessianApprox(const SVector<N>& x, double step) const;

    // the following methods are set as virtual since this class
    // does not force the definition of an analytical expression for both
    // gradient and hessian. They will be implemented by child classes.

    // a ScalarField implementing this method has an analytical expression for its gradient
    virtual std::function<SVector<N>(SVector<N>)> derive() const { return nullptr; };
    // a ScalarField implementing this method has an analytical expression for its hessian
    virtual std::function<SMatrix<N>(SVector<N>)> deriveTwice() const { return nullptr; };
  };

  // the following classes can be used to force particular regularity conditions
  // on the scalar field which might be required for some numerical methods.
  // i.e. Using a DifferentiableScalarfield as argument for a function forces to
  // define an analytical expression for the derivative of the field object passed as argument
  template <unsigned int N>
  class DifferentiableScalarField : public ScalarField<N> {
  protected:
    // gradient vector of scalar field f
    std::function<SVector<N>(SVector<N>)> df;
  
  public:
    // constructor
    DifferentiableScalarField(std::function<double(SVector<N>)> f_,      // base function
			      std::function<SVector<N>(SVector<N>)> df_  // gradient function
			      ) : ScalarField<N>(f_), df(df_) {};
  
    std::function<SVector<N>(SVector<N>)> derive() const override { return df; };
  };

  template <unsigned int N>
  class TwiceDifferentiableScalarField : public DifferentiableScalarField<N> {
  protected:
    // hessian matrix of scalar field f
    std::function<SMatrix<N>(SVector<N>)> ddf;
  
  public:
    // constructor
    TwiceDifferentiableScalarField(std::function<double(SVector<N>)> f_,      // base function
				   std::function<SVector<N>(SVector<N>)> df_, // gradient function
				   std::function<SMatrix<N>(SVector<N>)> ddf_ // hessian function
				   ) : DifferentiableScalarField<N>(f_, df_), ddf(ddf_) {};

    std::function<SMatrix<N>(SVector<N>)> deriveTwice() const override { return ddf; }
  };

#include "ScalarField.tpp"
}}
  
#endif // __SCALAR_FIELD_H__

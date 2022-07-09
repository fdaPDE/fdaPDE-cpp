#ifndef __SCALAR_FIELD_H__
#define __SCALAR_FIELD_H__

#include <functional>
#include <initializer_list>
#include "../Symbols.h"
#include "ScalarFieldExpressions.h"
#include "VectorField.h"

namespace fdaPDE{
namespace core{

  // a template class for handling general scalar fields. A field wrapped by this template doesn't guarantee any regularity condition.
  // N is the domain dimension. The wrapped field must be encoded in a lambda expression receiving an SVector<N> in input and returning a double
  // extend FieldExpr to allow for expression templates
  template <int N>
  class ScalarField : public FieldExpr<ScalarField<N>> {
  private:
    // approximation of first and second derivative using central differences
    double approxFirstDerivative (const SVector<N>& x, std::size_t i, double step) const;
    double approxSecondDerivative(const SVector<N>& x, std::size_t i, std::size_t j, double step) const;
    
  protected:
    // the function this class wraps
    std::function<double(SVector<N>)> f_;

  public:
    // constructor
    ScalarField(const std::function<double(SVector<N>)>& f) : f_(f) {};
    
    // preserve std::function syntax for evaluating a function at point, required for expression templates
    double operator()(const SVector<N>& x) const { return f_(x); };

    // approximation of gradient vector and hessian matrix without construction of a VectorField object
    SVector<N> approxGradient (const SVector<N>& x, double step) const;
    SMatrix<N> approxHessian  (const SVector<N>& x, double step) const;
    
    VectorField<N> derive(double step) const;
    virtual VectorField<N> derive() const; // returns exact gradient in DifferentiableScalarField
    
    std::function<SMatrix<N>(SVector<N>)> deriveTwice(double step) const;
    virtual std::function<SMatrix<N>(SVector<N>)> deriveTwice() const; // returns exact hessian in TwiceDifferentiableScalarField

    // discretize the scalar field over the given mesh with a discontinuous function which is constant at every mesh element.
    // The result of the discretization is a numeric vector where element i corresponds to the evaluation of the scalar field
    // at the midpoint of element i
    // template <unsigned int L, unsigned int K>
    // Eigen::Matrix<double, Eigen::Dynamic, 1> discretize(const Mesh<L, K>& mesh) const;
  };

  // the following classes can be used to force particular regularity conditions on the field which might be required for some numerical methods.
  // i.e. Using a DifferentiableScalarfield as argument for a function forces to define an analytical expression for the gradient vector
  template <int N>
  class DifferentiableScalarField : public ScalarField<N> {
  protected:
    // gradient vector of scalar field f
    VectorField<N> df_;
  
  public:
    // constructor
    DifferentiableScalarField(const std::function<double(SVector<N>)>& f,                         // base function
			      const std::initializer_list<std::function<double(SVector<N>)>>& df  // gradient function
			      ) : ScalarField<N>(f), df_(df) {};

    // allow the construction of the gradient VectorField from a single lambda
    DifferentiableScalarField(const std::function<double(SVector<N>)>& f,      // base function
			      const std::function<SVector<N>(SVector<N>)>& df  // gradient function
			      ) : ScalarField<N>(f), df_(VectorField<N>(df)) {};
    
    VectorField<N> derive() const override { return df_; };
  };

  template <int N>
  class TwiceDifferentiableScalarField : public DifferentiableScalarField<N> {
  protected:
    // hessian matrix of scalar field f
    std::function<SMatrix<N>(SVector<N>)> ddf_;
  
  public:
    // constructor
    TwiceDifferentiableScalarField(const std::function<double(SVector<N>)>& f,                         // base function
				   const std::initializer_list<std::function<double(SVector<N>)>>& df, // gradient function
				   const std::function<SMatrix<N>(SVector<N>)>& ddf                    // hessian function
				   ) : DifferentiableScalarField<N>(f, df), ddf_(ddf) {};

    // allow the construction of the gradient VectorField from a single lambda
    TwiceDifferentiableScalarField(const std::function<double(SVector<N>)>& f,       // base function
			           const std::function<SVector<N>(SVector<N>)>& df,  // gradient function
				   const std::function<SMatrix<N>(SVector<N>)>& ddf  // hessian function
				   ) : DifferentiableScalarField<N>(f, df), ddf_(ddf) {};
    
    std::function<SMatrix<N>(SVector<N>)> deriveTwice() const override { return ddf_; }
  };
  
#include "ScalarField.tpp"
}}
  
#endif // __SCALAR_FIELD_H__

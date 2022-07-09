# ScalarField

> :fontawesome-solid-file-code: /src/core/utils/ScalarField.h &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; :fontawesome-solid-sitemap: Extends: FieldExpr  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; :fontawesome-solid-gears: /test/core/utils/ScalarFieldTest.cpp

A template class for handling scalar fields \( f : \mathbb{R}^N \rightarrow \mathbb{R} \)

``` c++
template <int N> class ScalarField : public FieldExpr<ScalarField<N>>;
```

!!! note
	 A field wrapped by this class doesn't guarantee any regularity condition. You can wrap any scalar field which is just **evaluable** at any point. The definition of a formal derivative is not required. See [DifferentiableScalarField](DifferentiableScalarField.md) or [TwiceDifferentiableScalarField](TwiceDifferentiableScalarField.md) for inherited classes which force specific regularity conditions.

ScalarField allows you to wrap any C++ lambda encoding the analytical expression of your field and let you access to numerical approximations of both gradient and hessian as well as compose new scalar fields using a consistent arithmetic between field objects. You can even integrate a ScalarField over any triangulated domain. See [Mesh](../MESH/index.md) for details. 

ScalarField constitutes one of the most fundamental types in the fdaPDE architecture.


``` c++ linenums="1"
// define the field analytical expression: (e^(2x+y))/x
std::function<double(SVector<2>)> g = [](SVector<2> x) -> double { 
	return std::exp(2*x[0]+x[1])/x[0];
};

// to create a ScalarField object just wrap g as follow
ScalarField<2> field(g);

// or you can even omit the template argument thanks to C++17 template argument deduction
ScalarField field(g);

// scalar fields preserve the callable interface
std::cout << "evaluation of field at point" << std::endl;
SVector<2> p(4,1);
std::cout << field(p) << std::endl;
  
// you can access to gradient numerical approximation either directly
std::cout << "approximation of gradient at point" << std::endl;
std::cout << field.approxGradient(p, 0.01) << std::endl;

// or by building a VectorField callable object and evaluating it at p
VectorField<2> grad = field.derive();
std::cout << grad(p) << std::endl;

// you can even get an approximation of the hessian matrix at point
std::cout << "approximation of hessian at point" << std::endl;
std::cout << field.deriveTwice()(p) << std::endl;

// ScalarField can wrap also discontinuous functions: I_{x > 0, y > 0}
std::function<double(SVector<2>)> h = [](SVector<2> x) -> double { 
	if(x[0] > 0 && x[1] > 0)
		return 1;
	return 0;
};
ScalarField discontinuous_field(h);

// and allow you to combine different fields using standard arithmetic operators
auto newField = field/5 + 2*discontinuous_field - 4; 
```

!!! warning 
	The result of a field expression is not a ScalarField object. In particular a field expression is a callable object only and doesn't expose the same interface of ScalarField. See [Wrapping the result of a field expression in a ScalarField object](#wrapping-the-result-of-a-field-expression-in-a-scalarfield-object) for details on how obtain a ScalarField from a field expression.

## Field Expressions
Field expressions are supported by lazy evaluation via C++ expression templatates. The main goal is to let you in the position to write complex callables starting from simple ones in a formal way while preserving good performances. 

<table>
  <tr>
    <td style="text-align: center"> Supported arithmetic operators </td>
    <td style="text-align: center"> code example </td>
    <td style="text-align: center"> formal notation </td>
  </tr>
  <tr>
	<td style="text-align: center"> <tt> operator+ </tt> </td>
	<td style="text-align: center"> <tt> f1 + f2 </tt> </td>
	<td style="text-align: center"> \( f_1 + f_2 \) </td>
  </tr>
  <tr>
	<td style="text-align: center"> <tt> operator- </tt> </td>
	<td style="text-align: center"> <tt> f1 - f2 </tt> </td>
	<td style="text-align: center"> \( f_1 - f_2 \) </td>
  </tr>
  <tr>
  	<td style="text-align: center"> <tt> operator* </tt> </td>
	<td style="text-align: center"> <tt> f1 * f2 </tt> </td>
	<td style="text-align: center"> \( f_1 \cdot f_2 \) </td>
  </tr>
    <tr>
	<td style="text-align: center"> <tt> operator/ </tt> </td>
	<td style="text-align: center"> <tt> f1 / f2 </tt> </td>
	<td style="text-align: center"> \( f_1 / f_2 \) </td>
  </tr>
</table>

The above set of operators is available also between a ScalarField object and any integral or floating-point type.

``` c++ linenums="1"
// define the expression you want to wrap
std::function<double(SVector<2>)> ff1 = [](SVector<2> x) -> double { 
	return x.norm();
};

// wrap f1 in a ScalarField object
ScalarField<2> f1(ff1);

std::function<double(SVector<2>)> ff2 = [](SVector<2> x) -> double { 
	SMatrix<2> M;
	M << 1, 0.5,
	     0.5, 1;
	
	return (M*x).norm();
};
ScalarField f2(ff2);

// define field: norm(x) + norm(M*x)
auto f = f1 + f2;

// you can integrate this quantity over some domain
Mesh<2,2> domain(... domain data ...);
Integrator<2, 6> integrator;

std::cout << "integral of field over domain" << std::endl;
	std::cout << integrator.integrate(f, domain) << std::endl;
```

For more informations about mesh handling and integration see [Mesh](../MESH/index.md) and [Integrator]()

!!! info 
	An expression of kind \(f_1 + f_2\) in the previous example has internal type `FieldBinOp<ScalarField<2>, ScalarField<2>, plus<>>`. For this reason when using field expressions use always the `auto` keyword to let the compiler infer the correct type.

#### Wrapping the result of a field expression in a ScalarField object

Even if the result of a field expression is not a ScalarField there is an usefull workaround to wrap this one in a ScalarField object. The main point is that a field expression looses informations about the domain dimensions where the resulting field should be defined. Hence by providing this information, i.e. by wrapping the field expression in a C++ lambda function, you can obtain a valid ScalarField object.

``` c++ linenums="1"
// define f1 and f2 2D fields...
auto f = f1 + f2/f1;

// wrap f in a lambda expression
auto fWrapper = [f](SVector<2> x) -> double{ return f(x); };

// define a ScalarField object
ScalarField f_field(fWrapper);

// can now access to gradient and hessian approximations of the composed field f
```

## Methods

``` c++
ScalarField(const std::function<double(SVector<N>)>& f_);

```

!!! quote ""
	Constructor.

    | Args      | Description                                                 |
    |------------------------------------------------------------------------------------------|-------------------------------------------------------------|
    | `std::function<double(SVector<N>)> f_` | An `std::function` implementing the analytical expression of the field \( f \) taking an N dimensional point in input and returning a double in output                                                 |
	

``` c++
double operator()(const SVector<N>& x) const;

```

!!! quote ""
	Returns the evaluation of the field at point x. Preserves the `std::function` syntax.

    | Args                                                                                          | Description                                                 |
    |------------------------------------------------------------------------------------------|-------------------------------------------------------------|
    | `SVector<N>& x`  | The point in \( \mathbb{R}^N \) where to evaluate the field | 

``` c++
SVector<N> approxGradient (const SVector<N>& x, double step) const;
```
!!! quote ""
	Returns the numerical approximation of the gradient of \( f \) at point x. The partial derivatives of the field are approximated via central differences according to the following expression ( \( e_i \) denoting the normal unit vector along direction \( i \))
	$$ \frac{\partial f}{\partial x_i} \approx \frac{f(x + he_i) - f(x - he_i)}{2h} $$

    | Args                                                                                          | Description                                                 |
    |------------------------------------------------------------------------------------------|-------------------------------------------------------------|
    | `const SVector<N>& x`  | The point in \( \mathbb{R}^N \) where to approximate the gradient of the field |
	| `double step` | The value of \( h \) in the central difference formula |
	
	!!! warning 
		In principle the approximation of the partial derivatives require the field \( f \) to be just **evaluable** at point x. Anyway if the field is not differentiable at point x the approximation might cause problems.
		
``` c++
SMatrix<N> approxHessian(const SVector<N>& x, double step) const
```
!!! quote ""
	Returns the numerical approximation of the hessian matrix of \( f \) at point x. The partial derivatives of the field are approximated via central differences according to the following expressions ( \( e_i \) denoting the normal unit vector along direction \( i \))
	$$ \begin{aligned} 
	&\frac{\partial^2 f}{\partial x_i \partial x_i} &&\approx \frac{-f(x + 2he_i) + 16f(x + he_i) - 30f(x) + 16f(x - he_i) - f(x - 2he_i)}{12h^2} \\
	  &\frac{\partial^2 f}{\partial x_i \partial x_j} &&\approx \frac{f(x + he_i + he_j) - f(x + he_i - he_j) - 16f(x - he_i + he_j) + f(x - he_i - he_j)}{4h^2} 
	  \end{aligned}$$
	

    | Args                                                                                          | Description                                                 |
    |------------------------------------------------------------------------------------------|-------------------------------------------------------------|
    | `const SVector<N>& x`  | The point in \( \mathbb{R}^N \) where to approximate the hessian of the field |
	| `double step` | The value of \( h \) in the central difference formula |
	
	!!! warning 
		In principle the approximation of the partial derivatives require the field \( f \) to be just **evaluable** at point x. Anyway if the field is not differentiable at point x the approximation might cause problems.

``` c++
VectorField<N> derive(double step) const;
```
!!! quote ""
	Returns a [VectorField](VectorField.md) representing an approximation of \( \nabla f \). The obtained approximation along a given direction is obtained using the central difference formula as in `approxGradient()`.
	

    | Args                                                                                          | Description                                                 |
    |------------------------------------------------------------------------------------------|-------------------------------------------------------------|
	| `double step` | The value of \( h \) in the central difference formula |
	
	!!! note 
		A VectorField is a callable object. Evaluate the result of `gradient()` in a point \(x \in \mathbb{R}^n\) to obtain the gradient approximation.
		
		``` c++ linenums="1"
		// the following calls give the same result
		
		// evaluate the gradient of f in x once and store the result in grad1
		SVector<2> grad1 = f.approxGradient(x, 0.001);
		
		// produce a callable object
		VectorField<2> callableGradient = f.gradient(0.001);
		// evaluate the obtained callable in x
		SVector<2> grad2 = callableGradient(x);
		
		// you get grad1 == grad2
		```

``` c++
virtual VectorField<N> derive() const;
```
!!! quote ""
	Returns a [VectorField](VectorField.md) representing an approximation of \( \nabla f \). The obtained approximation along a given direction is obtained using the central difference formula as in `approxGradient()` with a fixed, not-controllable, value for the step size \( h \) equal to 0.001.
	
	!!! note 
	    This method is overridden in more specialized classes such [DifferentiableScalarField](DifferentiableScalarField.md) and [TwiceDifferentiableScalarField](TwiceDifferentiableScalarField.md) where the exact gradient is evaluated.
	
``` c++
std::function<SMatrix<N>(SVector<N>)> deriveTwice(double step) const;
```

!!! quote ""
	Returns a callable object representing an approximation of the hessian matrix \(H_f\). Approximation of second derivatives is obtained using the central difference formula as in `approxHessian()`.

    | Args                                                                                          | Description                                                 |
    |------------------------------------------------------------------------------------------|-------------------------------------------------------------|
	| `double step` | The value of \( h \) in the central difference formula |

``` c++
virtual std::function<SMatrix<N>(SVector<N>)> deriveTwice() const;
```

!!! quote ""
	Returns a callable object representing an approximation of the hessian matrix \(H_f\). Approximation of second derivatives is obtained using the central difference formula as in `approxHessian()` with a fixed, not-controllable, value for the step size \( h \) equal to 0.001.


	!!! note 
	    This method is overridden in more specialized classes like [TwiceDifferentiableScalarField](TwiceDifferentiableScalarField.md) where the exact hessian is evaluated.

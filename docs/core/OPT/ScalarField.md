# ScalarField

> :fontawesome-solid-file-code: core/OPT/ScalarField.h &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; :fontawesome-solid-sitemap: Extends: -

A template class for handling scalar fields \( f : \mathbb{R}^N \rightarrow \mathbb{R} \)

``` c++
template <unsigned int N> class ScalarField { ... };
```

!!! note
	 A field wrapped by this template doesn't guarantee any regularity condition. You can wrap any scalar field which is just **evaluable** at any point. The definition of a formal derivative is not required. See [DifferentiableScalarField](DifferentiableScalarField.md) or [TwiceDifferentiableScalarField](TwiceDifferentiableScalarField.md) for inherited classes which assume specific regularity conditions.

## Methods

``` c++
ScalarField(std::function<double(SVector<N>)> f_)

```

!!! quote ""
	Constructor.

    | Args      | Description                                                 |
    |------------------------------------------------------------------------------------------|-------------------------------------------------------------|
    | `std::function<double(SVector<N>)> f_` | An `std::function` implementing the analytical expression of the field \( f \) taking an N dimensional point in input and returning a double in output                                                 |

``` c++
inline double evaluateAtPoint(const Point<N>& x)

```

!!! quote ""
	Returns the evaluation of the field at point x.

    | Args                                                                                          | Description                                                 |
    |------------------------------------------------------------------------------------------|-------------------------------------------------------------|
    | `const SVector<N>& x`  | The point in \( \mathbb{R}^N \) where to evaluate the field | 


``` c++
inline double operator()(SVector<N>& x)

```

!!! quote ""
	Returns the evaluation of the field at point x. Preserves the `std::function` syntax.

    | Args                                                                                          | Description                                                 |
    |------------------------------------------------------------------------------------------|-------------------------------------------------------------|
    | `SVector<N>& x`  | The point in \( \mathbb{R}^N \) where to evaluate the field | 

``` c++
SVector<N> getGradientApprox (const SVector<N>& x, double step) const
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
SMatrix<N> getHessianApprox(const SVector<N>& x, double step) const
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

## Examples

``` c++ linenums="1"
// define the field analytical expression: 2*x^2 - 2*y^2*x
std::function<double(SVector<2>)> g = [](SVector<2> x) -> double { 
	return 2*std::pow(x[0],2) - 2*std::pow(x[1],2)*x[0]; 
	};

// create ScalarField
ScalarField<2> field(g);

std::cout << "evaluation of field at point" << std::endl;
std::cout << fun.evaluateAtPoint(SVector<2>(4,1)) << std::endl;
  
SVector<2> grad = fun.getGradientApprox(SVector<2>(2,1), 0.001);
    
std::cout << "approximation of gradient at point" << std::endl;
std::cout << grad << std::endl;

SMatrix<2> hessian = fun.getHessianApprox(SVector<2>(2,1), 0.001);

std::cout << "approximation of hessian at point" << std::endl;
std::cout << hessian << std::endl;
```

# DifferentiableScalarField

> :fontawesome-solid-file-code: core/OPT/ScalarField.h &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; :fontawesome-solid-sitemap: Extends: [ScalarField](ScalarField.md)

Template class used to represent a scalar field \( f : \mathbb{R}^N \rightarrow \mathbb{R} \) whose gradient function \( \nabla f : \mathbb{R}^N \rightarrow \mathbb{R}^N \) is known analitically at every point.

``` c++
template <unsigned int N> class DifferentiableScalarField : public ScalarField<N> { ... };
```

## Methods

``` c++
DifferentiableScalarField(std::function<double(SVector<N>)> f_, std::function<SVector<N>(SVector<N>)> df_)
```

!!! quote ""
	Constructor

    | Args      | Description                                                 |
    |------------------------------------------------------------------------------------------|-------------------------------------------------------------|
    | `std::function<double(SVector<N>)> f_` | An `std::function` implementing the analytical expression of the field \( f \) taking an N dimensional point in input and returning a double in output                                                 |
	| `std::function<SVector<N>(SVector<N>)> df_` | An `std::function` implementing the analytical expression of the vector field \( \nabla f \) taking an N dimensional point in input and returning an N dimensional point in output representing the gradient expression                                                 |

``` c++
std::function<SVector<N>(SVector<N>)> derive() const override;
```

!!! quote ""
	Returns the analytical expression of the field's gradient \( \nabla f : \mathbb{R}^N \rightarrow \mathbb{R}^N \) as a Callable object.

## Examples

Define a scalar field \( f(x,y) = 2x^2 - 2y^2x \) having exact gradient equal to \( \nabla f = [4x - 2y^2, -4yx]^T \)

``` c++ linenums="1"
// define the field analytical expression: 2*x^2 - 2*y^2*x
std::function<double(SVector<2>)> g = [](SVector<2> x) -> double { 
	return 2*std::pow(x[0],2) - 2*std::pow(x[1],2)*x[0]; 
	};

// define analytical expression of gradient field
std::function<SVector<2>(SVector<2>)> dg = [](SVector<2> x) -> SVector<2> { 
	return SVector<2>({4*x[0] - 2*std::pow(x[1],2), 
                       -4*x[1]*x[0]}); 
	};

// define differentiable field
DifferentiableScalarField<2> field(g, dg);

std::cout << "evaluation of field at point" << std::endl;
std::cout << field.evaluateAtPoint(SVector<2>(4,1)) << std::endl;

// get approximation of gradient at point
SVector<2> grad = field.getGradientApprox(SVector<2>(2,1), 0.001);
    
std::cout << "approximation of gradient at point" << std::endl;
std::cout << grad << std::endl;

// evaluate exact gradient at point
SVector<2> exactGrad = field.derive()(SVector<2>(2,1));

std::cout << "exact gradient at point" << std::endl;
std::cout << exactGrad << std::endl;
```


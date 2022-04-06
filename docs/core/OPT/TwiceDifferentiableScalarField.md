# TwiceDifferentiableScalarField

> :fontawesome-solid-file-code: core/OPT/ScalarField.h &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; :fontawesome-solid-sitemap: Extends: [DifferentiableScalarField](DifferentiableScalarField.md)

Template class used to represent a scalar field \( f : \mathbb{R}^N \rightarrow \mathbb{R} \) whose hessian function \( H(f) : \mathbb{R}^N \rightarrow \mathbb{R}^{N \times N} \) is known analitically at every point.

``` c++
template <unsigned int N> 
class TwiceDifferentiableScalarField : public DifferentiableScalarField<N> { ... };
```

## Methods

``` c++
TwiceDifferentiableScalarField(std::function<double(SVector<N>)> f_, 
	                           std::function<SVector<N>(SVector<N>)> df_, 
	                           std::function<SMatrix<N>(SVector<N>)> ddf_)
```

!!! quote ""
	Constructor

    | Args      | Description                                                 |
    |------------------------------------------------------------------------------------------|-------------------------------------------------------------|
    | `std::function<double(SVector<N>)> f_` | An `std::function` implementing the analytical expression of the field \( f \) taking an N dimensional point in input and returning a double in output                                                 |
	| `std::function<SVector<N>(SVector<N>)> df_` | An `std::function` implementing the analytical expression of the vector field \( \nabla f \) taking an N dimensional point in input and returning an N dimensional point in output representing the gradient expression                                                 |
	| `std::function<SMatrix<N>(SVector<N>)> ddf_` | An `std::function` implementing the analytical expression of the hessian function \( H(f) \) taking an N dimensional point in input and returning an N dimensional square matrix in output representing the hessian expression                                                 |

``` c++
std::function<SMatrix<N>(SVector<N>)> deriveTwice() const override;
```

!!! quote ""
	Returns the analytical expression of the field's hessian \( H(f) : \mathbb{R}^N \rightarrow \mathbb{R}^{N \times N} \) as a Callable object.

## Examples

Define a scalar field \( f(x,y) = 2x^2 + 2y^2 \) having exact gradient equal to \( \nabla f = [4x, 4y]^T \) and exact hessian matrix given by
$$ H(f) = \begin{bmatrix} 4 & 0 \\ 0 & 4 \end{bmatrix} $$

``` c++ linenums="1"
// define the field analytical expression: 2*x^2 + 2*y^2
std::function<double(SVector<2>)> g = [](SVector<2> x) -> double { 
	return 2*std::pow(x[0],2) + 2*std::pow(x[1],2); 
	};

// define analytical expression of gradient field
std::function<SVector<2>(SVector<2>)> dg = [](SVector<2> x) -> SVector<2> { 
	return SVector<2>({4*x[0], 
                       4*x[1]}); 
	};

// define analytical expression of hessian matrix
std::function<SMatrix<2>(SVector<2>)> ddg = [](SVector<2> x) -> SMatrix<2> { 
	return SMatrix<2>({{4, 0},
	                   {0, 4}});
	};

// define twice differentiable field
TwiceDifferentiableScalarField<2> field(g, dg, ddg);

std::cout << "evaluation of field at point" << std::endl;
std::cout << field.evaluateAtPoint(SVector<2>(4,1)) << std::endl;

// get approximation of hessian at point
SVector<2> hess = field.getHessianApprox(SVector<2>(2,1), 0.001);
    
std::cout << "approximation of gradient at point" << std::endl;
std::cout << hess << std::endl;

// evaluate exact hessian at point
SVector<2> exactHess = field.deriveTwice()(SVector<2>(2,1));

std::cout << "exact hessian at point" << std::endl;
std::cout << exactHess << std::endl;
```

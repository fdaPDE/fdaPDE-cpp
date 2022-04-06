# ExactNewtonOptimizer

> :fontawesome-solid-file-code: core/OPT/ExactNewtonOptimizer.h &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; :fontawesome-solid-sitemap: Extends: [IterativeOptimizer](Optimizer.md#IterativeOptimizer)
	
Template class to optimize a given [TwiceDifferentiableScalarField](TwiceDifferentiableScalarField.md) over \( \mathbb{R}^N \) using the Newton's iterative optimization method:
$$ x_{n+1} = x_{n} - \lambda H_f(x_n)^{-1} \cdot \nabla f(x_n) $$
where \( H_f(x_n) \in \mathbb{R}^{N \times N} \) denotes the Hessian matrix of the field evaluated at point \( x_n \).

``` c++
template <unsigned int N> class ExactNewtonOptimizer : public NewtonOptimizer<N> { ... };
```

!!! note
	The implementation of Newton's method provided by this class relies on the exact analytical expression of gradient and hessian of the field. See parent class [NewtonOptimizer](NewtonOptimizer.md) in case you want to use numerically approximations for these quantities or you have no analytical expressions for them.

!!! info
	At each iteration of the method to avoid the cost of matrix inversion, the quantity \( H_f(x_n)^{-1} \cdot \nabla f(x_n) \) is computed by solving the linear system \( H_f(x_n) z = \nabla f(x_n) \) in \( z \) using eigen's [QR decompostion with column-pivoting](https://eigen.tuxfamily.org/dox/classEigen_1_1ColPivHouseholderQR.html).

## Methods

``` c++
ExactNewtonOptimizer(double step_, const SVector<N>& x0_, unsigned int maxIteration_,
                     double tolerance_, const TwiceDifferentiableScalarField<N>& objective_)
```

!!! quote ""
	Constructor.

    | Args                                                                                          | Description                                                 |
    |------------------------------------------------------------------------------------------|-------------------------------------------------------------|
    | `double step_`             | The term \( \lambda \) in the iterative formulation of the Newton's method. |
    | `const SVector<N>& x0_` | The initial point from which the iterative method is started. |
	| `unsigned int maxIteration_` | The maximum number of iterations allowed. |
	| `double tolerance_` | The tolerance on the error of the obtained solution requested from the method. |
    | `TwiceDifferentiableScalarField<N>& objective_` | The objective function to optimize. |



``` c++
std::pair<SVector<N>, double> findMinimum() const override;
```

!!! quote ""
	Applies the optimization method to the objective passed as argument. Returns `std::pair<SVector<N>,double>` where the first element is the point in \( \mathbb{R}^N \) where the minimum is reached while the second one is the actual minimum value reached by the objective.
	
	!!! info
		Default stopping condition: the method stops either if the \( l^2 \) norm of the gradient of the objective \( \left\lVert \nabla f(x_n) \right\rVert \) reaches the required tolerance or if such tolerance is not reached before a `maxIteration` number of iterations.
	
	!!! tip
	    You can control the optimizer by exploiting the [IterativeOptimizer](Optimizer.md#IterativeOptimizer) interface.


## Examples

The following code finds the minimum of function \( g(x,y) = 2x^2 + x + 2y^2 \) using an exact Newton's method with \( \lambda = 0.01 \) starting from point \( (1,1) \).

``` c++ linenums="1"
// define target to optimize: 2*x^2 + x + 2*y^2
std::function<double(SVector<2>)> g = [](SVector<2> x) -> double { 
	return 2*std::pow(x[0],2) + 2*std::pow(x[1],2) + x[0]; 
  };

// define analytical expression of gradient and hessian of the field
std::function<SVector<2>(SVector<2>)> dg = [](SVector<2> x) -> SVector<2> { 
	return SVector<2>({4*x[0] + 1, 4*x[1]}); 
  };
  
std::function<SMatrix<2>(SVector<2>)> ddg = [](SVector<2> x) -> SMatrix<2> { 
	return SMatrix<2>({{4, 0},
                       {0, 4}});
  };

// define the objective to optimize
TwiceDifferentiableScalarField<2> objective(g, dg, ddg);

// perform newton optimization

double lambda = 0.01;                   // set learning rate
unsigned int max_iterations = 1000;     // set maximum number of iterations
double tolerance = 0.001;               // set tolerance

// define newton optimizer
ExactNewtonOptimizer<2> optNewton2DExact(lambda, SVector<2>(1,1), max_iterations, 
                                         tolerance, objective);
// optimize
std::pair<SVector<2>,double> min_g = optNewton2DExact.findMinimum();
```

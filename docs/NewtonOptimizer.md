# NewtonOptimizer

Template class to optimize a given [ScalarField](ScalarField.md) over \( \mathbb{R}^N \) using the Newton's iterative optimization method:
$$ x_{n+1} = x_{n} - \lambda H_f(x_n)^{-1} \cdot \nabla f(x_n) $$
where \( H_f(x_n) \in \mathbb{R}^{N \times N} \) denotes the Hessian matrix of the field evaluated at point \( x_n \).

``` c++
template <unsigned int N> class NewtonOptimizer{ ... };

```

!!! note
	The implementation of Newton's method provided by this class works by numerically approximating both gradient and hessian of the field at the given point. See [ExactNewtonOptimizer](ExactNewtonOptimizer.md) in case you want to use the analytical expression of these quantites during the iterative process

!!! info
	At each iteration of the method to avoid the cost of matrix inversion, the quantity \( H_f(x_n)^{-1} \cdot \nabla f(x_n) \) is computed by solving the linear system \( H_f(x_n) z = \nabla f(x_n) \) in \( z \) using eigen's [QR decompostion with column-pivoting](https://eigen.tuxfamily.org/dox/classEigen_1_1ColPivHouseholderQR.html).

## Methods

```
NewtonOptimizer(double step_)
```

!!! quote ""
	Constructor initializing some quantities in the internal representation of the class.

    | Args                                                                                          | Description                                                 |
    |------------------------------------------------------------------------------------------|-------------------------------------------------------------|
    | `double step_`             | The term \( \lambda \) in the iterative formulation of the Newton's method. |



``` c++
std::pair<SVector<N>, double> 
findMinimum(const SVector<N>& x0, unsigned int maxIteration, 
            double tolerance, ScalarField<N>& objective) const
```

!!! quote ""
	Applies the optimization method to the objective passed as argument. Returns `std::pair<SVector<N>,double>` where the first element is the point in \( \mathbb{R}^N \) where the minimum is reached while the second one is the actual minimum value reached by the objective.

    | Args      | Description                                                 |
    |------------------------------------------------------------------------------------------|-------------------------------------------------------------|
    | `const SVector<N>& x0` | The initial point from which the iterative method is started. |
	| `unsigned int maxIteration` | The maximum number of iterations allowed. |
	| `double tolerance` | The tolerance on the error of the obtained solution requested from the method. |
    | `ScalarField<N>& objective` | The objective function to optimize. |
	
	!!! note
		The method stops either if the \( l^2 \) norm of the gradient of the objective \( \left\lVert \nabla f(x_n) \right\rVert \) reaches the required tolerance or if such tolerance is not reached before a `maxIteration` number of iterations.

## Examples

The following code finds the minimum of function \( g(x,y) = 2x^2 + x + 2y^2 \) using Newton's method with \( \lambda = 0.01 \) starting from point \( (1,1) \).

``` c++ linenums="1"
// define target to optimize: 2*x^2 + x + 2*y^2
std::function<double(SVector<2>)> g = [](SVector<2> x) -> double { 
	return 2*std::pow(x[0],2) + 2*std::pow(x[1],2) + x[0]; 
  };

// wrap target in a ScalarField object
ScalarField<2> objective(g);

// perform newton optimization
  
// set learning rate
double lambda = 0.01;

// create optimizer
NewtonOptimizer<2> optNewton2D(lambda);

// find minimum of g
unsigned int max_iterations = 1000;
double tolerance = 0.001;
std::pair<SVector<2>,double> min_g = optNewton2D.findMinimum(SVector<2>(1,1), max_iterations, 
                                                             tolerance, objective);
```



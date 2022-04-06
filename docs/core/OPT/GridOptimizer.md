# GridOptimizer

> :fontawesome-solid-file-code: core/OPT/optimizers/Grid.h &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; :fontawesome-solid-sitemap: Extends: -

Template class to optimize a given [ScalarField](ScalarField.md) over an N-dimensional grid of equidistant nodes. 

``` c++
template <unsigned int N> class GridOptimizer{ ... };

```

The class simulates the construction of an N-dimensional rectangle \( [a_1, b_1] \times \ldots \times [a_N, b_N] \) splitted along each interval according to a given grid step (the grid step can possibly differ for each dimension). The search is performed using an exhaustive search over the N-dimensional grid. 

!!! warning
	Since this optimization algorithm applies just an exhaustive search without any heuristic the resulting time complexity is exponential in the dimension of the grid N. The method can be very slow for values of N bigger than 2.

## Methods
``` c++
GridOptimizer(std::array<std::pair<double,double>, N> domain, std::array<double, N> steps, 
	          const ScalarField<N>& objective_)

```

!!! quote ""
	Constructor. It sets the shape of the domain where to search for the optimum of the objective as well as the grid step for each dimension.

    | Args                                                                                          | Description                                                 |
    |------------------------------------------------------------------------------------------|-------------------------------------------------------------|
    | `std::array<std::pair<double,double>,N>` &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | An array of exactly N pairs where pair \( i \) represents infimum and superior limit of the interval \( [a_i, b_i] \) along the \( i \)-th dimension                                                 |
    | `std::array<double,N> steps`             | An array of exactly N doubles where the \( i \)-th element defines the grid step along the \( i \)-th dimension |
		| `const ScalarField<N>& objective_` | The objective function to optimize |

``` c++
std::pair<SVector<N>, double> findMinimum() override;
```

!!! quote ""
	Apply the search strategy to the objective passed as argument. Returns `std::pair<SVector<N>,double> ` where the first element is the point in the grid where the minimum is reached while the second one is the actual minimum value reached by the objective.

    !!! info
		During the search no grid is actually stored in memory making the call at least memory efficient.
	
## Examples

The following code finds the minimum of function \( g(x,y) = 2x^2 - 2y^2 \) over a grid of points defined in \( [1,5] \times [1,5] \subset \mathbb{R}^2 \) with grid step 0.01 for each dimension

``` c++ linenums="1"
// define target to optimize: 2*x^2 - 2*y^2
std::function<double(SVector<2>)> g = [](SVector<2> x) -> double { 
       return 2*std::pow(x[0],2) - 2*std::pow(x[1],2); 
};

// create a scalar field
ScalarField<2> objective(g);

// perform a 2D grid optimization

// set optimization domain: [1,5] x [1,5]
std::array<std::pair<double, double>, 2> domain2D = {
	std::pair<double, double>(1,5), 
	std::pair<double, double>(1,5)
};
  
// set grid step size
std::array<double,2> step2D = {0.01, 0.01};

// create optimizer
GridOptimizer<2> opt2D(domain2D, lambda2D, objective);
  
// find minimum of g
std::pair<array<double, 2>,double> min_g = opt2D.findMinimum();
```
